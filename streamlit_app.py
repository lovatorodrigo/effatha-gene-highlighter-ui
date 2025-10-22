import os, io, re, tempfile
import pandas as pd
import streamlit as st
from effatha_core import (
    Variant, TranscriptSeq, Highlight,
    extract_text_from_pdf_path, parse_report_text,
    fetch_transcript_gb, parse_gb_to_transcript,
    apply_hgvs_c_to_cds, wrap_seq_with_highlights, escape_html
)

# ---------------- UI ----------------
st.set_page_config(page_title="Effatha Gene Highlighter", layout="wide")
st.title("Effatha Gene Highlighter — Streamlit")
st.write("Suba um **PDF** do laudo (ou **TXT**), extraia variantes e visualize o **CDS** com mutações destacadas.")

# -------------- helpers de estado --------------
SS = st.session_state
def _init_state():
    SS.setdefault("uploaded_name", None)
    SS.setdefault("uploaded_bytes", None)
    SS.setdefault("parsed_text", None)
    SS.setdefault("variants", None)            # list[Variant]
    SS.setdefault("df_variants", None)         # pd.DataFrame (sem HTML)
    SS.setdefault("df_variants_html", None)    # pd.DataFrame (HGVS c. com <mark>)
    SS.setdefault("tx_by_gene", {})            # dict[str, TranscriptSeq]
    SS.setdefault("roi_left", 90)              # nt antes
    SS.setdefault("roi_right", 90)             # nt depois
    SS.setdefault("html_report", None)         # str
    SS.setdefault("xlsx_bytes", None)          # bytes
    SS.setdefault("processed", False)

_init_state()

# -------------- sidebar (form) --------------
with st.sidebar:
    st.header("Upload & Opções")

    with st.form("form_upload"):
        up = st.file_uploader("Laudo (PDF ou TXT)", type=["pdf","txt"], accept_multiple_files=False, key="uploader")

        # ROI configurável (antes/depois)
        left_default = SS.get("roi_left", 90)
        right_default = SS.get("roi_right", 90)
        col_l, col_r = st.columns(2)
        with col_l:
            roi_left = st.number_input("nt ANTES (à esquerda)", min_value=1, max_value=100000, value=left_default, step=1, key="roi_left_input")
        with col_r:
            roi_right = st.number_input("nt DEPOIS (à direita)", min_value=1, max_value=100000, value=right_default, step=1, key="roi_right_input")

        fetch_seq = st.checkbox(
            "Baixar transcrito (NCBI)", value=True,
            help="Usa eutils para obter NM_* e destacar mutações na CDS", key="fetch_flag"
        )

        # Preferir env; se não houver secrets, não quebra.
        try:
            default_email = st.secrets.get("NCBI_TOOL_EMAIL", "")
        except Exception:
            default_email = os.environ.get("NCBI_TOOL_EMAIL", "")
        ncbi_email = st.text_input("NCBI email (opcional)", value=default_email, key="ncbi_email")

        go = st.form_submit_button("Processar")

# -------------- caches --------------
@st.cache_data(show_spinner=False)
def _extract_text_from_pdf_bytes(pdf_bytes: bytes) -> str:
    with tempfile.NamedTemporaryFile(delete=False, suffix=".pdf") as tf:
        tf.write(pdf_bytes)
        temp_path = tf.name
    try:
        return extract_text_from_pdf_path(temp_path)
    finally:
        try: os.remove(temp_path)
        except Exception: pass

@st.cache_data(show_spinner=False)
def _parse_text(text: str):
    return parse_report_text(text)

@st.cache_data(show_spinner=True)
def _fetch_tx(nm: str, email: str | None):
    gb = fetch_transcript_gb(nm, email=email)
    return parse_gb_to_transcript(gb)

# -------------- helpers de ROI/arquivos --------------
_snv_re = re.compile(r"^c\.(\d+)[ACGT]>\w", re.IGNORECASE)
_del_re = re.compile(r"^c\.(\d+)_(\d+)del", re.IGNORECASE)
_ins_re = re.compile(r"^c\.(\d+)_(\d+)ins", re.IGNORECASE)
_dup_re = re.compile(r"^c\.(\d+)dup", re.IGNORECASE)

def _hgvs_anchor_wt(hgvs_c: str | None) -> int | None:
    """Âncora (posição CDS 1-based) no WT para centralizar ROI."""
    if not hgvs_c:
        return None
    m = _snv_re.match(hgvs_c)
    if m: return int(m.group(1))
    m = _del_re.match(hgvs_c)
    if m: return int(m.group(1))  # início do del
    m = _ins_re.match(hgvs_c)
    if m: return int(m.group(2))  # inserção após 'end'
    m = _dup_re.match(hgvs_c)
    if m: return int(m.group(1))  # duplicação da base em 'pos'
    return None

def _compute_roi_sequences(v: Variant, tx: TranscriptSeq | None, left: int, right: int) -> tuple[str, str]:
    """
    Retorna (ROI_WT, ROI_MUT). Strings sem quebras.
    Vazio para CNV/intrônica/splice ou quando não há transcrito.
    """
    if tx is None or not v.hgvs_c:
        return "", ""
    # intrônica/splice: sem posição CDS
    if "+" in v.hgvs_c or "-" in v.hgvs_c:
        return "", ""
    try:
        # MUTADO
        mut_seq, hlist = apply_hgvs_c_to_cds(tx.cds_seq, v.hgvs_c)
        if not hlist:
            return "", ""
        center_mut = hlist[0].start
        s_mut = max(1, center_mut - left)
        e_mut = min(len(mut_seq), center_mut + right)
        roi_mut = mut_seq[s_mut-1:e_mut]

        # WT usando âncora derivada do HGVS
        center_wt = _hgvs_anchor_wt(v.hgvs_c)
        if center_wt is None:
            center_wt = min(center_mut, len(tx.cds_seq))
        s_wt = max(1, center_wt - left)
        e_wt = min(len(tx.cds_seq), center_wt + right)
        roi_wt = tx.cds_seq[s_wt-1:e_wt]
        return roi_wt, roi_mut
    except Exception:
        return "", ""

def _roi_fasta_header(gene: str, tx: TranscriptSeq | None, hgvs_c: str, left: int, right: int) -> str:
    acc = tx.accession if tx else "NA"
    return f">{gene}|{acc}|{hgvs_c}|ROI-L{left}_R{right}"

def _seq_with_linebreaks(s: str, width: int = 60) -> str:
    return "\n".join(s[i:i+width] for i in range(0, len(s), width))

def _group_by_gene(variants: list[Variant]) -> dict[str, list[Variant]]:
    by_gene: dict[str, list[Variant]] = {}
    for v in variants:
        by_gene.setdefault(v.gene, []).append(v)
    return by_gene

# ---- NOVO: destacar nt na coluna HGVS c. (HTML) ----
_cell_snv = re.compile(r"^c\.(?P<pos>\d+(?:[+-]\d+)?)(?P<ref>[ACGT])>(?P<alt>[ACGT])$", re.IGNORECASE)
_cell_del = re.compile(r"^c\.(?P<s>\d+)_(?P<e>\d+)del(?P<seq>[ACGT]*)$", re.IGNORECASE)
_cell_ins = re.compile(r"^c\.(?P<s>\d+)_(?P<e>\d+)ins(?P<seq>[ACGT]+)$", re.IGNORECASE)
_cell_dup = re.compile(r"^c\.(?P<pos>\d+)dup(?P<seq>[ACGT]*)$", re.IGNORECASE)

def _mark_hgvs_c_html(hgvs: str | None) -> str:
    """Retorna HGVS c. com <mark> no trecho alterado, em HTML (escape seguro)."""
    if not hgvs:
        return "-"
    s = hgvs.strip()

    m = _cell_snv.match(s)  # inclui intrônica (pos com + ou -)
    if m:
        pos = escape_html(m.group("pos"))
        ref = escape_html(m.group("ref"))
        alt = escape_html(m.group("alt"))
        return f"c.{pos}<mark>{ref}</mark>&gt;<mark>{alt}</mark>"

    m = _cell_del.match(s)
    if m:
        a = escape_html(m.group("s")); b = escape_html(m.group("e"))
        seq = m.group("seq")
        if seq:
            return f"c.{a}_{b}del<mark>{escape_html(seq)}</mark>"
        return f"c.{a}_{b}<mark>del</mark>"

    m = _cell_ins.match(s)
    if m:
        a = escape_html(m.group("s")); b = escape_html(m.group("e"))
        seq = escape_html(m.group("seq"))
        return f"c.{a}_{b}ins<mark>{seq}</mark>"

    m = _cell_dup.match(s)
    if m:
        pos = escape_html(m.group("pos"))
        seq = m.group("seq")
        if seq:
            return f"c.{pos}dup<mark>{escape_html(seq)}</mark>"
        return f"c.{pos}<mark>dup</mark>"

    # fallback: só escapa
    return escape_html(s)

def _build_variants_df(variants: list[Variant],
                       tx_by_gene: dict[str, TranscriptSeq],
                       left: int, right: int) -> pd.DataFrame:
    rows = []
    for v in variants:
        tx = tx_by_gene.get(v.gene)
        roi_wt, roi_mut = _compute_roi_sequences(v, tx, left, right)
        rows.append({
            "gene": v.gene,
            "transcrito": v.transcript,
            "HGVS c.": v.hgvs_c,
            "HGVS p.": v.hgvs_p,
            "VAF%": v.vaf_pct,
            "tier": v.tier,
            "oncogenicidade": v.oncogenicity,
            "seção": v.section,
            "notas": ", ".join(v.notes),
            f"ROI WT (−{left}/+{right} nt)": roi_wt,
            f"ROI MUT (−{left}/+{right} nt)": roi_mut,
        })
    return pd.DataFrame(rows)

def _build_variants_df_html(df_plain: pd.DataFrame) -> pd.DataFrame:
    """Copia do DF com 'HGVS c.' contendo HTML com <mark> para uso no front (HTML) e relatório."""
    df = df_plain.copy()
    col = "HGVS c."
    if col in df.columns:
        df[col] = [ _mark_hgvs_c_html(x) for x in df[col].tolist() ]
    return df

def _build_html_report(df_html: pd.DataFrame,
                       grouped: dict[str, list[Variant]],
                       tx_by_gene: dict[str, TranscriptSeq],
                       left: int, right: int) -> str:
    css = """
    <style>
    body { font-family: -apple-system, Segoe UI, Roboto, Helvetica, Arial, sans-serif; margin: 24px; }
    h1 { margin-bottom: 0.2rem; }
    h2 { margin-top: 2rem; }
    table { border-collapse: collapse; width: 100%; font-size: 14px; }
    th, td { border: 1px solid #ddd; padding: 6px 8px; }
    th { background: #f4f4f4; text-align: left; }
    pre { background: #fafafa; padding: 10px; border: 1px solid #eee; overflow-x: auto; }
    .mut { margin-top: 0.5rem; margin-bottom: 1rem; }
    .caption { color:#666; font-size: 12px; margin: 0.3rem 0 0.8rem; }
    mark { background: #fffd80; }
    </style>
    """
    html = [f"<!DOCTYPE html><html><head><meta charset='utf-8'><title>Relatório — Effatha</title>{css}</head><body>"]
    html.append("<h1>Relatório — Effatha Gene Highlighter</h1>")
    html.append("<p>Resumo de variantes detectadas e CDS com mutações destacadas.</p>")
    html.append(df_html.to_html(index=False, escape=False))

    for gene, vs in grouped.items():
        html.append(f"<h2>{escape_html(gene)}</h2>")
        # Tabela por gene (sem destaque aqui, foco na coluna já destacada acima)
        gdf = pd.DataFrame([{
            "seção": v.section,
            "transcrito": v.transcript,
            "HGVS c.": _mark_hgvs_c_html(v.hgvs_c),
            "HGVS p.": v.hgvs_p,
            "VAF%": v.vaf_pct,
            "tier": v.tier,
            "oncogenicidade": v.oncogenicity,
            "notas": ", ".join(v.notes),
        } for v in vs])
        html.append(gdf.to_html(index=False, escape=False))

        tx = tx_by_gene.get(gene)
        if not tx:
            continue

        # WT combinado
        wt_highs: list[Highlight] = []
        for v in vs:
            if v.hgvs_c and v.hgvs_c.startswith("c.") and ("+" not in v.hgvs_c and "-" not in v.hgvs_c):
                try:
                    _, h = apply_hgvs_c_to_cds(tx.cds_seq, v.hgvs_c)
                    wt_highs.extend(h)
                except Exception:
                    pass
        html.append("<div class='mut'><strong>CDS (WT) — destaques combinados</strong></div>")
        html.append("<pre>" + wrap_seq_with_highlights(tx.cds_seq, wt_highs, line=90) + "</pre>")

        # Por variante (mutada + ROI MUT mostrado; WT/MUT completos vão na planilha/FASTA)
        for v in vs:
            if not v.hgvs_c:
                continue
            html.append(f"<div class='mut'><strong>{_mark_hgvs_c_html(v.hgvs_c)}</strong></div>")
            if "+" in v.hgvs_c or "-" in v.hgvs_c:
                html.append("<div class='caption'>Sem destaque CDS para variantes intrônicas/splice.</div>")
                continue
            try:
                mut_seq, hlist = apply_hgvs_c_to_cds(tx.cds_seq, v.hgvs_c)
                html.append("<pre>" + wrap_seq_with_highlights(mut_seq, hlist, line=90) + "</pre>")
                if hlist and hlist[0].kind in {"SNV","DEL","INS","DUP"}:
                    center = hlist[0].start
                    s = max(1, center - left); e = min(len(mut_seq), center + right)
                    roi_seq = mut_seq[s-1:e]
                    header = f">{gene}|{tx.accession}|{v.hgvs_c}|ROI-L{left}_R{right}"
                    fasta_block = f"{header}\n" + "\n".join(roi_seq[i:i+60] for i in range(0, len(roi_seq), 60))
                    html.append(f"<div class='caption'>ROI −{left}/+{right} nt (FASTA)</div>")
                    html.append("<pre>{}</pre>".format(escape_html(fasta_block)))
            except Exception as ex:
                html.append(f"<div class='caption'>Falha ao aplicar {escape_html(v.hgvs_c)}: {escape_html(str(ex))}</div>")
    html.append("</body></html>")
    return "\n".join(html)

def _build_xlsx_bytes(df_plain: pd.DataFrame) -> bytes:
    """Gera XLSX sem HTML na coluna HGVS (mantemos a versão 'limpa' para Excel)."""
    output = io.BytesIO()
    with pd.ExcelWriter(output, engine="xlsxwriter") as writer:
        sheet = "variantes"
        df_plain.to_excel(writer, sheet_name=sheet, index=False)
        wb = writer.book
        ws = writer.sheets[sheet]
        fmt_header = wb.add_format({"bold": True, "bg_color": "#F4F4F4", "border": 1})
        ws.set_row(0, 18, fmt_header)
        ws.freeze_panes(1, 0)
        ws.autofilter(0, 0, len(df_plain), len(df_plain.columns)-1)
        # larguras básicas
        for col_idx, col_name in enumerate(df_plain.columns):
            width = 12
            if col_name in ("transcrito", "HGVS c.", "HGVS p."): width = 20
            if "ROI" in col_name: width = 60
            if col_name == "oncogenicidade": width = 20
            if col_name == "notas": width = 24
            ws.set_column(col_idx, col_idx, width)
    return output.getvalue()

# -------------- processamento --------------
def _process_now(upload, left: int, right: int, fetch_flag: bool, email: str | None):
    # ler bytes e guardar
    SS.uploaded_name = upload.name
    SS.uploaded_bytes = upload.read()

    # extrair texto
    try:
        if upload.type == "application/pdf" or upload.name.lower().endswith(".pdf"):
            text = _extract_text_from_pdf_bytes(SS.uploaded_bytes)
        else:
            text = SS.uploaded_bytes.decode("utf-8", errors="ignore")
    except Exception as e:
        st.error(f"Falha ao ler arquivo: {e}")
        return

    SS.parsed_text = text

    # parsear variantes
    variants = _parse_text(text)
    if not variants:
        st.warning("Nenhuma variante detectada. Verifique o formato do laudo.")
        SS.variants = None
        SS.df_variants = None
        SS.df_variants_html = None
        SS.tx_by_gene = {}
        SS.html_report = None
        SS.xlsx_bytes = None
        SS.processed = False
        return

    SS.variants = variants
    SS.roi_left = int(left)
    SS.roi_right = int(right)

    # agrupar e buscar transcritos
    grouped = _group_by_gene(variants)
    tx_by_gene: dict[str, TranscriptSeq] = {}
    if fetch_flag:
        for gene, vs in grouped.items():
            tx = next((vv.transcript for vv in vs if vv.transcript), None)
            if tx:
                try:
                    tx_by_gene[gene] = _fetch_tx(tx, email or None)
                except Exception as e:
                    st.info(f"Não foi possível obter {tx}: {e}")
    SS.tx_by_gene = tx_by_gene

    # DF "limpo" e DF com HGVS marcado
    df_plain = _build_variants_df(variants, SS.tx_by_gene, SS.roi_left, SS.roi_right)
    df_html = _build_variants_df_html(df_plain)

    SS.df_variants = df_plain
    SS.df_variants_html = df_html

    # relatório HTML usa a versão com destaque; XLSX usa a versão limpa
    SS.html_report = _build_html_report(SS.df_variants_html, grouped, SS.tx_by_gene, SS.roi_left, SS.roi_right)
    SS.xlsx_bytes = _build_xlsx_bytes(SS.df_variants)

    SS.processed = True

# Acionar processamento
if go and up is not None:
    _process_now(up, roi_left, roi_right, st.session_state.fetch_flag, st.session_state.ncbi_email)

# -------------- renderização --------------
if SS.processed and SS.variants is not None:
    st.subheader("Variantes extraídas")

    # Tabela interativa (sem HTML)
    st.dataframe(SS.df_variants, use_container_width=True)

    # Tabela com destaque (HTML)
    st.markdown("**Tabela com destaque do nt (HGVS c.)**")
    st.markdown(SS.df_variants_html.to_html(index=False, escape=False), unsafe_allow_html=True)

    # downloads do relatório/xlsx
    col_a, col_b = st.columns(2)
    with col_a:
        st.download_button(
            "Baixar relatório (HTML com destaques)",
            data=SS.html_report or "",
            file_name="relatorio_effatha.html",
            mime="text/html",
            key="dl_html_report"
        )
    with col_b:
        st.download_button(
            "Baixar resumo (XLSX)",
            data=SS.xlsx_bytes or b"",
            file_name="variantes_effatha.xlsx",
            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
            key="dl_xlsx"
        )

    # por gene (visual)
    grouped = _group_by_gene(SS.variants)

    for gene, vs in grouped.items():
        st.markdown(f"### {gene}")
        tx = SS.tx_by_gene.get(gene)

        # Tabela do gene — com destaque na coluna HGVS c.
        gdf_plain = pd.DataFrame([{
            "seção": v.section,
            "transcrito": v.transcript,
            "HGVS c.": v.hgvs_c,
            "HGVS p.": v.hgvs_p,
            "VAF%": v.vaf_pct,
            "tier": v.tier,
            "oncogenicidade": v.oncogenicity,
            "notas": ", ".join(v.notes),
        } for v in vs])
        gdf_html = gdf_plain.copy()
        if "HGVS c." in gdf_html.columns:
            gdf_html["HGVS c."] = [ _mark_hgvs_c_html(x) for x in gdf_html["HGVS c."].tolist() ]

        st.dataframe(gdf_plain, use_container_width=True)
        st.markdown(gdf_html.to_html(index=False, escape=False), unsafe_allow_html=True)

        if tx:
            # WT combinado
            wt_highs: list[Highlight] = []
            for v in vs:
                if v.hgvs_c and v.hgvs_c.startswith("c.") and ("+" not in v.hgvs_c and "-" not in v.hgvs_c):
                    try:
                        _, h = apply_hgvs_c_to_cds(tx.cds_seq, v.hgvs_c)
                        wt_highs.extend(h)
                    except Exception:
                        pass
            st.markdown("**CDS (WT) — destaques combinados**")
            st.markdown("<pre>" + wrap_seq_with_highlights(tx.cds_seq, wt_highs, line=90) + "</pre>", unsafe_allow_html=True)

            # Por variante
            with st.expander("Ver por variante (mutada + ROI)", expanded=False):
                for idx, v in enumerate(vs):
                    if not v.hgvs_c:
                        continue
                    st.markdown(f"#### {_mark_hgvs_c_html(v.hgvs_c)}", unsafe_allow_html=True)
                    if "+" in v.hgvs_c or "-" in v.hgvs_c:
                        st.caption("Sem destaque CDS para variantes intrônicas/splice.")
                        continue
                    try:
                        mut_seq, hlist = apply_hgvs_c_to_cds(tx.cds_seq, v.hgvs_c)
                        st.markdown("<pre>" + wrap_seq_with_highlights(mut_seq, hlist, line=90) + "</pre>", unsafe_allow_html=True)

                        if hlist and hlist[0].kind in {"SNV","DEL","INS","DUP"}:
                            center = hlist[0].start
                            s = max(1, center - SS.roi_left); e = min(len(mut_seq), center + SS.roi_right)
                            roi_seq = mut_seq[s-1:e]

                            st.download_button(
                                label=f"Baixar ROI −{SS.roi_left}/+{SS.roi_right} nt ({v.hgvs_c})",
                                data=_roi_fasta_header(gene, tx, v.hgvs_c, SS.roi_left, SS.roi_right) + "\n" + _seq_with_linebreaks(roi_seq, 60),
                                file_name=f"{gene}_{tx.accession}_{re.sub(r'[^A-Za-z0-9_.-]', '_', v.hgvs_c)}_ROI-L{SS.roi_left}_R{SS.roi_right}.fasta",
                                mime="text/plain",
                                key=f"dl_roi_{gene}_{idx}"
                            )
                    except Exception as e:
                        st.warning(f"Falha ao aplicar {v.hgvs_c}: {e}")
        else:
            st.caption("Sem sequência de transcrito (ou fetch desativado). Apenas tabela exibida acima.")
else:
    st.info("Carregue um PDF/TXT e clique em **Processar**.")
