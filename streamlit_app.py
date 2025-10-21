import os, tempfile
import pandas as pd
import streamlit as st
from effatha_core import (
    Variant, TranscriptSeq,
    extract_text_from_pdf_path, parse_report_text,
    fetch_transcript_gb, parse_gb_to_transcript,
    apply_hgvs_c_to_cds, wrap_seq_with_highlights, escape_html
)

st.set_page_config(page_title="Effatha Gene Highlighter", layout="wide")
st.title("Effatha Gene Highlighter — Streamlit")
st.write("Suba um **PDF** do laudo (ou **TXT**), extraia variantes e visualize o **CDS** com mutações destacadas.")

with st.sidebar:
    st.header("Upload & Opções")
    up = st.file_uploader("Laudo (PDF ou TXT)", type=["pdf","txt"], accept_multiple_files=False)
    roi = st.slider("Janela ao redor da mutação (nt)", 20, 200, 90, 10)
    fetch_seq = st.checkbox(
        "Baixar transcrito (NCBI)", value=True,
        help="Usa eutils para obter NM_* e destacar mutações na CDS"
    )

    # Preferir variável de ambiente (Render). Só tentar st.secrets se houver, protegendo com try/except.
    def _default_ncbi_email() -> str:
        env = os.environ.get("NCBI_TOOL_EMAIL")
        if env:
            return env
        try:
            return st.secrets.get("NCBI_TOOL_EMAIL", "")
        except Exception:
            return ""

    ncbi_email = st.text_input("NCBI email (opcional)", value=_default_ncbi_email())
    go = st.button("Processar")

@st.cache_data(show_spinner=False)
def _extract_text_from_pdf_bytes(pdf_bytes: bytes) -> str:
    with tempfile.NamedTemporaryFile(delete=False, suffix=".pdf") as tf:
        tf.write(pdf_bytes)
        temp_path = tf.name
    try:
        return extract_text_from_pdf_path(temp_path)
    finally:
        try:
            os.remove(temp_path)
        except Exception:
            pass

@st.cache_data(show_spinner=False)
def _parse_text(text: str):
    return parse_report_text(text)

@st.cache_data(show_spinner=True)
def _fetch_tx(nm: str, email: str | None):
    gb = fetch_transcript_gb(nm, email=email)
    return parse_gb_to_transcript(gb)

if go and up is not None:
    try:
        if up.type == "application/pdf" or up.name.lower().endswith(".pdf"):
            text = _extract_text_from_pdf_bytes(up.read())
        else:
            text = up.read().decode("utf-8", errors="ignore")
    except Exception as e:
        st.error(f"Falha ao ler arquivo: {e}")
        st.stop()

    variants = _parse_text(text)
    if not variants:
        st.warning("Nenhuma variante detectada. Verifique o formato do laudo.")
        st.stop()

    # Tabela resumida
    df = pd.DataFrame([
        {
            "gene": v.gene,
            "transcrito": v.transcript,
            "HGVS c.": v.hgvs_c,
            "HGVS p.": v.hgvs_p,
            "VAF%": v.vaf_pct,
            "tier": v.tier,
            "oncogenicidade": v.oncogenicity,
            "seção": v.section,
        }
        for v in variants
    ])
    st.subheader("Variantes extraídas")
    st.dataframe(df, use_container_width=True)

    # Agrupar por gene
    by_gene = {}
    for v in variants:
        by_gene.setdefault(v.gene, []).append(v)

    for gene, vs in by_gene.items():
        st.markdown(f"### {gene}")
        tx = next((v.transcript for v in vs if v.transcript), None)
        tx_seq = None
        if fetch_seq and tx:
            try:
                tx_seq = _fetch_tx(tx, ncbi_email or None)
            except Exception as e:
                st.info(f"Não foi possível obter {tx}: {e}")

        # tabela do gene
        gdf = pd.DataFrame([
            {
                "seção": v.section,
                "transcrito": v.transcript,
                "HGVS c.": v.hgvs_c,
                "HGVS p.": v.hgvs_p,
                "VAF%": v.vaf_pct,
                "tier": v.tier,
                "oncogenicidade": v.oncogenicity,
                "notas": ", ".join(v.notes),
            } for v in vs
        ])
        st.dataframe(gdf, use_container_width=True)

        if tx_seq:
            # Destaques na WT
            wt_highs = []
            for v in vs:
                if v.hgvs_c and v.hgvs_c.startswith("c.") and ("+" not in v.hgvs_c and "-" not in v.hgvs_c):
                    try:
                        _, h = apply_hgvs_c_to_cds(tx_seq.cds_seq, v.hgvs_c)
                        wt_highs.extend(h)
                    except Exception:
                        pass
            st.markdown("**CDS (WT) — destaques combinados**")
            st.markdown("<pre>" + wrap_seq_with_highlights(tx_seq.cds_seq, wt_highs, line=90) + "</pre>", unsafe_allow_html=True)

            # Por variante
            with st.expander("Ver por variante (mutada + ROI)", expanded=False):
                for v in vs:
                    if not v.hgvs_c:
                        continue
                    st.markdown(f"#### {escape_html(v.hgvs_c)}", unsafe_allow_html=True)
                    if "+" in (v.hgvs_c or "") or "-" in (v.hgvs_c or ""):
                        st.caption("Sem destaque CDS para variantes intrônicas/splice.")
                        continue
                    try:
                        mut_seq, hlist = apply_hgvs_c_to_cds(tx_seq.cds_seq, v.hgvs_c)
                        st.markdown("<pre>" + wrap_seq_with_highlights(mut_seq, hlist, line=90) + "</pre>", unsafe_allow_html=True)
                        if hlist and hlist[0].kind in {"SNV","DEL","INS","DUP"}:
                            center = hlist[0].start
                            s = max(1, center - roi); e = min(len(mut_seq), center + roi)
                            roi_seq = mut_seq[s-1:e]
                            st.download_button(
                                label=f"Baixar ROI ±{roi} nt ({v.hgvs_c})",
                                data=f">{gene}|{tx_seq.accession}|{v.hgvs_c}|ROI{roi}\n" + "\n".join(roi_seq[i:i+60] for i in range(0, len(roi_seq), 60)),
                                file_name=f"{gene}_{tx_seq.accession}_{v.hgvs_c}_ROI{roi}.fasta",
                                mime="text/plain",
                            )
                    except Exception as e:
                        st.warning(f"Falha ao aplicar {v.hgvs_c}: {e}")
        else:
            st.caption("Sem sequência de transcrito (ou fetch desativado). Apenas tabela exibida acima.")
else:
    st.info("Carregue um PDF/TXT e clique em **Processar**.")
