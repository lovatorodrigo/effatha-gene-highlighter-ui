#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Effatha Gene Highlighter — v0.6

Novidades (v0.6):
- Bugfix: highlight correto para `c.Xdup` sem sequência explícita (marca a base inserida).
- Flag `--open`: abre o index.html no navegador ao final.
- Fallback de encoding ao ler TXT (UTF-8 → Latin-1).
- Relatórios extras:
  - variants_extracted.csv (tabela para Excel/Sheets)
  - report_all.html (um único HTML com todos os genes/variantes e destaques)
  - `--bundle` gera effatha_report.zip com tudo dentro (HTMLs, FASTAs, CSV, JSON).

Resumo do objetivo:
- Ler laudo (PDF ou TXT), extrair variantes (gene, NM_, HGVS c./p., VAF, categoria),
  baixar transcrito RefSeq (NM_*) no NCBI, obter CDS e gerar HTML com CDS completa
  e mutações destacadas + FASTAs WT/MUT/ROI por variante.
"""

from __future__ import annotations
import argparse
import re
import os
import json
import csv
import zipfile
import datetime as _dt
from dataclasses import dataclass, asdict
from typing import List, Dict, Optional, Tuple

# Importações preguiçosas para lidar com ambientes sem pdfminer
try:
    from pdfminer.high_level import extract_text as _pdfminer_extract_text  # type: ignore
except Exception:
    _pdfminer_extract_text = None  # type: ignore

try:
    import PyPDF2  # type: ignore
except Exception:
    PyPDF2 = None  # type: ignore

import requests  # type: ignore
from Bio import SeqIO  # type: ignore
from io import StringIO

# GUI (opcional) para seleção de arquivo (usada com --pick ou auto no PyCharm)
try:
    import tkinter as _tk  # type: ignore
    from tkinter import filedialog as _fd  # type: ignore
except Exception:
    _tk = None
    _fd = None

# -------------------------------
# Data structures
# -------------------------------
@dataclass
class Variant:
    gene: str
    transcript: Optional[str]
    hgvs_c: Optional[str]
    hgvs_p: Optional[str]
    vaf_pct: Optional[float]
    tier: Optional[str]
    oncogenicity: Optional[str]
    section: str  # 'oncogenic', 'cnv', 'vus'
    notes: List[str]

@dataclass
class TranscriptSeq:
    accession: str
    mrna_seq: str
    cds_range: Tuple[int, int]  # 1-based inclusive coordinates within mRNA
    cds_seq: str

@dataclass
class Highlight:
    start: int  # 1-based position in CDS
    end: int    # inclusive
    kind: str   # 'SNV'|'DEL'|'INS'|'DUP'|'OTHER'
    label: str  # e.g., c.725G>T

# -------------------------------
# Leitura de PDF (com fallback) e parser de TEXTO do laudo
# -------------------------------

def extract_text_from_pdf(pdf_path: str) -> str:
    """Extrai texto do PDF usando pdfminer.six se disponível; senão tenta PyPDF2.
    Lança um erro amigável se nenhuma opção estiver instalada.
    """
    if _pdfminer_extract_text is not None:
        return _pdfminer_extract_text(pdf_path)
    if PyPDF2 is not None:
        txt_parts: List[str] = []
        with open(pdf_path, "rb") as fh:
            reader = PyPDF2.PdfReader(fh)
            for page in reader.pages:
                t = page.extract_text() or ""
                txt_parts.append(t)
        return "\n".join(txt_parts)
    raise ModuleNotFoundError(
        "Nenhum extrator de PDF disponível. Instale 'pdfminer.six' OU 'PyPDF2', ou use --report_txt com texto do laudo."
    )

# --- Regexes para o parser de TEXTO ---
_gsym = r"[A-Z0-9]{2,}"
_nm = r"NM_\d+(?:\.\d+)?"
_c_hgvs_pat = r"c\.[0-9_]+(?:[+-]\d+)?(?:[ACGT]>[ACGT]|del(?:[ACGT]*)?|ins[ACGT]+|dup[ACGT]*)"
_p_hgvs_pat = (
    r"p\.[A-Za-z]{1}[a-z]{2}\d+(?:[A-Za-z]{1}[a-z]{2}|\*)?\b|"
    r"p\.[A-Za-z]{1}[a-z]{2}\d+[A-Za-z]{1}[a-z]{2}fs\*\d+|"
    r"p\.[A-Za-z]{1}[a-z]{2}=|p\.[A-Za-z]{1}[a-z]{2}\d+del"
)

_re_snvindel = re.compile(
    rf"(?P<gene>{_gsym})\s*\((?P<tx>{_nm})\)\s*(?P<c>{_c_hgvs_pat})\s*(?P<p>{_p_hgvs_pat})?\s*(?P<vaf>\d{{1,3}}[\.,]\d%)?",
    re.IGNORECASE,
)
_re_cnv = re.compile(rf"(?P<gene>{_gsym})\s*\((?P<tx>{_nm})\)\s*Amplifica[cç][aã]o", re.IGNORECASE)
_re_vus_block = re.compile(
    rf"(?P<gene>{_gsym})\s*\((?P<tx>{_nm})\)\s*(?P<c>{_c_hgvs_pat})\s*(?P<p>{_p_hgvs_pat})?\s*(?P<vaf>\d{{1,3}}[\.,]\d%)\s*VUS",
    re.IGNORECASE,
)

_re_onco_labels = {
    "oncogenic": re.compile(r"Oncog[eê]nica", re.IGNORECASE),
    "prob_oncogenic": re.compile(r"Provavelmente\s+Oncog[eê]nica", re.IGNORECASE),
    "tier2": re.compile(r"\(Tier\s*2\)", re.IGNORECASE),
    "tier3": re.compile(r"\(Tier\s*3\)", re.IGNORECASE),
}


def parse_report_text(text: str) -> List[Variant]:
    lines = [l.strip() for l in (text or "").splitlines() if l.strip()]
    blob = "\n".join(lines)
    variants: List[Variant] = []

    # SNV/Indel (Oncogênicas/Tier 2)
    for m in _re_snvindel.finditer(blob):
        vaf_raw = m.group("vaf")
        vaf = None
        if vaf_raw:
            vaf = float(vaf_raw.replace("%", "").replace(",", "."))
        window = blob[max(0, m.start()-80): m.end()+80]
        onc = "Oncogênica" if _re_onco_labels["oncogenic"].search(window) else (
              "Provavelmente oncogênica" if _re_onco_labels["prob_oncogenic"].search(window) else None)
        tier = "Tier 2" if _re_onco_labels["tier2"].search(window) else None
        variants.append(Variant(
            gene=m.group("gene"), transcript=m.group("tx"), hgvs_c=m.group("c"),
            hgvs_p=(m.group("p") or None), vaf_pct=vaf, tier=tier, oncogenicity=onc,
            section="oncogenic", notes=[],
        ))

    # CNV
    for m in _re_cnv.finditer(blob):
        variants.append(Variant(
            gene=m.group("gene"), transcript=m.group("tx"), hgvs_c=None, hgvs_p=None,
            vaf_pct=None, tier="Tier 2", oncogenicity="Oncogênica", section="cnv", notes=["amplificacao"],
        ))

    # VUS
    for m in _re_vus_block.finditer(blob):
        vaf = float(m.group("vaf").replace("%", "").replace(",", "."))
        variants.append(Variant(
            gene=m.group("gene"), transcript=m.group("tx"), hgvs_c=m.group("c"),
            hgvs_p=(m.group("p") or None), vaf_pct=vaf, tier="Tier 3", oncogenicity="VUS", section="vus", notes=[],
        ))

    # Dedup best-effort
    uniq: Dict[Tuple, Variant] = {}
    for v in variants:
        key = (v.gene, v.transcript, v.hgvs_c, v.section)
        if key not in uniq:
            uniq[key] = v
    return list(uniq.values())


def parse_pdf_variants(pdf_path: str) -> List[Variant]:
    text = extract_text_from_pdf(pdf_path)
    return parse_report_text(text)

# -------------------------------
# NCBI fetchers
# -------------------------------
ENTREZ_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


def fetch_transcript_gb(nm_accession: str, email: Optional[str] = None) -> str:
    params = {"db": "nuccore", "id": nm_accession, "rettype": "gb", "retmode": "text", "tool": "effatha_gene_highlighter"}
    if email:
        params["email"] = email
    r = requests.get(ENTREZ_EUTILS, params=params, timeout=30)
    r.raise_for_status()
    return r.text


def parse_gb_to_transcript(gb_text: str) -> TranscriptSeq:
    handle = StringIO(gb_text)
    rec = next(SeqIO.parse(handle, "genbank"))
    acc = rec.id
    cds_feature = None
    for feat in rec.features:
        if feat.type == "CDS":
            cds_feature = feat
            break
    if cds_feature is None:
        raise ValueError(f"Sem feature CDS em {acc}")
    cds_start = int(cds_feature.location.start) + 1
    cds_end = int(cds_feature.location.end)
    mrna_seq = str(rec.seq).upper()
    cds_seq = mrna_seq[cds_start-1:cds_end]
    return TranscriptSeq(accession=acc, mrna_seq=mrna_seq, cds_range=(cds_start, cds_end), cds_seq=cds_seq)

# -------------------------------
# HGVS c. applier (CDS space)
# -------------------------------

_re_snv = re.compile(r"c\.(?P<pos>\d+)(?P<ref>[ACGT])>(?P<alt>[ACGT])$")
_re_del = re.compile(r"c\.(?P<start>\d+)_(?P<end>\d+)del(?P<delseq>[ACGT]*)$")
_re_ins = re.compile(r"c\.(?P<start>\d+)_(?P<end>\d+)ins(?P<insseq>[ACGT]+)$")
_re_dup = re.compile(r"c\.(?P<pos>\d+)dup(?P<dupseq>[ACGT]*)$")


def apply_hgvs_c_to_cds(cds_seq: str, hgvs_c: str) -> Tuple[str, List[Highlight]]:
    s = cds_seq
    hs: List[Highlight] = []

    m = _re_snv.match(hgvs_c)
    if m:
        pos = int(m.group("pos"))
        ref, alt = m.group("ref"), m.group("alt")
        idx = pos - 1
        if idx < 0 or idx >= len(s):
            raise ValueError(f"Posição fora do CDS: {hgvs_c}")
        mutated = s[:idx] + alt + s[idx+1:]
        hs.append(Highlight(start=pos, end=pos, kind="SNV", label=hgvs_c))
        return mutated, hs

    m = _re_del.match(hgvs_c)
    if m:
        start = int(m.group("start"))
        end = int(m.group("end"))
        if start > end:
            start, end = end, start
        idx0, idx1 = start-1, end
        mutated = s[:idx0] + s[idx1:]
        hs.append(Highlight(start=start, end=end, kind="DEL", label=hgvs_c))
        return mutated, hs

    m = _re_ins.match(hgvs_c)
    if m:
        start = int(m.group("start"))
        end = int(m.group("end"))
        ins = m.group("insseq").upper()
        idx = end
        mutated = s[:idx] + ins + s[idx:]
        hs.append(Highlight(start=end+1, end=end+len(ins), kind="INS", label=hgvs_c))
        return mutated, hs

    m = _re_dup.match(hgvs_c)
    if m:
        pos = int(m.group("pos"))
        dupseq = m.group("dupseq")
        idx = pos
        if dupseq:
            mutated = s[:idx] + dupseq + s[idx:]
            end = pos + len(dupseq)
        else:
            base = s[pos-1]
            mutated = s[:idx] + base + s[idx:]
            end = pos + 1  # highlight cobre a base duplicada recém inserida
        hs.append(Highlight(start=pos+1, end=end, kind="DUP", label=hgvs_c))
        return mutated, hs

    return s, [Highlight(start=1, end=1, kind="OTHER", label=f"não aplicado: {hgvs_c}")]

# -------------------------------
# Render HTML com destaques
# -------------------------------

def wrap_seq_with_highlights(seq: str, highlights: List[Highlight], line=60) -> str:
    marks = []
    for h in highlights:
        marks.append((h.start, h.end, h))
    marks.sort(key=lambda x: x[0])

    html_parts = []
    i = 1
    for (start, end, h) in marks:
        if start < 1:
            start = 1
        if end > len(seq):
            end = len(seq)
        if i <= len(seq) and i < start:
            html_parts.append(escape_html(seq[i-1:start-1]))
        if start <= end:
            html_parts.append(f'<mark title="{escape_html(h.label)}">{escape_html(seq[start-1:end])}</mark>')
        i = end + 1
    if i <= len(seq):
        html_parts.append(escape_html(seq[i-1:]))

    merged = "".join(html_parts)
    out = []
    count = 0
    buf = []
    for ch in merged:
        buf.append(ch)
        if ch != '<':
            count += 1
        if count >= line and ch not in ('<', '>'):
            buf.append("<br>")
            count = 0
    out_html = "".join(buf)
    return out_html


def escape_html(s: str) -> str:
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")

# -------------------------------
# Orquestração / Relatórios
# -------------------------------

def ensure_dir(p: str):
    os.makedirs(p, exist_ok=True)


def write_fasta(path: str, header: str, seq: str):
    with open(path, "w", encoding="utf-8") as f:
        f.write(f">{header}\n")
        for i in range(0, len(seq), 60):
            f.write(seq[i:i+60] + "\n")


def resolve_outdir(user_out: Optional[str]) -> str:
    if user_out:
        return user_out
    env_out = os.environ.get("EFFATHA_GENE_OUT")
    if env_out:
        return env_out
    ts = _dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    return os.path.join("runs", f"highlight_{ts}")


def ensure_unique_outdir(base_out: str) -> str:
    if not os.path.exists(base_out):
        os.makedirs(base_out, exist_ok=True)
        return base_out
    try:
        if os.path.isdir(base_out) and len(os.listdir(base_out)) == 0:
            return base_out
    except FileNotFoundError:
        os.makedirs(base_out, exist_ok=True)
        return base_out
    ts = _dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    child = os.path.join(base_out, f"run_{ts}")
    os.makedirs(child, exist_ok=True)
    return child


def write_csv_summary(variants: List[Variant], outdir: str) -> str:
    path = os.path.join(outdir, "variants_extracted.csv")
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["gene", "transcript", "hgvs_c", "hgvs_p", "vaf_pct", "tier", "oncogenicity", "section", "notes"])
        for v in variants:
            w.writerow([v.gene, v.transcript or "", v.hgvs_c or "", v.hgvs_p or "",
                        v.vaf_pct if v.vaf_pct is not None else "",
                        v.tier or "", v.oncogenicity or "", v.section, "; ".join(v.notes)])
    return path


def build_gene_pages(variants: List[Variant], outdir: str, roi_window: int = 60, ncbi_email: Optional[str] = None) -> Dict[str, str]:
    """
    Constrói páginas por gene e retorna um dict {gene: html_string_completa} para gerar report_all.html depois.
    """
    by_gene: Dict[str, List[Variant]] = {}
    for v in variants:
        by_gene.setdefault(v.gene, []).append(v)

    index_items = []
    gene_html_map: Dict[str, str] = {}

    for gene, vs in by_gene.items():
        gene_dir = os.path.join(outdir, gene)
        ensure_dir(gene_dir)

        tx = next((v.transcript for v in vs if v.transcript), None)
        tx_seq: Optional[TranscriptSeq] = None
        if tx:
            try:
                gb = fetch_transcript_gb(tx, email=ncbi_email or os.environ.get("NCBI_TOOL_EMAIL"))
                tx_seq = parse_gb_to_transcript(gb)
                write_fasta(os.path.join(gene_dir, f"{tx}_CDS_WT.fasta"), f"{gene}|{tx}|CDS_WT", tx_seq.cds_seq)
            except Exception as e:
                tx_seq = None
                with open(os.path.join(gene_dir, "_errors.txt"), "a", encoding="utf-8") as ef:
                    ef.write(f"Falha ao obter {tx}: {e}\n")

        html = [
            "<html><head><meta charset='utf-8'><style>body{font-family:ui-monospace,monospace} mark{background:#fffd80} .chip{display:inline-block;padding:2px 8px;border-radius:12px;background:#eee;margin-right:6px}</style></head><body>",
            f"<h2 id='{escape_html(gene)}'>{escape_html(gene)}</h2>",
        ]
        html.append("<h3>Variantes</h3><table border='1' cellspacing='0' cellpadding='4'>")
        html.append("<tr><th>Seção</th><th>Transcrito</th><th>HGVS c.</th><th>HGVS p.</th><th>VAF%</th><th>Tier</th><th>Oncogenicidade</th><th>Notas</th></tr>")
        for v in vs:
            html.append("<tr>" + "".join([
                f"<td>{escape_html(v.section)}</td>",
                f"<td>{escape_html(v.transcript or '-') }</td>",
                f"<td>{escape_html(v.hgvs_c or '-') }</td>",
                f"<td>{escape_html(v.hgvs_p or '-') }</td>",
                f"<td>{'' if v.vaf_pct is None else f'{v.vaf_pct:.1f}'}</td>",
                f"<td>{escape_html(v.tier or '-') }</td>",
                f"<td>{escape_html(v.oncogenicity or '-') }</td>",
                f"<td>{escape_html(', '.join(v.notes))}</td>",
            ]) + "</tr>")
        html.append("</table>")

        if tx_seq:
            wt_highs: List[Highlight] = []
            for v in vs:
                if v.hgvs_c and re.match(r"^c\.\d", v.hgvs_c) and not re.search(r"[+-]", v.hgvs_c):
                    try:
                        _, hlist = apply_hgvs_c_to_cds(tx_seq.cds_seq, v.hgvs_c)
                        wt_highs.extend(hlist)
                    except Exception as e:
                        with open(os.path.join(gene_dir, "_errors.txt"), "a", encoding="utf-8") as ef:
                            ef.write(f"Falha destaque WT {gene} {v.hgvs_c}: {e}\n")
            html.append("<h3>CDS (WT) com destaques</h3><pre>")
            html.append(wrap_seq_with_highlights(tx_seq.cds_seq, wt_highs, line=90))
            html.append("</pre>")

            html.append("<h3>Por Variante</h3>")
            for v in vs:
                if not v.hgvs_c:
                    continue
                if re.search(r"[+-]", v.hgvs_c or ""):
                    html.append(f"<div class='chip'>Sem destaque CDS para {escape_html(v.hgvs_c)} (provável intrônica/splice)</div>")
                    continue
                try:
                    mut_seq, hlist = apply_hgvs_c_to_cds(tx_seq.cds_seq, v.hgvs_c)
                    base = re.sub(r"[^A-Za-z0-9_.-]", "_", f"{gene}_{tx_seq.accession}_{v.hgvs_c}")
                    write_fasta(os.path.join(gene_dir, base + "__CDS_MUT.fasta"), f"{gene}|{tx_seq.accession}|{v.hgvs_c}|CDS_MUT", mut_seq)
                    if hlist and hlist[0].kind in {"SNV", "DEL", "INS", "DUP"}:
                        h = hlist[0]
                        center = h.start
                        roi_start = max(1, center - roi_window)
                        roi_end = min(len(mut_seq), center + roi_window)
                        roi_seq = mut_seq[roi_start-1:roi_end]
                        write_fasta(
                            os.path.join(gene_dir, base + f"__ROI_{roi_window}nt.fasta"),
                            f"{gene}|{tx_seq.accession}|{v.hgvs_c}|ROI{roi_window}",
                            roi_seq,
                        )
                    html.append(f"<h4>{escape_html(v.hgvs_c)}</h4><pre>")
                    html.append(wrap_seq_with_highlights(mut_seq, hlist, line=90))
                    html.append("</pre>")
                except Exception as e:
                    with open(os.path.join(gene_dir, "_errors.txt"), "a", encoding="utf-8") as ef:
                        ef.write(f"Falha MUT {gene} {v.hgvs_c}: {e}\n")
        else:
            html.append("<p><em>Sem sequência de transcrito recuperada (ver _errors.txt)</em></p>")

        cnv_any = any(v.section == "cnv" for v in vs)
        if cnv_any:
            html.append("<h3>Eventos de CNV</h3><p>Amplificação reportada neste laudo. Não há destaque nucleotídico; considerar toda a região gênica.</p>")

        html.append("</body></html>")
        page_html = "\n".join(html)

        path_html = os.path.join(outdir, f"{gene}.html")
        with open(path_html, "w", encoding="utf-8") as f:
            f.write(page_html)

        index_items.append((gene, os.path.relpath(path_html, outdir)))
        gene_html_map[gene] = page_html

    idx = ["<html><head><meta charset='utf-8'><style>body{font-family:system-ui} a{display:block;margin:6px 0}</style></head><body>",
           "<h2>Effatha Gene Highlighter — Relatório</h2>"]
    for gene, rel in sorted(index_items):
        idx.append(f"<a href='{rel}'>{gene}</a>")
    idx.append("</body></html>")
    with open(os.path.join(outdir, "index.html"), "w", encoding="utf-8") as f:
        f.write("\n".join(idx))

    return gene_html_map


def write_master_report(gene_html_map: Dict[str, str], outdir: str) -> str:
    """
    Gera um único HTML (report_all.html) contendo todas as seções por gene.
    """
    parts = [
        "<html><head><meta charset='utf-8'>",
        "<style>body{font-family:system-ui} mark{background:#fffd80} pre{font-family:ui-monospace,monospace}</style>",
        "</head><body>",
        "<h1>Effatha Gene Highlighter — Relatório completo</h1>",
        "<p>Este arquivo consolida todos os genes/variantes com destaques e tabelas.</p>",
    ]
    for gene in sorted(gene_html_map.keys()):
        # extrai apenas o <body> da página do gene (ingênuo, mas suficiente)
        html = gene_html_map[gene]
        body = html
        m0 = re.search(r"<body[^>]*>(.*)</body>", html, flags=re.DOTALL | re.IGNORECASE)
        if m0:
            body = m0.group(1)
        parts.append(f"<hr><a id='{escape_html(gene)}'></a>")
        parts.append(body)
    parts.append("</body></html>")
    out = os.path.join(outdir, "report_all.html")
    with open(out, "w", encoding="utf-8") as f:
        f.write("\n".join(parts))
    return out


def make_zip_bundle(outdir: str, zip_name: str = "effatha_report.zip") -> str:
    """
    Empacota toda a pasta de saída num ZIP (exclui o próprio ZIP se já existir).
    """
    zip_path = os.path.join(outdir, zip_name)
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
        for dirpath, dirnames, filenames in os.walk(outdir):
            for fn in filenames:
                abspath = os.path.join(dirpath, fn)
                # Evitar incluir o zip dentro dele mesmo
                if os.path.abspath(abspath) == os.path.abspath(zip_path):
                    continue
                rel = os.path.relpath(abspath, outdir)
                zf.write(abspath, arcname=rel)
    return zip_path

# -------------------------------
# Resolução da fonte de entrada (para testes e para a IDE)
# -------------------------------

def _resolve_input_source(report: Optional[str], report_txt: Optional[str], pick_flag: bool) -> Tuple[str, str]:
    """Retorna (kind, path) onde kind é 'pdf' ou 'txt'. Pode abrir seletor se pick_flag True
    ou se estiver sob PyCharm (env PYCHARM_HOSTED=1) e Tkinter disponível.
    Também considera as variáveis de ambiente EFFATHA_DEFAULT_REPORT(_TXT).
    """
    # Args explícitos
    if report_txt:
        return ("txt", report_txt)
    if report:
        return ("pdf", report)

    # Defaults de ambiente
    env_txt = os.environ.get("EFFATHA_DEFAULT_REPORT_TXT")
    if env_txt and os.path.isfile(env_txt):
        return ("txt", env_txt)
    env_pdf = os.environ.get("EFFATHA_DEFAULT_REPORT")
    if env_pdf and os.path.isfile(env_pdf):
        return ("pdf", env_pdf)

    # Auto-picker na IDE se disponível
    pycharm = os.environ.get("PYCHARM_HOSTED") == "1"
    should_pick = pick_flag or pycharm
    if should_pick and _fd is not None:
        try:
            root = _tk.Tk(); root.withdraw()
            path = _fd.askopenfilename(title="Selecione o laudo (PDF ou TXT)", filetypes=[("PDF", "*.pdf"), ("Texto", "*.txt"), ("Todos", "*.*")])
            root.update(); root.destroy()
        except Exception:
            path = ""
        if not path:
            raise SystemExit("Nenhum arquivo selecionado. Forneça --report/--report_txt, use --pick ou defina EFFATHA_DEFAULT_REPORT(_TXT).")
        ext = os.path.splitext(path)[1].lower()
        return (("txt" if ext == ".txt" else "pdf"), path)

    # Nada resolvido
    raise SystemExit("Forneça --report (PDF) ou --report_txt (texto) — use --pick para selecionar ou defina EFFATHA_DEFAULT_REPORT(_TXT).")

# -------------------------------
# Self-tests (mínimos) — rodam com --selftest
# -------------------------------

def _selftest_apply_hgvs():
    cds = "ATGCCCAAAGGGTTT"
    mut, hs = apply_hgvs_c_to_cds(cds, "c.4C>T")
    assert mut[3] == 'T' and hs[0].kind == 'SNV'
    mut, hs = apply_hgvs_c_to_cds(cds, "c.2_4del")
    assert len(mut) == len(cds) - 3 and hs[0].kind == 'DEL'
    mut, hs = apply_hgvs_c_to_cds(cds, "c.5_5insAAA")
    assert mut[5:8] == 'AAA' and hs[0].kind == 'INS'
    mut, hs = apply_hgvs_c_to_cds(cds, "c.6dup")
    assert len(mut) == len(cds) + 1 and hs[0].kind == 'DUP'
    # novo: garantir que o highlight exista e cubra a base inserida (pos+1)
    assert hs[0].start == 7 and hs[0].end == 7


def _selftest_parser_text():
    fake = (
        "TP53 (NM_000546.6) c.725G>T p.Cys242Phe 76,7% (Tier 2) Provavelmente Oncogênica\n"
        "EGFR (NM_005228.5) Amplificação\n"
        "RB1 (NM_000321.3) c.228_255del p.Thr77Glufs*25 59,0% (Tier 2) Oncogênica\n"
        "PRSS8 (NM_002773.4) c.104-5G>A p.? 12,3% VUS\n"
    )
    vars = parse_report_text(fake)
    genes = sorted({v.gene for v in vars})
    assert set(genes) >= {"TP53", "EGFR", "RB1", "PRSS8"}
    tp53 = [v for v in vars if v.gene == "TP53"][0]
    assert tp53.transcript and tp53.hgvs_c == "c.725G>T"
    egfr = [v for v in vars if v.gene == "EGFR"][0]
    assert egfr.section == "cnv"
    prss8 = [v for v in vars if v.gene == "PRSS8"][0]
    assert prss8.section == "vus" and abs((prss8.vaf_pct or 0) - 12.3) < 1e-6


def _selftest_outdir_resolver():
    assert resolve_outdir("meu_out") == "meu_out"
    old_env = os.environ.get("EFFATHA_GENE_OUT")
    os.environ["EFFATHA_GENE_OUT"] = "out_env"
    assert resolve_outdir(None) == "out_env"
    if old_env is not None:
        os.environ["EFFATHA_GENE_OUT"] = old_env
    else:
        try:
            del os.environ["EFFATHA_GENE_OUT"]
        except KeyError:
            pass
    auto = resolve_outdir(None)
    assert auto.startswith(os.path.join("runs", "highlight_"))


def _selftest_unique_outdir():
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        base = os.path.join(tmp, "out")
        p1 = ensure_unique_outdir(base)
        assert p1 == base and os.path.isdir(p1)
        open(os.path.join(p1, "index.html"), "w").write("x")
        p2 = ensure_unique_outdir(base)
        assert p2 != p1 and p2.startswith(os.path.join(base, "run_"))
        assert os.path.isdir(p2)


def _selftest_wrap_highlight():
    seq = "ATGCCCAAAGGGTTT"
    _, hs = apply_hgvs_c_to_cds(seq, "c.4C>T")
    html = wrap_seq_with_highlights(seq, hs, line=10)
    assert "<mark" in html and "</mark>" in html


def _selftest_resolve_input_env():
    import tempfile
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as tf:
        tf.write("TP53 (NM_000546.6) c.725G>T\n")
        temp_path = tf.name
    old_txt = os.environ.get("EFFATHA_DEFAULT_REPORT_TXT")
    os.environ["EFFATHA_DEFAULT_REPORT_TXT"] = temp_path
    kind, path = _resolve_input_source(None, None, False)
    assert kind == "txt" and path == temp_path
    # cleanup env
    if old_txt is not None:
        os.environ["EFFATHA_DEFAULT_REPORT_TXT"] = old_txt
    else:
        try:
            del os.environ["EFFATHA_DEFAULT_REPORT_TXT"]
        except KeyError:
            pass


def run_selftests() -> bool:
    try:
        _selftest_apply_hgvs()
        _selftest_parser_text()
        _selftest_outdir_resolver()
        _selftest_unique_outdir()
        _selftest_wrap_highlight()
        _selftest_resolve_input_env()
        print("[SELFTEST] OK — todos os testes passaram.")
        return True
    except AssertionError as e:
        print("[SELFTEST] FALHOU:", e)
        return False

# -------------------------------
# CLI
# -------------------------------

def main():
    ap = argparse.ArgumentParser(description="Effatha Gene Highlighter")
    ap.add_argument("--report", help="Caminho para o PDF do laudo")
    ap.add_argument("--report_txt", help="Caminho para um TXT com o conteúdo do laudo (alternativa ao PDF)")
    ap.add_argument("--out", help="Diretório de saída (opcional). Se omitido, usa $EFFATHA_GENE_OUT ou runs/highlight_YYYYmmdd_HHMMSS")
    ap.add_argument("--window", type=int, default=60, help="Meia-janela (nt) para ROI ao redor da mutação")
    ap.add_argument("--ncbi_email", help="E-mail para identificar o cliente na NCBI eutils (ou use env NCBI_TOOL_EMAIL)")
    ap.add_argument("--pick", action="store_true", help="Abrir seletor de arquivo (PDF/TXT) se nenhum report for informado")
    ap.add_argument("--selftest", action="store_true", help="Executa testes internos e sai")
    ap.add_argument("--open", action="store_true", help="Abrir o index.html no navegador ao final")
    ap.add_argument("--bundle", action="store_true", help="Gerar ZIP (effatha_report.zip) com todos os arquivos de saída")
    args = ap.parse_args()

    if args.selftest:
        run_selftests()
        return

    base_out = resolve_outdir(args.out)
    outdir = ensure_unique_outdir(base_out)
    ensure_dir(outdir)
    print(f"[0/4] Diretório base: {base_out}")
    print(f"[0/4] Diretório desta execução: {outdir}")

    kind, src_path = _resolve_input_source(args.report, args.report_txt, args.pick)

    if kind == "txt":
        try:
            with open(src_path, "r", encoding="utf-8") as fh:
                text = fh.read()
        except UnicodeDecodeError:
            with open(src_path, "r", encoding="latin-1") as fh:
                text = fh.read()
        print("[1/4] Lendo TEXTO…")
        variants = parse_report_text(text)
    else:
        print("[1/4] Lendo PDF…")
        variants = parse_pdf_variants(src_path)

    with open(os.path.join(outdir, "variants_extracted.json"), "w", encoding="utf-8") as f:
        json.dump([asdict(v) for v in variants], f, ensure_ascii=False, indent=2)
    csv_path = write_csv_summary(variants, outdir)
    print(f"  → {len(variants)} variantes extraídas | CSV: {csv_path}")

    print("[2/4] Construindo páginas por gene…")
    gene_html_map = build_gene_pages(variants, outdir, roi_window=args.window, ncbi_email=(args.ncbi_email or os.environ.get("NCBI_TOOL_EMAIL")))

    print("[3/4] Gerando report_all.html…")
    master_path = write_master_report(gene_html_map, outdir)
    index_path = os.path.join(outdir, "index.html")
    print("       OK. Abra o index:", index_path)
    print("       OU o consolidado :", master_path)

    if args.open:
        try:
            import webbrowser, pathlib
            webbrowser.open(pathlib.Path(index_path).resolve().as_uri())
        except Exception as e:
            print(f"[WARN] Não foi possível abrir o navegador automaticamente: {e}")

    if args.bundle:
        zip_path = make_zip_bundle(outdir)
        print("       Bundle ZIP:", zip_path)

    print("[4/4] Pronto.")


if __name__ == "__main__":
    main()
