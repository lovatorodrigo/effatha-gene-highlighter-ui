from __future__ import annotations
import re, os
from dataclasses import dataclass
from typing import List, Dict, Optional, Tuple
from io import StringIO

# PDF extractors (lazy)
try:
    from pdfminer.high_level import extract_text as _pdfminer_extract_text  # type: ignore
except Exception:
    _pdfminer_extract_text = None
try:
    import PyPDF2  # type: ignore
except Exception:
    PyPDF2 = None

import requests  # type: ignore
from Bio import SeqIO  # type: ignore

# ----------------- data -----------------
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
    cds_range: Tuple[int, int]  # 1-based inclusive
    cds_seq: str

@dataclass
class Highlight:
    start: int
    end: int
    kind: str   # 'SNV'|'DEL'|'INS'|'DUP'|'OTHER'
    label: str

# ----------------- PDF/TXT parsing -----------------
def extract_text_from_pdf_path(pdf_path: str) -> str:
    """Extrai texto do PDF usando pdfminer.six se disponível; senão tenta PyPDF2."""
    if _pdfminer_extract_text is not None:
        return _pdfminer_extract_text(pdf_path)
    if PyPDF2 is not None:
        parts = []
        with open(pdf_path, 'rb') as fh:
            reader = PyPDF2.PdfReader(fh)
            for page in reader.pages:
                parts.append(page.extract_text() or '')
        return '\n'.join(parts)
    raise ModuleNotFoundError("Instale 'pdfminer.six' ou 'PyPDF2'.")

# Padrões de extração
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

# normalizador de símbolo gênico (remove ruídos como 'tumoralNotaEP300', '2KRAS', etc.)
def _normalize_gene_token(raw: str) -> str:
    """
    Limpa ruído comum de texto de PDF/OCR.
    - remove não alfanuméricos e põe em maiúsculas
    - corta dígitos iniciais (ex.: '2KRAS' -> 'KRAS')
    - captura o último padrão 'LETRAS(2+)[DIGITOS(0-3)]' no fim (ex.: 'TUMORALNOTAEP300' -> 'EP300')
    """
    if raw is None:
        return ""
    up = re.sub(r"[^A-Za-z0-9]", "", raw).upper()
    up = re.sub(r"^\d+", "", up)
    m = re.search(r"([A-Z]{2,}[0-9]{0,3})$", up)
    return m.group(1) if m else (up or raw)

def parse_report_text(text: str) -> List[Variant]:
    lines = [l.strip() for l in (text or "").splitlines() if l.strip()]
    blob = "\n".join(lines)
    variants: List[Variant] = []

    # SNV/Indel (Oncogênicas/Tier 2)
    for m in _re_snvindel.finditer(blob):
        vaf_raw = m.group("vaf")
        vaf = float(vaf_raw.replace('%','').replace(',','.')) if vaf_raw else None
        window = blob[max(0, m.start()-80): m.end()+80]
        onc = "Oncogênica" if _re_onco_labels["oncogenic"].search(window) else (
              "Provavelmente oncogênica" if _re_onco_labels["prob_oncogenic"].search(window) else None)
        tier = "Tier 2" if _re_onco_labels["tier2"].search(window) else None
        variants.append(Variant(
            gene=_normalize_gene_token(m.group("gene")),
            transcript=m.group("tx"),
            hgvs_c=m.group("c"),
            hgvs_p=(m.group("p") or None),
            vaf_pct=vaf,
            tier=tier,
            oncogenicity=onc,
            section="oncogenic",
            notes=[],
        ))

    # CNV
    for m in _re_cnv.finditer(blob):
        variants.append(Variant(
            gene=_normalize_gene_token(m.group("gene")),
            transcript=m.group("tx"),
            hgvs_c=None,
            hgvs_p=None,
            vaf_pct=None,
            tier="Tier 2",
            oncogenicity="Oncogênica",
            section="cnv",
            notes=["amplificacao"],
        ))

    # VUS
    for m in _re_vus_block.finditer(blob):
        vaf = float(m.group("vaf").replace('%','').replace(',','.'))
        variants.append(Variant(
            gene=_normalize_gene_token(m.group("gene")),
            transcript=m.group("tx"),
            hgvs_c=m.group("c"),
            hgvs_p=(m.group("p") or None),
            vaf_pct=vaf,
            tier="Tier 3",
            oncogenicity="VUS",
            section="vus",
            notes=[],
        ))

    # Dedup best-effort
    uniq: Dict[Tuple, Variant] = {}
    for v in variants:
        key = (v.gene, v.transcript, v.hgvs_c, v.section)
        if key not in uniq:
            uniq[key] = v
    return list(uniq.values())

# ----------------- NCBI -----------------
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

# ----------------- HGVS → CDS -----------------
_re_snv = re.compile(r"c\.(?P<pos>\d+)(?P<ref>[ACGT])>(?P<alt>[ACGT])$")
_re_del = re.compile(r"c\.(?P<start>\d+)_(?P<end>\d+)del(?P<delseq>[ACGT]*)$")
_re_ins = re.compile(r"c\.(?P<start>\d+)_(?P<end>\d+)ins(?P<insseq>[ACGT]+)$")
_re_dup = re.compile(r"c\.(?P<pos>\d+)dup(?P<dupseq>[ACGT]*)$")

def apply_hgvs_c_to_cds(cds_seq: str, hgvs_c: str):
    s = cds_seq
    hs: List[Highlight] = []

    m = _re_snv.match(hgvs_c)
    if m:
        pos = int(m.group("pos"))
        alt = m.group("alt")
        idx = pos - 1
        if idx < 0 or idx >= len(s):
            raise ValueError(f"Posição fora do CDS: {hgvs_c}")
        mutated = s[:idx] + alt + s[idx+1:]
        hs.append(Highlight(start=pos, end=pos, kind="SNV", label=hgvs_c))
        return mutated, hs

    m = _re_del.match(hgvs_c)
    if m:
        start = int(m.group("start")); end = int(m.group("end"))
        if start > end:
            start, end = end, start
        mutated = s[:start-1] + s[end:]
        hs.append(Highlight(start=start, end=end, kind="DEL", label=hgvs_c))
        return mutated, hs

    m = _re_ins.match(hgvs_c)
    if m:
        end = int(m.group("end"))
        ins = m.group("insseq").upper()
        mutated = s[:end] + ins + s[end:]
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
            end = pos + 1  # highlight cobre a base recém-inserida
        hs.append(Highlight(start=pos+1, end=end, kind="DUP", label=hgvs_c))
        return mutated, hs

    return s, [Highlight(start=1, end1 := 1, kind="OTHER", label=f"não aplicado: {hgvs_c}")]

# ----------------- HTML highlight -----------------
def escape_html(s: str) -> str:
    return s.replace('&','&amp;').replace('<','&lt;').replace('>','&gt;')

def wrap_seq_with_highlights(seq: str, highlights: List[Highlight], line: int = 90) -> str:
    """
    Renderiza a sequência com <mark> sem inserir <br> dentro de tags.
    Fazemos a quebra por largura apenas sobre o texto visível, emitindo <br> entre segmentos.
    """
    n = len(seq)
    if n == 0:
        return ""

    # Normaliza/clampa destaques
    norm: List[Tuple[int, int, str]] = []
    for h in (highlights or []):
        s = max(1, min(n, h.start))
        e = max(1, min(n, h.end))
        if s > e:
            s, e = e, s
        norm.append((s, e, h.label))

    # Pontos de corte dos segmentos
    cuts = {1, n + 1}
    for s0, e0, _ in norm:
        cuts.add(s0)
        cuts.add(e0 + 1)
    cuts = sorted(cuts)

    # Para cada segmento, se está marcado e quais labels carrega
    def labels_for(a: int, b: int) -> List[str]:
        labs = []
        for s0, e0, lab in norm:
            if not (e0 < a or s0 > b):
                if lab not in labs:
                    labs.append(lab)
        return labs

    segments: List[Tuple[int, int, bool, str]] = []
    for i in range(len(cuts) - 1):
        a, b = cuts[i], cuts[i + 1] - 1
        if a > b:
            continue
        labs = labels_for(a, b)
        segments.append((a, b, len(labs) > 0, ", ".join(labs)))

    # Emite HTML por linhas (line width), sem quebrar tags
    out: List[str] = []
    col = 0
    pos = 1  # 1-based

    while pos <= n:
        line_end = min(n, pos + (line - col) - 1)

        # Emite todos os segmentos que tocam esta janela
        for (a, b, marked, lab) in segments:
            if b < pos:
                continue
            if a > line_end:
                break
            s0 = max(a, pos)
            e0 = min(b, line_end)
            frag = escape_html(seq[s0 - 1:e0])
            if marked:
                out.append(f'<mark title="{escape_html(lab)}">{frag}</mark>')
            else:
                out.append(frag)

        col += (line_end - pos + 1)
        pos = line_end + 1
        if pos <= n:
            out.append("<br>")
            col = 0

    return "".join(out)
