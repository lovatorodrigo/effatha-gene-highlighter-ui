from effatha_core import parse_report_text

def test_parse_basic():
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
