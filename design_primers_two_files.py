#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Diseña primers discriminatorios usando DOS archivos FASTA separados:
  - target (Neandertal)
  - nontarget (Homo sapiens)
Coloca el 3′ del primer forward entre 901–960 (1-based) sobre el target.
Requiere: biopython (usa pairwise2; el warning de deprecación es normal).

Uso:
  python3 design_primers_two_files.py neander_target.fa human_nontarget.fa
o simplemente (si usas los nombres por defecto):
  python3 design_primers_two_files.py
"""

import sys, os, csv, re
from Bio import SeqIO, pairwise2
from Bio.Seq import Seq

# ----------------------- Parámetros -----------------------
FWD_3P_RANGE = (901, 960)       # rango deseado (1-based) para el 3′ del FWD
PRODUCT_RANGE = (80, 140)       # tamaño recomendado para ADN antiguo
PRIMER_MIN_LEN = 18
PRIMER_MAX_LEN = 24
TM_MIN, TM_MAX = 58.0, 62.0     # usaremos Tm (Wallace) aproximado
GC_MIN, GC_MAX = 40.0, 60.0
REQUIRE_BOTH_3P_SNP = False     # True: 3′ en SNP en ambos primers; False: al menos en FWD

MAX_RUN = 4
MAX_TERMINAL_GC_CLAMP = 3

# Alineación global (pairwise2)
ALIGN = dict(match=2, mismatch=-1, gap_open=-5, gap_ext=-1)

# ----------------------- Utilidades -----------------------
def clean_seq(s):
    s = s.upper().replace("\r","").replace("\n","")
    return re.sub(r"[^ACGTN]", "", s)

def read_first_fasta(path):
    recs = list(SeqIO.parse(path, "fasta"))
    if not recs:
        sys.exit(f"ERROR: No encontré secuencias en {path}")
    rec = recs[0]
    seq = clean_seq(str(rec.seq))
    if not seq:
        sys.exit(f"ERROR: Secuencia vacía tras limpieza en {path}")
    rec.seq = Seq(seq)
    return rec

def gc_percent(seq):
    s = seq.upper()
    atgc = sum(s.count(x) for x in "ATGC")
    if atgc == 0: return 0.0
    return 100.0*(s.count("G")+s.count("C"))/atgc

def tm_wallace(seq):
    s = seq.upper()
    a=s.count('A'); t=s.count('T'); g=s.count('G'); c=s.count('C')
    return 2*(a+t)+4*(g+c)

def has_bad_runs(seq, max_run=MAX_RUN):
    u = seq.upper()
    return any(base*(max_run+1) in u for base in "ATGC")

def terminal_gc_clamp_len(seq):
    cnt = 0
    for ch in reversed(seq.upper()):
        if ch in "GC": cnt += 1
        else: break
    return cnt

def screen_primer(seq):
    if not (PRIMER_MIN_LEN <= len(seq) <= PRIMER_MAX_LEN): return False, "len"
    tm = tm_wallace(seq)
    if not (TM_MIN <= tm <= TM_MAX): return False, "Tm"
    gc = gc_percent(seq)
    if not (GC_MIN <= gc <= GC_MAX): return False, "GC"
    if has_bad_runs(seq): return False, "runs"
    if terminal_gc_clamp_len(seq) > MAX_TERMINAL_GC_CLAMP: return False, "3'GCclamp"
    return True, "ok"

def global_align_strings(a, b):
    alns = pairwise2.align.globalms(a, b, ALIGN["match"], ALIGN["mismatch"], ALIGN["gap_open"], ALIGN["gap_ext"])
    at, ao, score, _, _ = alns[0]
    return at, ao

def build_maps(aln_t, aln_o):
    # maps: alignment index -> target pos (0b); and inverse
    alnidx_to_tpos = {}
    tpos_to_alnidx = {}
    t_i = 0
    for i, tc in enumerate(aln_t):
        if tc != '-':
            alnidx_to_tpos[i] = t_i
            tpos_to_alnidx[t_i] = i
            t_i += 1
    return alnidx_to_tpos, tpos_to_alnidx

def diagnostic_snp_idxs(aln_t, aln_o):
    snps = []
    for i,(tc,oc) in enumerate(zip(aln_t, aln_o)):
        if tc in "ATGC" and oc in "ATGC" and tc != oc:
            snps.append(i)
    return snps

def fwd_primer_from_endpos(tseq, end_pos_0b, length):
    start = end_pos_0b - (length - 1)
    if start < 0: return None
    return tseq[start:end_pos_0b+1]

def rev_primer_from_startpos(tseq, start_pos_0b, length):
    end = start_pos_0b + length
    if end > len(tseq): return None
    return str(Seq(tseq[start_pos_0b:end]).reverse_complement())

# ----------------------- Principal -----------------------
def main():
    # Archivos
    target_path = sys.argv[1] if len(sys.argv) > 1 else "neander_target.fa"
    nontarget_path = sys.argv[2] if len(sys.argv) > 2 else "human_nontarget.fa"

    target = read_first_fasta(target_path)
    other  = read_first_fasta(nontarget_path)
    tseq = str(target.seq)
    oseq = str(other.seq)

    # Alineación
    aln_t, aln_o = global_align_strings(tseq, oseq)
    alnidx_to_tpos, tpos_to_alnidx = build_maps(aln_t, aln_o)
    snp_idxs = diagnostic_snp_idxs(aln_t, aln_o)
    if not snp_idxs:
        sys.exit("No se detectaron SNPs diagnósticos entre target y non-target.")

    tlen = len(tseq)
    fwd_lo = max(1, FWD_3P_RANGE[0])
    fwd_hi = min(tlen, FWD_3P_RANGE[1])

    candidates = []

    for f_3p_pos_1b in range(fwd_lo, fwd_hi+1):
        f_3p_pos_0b = f_3p_pos_1b - 1
        aln_idx_f = tpos_to_alnidx.get(f_3p_pos_0b, None)
        if aln_idx_f is None or aln_idx_f not in snp_idxs:
            continue

        # Base en target vs non-target en ese SNP (para confirmar mismatch en 3′ del FWD)
        t_f = aln_t[aln_idx_f]; o_f = aln_o[aln_idx_f]
        if not (t_f in "ATGC" and o_f in "ATGC" and t_f != o_f):
            continue

        for aln_idx_r in snp_idxs if REQUIRE_BOTH_3P_SNP else list(range(len(aln_t))):
            # Si REQUIRE_BOTH_3P_SNP=False, permitimos que el REV no esté en SNP; si True, debe estar
            if aln_idx_r not in alnidx_to_tpos:
                continue
            r_3p_pos_0b = alnidx_to_tpos[aln_idx_r]
            if r_3p_pos_0b <= f_3p_pos_0b:
                continue

            product_len = r_3p_pos_0b - f_3p_pos_0b + 1
            if not (PRODUCT_RANGE[0] <= product_len <= PRODUCT_RANGE[1]):
                continue

            if REQUIRE_BOTH_3P_SNP:
                t_r = aln_t[aln_idx_r]; o_r = aln_o[aln_idx_r]
                if not (t_r in "ATGC" and o_r in "ATGC" and t_r != o_r):
                    continue

            # Probar longitudes
            for Lf in range(PRIMER_MIN_LEN, PRIMER_MAX_LEN+1):
                fwd_seq = fwd_primer_from_endpos(tseq, f_3p_pos_0b, Lf)
                if not fwd_seq: continue
                okf, whyf = screen_primer(fwd_seq)
                if not okf: continue

                for Lr in range(PRIMER_MIN_LEN, PRIMER_MAX_LEN+1):
                    rev_seq = rev_primer_from_startpos(tseq, r_3p_pos_0b, Lr)
                    if not rev_seq: continue
                    okr, whyr = screen_primer(rev_seq)
                    if not okr: continue

                    f_tm = tm_wallace(fwd_seq)
                    r_tm = tm_wallace(rev_seq)
                    f_gc = gc_percent(fwd_seq)
                    r_gc = gc_percent(rev_seq)
                    tm_diff = abs(f_tm - r_tm)

                    candidates.append({
                        "fwd_seq": fwd_seq, "rev_seq": rev_seq,
                        "f_tm": round(f_tm,1), "r_tm": round(r_tm,1),
                        "f_gc": round(f_gc,1), "r_gc": round(r_gc,1),
                        "f_3p_pos_1b": f_3p_pos_1b,
                        "r_3p_pos_1b": r_3p_pos_0b+1,
                        "product_len": product_len,
                        "tm_diff": round(tm_diff,1)
                    })

    if not candidates:
        sys.exit("No hay candidatos con estos filtros. Ajusta FWD_3P_RANGE / PRODUCT_RANGE o REQUIRE_BOTH_3P_SNP.")

    # Ordena por cercanía a 100 bp, ΔTm pequeño, mayor Tm media
    candidates.sort(key=lambda c: (abs(c["product_len"]-100), c["tm_diff"], -((c["f_tm"]+c["r_tm"])/2.0)))

    # Salidas
    os.makedirs("out", exist_ok=True)
    with open("out/primers_candidates.csv","w",newline="") as f:
        w=csv.writer(f)
        w.writerow(["pair_id","FWD_seq","REV_seq","FWD_Tm","REV_Tm","FWD_GC(%)","REV_GC(%)",
                    "FWD_3p_pos(1b)","REV_3p_pos(1b)","Amplicon_len_bp","Tm_diff"])
        for i,c in enumerate(candidates, start=1):
            w.writerow([f"pair_{i}", c["fwd_seq"], c["rev_seq"], c["f_tm"], c["r_tm"], c["f_gc"], c["r_gc"],
                        c["f_3p_pos_1b"], c["r_3p_pos_1b"], c["product_len"], c["tm_diff"]])

    with open("out/primers_candidates.fasta","w") as f:
        for i,c in enumerate(candidates[:50], start=1):
            f.write(f">pair_{i}_FWD\n{c['fwd_seq']}\n>pair_{i}_REV\n{c['rev_seq']}\n")

    # Top 1 para chequeo
    top = candidates[0]
    with open("out/check_primers.fa","w") as f:
        f.write(f">FWD_pair1\n{top['fwd_seq']}\n>REV_pair1\n{top['rev_seq']}\n")

    print(f"[OK] {len(candidates)} pares. Revisa ./out/primers_candidates.csv y ./out/primers_candidates.fasta")
    print(f"Top1: Amp={top['product_len']} bp | FWD(3′@{top['f_3p_pos_1b']}, Tm={top['f_tm']}, GC={top['f_gc']}%) / "
          f"REV(3′@{top['r_3p_pos_1b']}, Tm={top['r_tm']}, GC={top['r_gc']}%) | ΔTm={top['tm_diff']}")
    print("Par guardado en out/check_primers.fa")
    
if __name__ == "__main__":
    main()
