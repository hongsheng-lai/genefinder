#!/usr/bin/env python3
"""Compare the online (PAF projection) and offline (Liftoff) results for one gene.
Usage: compare.py ONLINE_DIR OFFLINE_DIR GENE   -> markdown on stdout"""
import sys
from collections import Counter


def primary(tsv):
    """hap -> first (= best) row."""
    rows = {}
    for l in list(open(tsv))[1:]:
        f = l.rstrip("\n").split("\t")
        rows.setdefault(f[1], f)
    return rows


def first_seqs(fa):
    seqs, hap, cur = {}, None, []
    for l in open(fa):
        if l.startswith(">"):
            if hap and hap not in seqs:
                seqs[hap] = "".join(cur)
            name = l[1:].split()[0].split("|")[1]
            hap, cur = "#".join(name.split("#")[:2]), []
            if " copy" in l:
                hap = None
        elif hap:
            cur.append(l.strip())
    if hap and hap not in seqs:
        seqs[hap] = "".join(cur)
    return seqs


def main():
    on_dir, off_dir, gene = sys.argv[1:4]
    on, off = primary(f"{on_dir}/{gene}.tsv"), primary(f"{off_dir}/{gene}.tsv")
    son, soff = first_seqs(f"{on_dir}/{gene}.fa"), first_seqs(f"{off_dir}/{gene}.fa")
    both = sorted(h for h in on if h in off and on[h][3] and off[h][3] and h != "CHM13#0")
    same_locus = [h for h in both if on[h][3] == off[h][3] and on[h][6] == off[h][6]
                  and int(on[h][4]) < int(off[h][5]) and int(off[h][4]) < int(on[h][5])]
    exact = [h for h in same_locus if on[h][4:6] == off[h][4:6]]
    ident = [h for h in same_locus if son.get(h) == soff.get(h)]
    d5 = Counter(abs(int(on[h][4]) - int(off[h][4])) for h in same_locus)
    d3 = Counter(abs(int(on[h][5]) - int(off[h][5])) for h in same_locus)
    dl = Counter(abs(len(son.get(h, "")) - len(soff.get(h, ""))) for h in same_locus)

    def dist(c):
        tot = sum(c.values()) or 1
        z = c.get(0, 0)
        le10 = sum(v for k, v in c.items() if k <= 10)
        return f"0 bp: {z}, <=10 bp: {le10}, max: {max(c) if c else 0} bp (n={tot})"

    print(f"### Online vs offline: {gene}\n")
    print("| | online | offline |\n|---|---|---|")
    for st in sorted(set(r[2] for r in on.values()) | set(r[2] for r in off.values())):
        print(f"| {st} | {sum(r[2] == st for r in on.values())} | {sum(r[2] == st for r in off.values())} |")
    print(f"\nHaplotypes found by both: {len(both)}; same contig, strand and overlapping interval: {len(same_locus)}")
    print(f"- identical start and end: {len(exact)}; identical sequence: {len(ident)}")
    print(f"- |start difference|: {dist(d5)}")
    print(f"- |end difference|: {dist(d3)}")
    print(f"- |length difference|: {dist(dl)}")
    diff_locus = [h for h in both if h not in same_locus]
    if diff_locus:
        print(f"- different locus ({len(diff_locus)}): " + ", ".join(
            f"{h} on={on[h][3]}:{on[h][4]}-{on[h][5]} off={off[h][3]}:{off[h][4]}-{off[h][5]}" for h in diff_locus[:10]))
    only_on = [h for h in on if on[h][3] and not off.get(h, ["", "", "", ""])[3]]
    only_off = [h for h in off if off[h][3] and not on.get(h, ["", "", "", ""])[3]]
    print(f"- found only online: {len(only_on)} {' '.join(sorted(only_on)[:12])}")
    print(f"- found only offline: {len(only_off)} {' '.join(sorted(only_off)[:12])}")


if __name__ == "__main__":
    main()
