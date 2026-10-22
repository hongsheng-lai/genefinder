#!/usr/bin/env python3
"""Offline unit tests (no network, no external tools): python3 -m unittest discover -s tests"""
import os
import subprocess
import sys
import tempfile
import unittest

SCRIPTS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts")
sys.path.insert(0, SCRIPTS)

from paf_project import project, project_many  # noqa: E402
from summarize_projection import chain_hits, union_len  # noqa: E402

# toy alignment: target 100..120 (20 bp), query 0..19; 5 match, 2 bp insertion, 5 match, 3 bp deletion, 7 match
CIGAR = "5=2I5=3D7="


class TestProjection(unittest.TestCase):
    def test_full_interval(self):
        self.assertEqual(project(CIGAR, 100, 0, 19, "+", 100, 120),
                         (0, 19, 100, 120, "5=2I5=3D7=", 17, 0, 2, 3, 2))

    def test_clip_ending_in_deletion(self):
        # interval end falls inside the deletion: projection stops at the last aligned base
        self.assertEqual(project(CIGAR, 100, 0, 19, "+", 103, 112),
                         (3, 12, 103, 110, "2=2I5=", 7, 0, 2, 0, 1))

    def test_reverse_strand(self):
        # first 5 target bases map to the last 5 query bases on a '-' alignment
        self.assertEqual(project(CIGAR, 100, 0, 19, "-", 100, 105)[:4], (14, 19, 100, 105))

    def test_interval_inside_deletion(self):
        self.assertIsNone(project(CIGAR, 100, 0, 19, "+", 111, 113))

    def test_many_equals_single(self):
        ivs = [(103, 112), (100, 120), (111, 113), (100, 105)]
        self.assertEqual(project_many(CIGAR, 100, 0, 19, "+", ivs),
                         [project(CIGAR, 100, 0, 19, "+", s, e) for s, e in ivs])


class TestStreamExtract(unittest.TestCase):
    def test_extract_and_revcomp(self):
        fasta = ">c1 desc\nACGTACGTAC\nGGGTTT\n>c2\nAAAA\n>c3\nCCCCGGGG\n"
        with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as fh:
            fh.write("c1:3-12\nc3:1-4\t-\n")
            regions = fh.name
        try:
            out = subprocess.run([sys.executable, os.path.join(SCRIPTS, "stream_extract.py"), "--regions", regions],
                                 input=fasta, capture_output=True, text=True, check=True).stdout
        finally:
            os.remove(regions)
        self.assertEqual(out, ">c1:3-12\nGTACGTACGG\n>c3:1-4/rc\nGGGG\n")


def piece(qname, qs, qe, ts, te, strand="+"):
    return dict(qname=qname, qlen=10 ** 6, qs=qs, qe=qe, strand=strand, ts=ts, te=te, m=te - ts, x=0, gaps=0)


class TestChaining(unittest.TestCase):
    def test_split_block_is_one_hit(self):
        hits = chain_hits([piece("h#1#c", 0, 500, 0, 500), piece("h#1#c", 800, 1300, 500, 1000)], 20000)
        self.assertEqual(len(hits), 1)
        self.assertEqual(hits[0]["n"], 2)
        self.assertEqual(union_len(hits[0]["tiv"]), 1000)

    def test_tandem_copy_is_separate_hit(self):
        # second piece re-covers the same CHM13 bases -> a second copy, not one longer hit
        hits = chain_hits([piece("h#1#c", 0, 1000, 0, 1000), piece("h#1#c", 3000, 4000, 0, 1000)], 20000)
        self.assertEqual(len(hits), 2)

    def test_far_pieces_split(self):
        hits = chain_hits([piece("h#1#c", 0, 500, 0, 500), piece("h#1#c", 90000, 90500, 500, 1000)], 20000)
        self.assertEqual(len(hits), 2)

    def test_contigs_and_strands_separate(self):
        hits = chain_hits([piece("h#1#a", 0, 500, 0, 500), piece("h#1#b", 0, 500, 500, 1000),
                           piece("h#1#a", 600, 900, 600, 900, "-")], 20000)
        self.assertEqual(len(hits), 3)


if __name__ == "__main__":
    unittest.main()
