#!/usr/bin/env python3
"""Unit tests, no network or external tools: python3 -m unittest discover -s tests"""
import io
import os
import struct
import sys
import unittest
import zlib

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from genefinder import build_paf_index, chain, project, status  # noqa: E402

# target 100..120, query 0..19: 5 match, 2 bp insertion, 5 match, 3 bp deletion, 7 match
CIGAR = "5=2I5=3D7="


class TestProject(unittest.TestCase):
    def test_full(self):
        self.assertEqual(project(CIGAR, 100, 0, 19, "+", 100, 120), (0, 19, 100, 120, 17, 0, 2))

    def test_end_inside_deletion(self):
        self.assertEqual(project(CIGAR, 100, 0, 19, "+", 103, 112), (3, 12, 103, 110, 7, 0, 1))

    def test_reverse(self):
        self.assertEqual(project(CIGAR, 100, 0, 19, "-", 100, 105)[:4], (14, 19, 100, 105))

    def test_inside_deletion_only(self):
        self.assertIsNone(project(CIGAR, 100, 0, 19, "+", 111, 113))


def piece(contig, q0, q1, t0, t1, strand="+"):
    return (contig, q0, q1, strand, t0, t1, t1 - t0, 0, 0)


class TestChain(unittest.TestCase):
    def test_split_block_one_hit(self):
        h = chain([piece("c", 0, 500, 0, 500), piece("c", 800, 1300, 500, 1000)], 1000)
        self.assertEqual(len(h), 1)
        self.assertEqual((h[0]["start"], h[0]["end"], h[0]["cov"]), (0, 1300, 1.0))

    def test_tandem_copy(self):
        self.assertEqual(len(chain([piece("c", 0, 1000, 0, 1000), piece("c", 3000, 4000, 0, 1000)], 1000)), 2)

    def test_far_apart(self):
        self.assertEqual(len(chain([piece("c", 0, 500, 0, 500), piece("c", 90000, 90500, 500, 1000)], 1000)), 2)

    def test_contig_strand(self):
        self.assertEqual(len(chain([piece("a", 0, 500, 0, 500), piece("b", 0, 500, 500, 1000),
                                    piece("a", 600, 900, 600, 900, "-")], 1000)), 3)


class TestStatus(unittest.TestCase):
    def test_flags(self):
        self.assertEqual(status([]), "missing")
        self.assertEqual(status([{"cov": 1.0}]), "ok")
        self.assertEqual(status([{"cov": 0.6}]), "partial")
        self.assertEqual(status([{"cov": 1.0}, {"cov": 0.9}]), "multi_copy")
        self.assertEqual(status([{"cov": 1.0}, {"cov": 0.1}]), "fragments")
        self.assertEqual(status([{"cov": 0.3}, {"cov": 0.2}]), "partial,fragments")


def bgzf(data):
    c = zlib.compressobj(6, zlib.DEFLATED, -15)
    body = c.compress(data) + c.flush()
    head = b"\x1f\x8b\x08\x04\0\0\0\0\0\xff\x06\0BC\x02\0" + struct.pack("<H", len(body) + 25)
    return head + body + struct.pack("<II", zlib.crc32(data), len(data))


class TestPafIndex(unittest.TestCase):
    def test_offsets(self):
        l1 = b"q1\t100\t0\t50\t+\tchrA\t1000\t10\t60\t50\t50\t60\tcg:Z:50M"
        l2 = b"q2\t100\t0\t50\t-\tchrB\t1000\t200\t250\t50\t50\t60\tcg:Z:50M"
        b1, b2 = bgzf(l1 + b"\n" + l2[:20]), bgzf(l2[20:] + b"\n")
        out = io.StringIO()
        self.assertEqual(build_paf_index(io.BytesIO(b1 + b2), out), 2)
        rows = [r.split("\t") for r in out.getvalue().splitlines()]
        self.assertEqual(rows[0], ["chrA", "10", "60", "0", str(len(b1)), "0", str(len(l1))])
        self.assertEqual(rows[1], ["chrB", "200", "250", "0", str(len(b1) + len(b2)), str(len(l1) + 1), str(len(l2))])


if __name__ == "__main__":
    unittest.main()
