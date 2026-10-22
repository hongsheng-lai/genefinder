#!/usr/bin/env python3
"""List objects in the public human-pangenomics bucket (anonymous HTTPS, paginated).

Usage: s3ls.py <prefix> [--delimiter /]
Prints: size<TAB>key   (or PRE<TAB>prefix for common prefixes when --delimiter given)
"""
import os
import sys
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

BUCKET = os.environ.get("HPRC_BUCKET", "https://human-pangenomics.s3.amazonaws.com").rstrip("/") + "/"
NS = "{http://s3.amazonaws.com/doc/2006-03-01/}"


def list_objects(prefix, delimiter=None):
    token = None
    while True:
        q = {"list-type": "2", "prefix": prefix}
        if delimiter:
            q["delimiter"] = delimiter
        if token:
            q["continuation-token"] = token
        with urllib.request.urlopen(BUCKET + "?" + urllib.parse.urlencode(q), timeout=120) as r:
            root = ET.fromstring(r.read())
        for c in root.findall(NS + "Contents"):
            yield int(c.find(NS + "Size").text), c.find(NS + "Key").text
        for p in root.findall(NS + "CommonPrefixes"):
            yield "PRE", p.find(NS + "Prefix").text
        if root.findtext(NS + "IsTruncated") != "true":
            break
        token = root.findtext(NS + "NextContinuationToken")


if __name__ == "__main__":
    args = sys.argv[1:]
    delim = None
    if "--delimiter" in args:
        i = args.index("--delimiter")
        delim = args[i + 1]
        del args[i:i + 2]
    try:
        for size, key in list_objects(args[0] if args else "", delim):
            print(f"{size}\t{key}")
    except BrokenPipeError:
        pass
