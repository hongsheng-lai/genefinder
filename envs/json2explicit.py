#!/usr/bin/env python3
"""Turn `micromamba create --dry-run --json` output into a conda @EXPLICIT lock file.

Usage: json2explicit.py solve.json > envs/genefinder.<platform>.lock.txt
"""
import json
import sys


def main():
    solve = json.load(open(sys.argv[1]))
    pkgs = solve["actions"]["FETCH"] or solve["actions"]["LINK"]
    platform = sorted({p["subdir"] for p in pkgs} - {"noarch"})
    print("# This file may be used to create an environment using:")
    print("# $ conda create --name <env> --file <this file>")
    print(f"# platform: {platform[0] if platform else 'noarch'}")
    print("@EXPLICIT")
    for p in sorted(pkgs, key=lambda p: p["name"]):
        print(f"{p['url']}#{p['md5']}")


if __name__ == "__main__":
    main()
