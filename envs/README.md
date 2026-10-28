# Conda environments

Two ways to get the tools genefinder needs, from loose to exact:

| File | What it pins | Use |
|---|---|---|
| [`../environment.yml`](../environment.yml) | the top-level tools (python 3.11, samtools/htslib 1.24, minimap2 2.31, impg 0.5.0, curl); conda solves the rest | `conda env create -f environment.yml` |
| `genefinder.<platform>.lock.txt` | every package, exact build, URL and md5 (conda `@EXPLICIT` format) | `conda create -n genefinder --file envs/genefinder.<platform>.lock.txt` |

Platforms: `linux-64`, `osx-64` (Intel Mac), `osx-arm64` (Apple silicon). The lock files were solved from
`environment.yml` (conda-forge + bioconda, strict channel priority) with micromamba 2.9.0 on 2026-10-28.
The `linux-64` lock is the environment used for the reference runs and tests on the Salz cluster. The macOS
locks are solver output and were not run.

An explicit lock file installs without a solver, so it gives the same environment every time, as long as the
packages stay on the channels. To activate the environment: `conda activate genefinder`. With micromamba or
mamba, use the same commands with `micromamba` or `mamba` in place of `conda`.

To regenerate the locks after changing `environment.yml` (needs only micromamba):

```bash
micromamba create -y -p ./tools/env -f environment.yml --channel-priority strict
micromamba env export -p ./tools/env --explicit --md5 > envs/genefinder.linux-64.lock.txt
# other platforms: solve without installing, then write the explicit list
micromamba create -n lock -f environment.yml --platform osx-arm64 --dry-run --json > solve.json
python3 envs/json2explicit.py solve.json > envs/genefinder.osx-arm64.lock.txt
```
