# Portable helpers for Linux and macOS (sourced by run.sh and scripts/*.sh; works with bash 3.2).

# GNU time with -v (resource logs). Linux: /usr/bin/time; macOS: `brew install gnu-time` provides gtime.
# BSD /usr/bin/time has no -v, so it is only used if it accepts -v. Empty = resource logging disabled.
# A preset GNU_TIME (even empty) overrides the detection.
if [[ -z ${GNU_TIME+x} ]]; then
  GNU_TIME=
  for _t in gtime /usr/bin/time; do
    if command -v "$_t" >/dev/null 2>&1 && "$_t" -v true >/dev/null 2>&1; then GNU_TIME=$(command -v "$_t"); break; fi
  done
  unset _t
fi

ncpu() { nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4; }

# apparent size of a directory in bytes (du -sb is GNU-only)
dir_bytes() { du -sk "$1" 2>/dev/null | awk '{print $1 * 1024}'; }

md5_of() { if command -v md5sum >/dev/null 2>&1; then md5sum "$1" | cut -d' ' -f1; else md5 -q "$1"; fi; }

os_name() {
  if [[ $(uname -s) == Darwin ]]; then echo "macOS $(sw_vers -productVersion 2>/dev/null)"
  else (. /etc/os-release 2>/dev/null; echo "${PRETTY_NAME:-$(uname -s)}"); fi
}

ram_gb() {
  if [[ $(uname -s) == Darwin ]]; then echo $(( $(sysctl -n hw.memsize) / 1073741824 ))
  else awk '/^MemTotal:/{printf "%d\n", $2 / 1048576}' /proc/meminfo; fi
}
