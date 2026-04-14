#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'USAGE'
Usage:
  tools/gen_xml_from_types.sh <types_file> <output_xml> [--parts N] [--logic-util U] [--relaxed-extra R]

Options:
  --parts N          Number of XML partitions/SLRs (default: 4)
  --logic-util U     Target utilization for LUT/FF (default: 0.80)
  --relaxed-extra R  Extra capacity ratio for DSP/IO/BRAM/MUX/CARRY/OTHER (default: 0.15)
USAGE
}

if [[ $# -lt 2 ]]; then
    usage
    exit 1
fi

types_file="$1"
output_xml="$2"
shift 2

parts=4
logic_util=0.80
relaxed_extra=0.15

while [[ $# -gt 0 ]]; do
    case "$1" in
        --parts)
            parts="$2"
            shift 2
            ;;
        --logic-util)
            logic_util="$2"
            shift 2
            ;;
        --relaxed-extra)
            relaxed_extra="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1" >&2
            usage
            exit 1
            ;;
    esac
done

if [[ ! -f "$types_file" ]]; then
    echo "Type file not found: $types_file" >&2
    exit 1
fi

awk -v parts="$parts" -v logic_util="$logic_util" -v relaxed_extra="$relaxed_extra" '
BEGIN {
    if (parts <= 0) {
        print "parts must be > 0" > "/dev/stderr";
        exit 2;
    }
    if (logic_util <= 0.0 || logic_util > 1.0) {
        print "logic-util must be in (0,1]" > "/dev/stderr";
        exit 2;
    }
    if (relaxed_extra < 0.0) {
        print "relaxed-extra must be >= 0" > "/dev/stderr";
        exit 2;
    }

    ntypes = 8;
    types[1] = "LUT";
    types[2] = "FF";
    types[3] = "DSP";
    types[4] = "IO";
    types[5] = "BRAM";
    types[6] = "MUX";
    types[7] = "CARRY";
    types[8] = "OTHER";
}
{
    t = $1;
    cnt[t]++;
}
END {
    for (slr = 0; slr < parts; ++slr) {
        printf("<SLR%d>\n", slr);
        for (i = 1; i <= ntypes; ++i) {
            t = types[i];
            c = cnt[t] + 0;
            cap = 0;
            if (c > 0) {
                base = (c / parts) / logic_util;
                if (t != "LUT" && t != "FF") {
                    base *= (1.0 + relaxed_extra);
                }
                cap = int(base);
                if (base > cap) {
                    cap += 1;
                }
            }
            printf("<%s><%d>\n", t, cap);
        }
        printf("</SLR%d>\n", slr);
    }
}
' "$types_file" > "$output_xml"

echo "Generated XML: $output_xml"
echo "  parts=$parts logic_util=$logic_util relaxed_extra=$relaxed_extra"
