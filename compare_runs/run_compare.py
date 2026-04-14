#!/usr/bin/env python3
import csv
import re
import subprocess
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BENCH_DIR = (ROOT / "../ConsMLP/ss_benchmarks").resolve()
PATOH = (ROOT / "../ConsMLP/other_tools/patoh").resolve()
CONVERTER = (ROOT / "../ConsMLP/other_tools/hgr2patoh.py").resolve()
LW = (ROOT / "build/ConsMLP_lw").resolve()

WORK = ROOT / "compare_runs"
CONV = WORK / "converted"
LOGS = WORK / "logs"
PARTS = WORK / "parts"
for d in [WORK, CONV, LOGS, PARTS]:
    d.mkdir(parents=True, exist_ok=True)

hgr_files = sorted(p for p in BENCH_DIR.glob("*.hgr") if p.is_file())
ks = [2, 3, 4, 5]

re_lw_cut = re.compile(r"Cut size:\s*([0-9]+)")
re_lw_imb = re.compile(r"Imbalance:\s*([0-9]+(?:\.[0-9]+)?)%")
re_pa_cut = re.compile(r"Cost:\s*([0-9]+)")
re_pa_imb = re.compile(r"Part Weights\s*:.*Max=\s*[-0-9.]+\s*\(([0-9eE+\-.]+)\)")

rows = []

if not LW.exists():
    raise SystemExit(f"Missing executable: {LW}")
if not PATOH.exists():
    raise SystemExit(f"Missing executable: {PATOH}")
if not CONVERTER.exists():
    raise SystemExit(f"Missing converter: {CONVERTER}")

for hgr in hgr_files:
    bench = hgr.stem
    u = CONV / f"{bench}.u"

    conv = subprocess.run(
        ["python3", str(CONVERTER), str(hgr), str(u)],
        capture_output=True,
        text=True,
    )
    if conv.returncode != 0:
        for k in ks:
            rows.append({
                "benchmark": bench,
                "k": k,
                "status": f"convert_fail:{conv.stderr.strip()[:120]}"
            })
        continue

    for k in ks:
        lw_log_path = LOGS / f"lw_{bench}_k{k}.log"
        pa_log_path = LOGS / f"patoh_{bench}_k{k}.log"
        lw_part = PARTS / f"lw_{bench}.part.{k}"

        # run ConsMLP_lw
        t0 = time.perf_counter()
        lw_run = subprocess.run(
            [
                str(LW), str(hgr),
                "-k", str(k),
                "-mode", "recursive",
                "-imbalance", "0.05",
                "-output", str(lw_part),
            ],
            capture_output=True,
            text=True,
        )
        lw_runtime = time.perf_counter() - t0
        lw_text = (lw_run.stdout or "") + "\n" + (lw_run.stderr or "")
        lw_log_path.write_text(lw_text)

        # run patoh
        t1 = time.perf_counter()
        pa_run = subprocess.run(
            [str(PATOH), str(u), str(k), "WI=1"],
            capture_output=True,
            text=True,
        )
        pa_runtime = time.perf_counter() - t1
        pa_text = (pa_run.stdout or "") + "\n" + (pa_run.stderr or "")
        pa_log_path.write_text(pa_text)

        row = {
            "benchmark": bench,
            "k": k,
            "lw_runtime_s": lw_runtime,
            "patoh_runtime_s": pa_runtime,
            "status": "ok",
        }

        if lw_run.returncode != 0:
            row["status"] = f"lw_fail:{lw_run.returncode}"
        if pa_run.returncode != 0:
            row["status"] = (row["status"] + "|") if row["status"] != "ok" else ""
            row["status"] += f"patoh_fail:{pa_run.returncode}"

        m = re_lw_cut.findall(lw_text)
        row["lw_cut"] = int(m[-1]) if m else None
        m = re_lw_imb.findall(lw_text)
        row["lw_imbalance_pct"] = float(m[-1]) if m else None

        m = re_pa_cut.findall(pa_text)
        row["patoh_cut"] = int(m[-1]) if m else None
        m = re_pa_imb.findall(pa_text)
        row["patoh_imbalance_pct"] = float(m[-1]) * 100.0 if m else None

        def rel(a, b):
            if a is None or b is None or b == 0:
                return None
            return (a - b) / b * 100.0

        row["runtime_diff_pct"] = rel(row["lw_runtime_s"], row["patoh_runtime_s"])
        row["cut_diff_pct"] = rel(row["lw_cut"], row["patoh_cut"])
        row["imbalance_diff_pct"] = rel(row["lw_imbalance_pct"], row["patoh_imbalance_pct"])

        if row["lw_imbalance_pct"] is not None and row["patoh_imbalance_pct"] is not None:
            row["imbalance_diff_pp"] = row["lw_imbalance_pct"] - row["patoh_imbalance_pct"]
        else:
            row["imbalance_diff_pp"] = None

        rows.append(row)

# write csv
csv_path = WORK / "results.csv"
fields = [
    "benchmark", "k", "status",
    "lw_runtime_s", "patoh_runtime_s",
    "lw_cut", "patoh_cut",
    "lw_imbalance_pct", "patoh_imbalance_pct",
    "runtime_diff_pct", "cut_diff_pct",
    "imbalance_diff_pct", "imbalance_diff_pp",
]
with csv_path.open("w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fields)
    w.writeheader()
    for r in rows:
        w.writerow(r)

# summary helpers
ok_rows = [r for r in rows if r.get("status") == "ok"]

def avg(key):
    vals = [r[key] for r in ok_rows if isinstance(r.get(key), (int, float))]
    return (sum(vals) / len(vals), len(vals)) if vals else (None, 0)

avg_runtime_diff, n_runtime = avg("runtime_diff_pct")
avg_cut_diff, n_cut = avg("cut_diff_pct")
avg_imb_diff, n_imb = avg("imbalance_diff_pct")
avg_imb_pp, n_imbpp = avg("imbalance_diff_pp")

lines = []
lines.append("# ConsMLP_lw vs PaToH 对比测试汇总")
lines.append("")
lines.append("- 测试数据目录: `../ConsMLP/ss_benchmarks`")
lines.append("- 工具1: `./build/ConsMLP_lw`，参数固定 `-mode recursive -imbalance 0.05`")
lines.append("- 工具2: `../ConsMLP/other_tools/patoh`（先将 `.hgr` 转 `.u`）")
lines.append("- 对比分区数: `k=2/3/4/5`")
lines.append("")
lines.append("## 逐项结果")
lines.append("")
lines.append("| benchmark | k | lw_runtime(s) | patoh_runtime(s) | lw_cut | patoh_cut | lw_imb(%) | patoh_imb(%) | runtime_diff% | cut_diff% | imbalance_diff% | imbalance_diff(pp) | status |")
lines.append("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|")

for r in rows:
    def f(v, nd=4):
        if v is None:
            return "N/A"
        if isinstance(v, int):
            return str(v)
        return f"{v:.{nd}f}"
    lines.append(
        f"| {r.get('benchmark','')} | {r.get('k','')} | {f(r.get('lw_runtime_s'))} | {f(r.get('patoh_runtime_s'))} | "
        f"{f(r.get('lw_cut'),0)} | {f(r.get('patoh_cut'),0)} | {f(r.get('lw_imbalance_pct'),3)} | {f(r.get('patoh_imbalance_pct'),3)} | "
        f"{f(r.get('runtime_diff_pct'),2)} | {f(r.get('cut_diff_pct'),2)} | {f(r.get('imbalance_diff_pct'),2)} | {f(r.get('imbalance_diff_pp'),3)} | {r.get('status','')} |"
    )

lines.append("")
lines.append("## 平均差值百分比（lw 相对 patoh）")
lines.append("")
lines.append("- 公式: `diff% = (lw - patoh) / patoh * 100%`")
lines.append(f"- Runtime 平均差值百分比: {avg_runtime_diff:.2f}% (样本数={n_runtime})" if avg_runtime_diff is not None else "- Runtime 平均差值百分比: N/A")
lines.append(f"- Cut 平均差值百分比: {avg_cut_diff:.2f}% (样本数={n_cut})" if avg_cut_diff is not None else "- Cut 平均差值百分比: N/A")
lines.append(f"- Imbalance 平均差值百分比: {avg_imb_diff:.2f}% (样本数={n_imb}, 仅统计 patoh imbalance>0 的样本)" if avg_imb_diff is not None else "- Imbalance 平均差值百分比: N/A（patoh imbalance 为 0 的样本无法按该公式计算）")
lines.append(f"- Imbalance 平均差值（百分点）: {avg_imb_pp:.3f} pp (样本数={n_imbpp})" if avg_imb_pp is not None else "- Imbalance 平均差值（百分点）: N/A")
lines.append("")
lines.append(f"- 明细 CSV: `{csv_path.relative_to(ROOT)}`")
lines.append(f"- 日志目录: `{LOGS.relative_to(ROOT)}`")

(ROOT / "analyze_sum").write_text("\n".join(lines) + "\n")
print(f"Done. rows={len(rows)}, ok={len(ok_rows)}")
print(f"Summary file: {ROOT / 'analyze_sum'}")
print(f"CSV file: {csv_path}")
