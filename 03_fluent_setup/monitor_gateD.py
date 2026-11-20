"""
monitor_gateD.py - live progress bars for the Gate D batch (reads the run logs; does not touch Fluent).

    .venv\\Scripts\\python.exe 03_fluent_setup\\monitor_gateD.py        # refreshes every 10 s, Ctrl-C to quit

Progress of one run = iterations done / about 1650 (50 warm-up + ~600 flow + ~1000 energy, as in Gate B);
it shows 100 % only when the run has written its result. Time estimates scale with cell count
(Gate B: 47 min for 449k cells) and are corrected by the actual speed of finished runs.
"""
import os
import re
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import run_gateD as rd  # noqa: E402

EXPECTED_ITS = 1650
MIN_PER_KCELL = 47.2 / 449.328          # Gate B speed, minutes per 1000 cells
CELLS = {"D3C": 107_712, "D3M": 239_540, "B3D": 449_328}
NAMES = {"D3C": "coarse", "D3M": "medium", "B3D": "fine"}
W = 30


def bar(frac, width=W):
    n = int(round(max(0.0, min(1.0, frac)) * width))
    return "[" + "#" * n + "-" * (width - n) + "]"


def run_state(re_, mesh):
    log = rd.RES / f"{rd.PREFIX}_Re{re_}_{mesh}_log.txt"
    done = rd.result_file(re_, mesh).exists()
    if not log.exists():
        return dict(state="waiting", frac=0.0, its=0, start=None, end=None, note="")
    txt = log.read_text(errors="replace")
    start = log.stat().st_ctime
    end = log.stat().st_mtime
    flow = [int(m) for m in re.findall(r"flow it (\d+)", txt)]
    settled = re.search(r"flow settled after (\d+) its", txt)
    energy = [int(m) for m in re.findall(r"iterations: (\d+)\s+T_out", txt)]
    its = 0
    if "Warm-up" in txt:
        its = 50
    if settled:
        its = 50 + int(settled.group(1)) + (max(energy) if energy else 0)
    elif flow:
        its = 50 + max(flow)
    fin = re.search(r"Finished: (\d+) iterations", txt)
    m = re.search(r"converged=(\w+), wall time ([\d.]+) min", txt)
    note = ""
    if "ABORTED" in txt or "FAILED" in txt or "FLOW STAGE ended" in txt:
        note = "FAILED - see log"
        state = "failed"
    elif done:
        state = "done"
    else:
        state = "running"
    if done:
        nu = re.search(r"Nu_window: ([\d.]+)", txt)
        note = f"Nu {float(nu.group(1)):.3f}" if nu else "done"
    frac = 1.0 if done else min(0.99, its / EXPECTED_ITS)
    return dict(state=state, frac=frac, its=int(fin.group(1)) if fin else its, start=start, end=end, note=note,
                wall=float(m.group(2)) if m else None)


def fmt_min(mins):
    if mins is None:
        return "   --  "
    h, mm = divmod(int(round(mins)), 60)
    return f"{h}h{mm:02d}m" if h else f"{mm:4d} min"


def main():
    os.system("")                         # enable ANSI on Windows consoles
    while True:
        now = time.time()
        rows, speed = [], []
        for mesh in rd.MESHES:
            for re_ in rd.RE_LIST:
                s = run_state(re_, mesh)
                rows.append((re_, mesh, s))
                if s["state"] == "done" and s.get("wall"):
                    speed.append(s["wall"] / (CELLS[mesh] / 1000))
        k = (sum(speed) / len(speed)) if speed else MIN_PER_KCELL
        est = {m: k * CELLS[m] / 1000 for m in rd.MESHES}
        total = sum(est[m] for _, m, _ in rows)
        left = 0.0
        out = ["\x1b[2J\x1b[H", f"Gate D mesh study   {time.strftime('%H:%M:%S')}   (refresh 10 s, Ctrl-C to close this view only)\n"]
        for re_, mesh, s in rows:
            if s["state"] == "running":
                elapsed = (now - s["start"]) / 60
                rem = max(0.0, est[mesh] * (1 - s["frac"]))
                if s["frac"] > 0.1:
                    rem = elapsed / s["frac"] * (1 - s["frac"])
                left += rem
                tail = f"{s['its']:5d} its  {elapsed:5.1f} min, ~{fmt_min(rem)} left"
            elif s["state"] == "done":
                tail = f"{s['note']}  ({fmt_min(s.get('wall'))})"
            elif s["state"] == "failed":
                tail = s["note"]
            else:
                left += est[mesh]
                tail = f"waiting (~{fmt_min(est[mesh])})"
            colour = {"done": "\x1b[32m", "running": "\x1b[33m", "failed": "\x1b[31m"}.get(s["state"], "\x1b[90m")
            out.append(f"{colour}Re {re_:>6,}  {NAMES[mesh]:6s} {mesh}  {bar(s['frac'])} {100 * s['frac']:5.1f} %  {tail}\x1b[0m")
        overall = 1 - left / total if total else 0
        n_done = sum(1 for *_, s in rows if s["state"] == "done")
        out.append(f"\nOVERALL  {bar(overall, 50)} {100 * overall:5.1f} %   {n_done}/{len(rows)} runs done, "
                   f"~{fmt_min(left)} left, finish about {time.strftime('%H:%M', time.localtime(now + 60 * left))}")
        print("\n".join(out), flush=True)
        time.sleep(10)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        pass
