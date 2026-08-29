#!/usr/bin/env python3
"""HANDOFF watcher — runs standalone, no agent session required.

Polls HANDOFF.md + the git log and appends a timestamped digest to
weekend_progress.md whenever something changes. Designed to be run from cron
(or a `while sleep` loop) on the PI's machine over a weekend, so that a
returning human sees a single chronological record of what moved.

  crontab:  */20 * * * * cd <repo> && python3 tools/watch_handoff.py >> tools/watch.log 2>&1

Checks each poll:
  * new HANDOFF entries (by header)
  * HANDOFF integrity (delegates to check_handoff.py — the append-only guard)
  * new commits
  * new/changed prod_* outputs in the sync tree
  * open [NEEDS DECISION] / [ACTION NEEDED] items still unanswered
State is kept in tools/.watch_state.json so each run reports only the delta.
"""
import json, os, re, subprocess, sys, hashlib
from datetime import datetime, timezone

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def _numpy_python():
    """Return an interpreter that can import numpy (the gates need it).
    Prefers the project conda env, falls back to whatever is running us."""
    cands = [os.environ.get("THERMOSPEC_PYTHON"),
             os.path.expanduser("~/.claude-science/conda/envs/thermospec/bin/python"),
             sys.executable]
    for c in cands:
        if c and os.path.exists(c):
            r = subprocess.run([c, "-c", "import numpy"], capture_output=True)
            if r.returncode == 0:
                return c
    return sys.executable
HANDOFF = os.path.join(REPO, "HANDOFF.md")
STATE = os.path.join(REPO, "tools", ".watch_state.json")
DIGEST = os.path.join(REPO, "weekend_progress.md")
SYNC = "/Users/phillipsm/Documents/Research/Publications/artemis-thermal-modeling/claude_session_sync"

HDR = re.compile(r'^## (\d{4}-\d{2}-\d{2}) — (CC|CS) → (?:CS|CC) — (.*)$', re.M)
FLAG = re.compile(r'\[(NEEDS DECISION|ACTION NEEDED)\]')

def git(*a):
    return subprocess.run(("git",)+a, cwd=REPO, capture_output=True, text=True).stdout.strip()

def load_state():
    if os.path.exists(STATE):
        try: return json.load(open(STATE))
        except Exception: pass
    return {"headers": [], "commit": "", "outputs": {}}

def scan_outputs():
    """Fingerprint prod_* / illum_* data products in the sync tree."""
    out = {}
    for root in (os.path.join(SYNC, "data"), os.path.join(SYNC, "figures")):
        if not os.path.isdir(root): continue
        for dirpath, _, files in os.walk(root):
            for fn in files:
                if not (fn.startswith(("prod_", "illum_", "big_psr_", "seasonal_")) or
                        fn.endswith((".png", ".npz"))):
                    continue
                p = os.path.join(dirpath, fn)
                try: st = os.stat(p)
                except OSError: continue
                out[os.path.relpath(p, SYNC)] = [st.st_size, int(st.st_mtime)]
    return out

def main():
    st = load_state()
    text = open(HANDOFF, encoding="utf-8").read()
    hdrs = ["## %s — %s → … — %s" % h for h in HDR.findall(text)]
    new_hdrs = [h for h in hdrs if h not in st["headers"]]

    commit = git("rev-parse", "--short", "HEAD")
    new_commits = []
    if st["commit"] and commit != st["commit"]:
        new_commits = git("log", "--oneline", f"{st['commit']}..HEAD").split("\n")
        new_commits = [c for c in new_commits if c]

    outs = scan_outputs()
    new_outs = [k for k, v in outs.items() if st["outputs"].get(k) != v]

    integ = subprocess.run([sys.executable, os.path.join(REPO, "tools", "check_handoff.py")],
                           cwd=REPO, capture_output=True, text=True)
    integ_ok = integ.returncode == 0
    integ_msg = integ.stdout.strip()

    # check_handoff.py compares against COMMITTED history, so it cannot see an
    # uncommitted clobber (variant (b) of the overwrite bug). The watcher's own
    # state file is a working-tree record, so use it as a second, independent
    # check: any header we saw on a previous poll must still be present.
    vanished = [h for h in st["headers"] if h not in hdrs]
    if vanished:
        integ_ok = False
        integ_msg = ((integ_msg + "\n") if integ_msg else "") + \
            "WATCHER: %d header(s) seen on a previous poll are now GONE (uncommitted overwrite):\n" % len(vanished) + \
            "\n".join("  " + v[:110] for v in vanished[:8])

    # open decision items: flagged entries whose author has not been replied to since
    blocks = re.split(r'(?m)^(?=## \d{4}-\d{2}-\d{2} — )', text)
    open_items = []
    for b in blocks:
        mh = HDR.search(b)
        if mh and FLAG.search(b):
            open_items.append((mh.group(1), mh.group(2), mh.group(3)[:72], FLAG.search(b).group(1)))
    open_items = open_items[:6]   # newest few only

    changed = bool(new_hdrs or new_commits or new_outs) or not integ_ok
    if changed:
        now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
        with open(DIGEST, "a", encoding="utf-8") as f:
            f.write(f"\n## {now}\n\n")
            if not integ_ok:
                f.write("**HANDOFF INTEGRITY FAILURE** (top-header overwrite bug):\n```\n"
                        + integ_msg + "\n```\n\n")
            if new_hdrs:
                f.write(f"**{len(new_hdrs)} new HANDOFF entr{'y' if len(new_hdrs)==1 else 'ies'}:**\n")
                for h in new_hdrs: f.write(f"- {h}\n")
                f.write("\n")
            if new_commits:
                f.write(f"**{len(new_commits)} new commit(s):**\n")
                for c in new_commits: f.write(f"- `{c}`\n")
                f.write("\n")
            if new_outs:
                f.write(f"**{len(new_outs)} new/changed data product(s):**\n")
                for k in sorted(new_outs)[:25]:
                    f.write(f"- `{k}` ({outs[k][0]/1e6:.2f} MB)\n")
                if len(new_outs) > 25: f.write(f"- …and {len(new_outs)-25} more\n")
                f.write("\n")
            # run the science gates on any new/changed run file
            new_runs = [k for k in new_outs if k.endswith(".npz") and
                        ("thermal" in k or "psr_floor" in k or "seasonal" in k)]
            if new_runs:
                paths = [os.path.join(SYNC, k) for k in new_runs[:6]]
                g = subprocess.run([_numpy_python(),
                                    os.path.join(REPO, "tools", "check_science_gates.py")] + paths,
                                   cwd=REPO, capture_output=True, text=True)
                verdict = "FAIL" if g.returncode else "PASS"
                f.write(f"**Science gates on {len(paths)} new run file(s): {verdict}**\n```\n"
                        + (g.stdout or g.stderr)[-2500:].strip() + "\n```\n\n")

            if open_items:
                f.write("**Open decision items:**\n")
                for dt, who, subj, kind in open_items:
                    f.write(f"- [{kind}] {dt} {who}: {subj}\n")
                f.write("\n")
        print(f"[{now}] digest updated: {len(new_hdrs)} entries, {len(new_commits)} commits, "
              f"{len(new_outs)} products, integrity={'OK' if integ_ok else 'FAIL'}")
    else:
        print(f"[{datetime.now(timezone.utc):%Y-%m-%d %H:%M UTC}] no change")

    json.dump({"headers": hdrs, "commit": commit, "outputs": outs}, open(STATE, "w"), indent=0)

if __name__ == "__main__":
    main()
