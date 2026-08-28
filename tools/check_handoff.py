#!/usr/bin/env python3
"""Guard against the HANDOFF top-header overwrite bug.

FAILURE MODE (observed 2026-08-28, 4 occurrences):
  An editor that anchors on the first '## ' header and REPLACES that line —
  instead of inserting a new entry BEFORE it — destroys the previous entry's
  header while leaving its body intact. The body is then silently attributed
  to the new entry's author.

  Two variants:
    (a) previous entry was committed -> git shows '-## …' (one line deleted)
    (b) previous entry was NOT yet committed -> git shows NO deletion at all,
        and the loss is invisible in history.

  Note a naive entry-COUNT check does not catch this: one header is removed
  and one added, so the count is unchanged. The correct invariant is that
  entry headers are APPEND-ONLY — a header that has ever existed in history
  must still exist in the file.

Checks:
  1. append-only: every '## ' header seen in any past revision is still present
  2. no orphaned body immediately after the preamble
Exit 1 on violation. Intended as a pre-commit hook or CI step.
"""
import re, subprocess, sys

HDR = re.compile(r'^## \d{4}-\d{2}-\d{2} — (?:CC|CS) → (?:CS|CC) — .*$', re.M)

def headers(text):
    return [h.strip() for h in HDR.findall(text)]

def git(*args):
    return subprocess.run(("git",) + args, capture_output=True, text=True).stdout

def historical_headers():
    """Union of entry headers across every revision that touched HANDOFF.md."""
    revs = git("log", "--format=%H", "--", "HANDOFF.md").split()
    seen = {}
    for rev in revs:
        blob = git("show", f"{rev}:HANDOFF.md")
        for h in headers(blob):
            seen.setdefault(h, rev[:9])
    return seen

def orphaned_preamble_body(text):
    lines = text.split("\n")
    try:
        start = next(i for i, l in enumerate(lines) if l.strip() == "---") + 1
    except StopIteration:
        return None
    for i in range(start, len(lines)):
        l = lines[i]
        if HDR.match(l):
            return None
        if l.strip() not in ("", "---"):
            return (i + 1, l[:70])
    return None

def main():
    cur_text = open("HANDOFF.md", encoding="utf-8").read()
    cur = set(headers(cur_text))
    hist = historical_headers()
    missing = [(h, rev) for h, rev in hist.items() if h not in cur]

    fail = False
    if missing:
        print(f"FAIL: {len(missing)} entry header(s) present in history are GONE from "
              f"HANDOFF.md. An entry was overwritten in place.")
        for h, rev in missing[:10]:
            print(f"  last seen in {rev}: {h[:96]}")
        print("  FIX: insert new entries BEFORE the top header; never replace that line.")
        fail = True

    orph = orphaned_preamble_body(cur_text)
    if orph:
        print(f"FAIL: orphaned body text at line {orph[0]} with no header above it: {orph[1]}")
        fail = True

    if not fail:
        print(f"OK: {len(cur)} entries, all {len(hist)} historical headers intact, "
              f"no orphaned bodies.")
    sys.exit(1 if fail else 0)

if __name__ == "__main__":
    main()
