#!/usr/bin/env python3
"""Live agenda for the CS <-> CC collaboration.

WHY THIS EXISTS
  HANDOFF.md is an append-only narrative: excellent as a record, useless as state.
  A header flagged [NEEDS DECISION] stays flagged forever, so counting flags tells
  you nothing about what is actually live. This file holds the state; HANDOFF.md
  keeps holding the reasoning.

CONTRACT
  Every agent, at the START of a working session on this project, runs:

      python3 tools/agenda.py --for CS     (or --for CC)

  Exit code 0 means nothing is waiting on you. Exit code 1 means something is.
  That makes it usable in a shell prompt, a cron digest, or a pre-flight check.

  When you finish an item:   python3 tools/agenda.py --close A3 --note "..."
  When you hand one over:    python3 tools/agenda.py --open --owner CC \
                                 --from CS --title "..." [--blocking]
"""
import argparse, json, os, sys, datetime

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PATH = os.path.join(ROOT, "agenda.json")
TODAY = datetime.date.today().isoformat()
OWNERS = ("CS", "CC", "PI")


def load():
    if not os.path.exists(PATH):
        return {"items": [], "next_id": 1}
    with open(PATH) as f:
        return json.load(f)


def save(d):
    with open(PATH, "w") as f:
        json.dump(d, f, indent=1)
        f.write("\n")


def fmt(it):
    age = (datetime.date.fromisoformat(TODAY) -
           datetime.date.fromisoformat(it["opened"])).days
    flag = "!" if it.get("blocking") else " "
    return (f"  {flag}{it['id']:>3}  {it['from']}->{it['owner']}  "
            f"{it['title']}\n        opened {it['opened']} ({age}d)"
            + (f" · {it['note']}" if it.get("note") else ""))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--for", dest="who", choices=OWNERS)
    ap.add_argument("--all", action="store_true", help="include closed items")
    ap.add_argument("--json", action="store_true")
    ap.add_argument("--open", action="store_true", help="add an item")
    ap.add_argument("--close", metavar="ID")
    ap.add_argument("--owner", choices=OWNERS)
    ap.add_argument("--from", dest="frm", choices=OWNERS)
    ap.add_argument("--title")
    ap.add_argument("--note", default="")
    ap.add_argument("--blocking", action="store_true",
                    help="the other side cannot proceed without it")
    a = ap.parse_args()
    d = load()

    if a.open:
        if not (a.owner and a.frm and a.title):
            sys.exit("--open needs --owner, --from and --title")
        it = {"id": f"A{d['next_id']}", "opened": TODAY, "owner": a.owner,
              "from": a.frm, "title": a.title, "status": "open",
              "blocking": a.blocking, "note": a.note}
        d["items"].append(it); d["next_id"] += 1; save(d)
        print("opened", it["id"]); return 0

    if a.close:
        for it in d["items"]:
            if it["id"] == a.close:
                if it["status"] == "done":
                    print(f"{a.close} was already closed on {it.get('closed')}"); return 0
                it["status"] = "done"; it["closed"] = TODAY
                if a.note: it["note"] = a.note
                save(d); print("closed", a.close); return 0
        sys.exit(f"no item {a.close}")

    items = [i for i in d["items"] if a.all or i["status"] == "open"]
    if a.who:
        items = [i for i in items if i["owner"] == a.who]
    if a.json:
        print(json.dumps(items, indent=1)); return 0

    if a.who:
        blocking = [i for i in items if i.get("blocking")]
        if not items:
            print(f"Nothing is waiting on {a.who}."); return 0
        print(f"{len(items)} open for {a.who}"
              + (f", {len(blocking)} blocking:" if blocking else ":"))
        for i in items: print(fmt(i))
        return 1

    for w in OWNERS:
        mine = [i for i in items if i["owner"] == w]
        if mine:
            print(f"\n{w} ({len(mine)}):")
            for i in mine: print(fmt(i))
    if not items: print("Agenda is clear.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
