#!/usr/bin/env python3
"""Replay recorded per-pair alignment times under different assignments and
schedules, with and without the FastGA-order constraint on self pairs.

Why: since src/fastga_order.c, a same-strand self pair (c,c) must start from
the end position left by its nearest earlier partner in FastGA's contig order
that reached the aligner, so it cannot start before those pairs are aligned,
wherever they are. This script estimates what that costs for load-balancing
schemes (proposal O3), from the per-pair records PGA_PAIR_TASKS=1 writes
(pairs.rank*.tsv: type, anchor, partner, seeds, seconds, nhits, alen, blen).

Usage:
  python3 scripts/simulate_pair_schedule.py <records dir> <stem>.gix <ranks>...

Model (all times are the recorded align_contigs seconds; everything else in
the alignment stage is ignored):
  task      one (type, anchor, partner) record
  group     all tasks of one unordered contig pair; a group runs on one rank
  deps      for self task (c,c,N): the N tasks (c,b) with b before c in
            FastGA's order, nearest first, up to the first with nhits > 0
            (for a pair of two contigs, a hit always reaches the aligner).
            If none has hits, the value comes from an earlier row: modelled
            as ready only once every rank has finished all other work.
  static    groups assigned before alignment: "length" is the current
            length-product rule; "oracle" assigns by the recorded cost
            (largest first, least-loaded rank), an upper bound for any cost
            model. "free" = each rank's total work; "ordered" = the order
            src/fastga_pipeline.c uses (earlier partners nearest first, later
            partners, ready self pairs, opposite strand, then self pairs set
            aside until their values arrive); "cheap" = the same, but among
            equally near earlier partners the one with fewer seeds first
            (seed counts are known before alignment); "+cost" = earlier
            partners that unlock the most expensive self pair first, from the
            recorded cost (what a cost model would aim for).
  dynamic   P workers take the largest ready group from a shared queue (a
            scheme that does not exist yet; seeds currently go to a fixed
            owner before alignment). "free" ignores deps; "ordered" makes a
            self group wait for its deps; "ordered+boost" also gives groups
            that unlock a self group that self group's cost as priority
            ("+cost").
  bound     max(total work / P, largest group): no schedule of whole groups
            can beat it.
"""
import collections
import csv
import glob
import heapq
import os
import struct
import sys


def read_perm(gix):
    with open(gix, "rb") as f:
        f.seek(4 * 4 + 0x1000000 * 8 + 4 * 3 + 8 + 4)
        n, = struct.unpack("<i", f.read(4))
        return list(struct.unpack("<%di" % n, f.read(4 * n)))


def load(recdir):
    tasks, clen = [], {}
    for path in sorted(glob.glob(os.path.join(recdir, "pairs.rank*.tsv"))):
        for r in csv.DictReader(open(path), delimiter="\t"):
            a, b = int(r["anchor"]), int(r["partner"])
            clen[a], clen[b] = int(r["alen"]), int(r["blen"])
            tasks.append((r["type"], a, b, float(r["seconds"]), int(r["nhits"]),
                          int(r["seeds"])))
    return tasks, clen


def length_assign(contigs, clen, P):
    """src/contig_assignment.c: pairs by length product, descending, (i,j)
    tie-break, each to the least-loaded rank (ties to the lower rank)."""
    cs = sorted(contigs)
    pairs = sorted(((clen.get(i, 0) * clen.get(j, 0), i, j)
                    for x, i in enumerate(cs) for j in cs[x:]),
                   key=lambda p: (-p[0], p[1], p[2]))
    heap = [(0, r) for r in range(P)]
    own = {}
    for w, i, j in pairs:
        load_, r = heapq.heappop(heap)
        own[(i, j)] = r
        heapq.heappush(heap, (load_ + w, r))
    return own


def oracle_assign(gcost, P):
    heap = [(0.0, r) for r in range(P)]
    own = {}
    for g, c in sorted(gcost.items(), key=lambda kv: (-kv[1], kv[0])):
        load_, r = heapq.heappop(heap)
        own[g] = r
        heapq.heappush(heap, (load_ + c, r))
    return own


def main():
    recdir, gix, ranks = sys.argv[1], sys.argv[2], [int(x) for x in sys.argv[3:]]
    perm = read_perm(gix)
    pos = {c: i for i, c in enumerate(perm)}
    tasks, clen = load(recdir)
    key = lambda a, b: (min(a, b), max(a, b))

    gcost = collections.defaultdict(float)
    gtasks = collections.defaultdict(list)
    for t in tasks:
        gcost[key(t[1], t[2])] += t[3]
        gtasks[key(t[1], t[2])].append(t)
    total = sum(gcost.values())
    biggest = max(gcost.values())

    # deps of each self N task, as (anchor, partner) N tasks
    earlier = collections.defaultdict(list)
    for ty, a, b, s, h, _ in tasks:
        if ty == "N" and a != b and pos[b] < pos[a]:
            earlier[a].append((pos[a] - pos[b], b, h))
    deps, rowcase, nself = {}, set(), 0
    for ty, a, b, s, h, _ in tasks:
        if ty == "N" and a == b:
            nself += 1
            chain = []
            for _, p, hits in sorted(earlier.get(a, [])):
                chain.append(p)
                if hits > 0:
                    break
            else:
                if pos[a] != 0:      # value from an earlier row (or unset)
                    rowcase.add(a)
            deps[a] = chain
    near = {}
    for a, lst in earlier.items():
        for k, (_, b, _) in enumerate(sorted(lst)):
            near[(a, b)] = k

    print(f"{recdir}: {len(tasks)} tasks, {len(gcost)} pairs, {nself} same-strand self pairs, "
          f"total {total:.0f} s, largest pair {biggest:.0f} s, "
          f"self pairs with an earlier-row value {len(rowcase)}")
    top = sorted(gcost.items(), key=lambda kv: -kv[1])[:5]
    print("  largest pairs: " + ", ".join(
        f"({i},{j}) {c:.0f} s{' self, deps ' + str(len(deps.get(i, []))) if i == j else ''}"
        for (i, j), c in top))

    # cost of the self pair a dep task unlocks (recorded), and seeds per task
    seeds_of = {(a, b): sd for ty, a, b, _, _, sd in tasks if ty == "N"}
    unlock_true = collections.defaultdict(float)
    for a, chain in deps.items():
        for b in chain:
            unlock_true[(a, b)] = max(unlock_true[(a, b)], gcost[key(a, a)])

    def static(own, boost=None, cheap_first=False):
        per = [dict(E=[], L=0.0, S=[], C=0.0) for _ in range(P)]
        for ty, a, b, s, h, _ in tasks:
            d = per[own[key(a, b)]]
            if ty != "N":
                d["C"] += s
            elif a == b:
                d["S"].append((a, s))
            elif pos[b] < pos[a]:
                first = -boost[(a, b)] if boost is not None else 0.0
                tie = seeds_of[(a, b)] if cheap_first else -clen.get(a, 0)
                d["E"].append((first, near[(a, b)], tie, a, b, s))
            else:
                d["L"] += s
        free = max(sum(e[5] for e in d["E"]) + d["L"] + sum(s for _, s in d["S"]) + d["C"]
                   for d in per)
        done = {}
        for d in per:
            t = 0.0
            for _, _, _, a, b, s in sorted(d["E"]):
                t += s
                done[(a, b)] = t
        # rank time before its final rounds, for the earlier-row case
        def ready_of(a):
            return max((done[(a, b)] for b in deps[a]), default=0.0)
        pre = []
        for d in per:
            t = sum(e[5] for e in d["E"]) + d["L"]
            later = []
            for a, s in d["S"]:
                if a in rowcase:
                    later.append((None, s))
                elif ready_of(a) <= t:
                    t += s
                else:
                    later.append((ready_of(a), s))
            t += d["C"]
            for r, s in sorted((x for x in later if x[0] is not None)):
                t = max(t, r) + s
            pre.append((t, [s for r, s in later if r is None]))
        rounds_start = max(t for t, _ in pre)
        ordered = max((rounds_start + sum(rs)) if rs else t for t, rs in pre)
        return free, ordered

    def dynamic(mode):
        # mode 0 free, 1 ordered, 2 ordered + boost
        unlock = collections.defaultdict(float)
        need = {}
        users = collections.defaultdict(list)      # dep group -> self groups
        for a, chain in deps.items():
            g = key(a, a)
            need[g] = len({key(a, b) for b in chain})
            for d in {key(a, b) for b in chain}:
                users[d].append(g)
                unlock[d] = max(unlock[d], gcost[g])
        rowg = {key(a, a) for a in rowcase} if mode else set()
        prio = {g: gcost[g] + (unlock[g] if mode == 2 else 0.0) for g in gcost}
        left = dict(need) if mode else {}
        ready = [(-prio[g], g) for g in gcost if g not in rowg and left.get(g, 0) == 0]
        heapq.heapify(ready)
        others = len(gcost) - len(rowg)              # groups that are not earlier-row self pairs
        running, free_workers, t, done_others = [], P, 0.0, 0
        while ready or running:
            while free_workers and ready:
                _, g = heapq.heappop(ready)
                heapq.heappush(running, (t + gcost[g], g))
                free_workers -= 1
            t, g = heapq.heappop(running)
            free_workers += 1
            if g not in rowg:
                done_others += 1
                if mode and done_others == others:   # earlier-row values: after everything else
                    for r in rowg:
                        heapq.heappush(ready, (-prio[r], r))
            if mode:
                for u in users.get(g, []):
                    left[u] -= 1
                    if left[u] == 0 and u not in rowg:
                        heapq.heappush(ready, (-prio[u], u))
        return t

    print(f"  {'P':>4} {'bound':>6} | {'length: free':>12} {'ordered':>7} {'cheap':>6} | "
          f"{'oracle: free':>12} {'ordered':>7} {'cheap':>6} {'+cost':>6} | "
          f"{'dynamic: free':>13} {'ordered':>7} {'+cost':>6}")
    for P in ranks:
        la, oa = length_assign(set(clen), clen, P), oracle_assign(gcost, P)
        lf, lo = static(la)
        _, lc = static(la, cheap_first=True)
        of, oo = static(oa)
        _, oh = static(oa, cheap_first=True)
        _, oc = static(oa, unlock_true)
        d0, d1, d2 = dynamic(0), dynamic(1), dynamic(2)
        print(f"  {P:>4} {max(total / P, biggest):>6.0f} | {lf:>12.0f} {lo:>7.0f} {lc:>6.0f} | "
              f"{of:>12.0f} {oo:>7.0f} {oh:>6.0f} {oc:>6.0f} | {d0:>13.0f} {d1:>7.0f} {d2:>6.0f}")


if __name__ == "__main__":
    main()
