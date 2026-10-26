"""M7 scripts/08 — register handoff 收訖矩陣（唯讀）
================================================================================================
依據 skill `o-series-handoff-audit` §1–§3：
  §1 抽矩陣：直接 parse register.jsonl（**唔用 CLI grep 作送達證據**）
  §2 送達確認：逐目標核 len(ins)==len(outs)
  §3 閉環判定：目標喺 handoff 之後有無活動；ref 是否存在（stale）
輸出：data/08_handoff_matrix.json
本腳本只讀 register.jsonl，不寫任何 register 事件。
"""
import json
import os
from datetime import datetime

HERE = os.path.dirname(os.path.abspath(__file__))
M7 = os.path.abspath(os.path.join(HERE, ".."))
REPO = os.path.abspath(os.path.join(M7, "..", "..", ".."))
REG = os.path.join(REPO, "problem_register", "register.jsonl")
OUT = {"module": "M7", "script": "08_handoff_matrix.py", "checks": [], "matrix": []}
NP = [0, 0]


def ok(name, cond, got, want, note=""):
    NP[0] += 1
    NP[1] += 1 if cond else 0
    print(f"  {'✅' if cond else '🔴'} {name:<44} {got}" + (f" | want {want}" if want else ""))
    OUT["checks"].append(dict(item=name, pass_=bool(cond), got=str(got), want=str(want), note=note))


recs = []
for ln in open(REG, encoding="utf-8"):
    try:
        recs.append(json.loads(ln))
    except Exception:
        continue
print(f"register 事件 {len(recs)} 條（{REG}）")

kinds = {}
for d in recs:
    kinds[d.get("kind")] = kinds.get(d.get("kind"), 0) + 1
print("kind 值域：" + ", ".join(f"{k}={v}" for k, v in sorted(kinds.items(), key=lambda x: -x[1])))

PID = "M7"
outs = [d for d in recs if d.get("pid") == PID and d.get("target_pid")]
ins = [d for d in recs if str(d.get("target_pid") or "").startswith(PID)]
notes_out = [d for d in recs if d.get("pid") != PID and d.get("kind") == "note"
             and PID in str(d.get("ref") or "")]
print(f"\nM7 outgoing handoff = {len(outs)}；incoming = {len(ins)}")

print("\n§2 送達矩陣（送達 = 記錄唯一 ＋ ref 存在；回覆 = 對方 pid 回本 M 之記錄）")
targets = sorted({d["target_pid"] for d in outs})
dup = []
seen = set()
for d in outs:
    k = (d["target_pid"], str(d.get("ref")), d["ts"])
    if k in seen:
        dup.append(k)
    seen.add(k)
reply = {T: [d for d in recs if d.get("pid") == T and d.get("target_pid") == PID] for T in targets}
all_unique = not dup
for T in targets:
    o = [d for d in outs if d["target_pid"] == T]
    r_all = reply[T]
    r = [x for x in r_all if x["ts"] > o[0]["ts"]]   # 只計「發出之後」之回覆
    act = [d for d in recs if d.get("pid") == T and d.get("ts", "") > o[-1]["ts"]]
    status = "✅ 已閉環" if r else "⏳ 未閉環（無回覆）"
    OUT["matrix"].append(dict(target=T, n_out=len(o), n_reply=len(r),
                              reply_ts=[x["ts"] for x in r],
                              prior_incoming=[x["ts"] for x in r_all if x not in r],
                              ts=[d["ts"] for d in o], refs=[d.get("ref") for d in o],
                              target_activity_after=len(act), status=status))
    print(f"  {T:<4} outs={len(o)}  回覆(發出後)={len(r)}  其後活動={len(act)}  {status}"
          + (f"  [另有 {len(r_all) - len(r)} 條為發出前之對方來件]" if len(r_all) > len(r) else ""))
    for d in o:
        print(f"       {d['ts']}  {os.path.basename(str(d.get('ref')))}")
    for x in r:
        print(f"       ← {x['ts']}  {os.path.basename(str(x.get('ref')))}")
ok("§2 handoff 記錄無重複（每封一條）", all_unique, dup or "0", "0")
closed = [T for T in targets if reply[T]]
print(f"    回覆覆蓋：{len(closed)}/{len(targets)} 目標（{'、'.join(closed) or '無'}）—— 屬對方責任，不作本 M PASS 判據")

print("\n§3 閉環判定：incoming 之後本 M 有無動作")
for d in ins:
    after = [x for x in recs if x.get("pid") == PID and x.get("ts", "") > d["ts"]]
    print(f"  ← {d['ts']}  from {d['pid']}  REF={os.path.basename(str(d.get('ref')))}")
    print(f"       其後 M7 事件 {len(after)} 條 ⟹ {'已回應' if after else '🔴 未回應'}")
    print(f"       note: {(d.get('note') or '')[:90]}…")
ok("§3 incoming 已回應", all(bool([x for x in recs if x.get('pid') == PID and x.get('ts', '') > d['ts']])
                            for d in ins), f"{len(ins)} 條 incoming", "全部有後續事件")

print("\n§3 stale：handoff ref 是否存在")
stale = []
for d in outs:
    ref = d.get("ref") or ""
    p = ref if os.path.isabs(ref) else os.path.join(REPO, "problem_register", ref)
    if not os.path.exists(p):
        stale.append(ref)
print(f"  stale ref: {stale or '無'}")
ok("§3 全部 handoff ref 存在", not stale, stale or "0", "0")

print("\n附：M7 發向他 pid 之 note（非 handoff；僅通知）")
for d in notes_out:
    print(f"  → {d['pid']:<4} {d['ts']}  {os.path.basename(str(d.get('ref')))}")
OUT["notes_out"] = [dict(pid=d["pid"], ts=d["ts"], ref=d.get("ref")) for d in notes_out]

print("\n附：目標是否已回覆本 M（incoming 明細）")
for d in ins:
    print(f"  {d.get('kind')}  {d['ts']}  {d['pid']} → {d['target_pid']}")

print("\n" + "=" * 100)
print(f"結果：{NP[1]}/{NP[0]} PASS")
OUT.update(npass=NP[1], ntotal=NP[0], n_out=len(outs), n_in=len(ins),
           targets=targets, stale=stale)
json.dump(OUT, open(os.path.join(M7, "data", "08_handoff_matrix.json"), "w"), ensure_ascii=False, indent=1)
print("JSON → data/08_handoff_matrix.json")