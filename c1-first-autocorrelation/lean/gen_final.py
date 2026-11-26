# Build Challenge.lean and Solution.lean for the proof board  Constant1a.c1a_le_frontier : C1a ≤ 1.5026709
import re,sys,fractions,hashlib
import gmpy2
BW=60; q=20; e=15; JP=8; n=2**e*JP
BOUND='1.5026709'
v0=[int(x) for x in re.split(r'[\s,]+',open('C1_262144_1.5026648295.txt').read().strip())]
assert len(v0)==n
mx=max(v0)
# heights rounded to q bits: w_i = round(v_i * (2^q - 1) / max v), exact integer rounding (half up)
w=[(2*x*(2**q-1)+mx)//(2*mx) for x in v0]
open('heights_q20.txt','w').write('\n'.join(map(str,w))+'\n')
nb=BW//4
assert BW%4==0
P=int.from_bytes(b''.join(x.to_bytes(8,'little') for x in w),'little') if BW==64 else sum(0 for _ in [])
if BW!=64:
    hexs=''.join('%0*x'%(nb,x) for x in reversed(w)); P=int(hexs,16)
Q=int(gmpy2.mpz(P)*gmpy2.mpz(P))
qh='%0*x'%(2*n*nb,Q)
digs=[int(qh[len(qh)-(m+1)*nb:len(qh)-m*nb],16) for m in range(2*n)]
T=max(digs); S=sum(w)
assert sum(digs)==S*S and S*S<2**BW
ratio=fractions.Fraction(2*n*T,S*S)
print('n',n,'T',T,'S',S,'S^2 bits',(S*S).bit_length(),'ratio',float(ratio))
assert ratio<=fractions.Fraction(BOUND)
s=[sum(w[j*2**e:(j+1)*2**e]) for j in range(JP)]
# compact encoding of P: low 14 bits of every height in 14-bit slots (L), high parts as a sparse list (X)
LW=14
low=[x&((1<<LW)-1) for x in w]; high=[x>>LW for x in w]
Lnum=int(''.join('%04x'%0 for _ in [])+'0',16)
Lnum=0
lb=bytearray()
bits=''.join(format(x,'014b') for x in reversed(low)); Lnum=int(bits,2)
exc=[(i,h) for i,h in enumerate(high) if h]
assert max(h for _,h in exc)<64
ent=[]; prev=0
for i,h in exc:                      # ascending cell index; gap to the previous exception (or to cell 0)
    assert 0<=i-prev<2**18
    ent.append((i-prev)*64+h); prev=i
ED=12; assert len(ent)<=2**ED
ent+= [0]*(2**ED-len(ent))
Xnum=0
for e_ in reversed(ent): Xnum=(Xnum<<24)|e_
# reference implementation of the kernel-side unpacking (must reproduce P)
def stage(k,X):
    b=2**k
    R=((1<<(n*BW))-1)//((1<<(2*b*BW))-1)
    mlo=((1<<(b*LW))-1)*R
    return (X&mlo)+((X&(mlo<<(b*LW)))<<(b*(BW-LW)))
Xs=Lnum
for k in reversed(range(18)): Xs=stage(k,Xs)
def gapSum(d,lst):
    if d==0: return (lst>>6)%262144
    return gapSum(d-1,lst%(1<<(24*2**(d-1))))+gapSum(d-1,lst>>(24*2**(d-1)))
def excB(d,lst):
    if d==0: return (lst%64)<<(BW*((lst>>6)%262144))
    lo=lst%(1<<(24*2**(d-1)))
    return excB(d-1,lo)+(excB(d-1,lst>>(24*2**(d-1)))<<(BW*gapSum(d-1,lo)))
acc=excB(ED,Xnum)
assert Xs+(acc<<LW)==P, 'unpacking mismatch'
print('exceptions',len(ent),'L hex chars',len('%x'%Lnum),'X hex chars',len('%x'%Xnum))
CH=65536; ph='%x'%Lnum; hc=CH//4
chunks=[ph[max(0,len(ph)-(i+1)*hc):len(ph)-i*hc] for i in range((len(ph)+hc-1)//hc)]
BLK=16; blocks=[chunks[i:i+BLK] for i in range(0,len(chunks),BLK)]
data=[]
for bi,bl in enumerate(blocks):
    ex='0'
    for c in reversed(bl): ex=f'(0x{c} + ({ex}) <<< {CH})'
    data.append(f'noncomputable def blk{bi} : ℕ := {ex}')
ex='0'
for bi in reversed(range(len(blocks))): ex=f'(blk{bi} + ({ex}) <<< {CH*BLK})'
data.append(f'/-- low 14 bits of the heights, in 14-bit slots -/\nnoncomputable def L : ℕ := {ex}')
data.append(f'/-- the heights of 15 bits or more: 24-bit entries `gap * 64 + high part`, in increasing cell order -/\nnoncomputable def X : ℕ := 0x{Xnum:x}')
data.append(f"""/-- one doubling step of the unpacking: in every block of `2 * 2^k` slots the upper half moves up -/
noncomputable def stage (k x : ℕ) : ℕ :=
  (x &&& (((1 <<< (2 ^ k * 14)) - 1) * (((1 <<< (2 ^ 18 * 60)) - 1) / ((1 <<< (2 ^ (k + 1) * 60)) - 1)))) +
    ((x &&& ((((1 <<< (2 ^ k * 14)) - 1) * (((1 <<< (2 ^ 18 * 60)) - 1) / ((1 <<< (2 ^ (k + 1) * 60)) - 1)))
      <<< (2 ^ k * 14))) <<< (2 ^ k * 46))

noncomputable def unpack : ℕ → ℕ → ℕ
  | 0, x => x
  | k + 1, x => unpack k (stage k x)

/-- sum of the gap fields of the `2^d` lowest entries -/
noncomputable def gapSum : ℕ → ℕ → ℕ
  | 0, lst => (lst >>> 6) % 262144
  | d + 1, lst => gapSum d (lst % (1 <<< (24 * 2 ^ d))) + gapSum d (lst >>> (24 * 2 ^ d))

/-- the high parts of the `2^d` lowest entries, placed in 60-bit slots relative to the cell before the first -/
noncomputable def excB : ℕ → ℕ → ℕ
  | 0, lst => (lst % 64) <<< (60 * ((lst >>> 6) % 262144))
  | d + 1, lst => excB d (lst % (1 <<< (24 * 2 ^ d))) +
      (excB d (lst >>> (24 * 2 ^ d))) <<< (60 * gapSum d (lst % (1 <<< (24 * 2 ^ d))))

/-- The packed step function: `P = ∑ i < 2^18, v i * 2^(60 i)`. -/
noncomputable def P : ℕ := unpack 18 L + (excB {ED} X) <<< 14""")
data.append(f'/-- The largest coefficient of the autoconvolution, `max_m ∑_(i+j=m) v i * v j`. -/\ndef T : ℕ := {T}')
data.append('def sv (j : ℕ) : ℕ := '+' else '.join([f'if j = {j} then {s[j]}' for j in range(JP-1)]+[str(s[JP-1])]))
kern=[f'theorem kP : P < 1 <<< (BW * 2 ^ {e} * {JP}) := by decide +kernel']
kern+=[f'theorem kS{j} : sumDigits {e} (slice {e} P {j}) = {s[j]} := by decide +kernel' for j in range(JP)]
kern+=[f'theorem kQ{j} : checkAll T {e} (slice {e} (P * P) {j}) = true := by decide +kernel' for j in range(2*JP)]
final=f'''
theorem cert_facts :
    (∑ i ∈ range (2 ^ {e} * {JP}), dig P i = {S}) ∧ ∀ m, conv (2 ^ {e} * {JP}) (fun i => dig P i) m ≤ T := by
  have hP : P < ((2 : ℕ) ^ BW) ^ (2 ^ {e} * {JP}) := by
    have h := kP
    rwa [Nat.one_shiftLeft, mul_assoc, pow_mul] at h
  have hS : ∀ j < {JP}, sumDigits {e} (slice {e} P j) = sv j := by
    intro j hj
    interval_cases j
'''+''.join(f'    · exact kS{j}\n' for j in range(JP))+f'''  have hQ : ∀ j < 2 * {JP}, checkAll T {e} (slice {e} (P * P) j) = true := by
    intro j hj
    interval_cases j
'''+''.join(f'    · exact kQ{j}\n' for j in range(2*JP))+f'''  have hsv : ∑ j ∈ range {JP}, sv j = {S} := by simp [Finset.sum_range_succ, sv]
  have hB : (∑ j ∈ range {JP}, sv j) ^ 2 < 2 ^ BW := by
    rw [hsv]; norm_num [BW]
  obtain ⟨hsum, hconv⟩ := cert_core {e} {JP} P T sv hP hS hQ hB
  exact ⟨hsum.trans hsv, hconv⟩
'''
core=open('dev/core_nat.txt').read()+open('dev/core_ana.txt').read()+open('dev/core_asm.txt').read()
src=open('formal-conjectures_1a.lean').read()
anchor='-- TODO: Formalise relationship to Sidon sets.'
assert anchor in src
doc=f'''/-- An upper bound from an explicit step function with $2^{{18}}$ equal steps on $[-1/4,1/4]$:
the record of the Commonplace frontier board `autocorrelation-c1-large` (exact value $1.50266482948\\ldots$),
with its heights rounded to $20$ bits, which gives $1.50267087879\\ldots$.
This is below the bound $1.5028628587$ asked for in `mem_Ico_c1a`. -/
@[category research solved, AMS 5 11 26]
theorem c1a_le_frontier : C1a ≤ {BOUND} := by
'''
challenge=src.replace(anchor,doc+'  sorry\n\n'+anchor)
proof=f'''  obtain ⟨hsum, hconv⟩ := C1aCert.cert_facts
  have hpos : 0 < ∑ i ∈ Finset.range (2 ^ {e} * {JP}), C1aCert.dig C1aCert.P i := by rw [hsum]; norm_num
  have h := C1aCert.C1a_le_of_cert (2 ^ {e} * {JP}) (by positivity) (fun i => C1aCert.dig C1aCert.P i)
    C1aCert.T hconv hpos
  refine h.trans ?_
  rw [hsum]
  rw [div_le_iff₀ (by positivity)]
  norm_num [C1aCert.T]
'''
# the helper development goes before the theorem; C1aCert refers to Constant1a.C1a, so it is placed after the definition
helpers='namespace C1aCert\n\nset_option maxRecDepth 100000\nset_option maxHeartbeats 0\n\nnoncomputable section\n'+core+'\n'+'\n'.join(data)+'\n\n'+'\n'.join(kern)+'\n'+final+'\nend\n\nend C1aCert\n\n'
helpers=helpers.replace('Constant1a.C1a','C1a')
solution=src.replace(anchor,helpers+doc+proof+'\n'+anchor)
open('Challenge.lean','w').write(challenge); open('Solution.lean','w').write(solution)
print('Challenge bytes',len(challenge.encode()),'Solution bytes',len(solution.encode()),'(limit 4194304)')
print('sha256 heights',hashlib.sha256(open('heights_q20.txt','rb').read()).hexdigest()[:16])
