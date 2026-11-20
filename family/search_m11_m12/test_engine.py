import sys, itertools, numpy as np, random
sys.argv = ['x']; exec(open('ls420.py').read().split("if __name__")[0])
# 1) row-check logic on the three known graphs (m = 13,14,15): 16 S = (4/9) H (f + m delta) must be 0/16-valued
W = '/mnt/jinstore/JinLab02/dmp131/srg/brouwer/work/'
for m, fn in [(13, W+'transfer525_20261116/results/455/connections_Z6_4.txt'), (14, W+'transfer525_20261116/results/490/connections_Z6_4.txt'), (15, W+'construct525_20261116/result/connections_Z6_4.txt')]:
    T = {a: set() for a in range(16)}
    for l in open(fn):
        z = tuple(int(x) for x in l.split()); b = tuple(x % 2 for x in z); v = tuple(x % 3 for x in z)
        T[int(''.join(map(str, b)), 2)].add(v)
    # here A index = binary number with standard dot product; build H for that labeling
    Hs = np.array([[(-1) ** bin(al & a).count('1') for a in range(16)] for al in range(16)])
    ok = True
    for b in range(40):
        f = np.array([sum(1 if sum(x*y for x, y in zip(P[b], t)) % 3 == 0 else 0 for t in T[a]) * 1.5 - len(T[a]) / 2 for a in range(16)])
        # character sum of symmetric set: sum cos = #(dot=0) - (#rest)/2 ; equivalently (3 n0 - n)/2
        y = f.copy(); y[0] += m
        S16 = (4 / 9) * (Hs @ y)
        ok &= bool(np.all(np.isclose(S16, 0) | np.isclose(S16, 16)))
    print('row check on known graph m=%d:' % m, ok)
# 2) plumbing: f_beta(h a) computed from reps via dual perms equals the direct character of phi(h) T_a
random.seed(1)
g1 = ((1,1,0,0),(0,1,0,0),(0,0,1,0),(0,0,0,1)); g2 = ((1,0,0,0),(0,1,0,0),(0,0,1,1),(0,0,0,1))
def apply(M, v): return tuple(sum(M[i][k]*v[k] for k in range(4)) % 3 for i in range(4))
T = {tuple(random.choice(V3)) for _ in range(20)}; mask = 0
for t in T:
    if any(t): mask |= 1 << pidx[norm(t)]
ch = chars(mask, False)
Ts = {p for j, p in enumerate(P) if mask >> j & 1}; Ts |= {tuple((-x) % 3 for x in p) for p in Ts}
dp = perm_dual(g1)
direct = [sum(1 if sum(x*y for x, y in zip(P[b], apply(g1, t))) % 3 == 0 else -0.5 for t in Ts) * 2 for b in range(40)]
print('plumbing (dual perm) ok:', np.allclose(ch[dp], direct))
print(ch[dp][:10], np.array(direct[:10]))
print(chars(mask,False)[:5], [sum(2 if sum(x*y for x,y in zip(P[b],t))%3==0 else -1 for t in Ts) for b in range(5)], len(Ts), bin(mask).count('1'))
