"""Verify the exported files without importing the generator or any packages."""
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
import sys


def read_graph(path):
    raw = path.read_bytes()
    values = [c - 63 for c in raw.rstrip(b'\n')]
    assert values[0] == 63 and all(0 <= t <= 63 for t in values)
    n = 4096 * values[1] + 64 * values[2] + values[3]
    nbits = n * (n - 1) // 2
    assert len(values) == 4 + (nbits + 5) // 6
    rows = [0] * n
    bit = 0
    for j in range(n):
        for i in range(j):
            if values[4 + bit // 6] & (32 >> (bit % 6)):
                rows[i] |= 1 << j
                rows[j] |= 1 << i
            bit += 1
    assert nbits % 6 == 0 or values[-1] % (1 << (-nbits % 6)) == 0
    return rows, sha256(raw).hexdigest()


def local_cliques(rows, vertex):
    # A triangle among the neighbours completes a 4-clique on vertex.
    remaining = rows[vertex]
    triangles = 0
    while remaining:
        x = remaining & -remaining
        remaining ^= x
        linked = remaining & rows[x.bit_length() - 1]
        while linked:
            y = linked & -linked
            linked ^= y
            triangles += (linked & rows[y.bit_length() - 1]).bit_count()
    return triangles


def verify(path):
    rows, digest = read_graph(path)
    n = len(rows)
    degrees = {r.bit_count() for r in rows}
    assert len(degrees) == 1
    k = degrees.pop()
    counts = [Counter(), Counter()]
    for j, row in enumerate(rows):
        for i in range(j):
            counts[(row >> i) & 1][(row & rows[i]).bit_count()] += 1
    assert all(len(c) == 1 for c in counts)
    mu, lam = [next(iter(c)) for c in counts]
    assert counts[0][mu] == n * (n - k - 1) // 2 and counts[1][lam] == n * k // 2
    # Check the four binary translations on the actual graph6 labelling.
    for t in (1, 2, 4, 8):
        lower = sum(1 << i for i in range(n) if (i & t) == 0)
        for i, r in enumerate(rows):
            assert ((r & lower) << t | (r >> t) & lower) == rows[i ^ t]
    local = [local_cliques(rows, v) for v in range(0, n, 16)]
    census = {value: 16 * local.count(value) for value in sorted(set(local))}
    return dict(file=path.name, sha256=digest, parameters=[n,k,lam,mu],
                adjacent_pairs=counts[1][lam], nonadjacent_pairs=counts[0][mu],
                four_cliques_through_vertex=census,
                representative_vertices={value: 16 * local.index(value) for value in census})


if __name__ == '__main__':
    for name in sys.argv[1:]:
        print(json.dumps(verify(Path(name))), flush=True)
