from Poly import R

def solve_NTRU_eq(f, g, q):
    F = R(f.a)
    G = R(g.a)

    return F, G

import random as rng
from enum import Enum
class NTRU_Methods(Enum):
    STANDART = 1
    TRANSPOSE = 2

def H(Data) -> R:
    return R([0])

def norm(p1, p2) -> float:
    return 0.0

def sample_random_bin_poly(N, d):
    poly = rng.sample([0, 1], N, counts=[N - d, d])
    print(poly)
    return R(poly)

def NTRU_init_parameters(N, q, d_f, d_g, B, delta, t=NTRU_Methods.TRANSPOSE):
    # Generate B private Lattice bases and one public Lattice basis
    lattice_f = []
    lattice_f_shtrih = []
    lattice_h = []
    for _ in range(B, -1, -1):
        f, g = [sample_random_bin_poly(N, d_f), sample_random_bin_poly(N, d_g)]
        F, G = solve_NTRU_eq(f, g, q)
        if t is NTRU_Methods.STANDART:
            lattice_f.append(f)
            lattice_f_shtrih.append(F)
        elif t is NTRU_Methods.TRANSPOSE:
            lattice_f.append(f)
            lattice_f_shtrih.append(g)
        else:
            raise RuntimeError(f"Incorrect NTRU method {t}!")
        
        h = (lattice_f[-1].get_inv(q) * lattice_f_shtrih[-1]) % q
        lattice_h.append(h)
    
    pub_out = [N, q, d_f, d_g, B, t, delta, lattice_h[-1]]
    private_out = [lattice_f, lattice_f_shtrih, lattice_h]
    return pub_out, private_out

def NTRU_sign(Data: bytes, private_key, pub_key):
    f, f_shtrih, h = private_key
    N, q, _, _, B, _, delta = pub_key[:-1]
    r = 0

    while True:    
        s = R([0], N)
        m = m_0 = H(Data + r.to_bytes(16, 'little'))

        # Perturb the point using the private lattices
        for i in range(0, B):
            x = (m * f_shtrih[i]).float_scalar_mult(-1/q)
            y = (m * f[i]).float_scalar_mult(1/q)
            s_i = (x * f[i] + y * f_shtrih[i])
            m = s_i * (h[i] - h[i + 1]) % q
            s = s + s_i

        # Sign the perturbed point using the public lattice
        x = (m * f_shtrih[-1]).float_scalar_mult(-1/q)
        y = (m * f[-1]).float_scalar_mult(1/q)
        s_0 = x * f[-1] + y * f_shtrih[-1]
        s = s + s_0

        # Check the signature
        b = norm(s, (s * h - m_0) % q)
        if b < delta:
            break
        else:
            r += 1
    
    return [Data, r, s]

def NTRU_verify(signature_pack, pub_key):
    Data, r, s = signature_pack
    N, q, _, _, B, _, delta, h = pub_key

    m = H(Data + r.to_bytes(16, 'little'))
    
    b = norm(s, (s * h - m) % q)
    return b < delta