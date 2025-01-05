from Poly import *
import random as rng
from enum import Enum
class NTRU_Methods(Enum):
    STANDART = 1
    TRANSPOSE = 2

def sample_random_bin_poly(N, d):
    poly = rng.sample([0, 1], N, counts=[N - d, d])
    return R(N, poly)

def sample_oboroten_bin_poly(N, d, p, e):
    i = 0
    while True:
        if i % 10 == 0:
            print(f'try: {i}')
        i += 1
        poly = sample_random_bin_poly(N, d)
        if poly.get_inv_q(p, e) != R(N, np.zeros(N)):
            return poly


def NTRU_init_parameters(N, q, d_f, d_g, B, delta, t=NTRU_Methods.TRANSPOSE):
    # Generate B private Lattice bases and one public Lattice basis
    lattice_f = []
    lattice_f_shtrih = []
    lattice_h = []
    for _ in range(B, -1, -1):
        f, g = [sample_oboroten_bin_poly(N, d_f, q[0], q[1]), sample_oboroten_bin_poly(N, d_g, q[0], q[1])]
        F, G = solve_NTRU_eq(f, g, (q[0]**q[1]))
        if t is NTRU_Methods.STANDART:
            lattice_f.append(f)
            lattice_f_shtrih.append(F)
        elif t is NTRU_Methods.TRANSPOSE:
            lattice_f.append(f)
            lattice_f_shtrih.append(g)
        else:
            raise RuntimeError(f"Incorrect NTRU method {t}!")
        
        h = (lattice_f[-1].get_inv_q(q[0], q[1]) * lattice_f_shtrih[-1]) % (q[0]**q[1])
        lattice_h.append(h)
    
    pub_out = [N, (q[0]**q[1]), d_f, d_g, B, t, delta, lattice_h[-1]]
    private_out = [lattice_f, lattice_f_shtrih, lattice_h]
    return pub_out, private_out

def NTRU_sign(Data: bytes, private_key, pub_key):
    f, f_shtrih, h = private_key
    N, q, _, _, B, _, delta = pub_key[:-1]
    r = 0

    while True:    
        s = R(N, [0]*N)
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
        b = norm(s, (s * h[-1] - m_0) % q)
        print(f"b = {b}")
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


