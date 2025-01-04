import numpy as np
import scipy.signal as scs
from Crypto.Hash import SHA1 

class NPoly:
    def __init__(self, data):
        if isinstance(data, np.ndarray) or isinstance(data, list):
            msp = np.where(np.flip(data) != 0)[0]
            if len(msp) == 0:
                self.n = 0
                self.a = np.zeros(1, dtype=np.int32)
            else:
                self.n = len(data) - msp[0] - 1
                self.a = np.round(np.array(data[0:self.n+1]))
        elif isinstance(data, dict):
            self.n = max(data.keys())
            self.a = np.zeros(self.n+1, dtype=np.int32)
            for (p, c) in data.items():
                if not (int(p) == p and p >= 0 and int(c) == c):
                    raise RuntimeError("Incorrect dict: (power, coefficent) pairs")

                self.a[p] = c

        else:
            raise "p!Ш0L n@x1y"
        
    def __str__(self):
        s = ""
        for p in range(self.n, -1, -1):
            if self.a[p] == 0:
                continue

            s = s + (f"{self.a[p]}x^{p} + " if p != 0 else f"{self.a[p]}")
                 
        return s

    def __getitem__(self, key):
        return self.a[key]

    def __eq__(self, other):
        return (np.array_equal(self.a, other.a) and self.n == other.n)

    def __mod__(self, q: int):
        if not (q > 1):
            raise RuntimeError("Incorrect modulus")
        
        return NPoly(np.round(self.a % q).astype(np.int32))

    def __mul__(self, other):
        res = scs.fftconvolve(self.a, other.a)
    
        return NPoly(np.round(np.real(res)).astype(np.int32))
    
    def __add__(self, other):
        if self.n < other.n:
            res = other.a
            res[:len(self.a)] += self.a
        else:
            res = self.a
            res[:len(other.a)] += other.a

        return NPoly(np.round(res).astype(np.int32))

    def __sub__(self, other):
        return self + NPoly(-other.a)

    def div_mod(self, other, mod):
        # a = q*b + r
        a = self % mod
        b = other % mod
        q = NPoly({0: 0})

        while (a.n >= b.n) and (a != NPoly({0: 0})):
            qi = NPoly({a.n - b.n :a[a.n]*pow(int(b[b.n]), -1, mod)}) % mod
            a = (a - qi*b) % mod
            q = q + qi

        return q, a


class R:
    def __init__(self, n: int, data):
        if not (int(n) == n and n > 1):
            raise "p!Ш0L n@x1y"

        self.n = n

        if isinstance(data, np.ndarray) or isinstance(data, list):
            if len(data) > n: 
                raise RuntimeError("data is too big for that modulus")

            self.a = np.zeros(self.n, dtype=np.int32)
            self.a[:len(data)] = np.round(np.array(data))
        elif isinstance(data, dict):
            m = max(data.keys())
            if m >= n: 
                raise RuntimeError("max power is too big for that modulus")

            self.a = np.zeros(self.n, dtype=np.int32)
            for (p, c) in data.items():
                if not (int(p) == p and p >= 0 and int(c) == c):
                    raise RuntimeError("Incorrect dict: (power, coefficent) pairs")

                self.a[p] = c

        else:
            raise "p!Ш0L n@x1y"
    

    def __str__(self):
        return np.array2string(self.a)

    def __getitem__(self, key):
        return self.a[key]

    def __mul__(self, other):
        if not (self.n == other.n):
            raise RuntimeError("Can't multiply polynomials from different rings")
        
        n = self.n

        res = np.zeros(n)
        other_a = np.flip(other.a)
        for k in range(n):
            res[k] = self.a.dot(np.roll(other_a, -n+k+1))

        return R(n, res)

    def _msp(self):
        return self.n - np.argmax(np.flip(self.a) != 0) - 1

    def get_inv(self, mod: int):
        b = NPoly(self.a) % mod
        a = NPoly({self.n: 1, 0: -1}) % mod

        r = [a, b]
        u_a = [NPoly([1]), NPoly([0])]
        v_b = [NPoly([0]), NPoly([1])]


        while (r[-1] != NPoly({0: 0})) and (r[-1] != NPoly({0: 1})):
            q, ri = r[-2].div_mod(r[-1], mod)
            r.append(ri % mod)
            u_a.append((u_a[-2] - q*u_a[-1]) % mod)
            v_b.append((v_b[-2] - q*v_b[-1]) % mod)

            print(r[-1])


        if r[-1] == NPoly({0: 1}):
            coeffs = np.zeros(self.n)
            v = v_b[-1]
            if v.n >= self.n:
                _, v = v.div_mod(a, mod)
            coeffs[:len(v.a)] = v.a
            return R(self.n, coeffs)
        else:
            return R(self.n, np.zeros(self.n))

    def float_scalar_mult(self, s):
        return R(self.n, np.round(self.a * s))

    def _norm(self):
        return np.sqrt((self.a.dot(self.a) - (sum(self.a)**2 / self.n)))

    def __mod__(self, q: int):
        if not (q > 1):
            raise RuntimeError("Incorrect modulus")
        
        return R(self.n, self.a % q)

    def __add__(self, other):
        if not (self.n == other.n):
            raise RuntimeError("Can't multiply polynomials from different rings")
        
        return R(self.n, self.a + other.a)
        

def H(m):
    d1 = SHA1.new(data=m).digest()

    df = SHA1.new(d1 + bytes([0]))
    for i in range(1, 13):
        df += SHA1.new(d1 + bytes([i]))

    a = []
    for b in df:
        a.append(b & 127)

    return R(251, a)
    

def norm(p1, p2):
    return np.sqrt(p1._norm()**2 + p2._norm()**2)


def solve_NTRU_eq(f, g, q):
    F = R(f.a)
    G = R(g.a)

    return F, G


# a1 = np.arange(4)
# a2 = np.arange(4)
# n = 4
# for k in range(4):
#     print(a1 + np.roll(np.flip(a2), -n+k+1))

a = R(6,np.array([1, 5, 1, 1, 5, 6]))
b = NPoly({0:1, 3:5, 4:1})

print(a.get_inv(11))

print((a * a.get_inv(11)) % 11)
# print(f'a =  {a}')
# print(f'b =  {b}\n')

# q, r = a.div_mod(b, 11)

# print('a = q*b + r')
# print(f'q =  {q}')
# print(f'r =  {r}')

# print('\n Check')

# print(f'q*b + r =  {(q*b + r) % 11}')
