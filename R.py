import numpy as np
import scipy.signal as scs


class NPoly:
    def __init__(self, data):
        if isinstance(data, np.ndarray):
            msp = np.where(np.flip(data) != 0)[0]
            if len(msp) == 0:
                self.n = 0
                self.a = np.zeros(1)
            else:
                self.n = len(data) - msp[0] - 1
                self.a = data[0:self.n+1]
        else:
            raise "p!Ш0L n@x1y"
        
    def __str__(self):
        return np.array2string(self.a)

    def __mul__(self, other):
        res = scs.fftconvolve(self.a, other.a)
    
        return NPoly(np.round(res))
    
    def __add__(self, other):
        return NPoly(self.a + other.a)


class R:
    def __init__(self, data):
        if isinstance(data, np.ndarray):
            self.a = data
            self.n = len(data)
        else:
            raise "p!Ш0L n@x1y"
    
    def __str__(self):
        return np.array2string(self.a)

    def __mul__(self, other):
        if not (self.n == other.n):
            raise RuntimeError("Can't multiply polynomials from different rings")
        
        n = self.n

        res = np.zeros(n)
        other_a = np.flip(other.a)
        for k in range(n):
            res[k] = self.a.dot(np.roll(other_a, -n+k+1))

        return R(res)
    
    def _msp(self):
        return self.n - np.argmax(np.flip(self.a) != 0) - 1

    def get_inv(self, q: int):
        b = NPoly(self.a)

        a = np.zeros(self.n+1)
        a[0] = -1
        a[-1] = 1
        a = NPoly(b)


        u_a = 1
        v_a = 0

        u_b = 0
        v_b = 1

        
        return self

    def __mod__(self, q: int):
        if not (q > 1):
            raise RuntimeError("Incorrect modulus")
        
        return R(self.a % q)

    def __add__(self, other):
        if not (self.q == other.q and self.n == other.n):
            raise RuntimeError("Can't multiply polynomials from different rings")
        
        return R(self.q, self.a + other.a)
        



def solve_NTRU_eq(f, g, q):
    F = R(f.a)
    G = R(g.a)

    return F, G


# a1 = np.arange(4)
# a2 = np.arange(4)
# n = 4
# for k in range(4):
#     print(a1 + np.roll(np.flip(a2), -n+k+1))


a = NPoly(np.array([1, -3, 0, 0, 0]))
b = NPoly(np.array([1, 0, 0, 5]))

print(a)
print(b)

c = a * b

print(c)