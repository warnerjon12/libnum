#-*- coding:utf-8 -*-

"""
Some factorization methods are listed here
"""

import math
import random

from functools import reduce
from .compat import xrange
from .primes import primes, prime_test
from .common import _gcd, gcd, nroot, extract_prime_power

__all__ = "factorize unfactorize".split()


_PRIMES_P1 = primes(100)


def rho_pollard_reduce(n, f):
    # use Pollard's (p-1) method to narrow down search
    a = random.randint(2, n - 2)
    for p in _PRIMES_P1:
        a = pow(a, p, n)

    b = a
    Q = 1
    while True:
        a = f(a)
        b = f(f(b))
        Q = (Q * (b - a)) % n

        g = gcd(Q, n)
        if g == n:
            a = b = random.randint(2, n - 2)
            Q = 1
        elif g != 1:
            return g


def brent_reduce(n, x0=0, m=1, f=None):
    # Brent's improved Pollard-Rho implementation
    if f is None:
        f = lambda x: (x ** 2 + 3) % n

    y, r, q = x0, 1, 1
    while True:
        x = y
        for i in xrange(r):
            y = f(y)
        k = 0
        while True:
            ys = y
            for i in xrange(min(m, r - k)):
                y = f(y)
                q = (q * (x - y)) % n
            G = _gcd(q, n)
            k += m
            if k >= r or G > 1:
                break
        r *= 2
        if G > 1:
            break
    if G == n:
        while True:
            ys = f(ys)
            G = _gcd(abs(x - ys), n)
            if G > 1:
                break
    return G


def factorize(n, td_primes=primes(100)):
    """
    Return {p: e for p**e in prime_powers(n)}.
    """
    prime_factors = {}

    if n < 0:
        n = -n
        prime_factors[-1] = 1

    if n == 0:
        prime_factors[0] = 1
        return prime_factors

    for p in td_primes:
        e, n = extract_prime_power(n, p)
        if e > 0:
            prime_factors[p] = e

    def combine(d1, d2):
        for k in d2:
            d1[k] = d1.get(k, 0) + d2[k]
        return d1

    def brent_factorize(n):
        if n == 1:
            return {}
        elif prime_test(n):
            return {n: 1}
        for c in xrange(1, min(0x7fffffff, n)):
            p = brent_reduce(n, f=lambda x: (x ** 2 + c) % n)
            if p != n:
                if prime_test(p):
                    e, n = extract_prime_power(n, p)
                    return combine({p: e}, brent_factorize(n))
                else:
                    bfactors = brent_factorize(p)
                    for p in bfactors:
                        e, n = extract_prime_power(n, p)
                        bfactors[p] = e
                    return combine(bfactors, brent_factorize(n))
        raise ValueError("Failed to factorize %d" % n)

    return combine(prime_factors, brent_factorize(n))


def unfactorize(factors):
    return reduce(lambda acc, p_e: acc * (p_e[0] ** p_e[1]),
                  factors.items(), 1)


def is_power(n):
    limit = int(math.log(n, 2))
    for power in xrange(limit, 1, -1):
        p = nroot(n, power)
        if pow(p, power) == n:
            return p, power
    return False
