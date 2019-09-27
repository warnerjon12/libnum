#-*- coding:utf-8 -*-

from libnum.compat import xrange
from libnum.primes import *
from libnum.tests.primes10000 import primes10000
from utcompat import *


def test_primes():
    p = primes(primes10000[-1])
    assertEqual(len(p), len(primes10000))
    assertEqual(p[0], 2)
    assert all(p[i] == pt for i, pt in enumerate(primes10000))

    assertEqual(list(primes(-1)), [])
    assertEqual(list(primes(1)), [])
    assertEqual(list(primes(2)), [2])
    assertRaises(TypeError, primes, "")
    assertRaises(TypeError, primes, 1000000, "fake")


def test_genprime():
    for size in (2, 10, 64, 128, 129, 256):
        for ntry in xrange(10):
            p = generate_prime(size, k=25)
            assertEqual(len_in_bits(p), size)

            assertTrue(prime_test_fermat(p, k=25))
            assertTrue(prime_test_solovay_strassen(p, k=25))
            assertTrue(prime_test_miller_rabin(p, k=25))
            assertTrue(prime_test(p, k=25))

    assertRaises(ValueError, generate_prime, 1)
    assertRaises(TypeError, generate_prime, "")


def test_genprime_str():
    begin = "preved medved \xde\xad\xbe\xef\x00\x00\x00\x00"
    n = generate_prime_from_string(begin)
    s = n2s(n)
    assertTrue(s.startswith(begin))
    assertTrue(prime_test(n, 25))

    assertRaises(TypeError, generate_prime_from_string, 31337)
    assertRaises(ValueError, generate_prime_from_string, "test", 8)
    assertRaises(ValueError, generate_prime_from_string, "test", -8)


def do_test_prime_test(func):
    for p in (3, 1993, 17333, 1582541, 459167430810992879232575962113190418519):
        assertTrue(func(p, 50))

    for p in primes(1000):
        assertTrue(func(p, 50))

    for not_p in (4, 1994, 1995, 16231845893292108971):
        assertFalse(func(not_p, 50))

    assertRaises(TypeError, func, "test")


def test_fermatest():
    return do_test_prime_test(prime_test_fermat)


def test_solovay():
    return do_test_prime_test(prime_test_solovay_strassen)


def test_miller():
    return do_test_prime_test(prime_test_miller_rabin)
