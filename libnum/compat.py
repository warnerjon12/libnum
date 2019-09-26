import math

try:
    basestring = basestring
    xrange = xrange
    long = long

    def math_log2(x):
        return math.log(x) / math.log(2)
except NameError:
    basestring = str
    xrange = range
    long = int
    math_log2 = math.log2
