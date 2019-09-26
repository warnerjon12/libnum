import math

try:
    import __builtins__
    basestring = basestring
    xrange = xrange
    long = long
    builtins = __builtins__

    def math_log2(x):
        return math.log(x) / math.log(2)
except NameError:
    import builtins
    basestring = str
    xrange = range
    long = int
    builtins = builtins
    math_log2 = math.log2
