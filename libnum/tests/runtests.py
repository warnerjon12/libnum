#!/usr/bin/env python

"""
python runtests.py -py
  Use py.test to run tests (more useful for debugging)

python runtests.py -coverage
  Generate test coverage report. Statistics are written to /tmp

python runtests.py -profile
  Generate profile stats (this is much slower)

python runtests.py -local
  Insert "../.." at the beginning of sys.path to use local libnum

python runtests.py -verbose
  Do not suppress printing from within tests

Additional arguments are used to filter the tests to run. Only files that have
one of the arguments in their name are executed.

"""

import sys
import os
import traceback

profile = False
if "-profile" in sys.argv:
    sys.argv.remove("-profile")
    profile = True

coverage = False
if "-coverage" in sys.argv:
    sys.argv.remove("-coverage")
    coverage = True

if "-local" in sys.argv:
    sys.argv.remove("-local")
    importdir = os.path.abspath(os.path.join(os.path.dirname(sys.argv[0]),
                                             "../.."))
else:
    importdir = ""

verbose = False
if "-verbose" in sys.argv:
    sys.argv.remove("-verbose")
    verbose = True

# TODO: add a flag for this
testdir = ""

def testit(importdir="", testdir=""):
    """Run all tests in testdir while importing from importdir."""
    if importdir:
        sys.path.insert(1, importdir)
    if testdir:
        sys.path.insert(1, testdir)
    import os.path
    import libnum
    print("libnum imported from %s" % os.path.dirname(libnum.__file__))
    print("libnum version: %s" % libnum.__version__)
    print("Python version: %s" % sys.version)
    print("")
    if "-py" in sys.argv:
        sys.argv.remove("-py")
        import py
        py.test.cmdline.main()
    else:
        import glob
        from timeit import default_timer as clock
        modules = []
        args = sys.argv[1:]
        # search for tests in directory of this file if not otherwise specified
        if not testdir:
            pattern = os.path.dirname(sys.argv[0])
        else:
            pattern = testdir
        if pattern:
            pattern += '/'
        pattern += "test*.py"
        # look for tests (respecting specified filter)
        for f in glob.glob(pattern):
            name = os.path.splitext(os.path.basename(f))[0]
            # If run as a script, only run tests given as args, if any are given
            if args and __name__ == "__main__":
                ok = False
                for arg in args:
                    if arg in name:
                        ok = True
                        break
                if not ok:
                    continue
            module = __import__(name)
            priority = module.__dict__.get("priority", 100)
            if priority == 666:
                modules = [[priority, name, module]]
                break
            modules.append([priority, name, module])
        # execute tests
        modules.sort()
        tstart = clock()
        for priority, name, module in modules:
            print(name)
            for f in sorted(module.__dict__.keys()):
                if f.startswith("test_"):
#                     if coverage and ("numpy" in f):
#                         continue
                    sys.stdout.write("    " + f[5:].ljust(25) + " ")
                    t1 = clock()

                    if not verbose:
                        from nostdout import nostdout
                    else:
                        from nostdout import stdout as nostdout

                    try:
                        with nostdout():
                            module.__dict__[f]()
                    except:
                        etype, evalue, trb = sys.exc_info()
                        if etype in (KeyboardInterrupt, SystemExit):
                            raise
                        print("")
                        print("TEST FAILED!")
                        print("")
                        traceback.print_exc()
                    t2 = clock()
                    print("ok " + "       " + ("%.7f" % (t2-t1)) + " s")
        tend = clock()
        print("")
        print("finished tests in " + ("%.2f" % (tend-tstart)) + " seconds")
        # clean sys.path
        if importdir:
            sys.path.remove(importdir)
        if testdir:
            sys.path.remove(testdir)

if __name__ == "__main__":
    if profile:
        import cProfile
        cProfile.run("testit('%s', '%s')" % (importdir, testdir), sort=1)
    elif coverage:
        import trace
        tracer = trace.Trace(ignoredirs=[sys.prefix, sys.exec_prefix],
            trace=0, count=1)
        tracer.run("testit(importdir, testdir)")
        r = tracer.results()
        r.write_results(show_missing=True, summary=True, coverdir="/tmp")
    else:
        testit(importdir, testdir)
