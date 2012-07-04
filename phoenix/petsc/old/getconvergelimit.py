#!/usr/bin/env python

import commands
from numpy import array, double, int
from getopt import getopt, GetoptError
from sys import argv, exit

usage="""
GETCONVERGELIMIT.PY  Do bisection to find convergence limit on implicit time
step for porous.c.

Example:
  ./getconvergelimit.py --left=30 --right=3010
"""

def usagefailure(message):
    print message
    print
    print usage
    exit(2)

def evalcase(tenddays):
    """runs porous.c with -por_tend_days <tenddays>"""

    # change this to "mpiexec -n 8" or similar to run on multiple processes
    mpido=""
    theopts = " -da_grid_x 21 -da_grid_y 21 -por_dosteps 1 -por_tend_days "

    command = mpido + " ./porous " + theopts + str(tenddays)
    print "  doing:"
    print "    " + command
    try:
      (status,output) = commands.getstatusoutput(command)
    except KeyboardInterrupt:
      exit(2)
    print "status = " + str(status)
    if status == 0:
      return +1.0
    else:
      return -1.0
    

if __name__ == "__main__":

    try:
      opts, args = getopt(argv[1:], "", 
                          ["left=", "right="])
    except GetoptError:
      usagefailure('ERROR: INCORRECT COMMAND LINE ARGUMENTS FOR getconvergelimit.py')
    lefttenddays = 30.0
    righttenddays = 3010.0
    for (opt, optarg) in opts:
        if opt in ("--left"):
            lefttenddays = float(optarg)
        if opt in ("--right"):
            righttenddays = float(optarg)
        if opt in ("--help", "--usage"):
            print usage
            exit(0)

    # line search:  combines bisection with false position at end
    F_left = evalcase(lefttenddays)
    F_right = evalcase(righttenddays)

    print F_left
    print F_right

    if F_left * F_right > 0.0:
      print "************ NO BRACKET ***********"
    else:
      count = 2
      while True:
        c = 0.5 * (lefttenddays + righttenddays)
        result = evalcase(c)
        count = count + 1
        if count >= 15:
          break
        if F_left * result > 0.0:
          lefttenddays = c
          F_left = result
        else:
          righttenddays = c
          F_right = result

