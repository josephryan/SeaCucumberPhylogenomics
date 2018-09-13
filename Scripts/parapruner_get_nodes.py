#!/usr/bin/python

import sys, getopt
from ete3 import Tree

sys.tracebacklimit = 0

def main(argv):
   inputfile = ''
   try:
      opts, args = getopt.getopt(argv,"ht:",["tree="])
   except getopt.GetoptError:
      print 'test.py -t <treefile>'
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print 'test.py -t <tree>'
         sys.exit()
      elif opt in ("-t", "--tree"):
         tree = arg

   t = Tree(tree)
   for leaf in t:
      print leaf.name

if __name__ == "__main__":
   main(sys.argv[1:])


