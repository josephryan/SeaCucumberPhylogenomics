#!/usr/bin/python

# version 0.02

import sys, getopt
from ete3 import Tree

def main(argv):
   inputfile = ''
   try:
      opts, args = getopt.getopt(argv,"ht:p:o:",["tree=","prune=","outtree="])
   except getopt.GetoptError:
      print 'prune.py -t <treefile> -p <comma_sep_leaves2prune> -o <outtree>'
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print 'prune.py -t <tree> -p <comma_sep_leaves2prune> -o <outtree>'
         sys.exit()
      elif opt in ("-t", "--tree"):
         tree = arg
      elif opt in ("-p", "--prune"):
         prune = arg
      elif opt in ("-o", "--outtree"):
         outtree = arg

   t = Tree(tree)
   list = prune.split(',')
   t.prune(list)
   print t.write(outfile=outtree)

if __name__ == "__main__":
   main(sys.argv[1:])


