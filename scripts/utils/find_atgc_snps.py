#!/usr/bin/python3
from sys import argv

bimfile = argv[-1]
assert bimfile.endswith('.bim')

with open(bimfile, 'r') as file:
    for line in file:
        line = line.strip().split()
        alleles = line[-2:]
        if sorted(alleles) in (['C','G'], ['A','T']):
            print(line[1])
            