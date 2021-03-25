#!/usr/bin/env python
import sys

for line in sys.stdin:
    line = line.strip()
    if line[0] == '#':
        print(line)
    row = line.split('\t')
    count = sum(x.split(':')[0] != './.' for x in row[9:])
    print(line + '\t' + str(count))
