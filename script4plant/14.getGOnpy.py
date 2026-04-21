#!/usr/bin/python3

import numpy as np
with open('go_term.txt') as f:
   content = [ rec.strip() for rec in f.readlines() ]

go_dict = dict()
for rec in content:
    t,term,gene = rec.split('\t')
    if t not in go_dict:
        go_dict[t] = dict()
    if term not in go_dict[t]:
        go_dict[t][term] = list()
    go_dict[t][term].append(gene)

np.save('go_dict.npy', go_dict)
