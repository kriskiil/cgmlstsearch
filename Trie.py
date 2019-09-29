#!/usr/bin/env python3

import numpy as np

class Tries:
    def __init__(self,n,iterable,Capacity,seqs):
        self.tries = list()
        for i in range(n):
            self.tries.append(Trie(iterable,Capacity,seqs))
    def search(self,item):
        hits = set()
        for trie in self.tries:
            hits.update(trie.search(item))
        return hits
    def __len__(self):
        return len(self.tries[0])

class Trie:
    seqs = None
    def __init__(self,iterable,Capacity,seqs=None):
        if seqs is not None:
            Trie.seqs = seqs
        self.d = None
        self.l = list()
        self.C = Capacity
        self.i = np.random.randint(low=0,high=len(Trie.seqs[0])-1)
        self.length = 0
        for item in iterable:
            self.add(item)

    def add(self,n):
        val = self._val(n)
        self.length += 1
        if self.d is None:
            self.l.append(n)
            if len(self.l) > self.C:
                self.d = dict()
                for seq in self.l:
                    val = self._val(seq)
                    self.d.setdefault(val[self.i],Trie([],self.C)).add(seq)
                del self.l
        else:
            self.d.setdefault(val[self.i],Trie([],self.C)).add(n)

    def search(self,item):
        """Returns up to C matches"""
        if self.d is None:
            return self.l
        else:
            return self.d[item[self.i]].search(item)

    def _val(self,item):
        return Trie.seqs[item]
    
    def __len__(self):
        return self.length
