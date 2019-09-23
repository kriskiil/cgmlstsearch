#!/usr/bin/env python3

class Trie:
    def __init__(self,iterable,C,i=0,key=None):
        self.d = None
        self.l = list()
        self.C = C
        self.i = i
        self.length = 0
        self.key = key
        for item in iterable:
            self.add(item)

    def add(self,item):
        val = self._val(item)
        self.length += 1
        if self.d is None:
            self.l.append(item)
            if len(self.l) > self.C and self.i<len(val)-1:
                self.d = dict()
                for seq in self.l:
                    val = self._val(seq)
                    self.d.setdefault(val[self.i],Trie([],self.C,self.i+1,self.key)).add(seq)
                del self.l
        else:
            self.d.setdefault(val[self.i],Trie([],self.C,self.i+1,self.key)).add(item)

    def search(self,item):
        """Returns up to C matches"""
        #print(self.i,self.length)
        if self.d is None:
            return self.l
        else:
            return self.d[self._val(item)[self.i]].search(item)

    def _val(self,item):
        if self.key is not None:
            return self.key(item)
        else:
            return item

    def __len__(self):
        return self.length
