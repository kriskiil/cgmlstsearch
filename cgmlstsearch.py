#!/usr/bin/env python3

import argparse
import os.path
import numpy as np
import Trie
import pickle
from copy import copy
#import cProfile

## dtype
dtype = np.uint8

def readargs():
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--naive',action='store_true')
    group.add_argument('--trie',action='store_true')
    parser.add_argument('--ntypes', type=int, default=100)
    parser.add_argument('--nseqs', type=int,default=200000)
    parser.add_argument('--schemalength', type=int, default=3500)
    parser.add_argument('--seed',type=int, default=42)
    parser.add_argument('--diversity', type=int, default=100)
    parser.add_argument('--seqs', type=str, default = ".seqs.npy")
    parser.add_argument('--index',type=str,default = ".seqs.idx")
    parser.add_argument('--distance',type=int, default = 10)
    parser.add_argument('--create-seqs', action='store_true', help="Force creation of seqs and indexes")
    parser.add_argument('--create-index', action='store_true', help="Force creation of indices")
    return parser.parse_args()

def create_seqs(nseqs,l,diversity,seed):
    """Create random seqs"""
    allele_maxima = np.ones(l)
    originalseq = np.ones(l)
    #originalseq = np.random.randint(low=1,high=diversity,size=l)
    seqs = np.ndarray((nseqs,l),dtype=dtype)
    seqs[0] = originalseq
    #distances = np.random.poisson(diversity,size=nseqs-1)
    distances = np.random.binomial(l,diversity/l,size=nseqs-1)
    n = 0
    for i in distances:
        n += 1
        sites = np.random.choice(l,size=i, replace=False)
        originalseq = seqs[np.random.choice(n)]
        newseq = copy(originalseq)
        allele_maxima[sites] += 1
        newseq[sites] = allele_maxima[sites]
        #newseq[sites] = np.random.randint(low=1,high=diversity,size=i)
        seqs[n] = newseq
    return seqs

def search_seqs(seqs,query,maxdist):
    ## Naive search
    hits=list()
    for s in seqs:
        result = compare(s,query,maxdist)
        if result is not None:
            hits.append(result)
    return hits

def compare(s,query,maxdist):
    assert(len(s)==len(query))
    d = 0
    for a,b in zip(s,query):
        if a!=b:
            d += 1
            if d>maxdist:
                return
    return s

def index_trie(seqs,indexpath):
    if os.path.exists(indexpath) and not args.create_index and not args.create_seqs:
        index = pickle.load(open(indexpath,'rb'))
    else:
        index = Trie.Tries(10,range(len(seqs)),200,seqs)
        pickle.dump(index,open(indexpath,'wb'))
    assert(len(seqs)==len(index))
    return index

def search_trie(index,seqs,query,maxdist):
    idx = index.search(query)
    print(len(idx))
    hits=list()
    for i in idx:
        result = compare(seqs[i],query,maxdist)
        if result is not None:
            hits.append(result)
    return hits

def dist(seq,query):
    assert(len(seq) == len(query))
    d = 0
    for a,b in zip(seq,query):
        if a!=b:
            d+=1
    return d


if __name__=="__main__":
    args = readargs()
    if os.path.exists(args.seqs) and args.create_seqs == False:
        seqs = np.memmap(args.seqs,mode='r+',dtype=dtype,shape=(args.nseqs,args.schemalength))
    else:
        seqs = create_seqs(args.nseqs,args.schemalength,
                           args.diversity,args.seed)
        mm = np.memmap(args.seqs,dtype=dtype,mode='w+',shape=(args.nseqs,args.schemalength))
        mm[:] = seqs[:]
    #s = seqs[np.random.choice(args.nseqs)]
    s = seqs[10]
    hits=[]
    if args.naive:
        #cProfile.run("search_seqs(seqs,s,args.distance)")
        hits = search_seqs(seqs,s,args.distance)
    elif args.trie:
        index = index_trie(seqs,args.index)
        #cProfile.run("search_trie(index,seqs,s,args.distance)")
        hits = search_trie(index,seqs,s,args.distance)
    print(len(hits))

