#!/usr/bin/env python3

import random
import argparse
import os.path
import sys
import pickle
import array
import Trie
from scipy.stats import binom
from itertools import compress
from collections import Counter
from copy import copy
import timeit
import cProfile


def readargs():
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--naive',action='store_true')
    group.add_argument('--triangle',action='store_true')
    group.add_argument('--triesearch',action='store_true')
    group.add_argument('--heuristic',action='store_true')
    parser.add_argument('--ntypes', type=int, default=100)
    parser.add_argument('--nseqs', type=int,default=200000)
    parser.add_argument('--schemalength', type=int, default=3500)
    parser.add_argument('--seed',type=int, default=12)
    parser.add_argument('--diversity', type=int, default=1000)
    parser.add_argument('--seqs', type=str, default = ".seqs")
    parser.add_argument('--index',type=str,default = ".seqs.idx")
    parser.add_argument('--distance',type=int, default = 10)
    parser.add_argument('--create_seqs', action='store_true', help="Force creation of seqs and indexes")
    parser.add_argument('--createindex', action='store_true', help="Force creation of indices")
    return parser.parse_args()

def create_seqs(ntypes,nseqs,l,diversity,seed):
    """Create random seqs"""
    seqs = list()
    for i in range(ntypes):
        seqs.append(array.array('I',[random.randint(0,1000) for j in range(l)]))
    for i in range(nseqs-ntypes):
        s = copy(random.sample(seqs,1)[0])
        assert len(s)==l
        d = random.randint(0,200)
        positions = random.sample(range(l),d)
        for pos in positions:
            s[pos] = random.randint(0,l)
        seqs.append(s)
    return seqs


def search_seqs(seqs,query,maxdist):
    ## Naive search
    hits=list()
    for k,s in seqs:
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
    
def compare_heuristic(s,query,maxdist,softrange):
    d = 0
    i=0
    p=0
    checkpoint,low,high = softrange[p]
    p+=1
    for a,b in zip(s,query):
        i+=1
        if a!=b:
            d += 1
            if d>maxdist:
                return
        if i%checkpoint == 0:
            if d < low:
                return s
            elif d>high:
                return
            if len(softrange)>p:
                checkpoint,low,high = softrange[p]
                p+=1
    return s

def search_seqs_heuristic(seqs,query,maxdist):
    ## Heuristic search
    softrange = []
    for i in [10,100,1000]:
        b = binom(i,maxdist/len(query))
        softrange.append((i,b.ppf([0.001]),b.ppf([0.999])))
    hits=list()
    for k,s in seqs:
        result = compare_heuristic(s,query,maxdist,softrange)
        if result is not None:
            hits.append(result)
    return hits

def get_index_triangle(seqs,indexpath):
    if os.path.exists(indexpath):
        base, index = pickle.load(open(indexpath,'rb'))
    else:
        base = seqs[0][1]
        index = array.array('I')
        print(seqs[0][1][:10])
        for seq in seqs:
            index.append(dist(base,seq[1]))
        pickle.dump((base,index),open(indexpath,'wb'))
    assert(len(seqs)==len(index))
    return base,index

def index_seqs(seqs,indexpath,capacity):
    if os.path.exists(indexpath) and not args.createindex:
        index = pickle.load(open(indexpath,'rb'))
    else:
        index = dict()
        for seq in seqs:
            for i in range(args.schemalength):
                try:
                    index[i][seq[1][i]].add(seq[0])
                except KeyError:
                    index[i][seq[1][i]] = {seq[0]}
        pickle.dump(index,open(indexpath,'wb'))
    return index

def key(x):
    return x[1]

def index_trie(seqs,indexpath):
    if os.path.exists(indexpath) and not args.createindex:
        index = pickle.load(open(indexpath,'rb'))
    else:
        index = Trie.Trie(seqs,100,key=key)
        pickle.dump(index,open(indexpath,'wb'))
    assert(len(seqs)==len(index))
    return index

def search_trie(index,query,maxdist):
    seqs = index.search(query)
    #print(len(seqs),key(seqs[0])[:10],key(query)[:10])
    return search_seqs(seqs,key(query),maxdist)

def dist(seq,query):
    assert(len(seq) == len(query))
    d = 0
    for a,b in zip(seq,query):
        if a!=b:
            d+=1
    return d

def search_seqs_triangle(seqs,query,maxdist,indexpath):
    ## Triangle inequality search -- not really helpful

    ## Create distance index
    base,index = get_index_triangle(seqs,indexpath)
    dist_to_base = dist(base,query)
    subset = compress(seqs,[i>dist_to_base-maxdist and i <= dist_to_base+maxdist for i in index])
    ## Search using triangle inequality and index
    hits=list()
    for k,s in subset:
        result = compare(s,query,maxdist)
        if result is not None:
            hits.append(result)        
    return hits

if __name__=="__main__":
    args = readargs()
    random.seed(args.seed)
    if os.path.exists(args.seqs) and args.create_seqs == False:
        seqs = pickle.load(open(args.seqs,'rb'))
    else:
        seqs = create_seqs(args.ntypes,args.nseqs,args.schemalength,
                           args.diversity,args.seed)
        z = zip(range(args.nseqs),seqs)
        seqs = sorted(z,key=lambda x:x[1])
        #seqs = z
        pickle.dump(seqs,open(args.seqs,'wb'))
    s = random.sample(seqs,1)[0]
    hits=[]
    if args.naive:
        cProfile.run("search_seqs(seqs,s[1],args.distance)")
        hits = search_seqs(seqs,s[1],args.distance)
    elif args.heuristic:
        hits = search_seqs_heuristic(seqs,s[1],args.distance)
        cProfile.run("search_seqs_heuristic(seqs,s[1],args.distance)")
    elif args.triangle:
        hits = search_seqs_triangle(seqs,s[1],args.distance,args.index)
    elif args.triesearch:
        index = index_trie(seqs,args.index)
        #print(timeit.timeit("search_trie(index,s,args.distance)", setup="from __main__ import search_trie",number=100))
        cProfile.run("search_trie(index,s,args.distance)")
        hits = search_trie(index,s,args.distance)
    print(len(hits))


