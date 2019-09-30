# cgmlstsearch

cgMLSTsearch is a proof of concept implementation for quickly searching a cgMLST profile against a database. 
So far there is no handling of missing data.

The implementation drew inspiration from this paper:
https://arxiv.org/ftp/arxiv/papers/1505/1505.03090.pdf

# Usage

--create-sequences Creates a new set of sequences

--create-index Recreates the index

--distance [int] the distance from the query within which matches are returned.

--naive Naive search that compares the query to each database sequence, and drops off when the search distance is exceeded.

--trie Search that uses a forest of random order tries to limit the search space.

--heuristic This option changes the comparison to one that uses a binomial confidence interval to end the comparison early, 
for very similar or disimilar profiles. Mainly useful if comparisons are slow due to high distance (=slow dropoff and more comparisons)

# Performance

On a laptop cpu, actually running the program takes less than 2 seconds to search 200.000 sequences. The actual search (with indexing) takes 89 milliseconds according to cProfile.

The search is parallelizable but this is hardly necessary with current database sizes.

# Caveats

In the current implementation of indexed search, hits will be missed when using a large search radius.
