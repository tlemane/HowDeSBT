# Tutorial, Creating a Bloom Filter Tree

(1) Estimate the best bloom filter size.

_This section is not yet written. For now, we magically assume the size is
500,000._

(2) Convert the fasta files to bloom filter bit vectors.

_Difference vs bloomtree-allsome: here you have to set the minimum abundance.
The default was 3 in bloomtree-allsome, in sabutan the default is 1._

```bash  
bf_size=500000
min_abundance=3
ls experiment*.fa \
  | while read f ; do
      bf=`echo ${f} | sed "s/\.fa/.bf/"`
      echo "=== converting ${f} to ${bf} ==="
      sabutan makebf --min=${min_abundance} K=20 --bits=${bf_size} ${f}
      done
```

The result of this step is a bloom filter for each of the fasta files, named
expriment1.bf, expriment2.bf, etc.

(2-ALT) Convert kmer files to bloom filter bit vectors.

Kmer files can be used as input instead of fasta or fastq files. In this case,
each line of the sequence input files is a single kmer, as the first field in
the line. Any additional fields on the line are ignored. For example, with
--k=20 this might be
```bash  
ATGACCAGATATGTACTTGC
TCTGCGAACCCAGACTTGGT
CAAGACCTATGAGTAGAACG
```

_Note that any minimum abundance is expected to have been applied when the kmer
list was created._

```bash  
bf_size=500000
ls experiment*.kmers \
  | while read f ; do
      bf=`echo ${f} | sed "s/\.kmers/.bf/"`
      echo "=== converting ${f} to ${bf} ==="
      sabutan makebf --kmersin K=20 --bits=${bf_size} ${f}
      done
```

(3) Create a tree topology.

_Difference vs bloomtree-allsome: bfcluster/sbuild did both the clustering
*and* built all the uncompressed nodes. Here we just cluster and leave the
building of nodes for later steps._

_Also, ../bfcluster/sbuild hardwired the number of bits to 500K (used for
determining Hamming distance between filters). In this example I'm using 50K,
but to reproduce what you got with bloomtree-allsome you'd want 500K._


```bash  
cluster_bits=50000
ls experiment*.bf > leafnames
sabutan cluster --list=leafnames --bits=${cluster_bits} \
  --tree=example.sbt --node=node{node}.bf
```

The result of this step is a tree topology file, example.sbt. Note that no
bloom filters are actually created in this step.

(4) Build the "determined,brief" tree, compressed as RRR.

_Note that in earlier versions of sabutan this was performed in two steps, but
now RRR can be created directly without having to build the uncompressed tree._

```bash  
sabutan build --determined,brief --rrr --tree=example.sbt \
  --outtree=detbrief.rrr.sbt
```

The result of this step is bloom filter files for the leaves and internal nodes,
in "determined,brief" format and compressed with RRR. And a new topology file
named detbrief.rrr.sbt. The bloom filter files are named experiment1.detbrief.rrr.bf,
..., node1.detbrief.rrr.bf, etc.

(5) Run a batch of queries.

_Note that our example has queries as fasta, where the bloomtree-allsome example
just had queries as separate lines of the input file. Sabutan can accept that
format (if there's no fasta header, it assumes queries are one per line). But
if fasta input is used, the output identifies the queries by their fasta
headers._

```bash  
sabutan query --tree=detbrief.rrr.sbt \
    queries1.fa=0.5 queries2.fa=0.7 \
  > queryresults
```

The resulting queryresults should be identical to queryresults.expected.

(6) Ordering query results by how good they are.

_Note that this describes a short term solution that we plan to improve upon.
It requires that the query files are real fasta files with sequence headers._

```bash  
../scripts/order_query_results.sh \
  queries1.fa queries2.fa \
  queryresults \
  queryresults.kmers
```

The resulting queryresults.kmers should be identical to
queryresults.kmers.expected, which is also shown below. The first column shows
the name of a query sequence, as appeared within a fasta file. The second
column gives the name of an experiment which the query 'hit'. The third column
shows the fraction of the query's kmers that are present in the experiment
(more correctly, in the the experiment's bloom filter). The fourth column gives
that fraction as a a number. Note that queries are sorted by name, which may
not match their order in the sabutan output file.

```bash  
EXAMPLE_QUERY1 experiment8  481/481 1.000000
EXAMPLE_QUERY1 experiment10 355/481 0.738046
EXAMPLE_QUERY3 experiment2  481/481 1.000000
EXAMPLE_QUERY4 experiment1  387/481 0.804574
EXAMPLE_QUERY5 experiment3  481/481 1.000000
EXAMPLE_QUERY5 experiment6  481/481 1.000000
```
