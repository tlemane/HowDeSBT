# Tutorial, Creating a Bloom Filter Tree

(1) Estimate the best bloom filter size.
_This section is not yet written.  For now, we magically assume the size is
500,000._

(2) Convert the fasta files to bloom filter bit vectors.

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

Kmer files can be used as input instead of fasta or fastq files.  In this case,
each line of the sequence input files is a single kmer, as the first field in
the line. Any additional fields on the line are ignored.  For example, with
--k=20 this might be
```bash  
ATGACCAGATATGTACTTGC
TCTGCGAACCCAGACTTGGT
CAAGACCTATGAGTAGAACG
```

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

```bash  
cluster_bits=50000
ls experiment*.bf > leafnames
sabutan cluster --list=leafnames --bits=${cluster_bits} \
  --tree=example.sbt --node=node{node}.bf
```

The result of this step is a tree topology file, example.sbt.  Note that no
bloom filters are actually created in this step.

(4) Build the "determined,brief" tree.

```bash  
sabutan build --determined,brief --tree=example.sbt --outtree=detbrief.sbt
```

The result of this step is bloom filter files for the leaves and internal nodes,
in "determined,brief" format, and a new topology file named detbrief.sbt.  The
bloom filter files are named experiment1.detbrief.bf, ..., node1.detbrief.bf,
etc.

(4b) Convert the "determined,brief" tree to RRR.

Note that eventually this step will be combined with step 4.

```bash  
cat detbrief.sbt \
  | sed "s/\.detbrief\.bf/.detbrief.rrr.bf/" \
  > detbrief.rrr.sbt

cat detbrief.rrr.sbt \
  | tr -d "*" \
  | sed "s/\.detbrief\.rrr\.bf//" \
  | while read node ; do
      nodeNum=$((nodeNum+1))
      echo "=== RRR-compressing ${node} ==="
      sabutan compressbf ${node}.detbrief.bf --rrr
      done
```

The result of this step is bloom filter files for the leaves and internal nodes,
in "determined,brief" format and compressed with RRR.  And a new topology file
named detbrief.rrr.sbt.  The bloom filter files are named experiment1.detbrief.rrr.bf,
..., node1.detbrief.rrr.bf,
etc.

(5) Run a batch of queries.
_Note that our example has queries as fasta, where the bloomtree-allsome example
just had queries as separate lines of the input file.  Sabutan can accept that
format (if there's no fasta header, it assumes queries are one per line).  But
if fasta input is used, the output identifies the queries by their fasta
headers._

```bash  
${sabutan} query --tree=detbrief.rrr.sbt \
    queries.fa --threshold=0.5 \
  > queryresults
```
