import os
import pysam
import shutil
from bx.intervals.cluster import ClusterTree
from collections import defaultdict, Counter
import numpy as np
import networkx as nx
from itertools import groupby, combinations
import community

from .step import StepChunk
from ..mlib import util

class BinMetaReadsStep(StepChunk):

  @staticmethod
  def get_steps(options):
    yield BinMetaReadsStep(options)

  def __init__(
    self,
    options,
  ):
    self.options = options
    util.mkdir_p(self.outdir)

  def __str__(self):
    return self.__class__.__name__

  @property
  def outdir(self):
    return os.path.dirname(self.options.groups_pickle_path)

  def outpaths(self, final=False):
    paths = {}
    paths['bins.p'] = self.options.bins_pickle_path
    paths['bins2.p'] = self.options.bins2_pickle_path
    paths['shit.p'] = 'shit'
    return paths

  def run(self):
    self.logger.log('determine bins from seed contigs')

    ctg_size_map = util.get_fasta_sizes(self.options.ctgfasta_path)
    
    hmpfasta_path = '/scratch/PI/serafim/abishara/reference/hmp/all_seqs.fa'
    ctghmp_bam_path = '/scratch/users/abishara/scratch/scratch.metagenome/align-contig-hmp.sorted.bam'

    # find hmp seeds
    hmp_seedctg_set = set([
      'BACT_852|gi|145218786|ref|NZ_AAXE02000112.1|',
    ])
    self.logger.log('find hmp seed contigs')
    hmp_seedctg_set = get_hmp_seeds(
      hmpfasta_path,
      ctghmp_bam_path,
    )
    hmpctg_size_map = util.get_fasta_sizes(hmpfasta_path)
    hmp_bases = sum(map(lambda(c): hmpctg_size_map[c], hmp_seedctg_set))
    self.logger.log('  {} seeds covering {} bases'.format(
      len(hmp_seedctg_set),
      hmp_bases,
    ))

    self.logger.log('find idba0 contigs mapping to hmp seeds')
    hmp_ctgset_map = get_hmpmapped_ctgs(
      hmp_seedctg_set,
      ctghmp_bam_path,
    )
    filtqname_set = set()
    for hmp_ctg in hmp_ctgset_map:
      filtqname_set |= hmp_ctgset_map[hmp_ctg]

    util.write_pickle(
      os.path.join(self.options.working_dir, 'hmp-ctgset.p'),
      filtqname_set,
    )
    die

    self.logger.log('  {} idba0 contigs mapping to seeds'.format(len(filtqname_set)))

    tmp_hits_path = '/home/abishara/tmp/tmp_bcode-ctg-hits.txt'
    hits_path      = self.options.bcode_ctg_hits_path + '_bcode-ctg-hits.txt'
    bcode_idx_path = self.options.bcode_ctg_hits_path + '_bcode-idx.p'
    ctg_idx_path   = self.options.bcode_ctg_hits_path + '_ctg-idx.p'

    bcode_idx_map = util.load_pickle(bcode_idx_path)
    idx_bcode_map = {v: k for k, v in bcode_idx_map.items()}
    ctg_idx_map = util.load_pickle(ctg_idx_path)
    idx_ctg_map = {v: k for k, v in ctg_idx_map.items()}

    # some contigs may have no barcode hits...
    filtqname_set = set(filter(
      lambda(qname): qname in ctg_idx_map,
      filtqname_set,
    ))
    self.logger.log('  {} idba0 contigs mapping to seeds with barcode hits'.format(len(filtqname_set)))
    filt_ctgidx_set = set(map(
      lambda(qname): ctg_idx_map[qname],
      filtqname_set,
    ))

    self.logger.log('parsing barcode-contig hits to create graph')
    G = nx.Graph()
    ctg_bcodes_map = defaultdict(Counter)
    total_reads = 0
    #fout = open(tmp_hits_path, 'w')
    with open(tmp_hits_path) as f:
    #with open(hits_path) as f:
      for k, lines in groupby(f, lambda(l): l.split('\t')[0]):
        hits_set = set()
        for line in lines:
          bcode_idx, ctg_idx, cnt = map(int, line.split('\t'))
          # skip irrelevant contig hits
          if ctg_idx not in filt_ctgidx_set:
            continue
          total_reads += cnt
          G.add_node(ctg_idx)
          hits_set.add(ctg_idx)
          ctg_bcodes_map[ctg_idx][bcode_idx] += cnt
          #fout.write(line)
        for u, v in combinations(hits_set, 2):
          if G.has_edge(u,v):
            G[u][v]['h'] += 1
          else:
            G.add_edge(u,v)
            G[u][v]['h'] = 1
    #fout.close()
    
    print 'total reads', total_reads
    self.logger.log('  - done, {} nodes, {} edges'.format(
      G.number_of_nodes(),
      G.number_of_edges(),
    ))

    self.logger.log('filter for edges with >= 10 barcode hits')
    to_remove = []
    for u, v in G.edges_iter():
      if G[u][v]['h'] < 10:
        to_remove.append((u,v))
    G.remove_edges_from(to_remove)

    to_remove = []
    for v, d in G.degree_iter():
      if d == 0:
        to_remove.append(v)
    G.remove_nodes_from(to_remove)

    self.logger.log('  - done, {} nodes, {} edges'.format(
      G.number_of_nodes(),
      G.number_of_edges(),
    ))

    self.logger.log('find cliques')
    cliques = sorted(filter(
      lambda(c): len(c) >= 4,
      nx.find_cliques(G),
    ), reverse=True, key=len)

    seen_set = set()
    n_cliques = []
    for c in cliques:
      #if not set(c).issubset(seen_set):
      #  n_cliques.append(c)
      if len(set(c)) >= len(set(c) & seen_set) + 2:
        n_cliques.append(c)
      #print 'size {}, seen {}'.format(
      #  len(c),
      #  len(set(c) & seen_set),
      #), sorted(c)

      seen_set |= set(c)
    print 'new cliques', len(n_cliques)
    cliques = n_cliques


    clique_ctgidx_set = set([n for c in cliques for n in c])
    self.logger.log('  {} cliques of size >= 4'.format(len(cliques)))
    self.logger.log('  {} contigs out of {} included in cliques'.format(
      len(clique_ctgidx_set),
      len(filt_ctgidx_set),
    ))

    def get_bcode_keys(c):
      bcode_read_counts = Counter()
      bcode_ctg_counts = Counter()
      for ctg in c:
        bcode_read_counts.update(ctg_bcodes_map[ctg])
        bcode_ctg_counts.update(ctg_bcodes_map[ctg].keys())
     
      bcode_keys = []
      for bcode, ctg_counts in bcode_ctg_counts.most_common():
        if ctg_counts <= 1:
          continue
        read_counts = bcode_read_counts[bcode]
        bcode_keys.append((bcode, ctg_counts, read_counts))

      bcode_keys.sort(key=lambda(x): x[1:], reverse=True)

      return bcode_keys
        
    # create a bin for each clique
    bins = []
    bins2 = []
    skipped = 0
    all_bcode_set = set()
    for i, c in enumerate(cliques):
      bcode_keys = get_bcode_keys(c)
      total_reads = sum(map(lambda(x): x[2], bcode_keys))
      # compute coverage and downsample barcodes if necessary
      size = sum(map(lambda(ctg): ctg_size_map[idx_ctg_map[ctg]], c))
      min_ctg_size = min(map(lambda(ctg): ctg_size_map[idx_ctg_map[ctg]], c))
      if min_ctg_size > 4000:
        continue
      cov = 95. * total_reads / size

      bcode_set = set(map(lambda(x): idx_bcode_map[x[0]], bcode_keys))
      ctg_set = set(map(lambda(x): idx_ctg_map[x], c))
      #print 'binid      ', 'bin.{}'.format(i)
      #print '  - size       ', size
      #print '  - total_reads', total_reads
      #print '  - ctg_set', ctg_set
      #print '  - num barcodes', len(bcode_set)
      #print '  - cov', cov
      #die
      if cov <= 10000 and cov > 10:
        binid = 'bin.{}'.format(i)
        downsample_rate = 200. / cov
        downsample_rate = 800. / len(bcode_set)
        bcode_set = set(map(lambda(x): idx_bcode_map[x[0]], bcode_keys))
        all_bcode_set |= bcode_set
        ctg_set = set(map(lambda(x): idx_ctg_map[x], c))
        ds_bcode_set = None
        #print '  - num barcodes', len(bcode_set), cov
        if downsample_rate < 1:
          e = int(downsample_rate * len(bcode_set))
          ds_bcode_set = set(map(lambda(x): idx_bcode_map[x[0]], bcode_keys[:e]))
          ds_total_reads = sum(map(lambda(x): x[2], bcode_keys[:e]))
          ds_cov = 95. * ds_total_reads / size
          #print '    - ds cov', ds_cov
          
        bins.append((binid, bcode_set))
        bins2.append((binid, ds_bcode_set, ctg_set))
    bins.append(('bin.union', all_bcode_set))
      
    self.logger.log('created {} bins from seeds'.format(len(bins)))
    self.logger.log('  - {} skipped for not containing small contig'.format(skipped))
    util.write_pickle(self.options.bins_pickle_path, bins)
    util.write_pickle(self.options.bins2_pickle_path, bins2)

#--------------------------------------------------------------------------
# helpers
#--------------------------------------------------------------------------
def get_hmpmapped_ctgs(
  hmp_seed_set,
  ctghmp_bam_path,
):

  fin = pysam.Samfile(ctghmp_bam_path, 'rb')
  hmp_ctgset_map = defaultdict(set)
  for i, read in enumerate(fin):
    if (
      read.is_unmapped or
      read.query_alignment_length < 0.8 * read.query_length
    ):
      continue

    hmp_ctg = fin.get_reference_name(read.tid)
    if hmp_ctg in hmp_seed_set:
      hmp_ctgset_map[hmp_ctg].add(read.qname)
  fin.close()
  return hmp_ctgset_map


def get_hmp_seeds(
  hmpfasta_path,
  ctghmp_bam_path,
):

  fin = pysam.Samfile(ctghmp_bam_path, 'rb')
  regions_map = defaultdict(lambda: ClusterTree(1,1))
  for i, read in enumerate(fin):
    if read.is_unmapped:
      continue
    hmp_ctg = fin.get_reference_name(read.tid)
    regions_map[hmp_ctg].insert(read.pos, read.aend, i)
  fin.close()

  hmpctg_size_map = util.get_fasta_sizes(hmpfasta_path)
  fasta = pysam.FastaFile(hmpfasta_path)
  cand_hmpctg_set = set()
  for hmp_ctg in regions_map.keys():
    bpcov = 0
    for b, e, _ in regions_map[hmp_ctg].getregions():
      bpcov += (e - b)

    size = hmpctg_size_map[hmp_ctg]
    cov = 1. * bpcov / size
    if size > 50000 and cov > 0.9:
      cand_hmpctg_set.add(hmp_ctg)

  return cand_hmpctg_set

