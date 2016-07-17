import os
import pysam
import subprocess
from collections import defaultdict, Counter
import glob
import shutil
from itertools import izip
import random
import numpy as np
from bx.intervals.cluster import ClusterTree

from ..assembler_tools.architect import architect
from ..assembler_tools.architect import local_assembly

from .step import StepChunk
from ..mlib import util
from ..mlib.fq_idx import FastqIndex

# NOTE must be in path
idbabin_path = 'idba_ud'

wd = os.path.dirname(os.path.abspath(__file__))
architect_scripts_path = os.path.join(
  wd,
  '..',
  'assembler_tools/architect/scripts/',
)
edgesbin_path = os.path.join(architect_scripts_path, 'pe-connections.py')
containmentbin_path = os.path.join(architect_scripts_path, 'bam_to_containment.py')

MIN_SEED_SIZE = 500
MIN_SEED_SIZE = 400

class AssembleMetaBinnedStep(StepChunk):

  @staticmethod
  def get_steps(options):
    bins = util.load_pickle(options.bins_pickle_path)
    bins2 = util.load_pickle(options.bins2_pickle_path)
    assert len(bins[:-1]) == len(bins2)
    
    for i, ((binid, bcode_set), (_, bcode_counts)) in \
        enumerate(izip(bins[:-1], bins2)):
      #if binid != 'bin.contig-100_5761':
      #  continue
      yield AssembleMetaBinnedStep(options, binid, bcode_set, bcode_counts)

  def __init__(
    self,
    options,
    binid,
    bcode_set,
    bcode_counts,
  ):
    self.options = options
    self.binid = binid
    self.bcode_set = bcode_set
    self.bcode_counts = bcode_counts
    util.mkdir_p(self.outdir)

  def __str__(self):
    return '{}.{}'.format(self.__class__.__name__, self.binid)

  @property
  def outdir(self):
    return self.options.get_bin_dir(self.binid, final=True)

  def outpaths(self, final=False):
    paths = {}
    paths['local-asm.p'] = os.path.join(self.outdir, 'local-asm.p')
    paths['local-asm-merged.fa'] = os.path.join(self.outdir, 'local-asm-merged.fa')
    #paths['shit'] = 'shit'
    return paths

  def clean_working(self):
    bin_path = self.options.get_bin_dir(self.binid)
    self.logger.log('removing bin directory {}'.format(bin_path))
    shutil.rmtree(bin_path)
    return 

  def run(self):

    self.logger.log('assembling barcoded reads for this bin')
    root_ctg = self.binid[4:]

    ctg_size_map = util.get_fasta_sizes(self.options.ctgfasta_path)
    seed_ctgs = set(filter(
      lambda(c): ctg_size_map[c] >= MIN_SEED_SIZE,
      ctg_size_map.keys(),
    ))

    fqdir_path = self.options.get_bin_fq_dir(self.binid)
    asmrootdir_path = self.options.get_bin_asm_dir(self.binid)
    util.mkdir_p(asmrootdir_path)
    util.mkdir_p(fqdir_path)

    idx_bcode_map = {v: k for k, v in self.options.bcode_idx_map.items()}
    bcode_counts = dict(map(
      lambda(i, c): (idx_bcode_map[i], c),
      self.bcode_counts.items(),
    ))

    self.logger.log('determing local assemblies')
    local_asms = self.get_local_asms(
      seed_ctgs,
      root_ctg, 
      self.bcode_set,
      bcode_counts,
      self.options.reads_ctg_bam_path,
    )
    util.write_pickle(
      os.path.join(asmrootdir_path, 'local-asm.p'),
      local_asms,
    )
    self.logger.log('found {} local assemblies'.format(len(local_asms)))

    rootfq_path = self.options.longranger_fqs_path
    tenxfq_paths = list(glob.glob(rootfq_path + '/chnk*/files/*fastq'))
    for i, (n_ctg, bcode_set, ds_bcode_set) in enumerate(local_asms):
      ## FIXME remove
      #continue
      self.logger.log('assembling with neighbor {}'.format(n_ctg))
      self.logger.log('  - {} orig barcodes'.format(len(bcode_set)))
      self.logger.log('  - {} downsampled barcodes'.format(len(ds_bcode_set)))

      lrhintsfa_path = os.path.join(fqdir_path, 'lr-hints.{}.fa'.format(i))
      readsfa_path = os.path.join(fqdir_path, 'reads.{}.fa'.format(i))
      asmdir_path = os.path.join(asmrootdir_path, 'local-asm.{}'.format(i))

      # filter input reads
      get_bcode_reads(
        tenxfq_paths,
        readsfa_path,
        ds_bcode_set,
      )
      # create long read hints
      get_lrhints_fa(
        [root_ctg, n_ctg],
        self.options.ctgfasta_path,
        lrhintsfa_path,
      )

      # initial idba assembly
      self.logger.log('  - performing idba assembly')
      self.do_idba_assemble(
        readsfa_path,
        asmdir_path,
        lrhintsfa_path,
      )

    # merge output contigs from local assemblies
    self.logger.log('merge long output contigs from local assemblies')
    mergedasm_path = os.path.join(asmrootdir_path, 'local-asm-merged.fa')
    total_asm_contigs = 0
    total_asm_bp = 0
    with open(mergedasm_path, 'w') as fout:
      for i, (n_ctg, _, _) in enumerate(local_asms):
        asmdir_path = os.path.join(asmrootdir_path, 'local-asm.{}'.format(i))
        contig_path = os.path.join(asmdir_path, 'contig.fa')
        try:
          fasta = pysam.FastaFile(contig_path)
        except Exception as e:
          self.logger.log('could not open fasta {}'.format(contig_path))
          self.logger.log(str(e))
          self.logger.log('skipping')
          continue
        for contig in sorted(
          fasta.references,
          key=lambda(c): fasta.get_reference_length(c),
          reverse=True,
        ):
          seq = str(fasta.fetch(contig).upper())
          if len(seq) < 2000:
            break
          total_asm_contigs += 1
          total_asm_bp += len(seq)
          fout.write('>{}.{}${}.{}\n'.format(root_ctg, n_ctg, contig, i))
          fout.write(str(seq) + '\n')

    self.logger.log('  - {} contigs covering {} bases'.format(
      total_asm_contigs,
      total_asm_bp))
      
    def copyfinal(src, dest):
      shutil.copyfile(
        os.path.join(self.options.get_bin_asm_dir(self.binid), src),
        os.path.join(self.options.get_bin_dir(self.binid, final=True), dest),
      )

    self.logger.log('copying deliverables to final')
    copyfinal('local-asm-merged.fa', 'local-asm-merged.fa')
    copyfinal('local-asm.p', 'local-asm.p')

    self.logger.log('done')

#--------------------------------------------------------------------------
# helpers
#--------------------------------------------------------------------------
  def do_idba_assemble(
    self,
    readsfa_path,
    asmdir_path,
    lrhintsfa_path,
  ):
    #cmd = '{} -r {} -l {} -o {}'.format(
    cmd = '{} --maxk 60 -r {} -l {} -o {}'.format(
      idbabin_path,
      readsfa_path,
      lrhintsfa_path,
      asmdir_path,
    )
    cmd = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE)
    retcode = cmd.wait()
    print "::::", retcode
    contigs_path = os.path.join(asmdir_path, 'contig.fa')
    if not os.path.exists(contigs_path):
      if retcode != 0 and  "invalid insert distance" in cmd.stderr.read():
        self.logger.log("idba_ud internal error; this is probably a bug in idba_ud, so we have to bail here")
        open(contigs_path, "w")
      else:
        error_message = "assembly failed to produce contig.fa"
        self.logger.log(error_message)
        raise Exception()
    elif retcode != 0:
      self.logger.log(
        "something funky happened while running idba_ud (got error code "
        "{}) but there's not much we can do so continuing".format(retcode))


  def get_local_asms(
    self,
    seed_ctgs,
    root_ctg, 
    root_bcode_set,
    bcode_counts,
    bam_path,
  ):
    SMALL_CTG_SIZE = 1000
    COV_BIN_SIZE = 100
  
    fhandle = pysam.Samfile(bam_path, 'rb')
    ctg_size_map = util.get_fasta_sizes(self.options.ctgfasta_path)

    def get_barcode(read):
      filt_list = filter(lambda(k, v): k == 'BX', read.tags)
      if filt_list == []: 
        return None
      else:
        k, v = filt_list[0]
        return v
  
    def get_ctg_links(ctg):
      bcodes = set()
      bcode_counts = Counter()
      link_ctg_counts = Counter()
      link_ctg_qnames = defaultdict(set)
      link_regions_map = defaultdict(lambda: ClusterTree(200, 1))
      cov_bin_counts = Counter()
      for i, read in enumerate(fhandle.fetch(ctg)):
        if read.is_unmapped:
          continue
        # skip low quality links for now
        if read.mapq < 10:
          continue
        bcode = get_barcode(read)
        if bcode == None:
          continue
        bcodes.add(bcode)
        bcode_counts[bcode] += 1
        p_ctg = fhandle.getrname(read.next_reference_id)    
        link_ctg_counts[p_ctg] += 1
        link_ctg_qnames[p_ctg].add(read.qname)
        link_regions_map[p_ctg].insert(read.pos, read.aend, i)
  
        # update coverage profile
        bin_idx = read.pos / COV_BIN_SIZE
        cov_bin_counts[bin_idx] += 1
        
      # filter for only links in which the paired end reads are clustered
      # near each other
      #print 'initial cand links', len(link_ctg_counts)
      for p_ctg, regions in link_regions_map.items():
        regions = link_regions_map[p_ctg]
        max_hits = max(map(lambda(x): len(x[2]), regions.getregions()))
        if max_hits < 3:
          del link_ctg_counts[p_ctg]
          del link_ctg_qnames[p_ctg]
      #print 'post cand links', len(link_ctg_counts)

      return (
        bcodes,
        bcode_counts,
        link_ctg_counts,
        link_ctg_qnames,
        cov_bin_counts,
      )

    def is_whack_coverage(cov_bin_counts):
      max_cov = max(cov_bin_counts.values())
      med_cov = np.median(cov_bin_counts.values())
      return max_cov > 20 * med_cov

    ( _, _, 
      link_ctg_counts,
      link_ctg_qnames,
      cov_bin_counts,
    ) = get_ctg_links(root_ctg)
  
    # don't locally assemble any small contigs with out of whack coverage
    if ctg_size_map[root_ctg] < SMALL_CTG_SIZE and is_whack_coverage(cov_bin_counts):
      max_cov = max(cov_bin_counts.values())
      med_cov = np.median(cov_bin_counts.values())
      self.logger.log('max coverage {} order of magnitude higher than median {}'.format(
        max_cov, med_cov))
      return []
  
    link_cand_ctgs = filter(
      lambda(c): (c != root_ctg and  link_ctg_counts[c] >= 3),
      link_ctg_counts,
    )

    self.logger.log('{} initial link candidates to examine'.format(
      len(link_cand_ctgs)))
    # examine recipricol links 
    link_ctgs = []
    # <ctg> : set(bcodes)
    link_ctg_bcodes = {}
    # <ctg> : Counter(bcodes)
    link_ctg_bcode_counts = {}
    for i, n_ctg in enumerate(link_cand_ctgs):
      print 'n_ctg', n_ctg
      (
        link_ctg_bcodes[n_ctg],
        link_ctg_bcode_counts[n_ctg],
        _link_ctg_counts,
        _link_ctg_qnames,
        _cov_bin_counts,
      ) = get_ctg_links(n_ctg)
      # link is not recripricol
      if _link_ctg_counts[root_ctg] < 3:
        continue
      # link contig has whacky coverage profile
      # FIXME consider removing this, may not be necessary
      if (
        ctg_size_map[n_ctg] < SMALL_CTG_SIZE and
        is_whack_coverage(_cov_bin_counts)
      ):
        continue
      # not enough barcode overlap
      if len(link_ctg_bcodes[n_ctg] & root_bcode_set) < 15:
        continue
      
      print '  - added'
      link_ctgs.append(n_ctg)

    self.logger.log('  - {} pass reciprocal filtering'.format(
      len(link_ctgs)))

    # only return assemblies that aren't going to be duplicated in other
    # bins
    n_ctgs = filter(
      lambda(c): (
        c != root_ctg and 
        # only perform local asm in this bin if neighbor contig is
        # lexicographically smaller
        (c < root_ctg or c not in seed_ctgs)
      ),
      link_ctgs,
    )
    ## FIXME remove
    #n_ctgs = link_ctgs

    local_asms = []
    for n_ctg in n_ctgs:
      assert link_ctg_counts[n_ctg] >= 3
      link_bcode_counts = link_ctg_bcode_counts[n_ctg]
      i_bcode_set = link_ctg_bcodes[n_ctg] & root_bcode_set
      ds_bcode_set = i_bcode_set
      # downsample to barcodes with the most reads
      if len(i_bcode_set) > 400:
        ds_bcode_set = set(sorted(
          i_bcode_set,
          reverse=True,
          key=lambda(b): bcode_counts[b] + link_bcode_counts[b],
        )[:400])
      local_asms.append(
        (n_ctg, i_bcode_set, ds_bcode_set)
      )

    # do local reassembly of the root contig
    ds_bcode_set = root_bcode_set
    if len(root_bcode_set) > 400:
      ds_bcode_set = set(random.sample(root_bcode_set, 400))
  
    local_asms.append(
      (root_ctg, root_bcode_set, ds_bcode_set)
    )
  
    return local_asms
  
    # DEBUG
    #--------------------------------------------------------------------------
    # for each local assembly, dump the paired end reads
    ctgs = list(n_ctgs)
    ctgs.append(root_ctg)
  
    fouts_map = {}
    nctg_fout_map = {}
    paths = []
    for n_ctg in n_ctgs:
      path = os.path.join(self.outdir, '{}.bam'.format(n_ctg))
      paths.append(path)
      nctg_fout_map[n_ctg] = pysam.Samfile(path, 'wb', template=fhandle)
    path = os.path.join(self.outdir, 'all.bam')
    paths.append(path)
    nctg_fout_map['all'] = pysam.Samfile(path, 'wb', template=fhandle)
  
    qname_nctgs_map = defaultdict(set)
    for n_ctg, qname_set in link_ctg_qnames.items():
      if n_ctg == root_ctg or n_ctg not in n_ctgs:
        continue
      for qname in qname_set:
        qname_nctgs_map[qname].add(n_ctg)
  
    for ctg in ctgs:
      for read in fhandle.fetch(ctg):
        if read.qname not in qname_nctgs_map:
          continue
        for n_ctg in qname_nctgs_map[read.qname]:
          # FIXME remove
          if n_ctg not in nctg_fout_map:
            continue
          nctg_fout_map[n_ctg].write(read)
        nctg_fout_map['all'].write(read)
  
    # close
    fhandle.close()
    for fout_handle in nctg_fout_map.values():
      fout_handle.close()
  
    # do local reassembly of the root contig
    ds_bcode_set = root_bcode_set
    if len(root_bcode_set) > 400:
      ds_bcode_set = set(random.sample(root_bcode_set, 400))
  
    # FIXME TODO
    # for large contigs only locally reassemble the end 10kb
    local_asms.append(
      (root_ctg, root_bcode_set, ds_bcode_set)
    )
  
    return local_asms
  
def get_lrhints_fa(
  ctgs,
  infa_path,
  outfa_path,
):
  with open(outfa_path, 'w') as fout:
    ctg_fasta = pysam.FastaFile(infa_path)
    for ctg in ctgs:
      seq = str(ctg_fasta.fetch(ctg).upper())
      for _ in xrange(10):
        fout.write('>{}\n{}\n'.format(ctg, seq))

def get_bcode_reads(infq_paths, outfa_path, bcode_set):
  seen_set = set()
  with open(outfa_path, 'w') as fout:
    for fq_path in infq_paths:
      with FastqIndex(fq_path) as idx:
        for bcode in bcode_set & idx.bcode_set:
          for e in idx.get_reads(bcode):
            # tag qname with barcode
            nqname = '{}${}'.format(bcode, e.qname)
            fout.write('>{}\n'.format(nqname))
            fout.write('{}\n'.format(e.seq1))
            fout.write('>{}\n'.format(nqname))
            fout.write('{}\n'.format(e.seq2))
            seen_set.add(e.bcode)
  assert seen_set.issubset(bcode_set), 'not all barcodes loaded'
  return

#def filter_alignments(in_path, out_path):
#  def get_clip_info(read):
#    (c, l) = read.cigar[0]
#    lclip = (c in [4,5])
#    (c, l) = read.cigar[-1]
#    rclip = (c in [4,5])
#  
#    return (lclip, rclip)
  
