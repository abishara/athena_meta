Athena
--------

Athena is a read cloud assembler for metagenomes.


Installation
============

To run Athena through a Docker image with Athena and its prerequisities
already installed, please skip to section Docker (and example dataset).

To install Athena in your native environment, the following prerequisites
must be installed:

* `idba_ud <https://github.com/abishara/idba/releases/tag/1.1.3a1>`_ --
please use **this** version, which is modified both to handle longer
short-read lengths and to locally assemble subsampled barcoded reads
clouds.  Ensure all compiled binaries, including ``idba_subasm``, are in
your ``$PATH``
* `samtools and htslib <http://www.htslib.org/download/>`_ -- version 1.3 or later of ``samtools`` must all be in your ``$PATH``
* `bwa-mem <https://github.com/lh3/bwa/releases>`_
* `canu <https://github.com/marbl/canu>`_ -- version 1.3 or later

We recommend setting up a `virtualenv
<http://docs.python-guide.org/en/latest/dev/virtualenvs/>`_ prior to
installing Athena (or using `virtualenvwrapper
<http://www.simononsoftware.com/virtualenv-tutorial-part-2/>`_):

.. code-block:: bash

    sudo pip install virtualenv
    virtualenv athena_meta

Then, to install 

.. code-block:: bash

    cd /path/to/athena_meta
    pip install -r requirements.txt
    pip install .

To test that Athena is installed correctly, you can simply run
``athena-meta`` from the commandline, which should show help text without
any error messages.


Running Athena
================

Overview:

1. Generate input seed contigs for Athena with Spades/idba_ud.  Align barcoded input reads to seed contigs with BWA.
2. Setup a ``config.json`` file, which specifies inputs to Athena and your compute (eg cluster) setup
3. Run Athena


1. Generate inputs
"""""""""""""""""""""""""""""""""""

Input read clouds must be specified as an uncompressed paired-end
interleaved FASTQ, with the following tag information as in the example
read pair below:

.. code-block:: txt

  @NS500418:354:H27G3BGXY:3:12612:25572:11380 RG:Z:readgroup1 BC:Z:GCCAATTCAAGTTT-1
  TTCCATGTGGAAGTAGTTGTATTTGACGTAGCCCGCCATACCGTTTTCTGACATGAAGCGGTAATTCTCCTCAGAACCGTAGCCGGATACGGCCACCACCGTATGGGCCAACCTGTCATATCTGCTTGAGAAGGATTG
  +
  EEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE/EEEEEEEEEEEEEEEEEEE6EEAEEEEAEEEEEEEEEEEEEEAEEEEEEEEEEEEEAEEEEEE
  @NS500418:354:H27G3BGXY:3:12612:25572:11380 RG:Z:readgroup1 BC:Z:GCCAATTCAAGTTT-1
  CACGTGGTCTGGCGGGTCTCGCGCCACCTCTGGTTCGCCGTGGCCCTAACGGACAAGGACGCTACTTTCATGAGAATGAAGGAGGATGCCATGCGTAACGGCCAGACAAAGCCCGGTTACAACCTCCAGAACGGCACCGAGAACCAGA
  +
  AAA<A/AE///<///<EEEE<AE/AEA/EAE</EAAEEE/EEAE6AE/<E<<EE<AA<AEEAA/E6/6AEA/</EEA/A/AAAEE</<EEA6<<AA<<EEEEEA//EA<<AE<EA/66<EA/EE6<A////A/AA6EA/66/6AA/A6 

For each barcoded read, there must be a tag (either BC or BX, but not
both) specifying the barcode.  Note that the query name line for each read
can have mulitple space or tab delimited tags. 
 
Run Spades or idba_ud out of the box to assemble your input barcoded read
clouds into seed contigs.  An example using Spades:

.. code-block:: bash

   spades.py --12 /path/to/reads  -o /path/to/spades/out

Create a BWA index for the ouput short-read draft assembly.
Then run BWA MEM specifying -C, to pass the FASTQ tags through to the BAM, and
-p, to specify that the paired-end FASTQ is interleaved:

.. code-block:: bash

   bwa index /path/to/spades/out/contigs.fasta
   bwa mem -C -p /path/to/spades/out/contigs.fasta /path/to/reads | samtools sort -o align-reads.spades-contigs.bam -

2. Setup a configuration file
"""""""""""""""""""""""""""""

See the examples directory for an example ``config.json`` file.

The configuration file is in the `JSON <http://www.json.org>`_ format, and contains the following three parts:

1. input barcoded reads (FASTQ).  Must be uncompressed interleaved paired end reads, which specify barcodes with the BC tag as specified above.
2. input seed contigs (FASTA).  Must be a path to a BWA built index.
3. BWA alignments of barcoded input reads to input seeds (BAM)
4. compute cluster settings

**Data Inputs** The following paths must be defined:

* ``input_fqs``: path to input uncompressed interleaved paired-end FASTQ (must specify barcodes with the BC tag as specified above.)
* ``ctgfasta_path``: path to input seed contigs (must be BWA indexed)
* ``reads_ctg_bam_path``: path to BAM of input reads BWA aligned to input seed contigs (alignments must have BC tag with barcode information)

**Compute cluster settings** This defines the compute environment being
used to perform assembly.  Athena manages the environment using
`ipython-cluster-helper
<https://github.com/roryk/ipython-cluster-helper>`_, however, support for
Canu OLC to use a cluster is still under develoopment.  For now, we encourage
running Athena on a large multicore machine.  

A multiprocessing setup looks like this:

.. code-block:: json

  "cluster_settings": {
    "cluster_type": "multiprocessing",
    "processes": 8
  }

Where ``processes`` specifies the maximum number of separate jobs (1
processor per job) to allow in flight.  Each job can use up to 4G of
memory, so be sure not oversubscribe the host machine.

To use a compute cluster (not yet fully supported), a setup looks like this:

.. code-block:: json

  "cluster_settings": {
    "cluster_type": "IPCluster",
    "processes": 128,
    "cluster_options": {
      "scheduler": "slurm",
      "queue": "normal",
      "extra_params": {"mem":16}
    }
  }

``scheduler`` may be any of the clusters supported by
`ipython-cluster-helper`. Currently, these are
Platform LSF ("lsf"), Sun Grid Engine ("sge"), Torque ("torque"), and
SLURM ("slurm").  The Canu OLC step will run as a single job with a single
node.

3. Run Athena
"""""""""""""""

To run Athena, use the ``athena-meta /path/to/config.json`` command. 

Note that the ``athena-meta`` command will continue running until all
steps have completed. The ``athena-meta`` command itself is lightweight,
and so can be run from a head node if the configuration is setup to use a
cluster.  If running on a local machine in multiprocessing mode, please be
aware that some subassembly problems can require up to 4G of memory.
Adjust the number of ``processes`` to prevent oversubscription of the
machine.

The output assembled contigs will be placed in a subdirectory of the one
``config.json`` resides in (in this case
``/path/to/results/olc/athena.asm.fa``.) Logging output for each step will
also be in the subdirectory ``logs`` (in this case ``/path/to/logs``),
which can be used to debug in event of an error.

Docker (and example dataset)
============================

A docker image is available for Athena.  To download and run
``athena-meta`` on the example read clouds (~108MB), you can run the
following commands:

.. code-block:: bash
    
    # use 'curl -O' if you're on a mac without wget
    wget https://storage.googleapis.com/gbsc-gcp-lab-bhatt-public/readclouds-meta-asm-example.tar.gz
    tar -xzf readclouds-meta-asm-example.tar.gz

Assuming `docker <https://docs.docker.com/engine/installation/>`_ is
installed, the following command can be used to assemble the example read
clouds from within docker (make sure you are in the same directory where
you downloaded and extracted readclouds-meta-asm-example.tar.gz):

.. code-block:: bash

    docker run -v `pwd`:/data -w /data/readclouds-meta-asm-example abishara/athena-meta-docker athena-meta config.json

This requires ~16GB of memory to run (for OLC assembly) and will take ~20
minutes to complete. If you are running docker for Mac, please make sure
that your docker client has access to at least 16GB of memory (you may
need to set in Preferences).

The output can be found in native host directory of
``readclouds-meta-asm-example``.

Troubleshooting
===============

The ``athena-meta`` command may be run multiple times to resume the pipeline.

If you are having trouble installing or running Athena, the docker file
(see above) may help you diagnose the issue.

If an error arises, the output from ``athena-meta`` or the log files may
be informative.

**ShortSequence: Sequence is too long.** If you get this error during
assembly, please make sure you are using `the right fork of idba_ud
<https://github.com/abishara/idba/releases/tag/1.1.3a1>`_.

Please submit issues on the `github page for Athena
<https://github.com/abishara/athena_meta/issues>`_.

