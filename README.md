# Athena

Athena is a read cloud assembler for metagenomes.

**Recent updates**

* **v1.3 release**: Updates to command-line arguments and logging.

## Installation

To run Athena through a Docker image with Athena and its prerequisities
already installed, please skip to section Docker (and example dataset).

To install Athena in your native environment, the following prerequisites
must be installed:

* python 2.7 on Mac and Linux; Athena is not compatible with python 3.x
* [idba_ud](https://github.com/abishara/idba/releases/tag/1.1.3a1) --
  please use **this** version, which is modified both to handle longer
  short-read lengths and to locally assemble subsampled barcoded reads
  clouds.  Ensure all compiled binaries, including `idba_subasm`, are in
  your `$PATH`
* [samtools and htslib](http://www.htslib.org/download/) -- version 1.3 or
  later of `samtools` must all be in your `$PATH`
* [bwa-mem](https://github.com/lh3/bwa/releases)
* [flye](https://github.com/fenderglass/Flye) -- version 2.3.1 

We recommend setting up a [virtualenv](http://docs.python-guide.org/en/latest/dev/virtualenvs/) prior to
installing Athena (or using [virtualenvwrapper](http://www.simononsoftware.com/virtualenv-tutorial-part-2/)):

```bash
sudo pip install virtualenv
virtualenv athena_meta
```

Then, to install 

```bash
cd /path/to/athena_meta
pip install .
```

To test that Athena is installed correctly, you can simply run
`athena-meta` from the commandline, which should show help text without
any error messages.


## Running Athena

Overview:

1. Generate input seed contigs for Athena with metaspades/idba_ud.  Align
   barcoded input reads to seed contigs with `bwa`.
2. Setup a `config.json` file, which specifies inputs to Athena
3. Run Athena


#### 1. Generate inputs

Input read clouds must be specified as an uncompressed paired-end
interleaved FASTQ, with the following tag information as in the example
read pair below:

```text
@NS500418:354:H27G3BGXY:3:12612:25572:11380	RG:Z:rg-1	BC:Z:GCCAATTCAAGTTT-1
TTCCATGTGGAAGTAGTTGTATTTGACGTAGCCCGCCATACCGTTTTCTGACATGAAGCGGTAATTCTCCTCAGAACCGTAGCCGGATACGGCCACCACCGTATGGGCCAACCTGTCATATCTGCTTGAGAAGGATTG
+
EEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE/EEEEEEEEEEEEEEEEEEE6EEAEEEEAEEEEEEEEEEEEEEAEEEEEEEEEEEEEAEEEEEE
@NS500418:354:H27G3BGXY:3:12612:25572:11380	RG:Z:rg-1	BC:Z:GCCAATTCAAGTTT-1
CACGTGGTCTGGCGGGTCTCGCGCCACCTCTGGTTCGCCGTGGCCCTAACGGACAAGGACGCTACTTTCATGAGAATGAAGGAGGATGCCATGCGTAACGGCCAGACAAAGCCCGGTTACAACCTCCAGAACGGCACCGAGAACCAGA
+
EEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE6EEEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEEE6EEEEEEEEEEEEAEEE
```

The barcoded interleaved FASTQ must satisfy the following:

* For each barcoded read, there must be a tag (either BC or BX, but not
  both) specifying the barcode.  Each barcode must also end with a '-'
  followed by an integer sample identifier.
* The query name line for each read can have mulitple tags, but these
  **must** be tab-delimited to be compatible with `bwa mem` specifying `-C`.
* The input FASTQ file **must** be barcode-sorted such that all reads with
  the same attached barcode appear in a contiguous block.
 
Run `metaspades` or `idba_ud` out of the box to assemble your input barcoded read
clouds into seed contigs.  An example using metaspades:


```bash
metaspades.py --12 /path/to/reads  -o /path/to/metaspades/out
```

Create a `bwa` index for the ouput short-read draft assembly.  Then run
`bwa mem` specifying `-C`, to pass the FASTQ tags through to the BAM, and
`-p`, to specify that the paired-end FASTQ is interleaved:

```bash
bwa index /path/to/metaspades/out/contigs.fasta
bwa mem -C -p /path/to/metaspades/out/contigs.fasta /path/to/reads | samtools sort -o align-reads.metaspades-contigs.bam -
samtools index align-reads.metaspades-contigs.bam
```

Note that the resulting BAM must be position sorted and indexed.

#### 2. Setup a configuration file

The configuration file is in the [JSON](http://www.json.org) format, and
contains the following parts:

1. `input_fqs`: path to barcoded reads (FASTQ). Must be uncompressed interleaved
   paired end reads, which specify barcodes with the BC tag as specified
   above.
2. `ctgfasta_path`: path to seed contigs (FASTA), which must be `bwa` indexed
3. `reads_ctg_bam_path`:  alignments of barcoded input reads to seed
   contigs (alignments **must** have BC tag with barcode information per
   read).
4. (optional) `cluster_settings`: cluster compute environment to be used
   to perform assembly if a batch queueing submission system is available.
   Athena manages the environment using
   [ipython-cluster-helper](https://github.com/roryk/ipython-cluster-helper)

A minimal `config.json` file contains the following:
```json
{
    "input_fqs": "/path/to/fq",
    "ctgfasta_path": "/path/to/seeds.fa",
    "reads_ctg_bam_path": "/path/to/reads_2_seeds.bam"
}
```

An example `cluster_settings` entry specifying a compute cluster
contains the following:

```json
{
    "cluster_settings": {
        "cluster_type": "IPCluster",
        "processes": 128,
        "cluster_options": {
            "scheduler": "slurm",
            "queue": "normal",
            "extra_params": {"mem":16}
        }
    }
}
```

`scheduler` may be any of the clusters supported by
`ipython-cluster-helper`. Currently, these are Platform LSF ("lsf"), Sun
Grid Engine ("sge"), Torque ("torque"), and SLURM ("slurm").
`processes` specifies the size of the job array to be used.

#### 3. Run Athena

To check all prerequisites are installed, run ``athena-meta --check_prereqs``.

To run a tiny test assembly to check that Athena is properly setup, run ``athena-meta --test``.

To run Athena on an input dataset, run ``athena-meta --config /path/to/config.json``.

Note that the `athena-meta` command will continue running until all steps
have completed.   `athena-meta` runs locally with a single thread by
default, but can be run using multiple threads by specifying `--threads`.
Please be aware that each thread can required up to 4Gb of memory during
the subassembly step and so the number of threads should be adjusted
accordingly.  If the config file provided specifies a cluster environment,
the `athena-meta` command itself can be run from a head node as it is
itself a lightweight process.

The output assembled contigs will be placed in a subdirectory of the one
`config.json` resides in (in this case
`/path/to/results/olc/athena.asm.fa`.) Logging output for each step will
also be in the subdirectory `logs` (in this case `/path/to/logs`), which
can be used to debug in event of an error.

## Docker (and example dataset)

A docker image is available for Athena.  To download and run
``athena-meta`` on the example read clouds (~46MB), you can run the
following commands:

```bash
# use 'curl -O' if you're on a mac without wget
wget https://storage.googleapis.com/gbsc-gcp-lab-bhatt-public/readclouds-l-gasseri-example.tar.gz

tar -xzf readclouds-l-gasseri-example.tar.gz
```

Assuming [docker](https://docs.docker.com/engine/installation/) is
installed, the following command can be used to assemble the example read
clouds from within docker (make sure you are in the same directory where
you downloaded and extracted `readclouds-meta-asm-example.tar.gz`):


```bash
docker run -v `pwd`:/data -w /data/readclouds-l-gasseri-example abishara/athena-meta-flye-docker athena-meta --config config.json
```

This requires ~16GB of memory to run (for overlap assembly) and will take ~20
minutes to complete. If you are running docker for Mac, please make sure
that your docker client has access to at least 16GB of memory (you may
need to set in Preferences).

The output can be found in native host directory of
`readclouds-meta-asm-example`.

## Citing Athena

Please cite the following publication:

* A. Bishara and E. Moss, et al.  High-quality genome sequences of
  uncultured microbes by assembly of read clouds. *Nature Biotechnology
  2018* (https://doi.org/10.1038/nbt.4266).

## Troubleshooting

The `athena-meta` command may be run multiple times to resume from the
last step successfully completed.

If you are having trouble installing or running Athena, the docker file
(see above) may help you diagnose the issue.

If an error arises, the output from `athena-meta` or the log files may
be informative.

**ShortSequence: Sequence is too long.** If you get this error during
assembly, please make sure you are using [the right fork of
idba_ud](https://github.com/abishara/idba/releases/tag/1.1.3a1).

Please submit issues on the [github page for
Athena](https://github.com/abishara/athena_meta/issues).

