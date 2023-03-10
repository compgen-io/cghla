import sys
import time
import subprocess
import gzip
import os
import contextlib
import typing
from collections import namedtuple

import cghla.flanking
import cghla.extract
import cghla.score
import cghla.samfile
import cghla.motif_similarity


@contextlib.contextmanager
def open_file(fname, mode) -> typing.IO[str]:
    if fname == '-':
        # sys.stderr.write('using stdin\n')
        yield sys.stdin
        return

    f = open(fname, 'rb')
    magic = f.read(2)
    f.close()

    if magic[0] == 0x1F and magic[1] == 0x8B:
        # sys.stderr.write('gzipped input file\n')
        f = gzip.open(fname, mode)
    else:
        # sys.stderr.write('uncompressed input file\n')
        f = open(fname, mode)

    yield f
    f.close()
    return


class HLAallele(object):
    def __init__(self, accn, allele):
        self._accn = accn
        self._allele = allele
        self._gene = allele.split('*')[0]
        self._org = accn.split(':')[0]

    @staticmethod
    def read_fasta(fname: str) -> list['HLAallele']:
        ''' 
        Reads the allele accn and allele name/group from a reference HLA FASTA file
        '''

        alleles: list[HLAallele] = []
        with open_file(fname, 'rt') as f:
            for line in f:
                if not line.strip() or line[0] != '>':
                    continue

                spl = line.strip('\n').split()  # on any whitespace
                accn = spl[0][1:]               # remove '>' prefix
                allele = spl[1]
                gene = allele.split('*')[0]

                alleles.append(HLAallele(accn, allele))

        return alleles


class HLAalleles(object):
    def __init__(self, fname):
        self.allele_gene = {}
        self.allele_accn = {}
        self.genes = {}
        for allele in HLAallele.read_fasta(fname):
            self.allele_gene[allele._accn] = allele._gene
            self.allele_accn[allele._accn] = allele._allele

            if not allele._gene in self.genes:
                self.genes[allele._gene] = []
            self.genes[allele._gene].append(allele._accn)


    def sorted_alleles(self, gene=None):
        if gene:
            tmp = []
            for accn in self.genes[gene]:
                allele = self.allele_accn[accn]
                tmp.append((allele, accn))
            
            out = []
            for e in sorted(tmp):
                out.append(e[1])
            
            return out
        
        else:
            tmp = []
            for accn in self.allele_accn:
                allele = self.allele_accn[accn]
                tmp.append((allele, accn))
            
            out = []
            for e in sorted(tmp):
                out.append(e[1])
            
            return out

    @staticmethod
    def to_4digit(allele):
        spl = allele.split(':')
        return '%s:%s' % (spl[0], spl[1])


class Logger(object):
    def __init__(self, fobj=sys.stderr, min_sec=1):
        self.__fobj = fobj
        self.__last = ''
        self.__last_time = 0
        self.__min_sec = min_sec
    
    def write(self, msg=None, func=None, /, *args, **kwargs):
        curtime = time.time() 
        if not msg or not '\n' in msg:
            if curtime - self.__last_time < self.__min_sec:
                return

        if not msg and func:
            msg = func(*args, **kwargs)
        elif func:
            msg = msg + func(*args, **kwargs)

        if self.__last:
            self.__fobj.write('\r')
            for v in range(len(self.__last)):
                self.__fobj.write(' ')
            self.__fobj.write('\r')

        self.__fobj.write(msg)
        self.__fobj.flush()
        self.__last = msg
        self.__last_time = curtime

    def clear(self):
        self.write('')


def write_fasta(name, seq, wrap=-1):
    sys.stdout.write('>%s\n' % name)
    if wrap > 0:
        while seq:
            sys.stdout.write('%s\n' % (seq[:wrap]))
            seq = seq[wrap:]
    else:
        sys.stdout.write(seq)


def get_faidx_subseq(fasta, chrom, start, end) -> str: 
    '''
    Note: start is 0-based!!!
    '''

    proc = subprocess.run(['samtools', 'faidx', fasta, '%s:%s-%s' % (chrom, start+1, end)], stdout=subprocess.PIPE)
    retval = proc.stdout.decode('utf-8')
    seq = ""
    for i, line in enumerate(retval.split('\n')):
        if i > 0:
            seq += line.strip()

    return seq


def revcomp(seq):
    ret = ''
    for b in seq.upper()[::-1]:
        if b == 'A':
            ret += 'T'
        elif b == 'T':
            ret += 'A'
        elif b == 'C':
            ret += 'G'
        elif b == 'G':
            ret += 'C'
        else:
            ret += 'N'

    return ret


def hla_align_ref(ref, hla, threads=1):
    '''
hla_align_ref - Aligns the HLA FASTA file to a genome reference.

This will write a SAM file to stdout. This file is used for 
determining the appropriate flanking regions for each HLA loci.

The genome reference FASTA must be indexed with bwa.

Requires: bwa

Arguments:
--ref genome.fa (bwa indexed)
--hla hla.fa

--threads N (optional)
'''

    args = ['bwa', 'mem']
    
    if threads > 1:
        args.append('-t')
        args.append(str(threads))

    args.append(ref)
    args.append(hla)    

    subprocess.call(args, stdout=sys.stdout)


def mask_ref(ref, bed):
    '''
mask_ref - Using the HLA BED file coordinates, mask out the HLA regions from genome FASTA file.

This will write a new FASTA file to stdout with the HLA regions masked out.

Requires: ngsutilsj

Arguments:
--ref genome.fa
--bed hla.bed
'''

    args = ['ngsutilsj', 'fasta-mask', '--bed', bed, ref]
    
    subprocess.call(args, stdout=sys.stdout)


def hla_flanking_fasta(hla, sam, ref, flanking=1000):
    '''
hla_flanking_fasta - Given a SAM alignment, reference HLA FASTA, and a reference genome, 

Return a new HLA FASTA file with X bp of extra flanking sequence (default 1kb).
This will write a new FASTA file to stdout.

Requires: samtools

Arguments:
--ref genome.fa (must be faidx indexed)
--hla hla.fa
--sam hla.genome.sam (gzip optional)
--flanking (default: 1000)
'''

    cghla.flanking.hla_flanking_fasta(hla, sam, ref, flanking)


def hla_to_bed(hla, sam):
    '''
hla_to_bed - Given a SAM alignment, reference HLA FASTA, 
find the bounding region that contains each HLA gene (-A, -B, -C, etc)...

This will write a new BED file to stdout.

Arguments:
--hla hla.fa
--sam hla.genome.sam (gzip optional)
'''

    cghla.flanking.hla_to_bed(hla, sam)


def extract_reads(bam, bed, out, bam2=None, unmapped=None):
    '''
extract_reads - Extract reads from an aligned BAM file (and optionally a failed/unmapped BAM file)

This will write two new FASTQ files.

Requires: samtools, ngsutilsj

Arguments:
--bam aligned BAM file
--bed HLA BED file
--bam2 failed/unmapped BAM file to find paired reads (optional, see: ngsutilsj bam-extract --bam2)
--unmapped failed/unmapped BAM file (optional)
--out Base name for output FASTQ files (will write {out}_R1.fastq.gz, {out}_R2.fastq.gz)
'''

    cghla.extract.extract_reads(bam, bed, out, bam2, unmapped)


def score_pairs(sam, hla, min_as=0):
    '''
score_pairs - Score alignments for potential HLA allele pairs. 

This will write a new score matrix file to stdout.

Arguments:
--sam Reads aligned to HLA (flanked) FASTA
--hla HLA FASTA file
--min_as minimum alignment score (default 0)
'''

    cghla.score.score_sam(sam, hla, min_as)



def predict(scores, hla, motifs, thres=0.95):
    '''
predict - Baseds on the pair-wise allele scores, determine which pair of alleles are the most likely genotype. 

This will write a new score matrix file to stdout.

Arguments:
--scores scores.txt.gz    Pair-wise scores file (you can include more than one score file)
--hla file.fasta          HLA FASTA file
--motifs file.csv         MHCFlurry common motifs file (mhcflurry.allele_sequences.csv)
--thres val               report secondary genotypes that fall within this threshold of the best score (0-1.0, default 0.95)
'''

    cghla.score.predict(scores, hla, motifs, thres)


def align_to_hla(ref, fq1, fq2=None, threads=1):
    '''
align_to_hla - Align reads to the HLA reference FASTA file. This can be either
a (flanked) WGS FASTA file (_gen.fasta) or the mRNA FASTA (_nuc.fasta). The 
FASTA file can contain multiple genes (HLA-A, -B, -C). Both WGS and RNAseq 
inputs will be aligned using BWA MEM (with -a argument).

This will write a SAM file to stdout. Unmapped reads will be silently removed.

Requires: bwa

Arguments:
--ref Reference FASTA file to align reads to (bwa indexed)
--fq1 FASTQ file to align (R1)
--fq2 FASTQ file to algin (R2, optional)

--threads N (optional)
'''

    args = ['bwa', 'mem', '-a']
    
    if threads > 1:
        args.append('-t')
        args.append(str(threads))

    args.append(ref)
    args.append(fq1)
    
    if fq2:
        args.append(fq2)

    proc = subprocess.Popen(args, bufsize=1, text=True, stdout=subprocess.PIPE, stderr=None)

    for line in iter(proc.stdout.readline, ''):
        if not line.strip() or line[0] == '@':
            sys.stdout.write(line)
            continue

        cols = line.strip('\n').split('\t')
        flags = cghla.samfile.parse_flags(int(cols[1]))

        if flags.unmapped:
            # remove any unmapped reads (there may be a lot)
            continue

        sys.stdout.write(line)

    proc.stdout.close()
    proc.wait()


def similarity(motifs, alleles=None, f=None, length=9, cutoff_fraction=0.01, org=None):
    '''
similarity - Calulate the cosine similarity between MHC alleles' peptide binding motifs

If two allele are given, the cosine similarity between these two are calulated. Otherwise, all
comparisons are made between all pairs of alleles present in the motifs file.

Arguments:
--motifs file.csv    MHCFlurry frequency matrix file (mhcflurry.ba.frequency_matrices.csv.gz)

Optional arguments:

--alleles a1,a2,...   Comma separated list of alleles to compare
--f file.txt          File with one allele to compare per line

--org org             Add this org to allele names from --f (ex: --org HLA)

--length              Length of peptides to use (default: 9)
--cutoff_fraction     Cut-off fraction to use (default 0.01)

'''
    valid = None
    if alleles:
        valid = alleles.split(',')
    
    if f:
        valid = []
        with open(f, 'rt') as f:
            for line in f:
                if org:
                    valid.append('%s-%s' % (org, line.strip()))
                else:
                    valid.append(line.strip())
                

    cghla.motif_similarity.main(motifs, length, cutoff_fraction, valid)
