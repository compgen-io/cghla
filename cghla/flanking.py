import sys
import cghla
import cghla.samfile

def hla_flanking_fasta(hla, sam, ref, flanking=1000):
    '''
hla_flanking_fasta - Given a SAM alignment, reference HLA FASTA, and a reference genome, 

Return a new HLA FASTA file with X bp of extra flanking sequence (default 1kb).
This will write a new FASTA file to stdout. All flanking sequences will start at the
same position, regardless of where the specific allele has mapped. This means that each
flanking sequence will be comparable across all alleles.

Requires: samtools

Arguments:
--ref genome.fa (must be faidx indexed)
--hla hla.fa
--sam hla.genome.sam (gzip optional)
--flanking (default: 1000)
'''

    alleles = cghla.HLAalleles(hla)

    gene_windows = find_allele_gene_window(hla, sam)

    gene_chrom = {}
    for gene in gene_windows:
        gene_chrom[(gene, gene_windows[gene][0])] = (gene_windows[gene][1], gene_windows[gene][2])


    with cghla.open_file(sam, 'rt') as f:
        for line in f:
            if not line.strip() or line[0] == '@':
                continue

            cols = line.strip('\n').split('\t')

            qname = cols[0]
            flags = cghla.samfile.parse_flags(int(cols[1]))
            chrom = cols[2]
            pos = int(cols[3])
            cigar = cols[5]
            seq = cols[9]

            if flags.unmapped:
                continue

            refend = pos + cghla.samfile.cigar_to_reflen(cigar)

            gene = alleles.allele_gene[qname]

            gene_flank_start = gene_chrom[(gene, chrom)][0] - flanking - 1
            gene_flank_end = gene_chrom[(gene, chrom)][1] + flanking

            left = cghla.get_faidx_subseq(ref, chrom, gene_flank_start, pos-1)
            right = cghla.get_faidx_subseq(ref, chrom, refend-1, gene_flank_end)
            
            # we are using this seq, so we don't need to worry about +/- strandedness
            # seq in SAM is always +, left/right from FASTA is also +
            flanking_seq = '%s%s%s' % (left, seq, right)

            cghla.write_fasta('%s %s %s:%s-%s' % (qname, alleles.allele_accn[qname], chrom, pos, refend), '%s%s%s' % (left, seq, right), 60)


def find_allele_gene_window(hla, sam):
    allele_genes = cghla.HLAalleles(hla)

    starts = {}
    ends = {}

    with cghla.open_file(sam, 'rt') as f:
        for line in f:
            if not line.strip() or line[0] == '@':
                continue

            cols = line.strip('\n').split('\t')

            qname = cols[0]
            flags = cghla.samfile.parse_flags(int(cols[1]))
            chrom = cols[2]
            pos = int(cols[3])
            cigar = cols[5]

            if flags.unmapped:
                continue

            refend = pos + cghla.samfile.cigar_to_reflen(cigar)
            gene = allele_genes.allele_gene[qname]

            if not (gene, chrom) in starts:
                starts[(gene, chrom)] = pos
                ends[(gene, chrom)] = refend
            else:
                if pos < starts[(gene, chrom)]:
                    starts[(gene, chrom)] = pos
                if refend > ends[(gene, chrom)]:
                    ends[(gene, chrom)] = refend

    ret = {}
    for (gene, chrom) in starts:
        ret[gene] = (chrom, starts[(gene, chrom)], ends[(gene, chrom)])
    
    return ret



def hla_to_bed(hla, sam):
    '''
hla_to_bed - Given a SAM alignment, reference HLA FASTA, 
find the bounding region that contains each HLA gene (-A, -B, -C, etc)...

This will write a new BED file to stdout.

Arguments:
--hla hla.fa
--sam hla.genome.sam (gzip optional)
'''

    gene_windows = find_allele_gene_window(hla, sam)

    for gene in gene_windows:
        # zero-based starts
        sys.stdout.write('%s\t%s\t%s\tHLA-%s\n' % (gene_windows[gene][0], gene_windows[gene][1]-1, gene_windows[gene][2], gene))

