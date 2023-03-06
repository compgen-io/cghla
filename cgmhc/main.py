import cgmhc
import sys

def usage():
    sys.stderr.write("""Usage: cgmhc cmd {args}

Valid commands (in order):
    hla_align_ref         Align HLA reference alleles FASTA to a genome reference (bwa)
    hla_flanking_fasta    Add flanking reference sequence to HLA alleles
    hla_to_bed            Calculate BED regions covering HLA alleles
    mask_ref              Mask out the HLA regions from a reference FASTA
    extract_reads         Extract sample reads from an aligned BAM file
    align_to_hla          Align extracted reads to HLA (flanking) FASTA file
    score_pairs           Generate scores for HLA allele pairs
    predict               Find the most likely allele pairs

""")

def help(func):
    sys.stderr.write('Usage: cgmhc %s {args}\n' % func.__name__) 
    sys.stderr.write('%s\n' % func.__doc__)
    sys.exit(1)

def main(argv):
    cmd = None
    arg_kv = {}

    last = None
    for arg in argv:
        if not cmd and arg[0] != '-':
            cmd = arg
        elif last:
            try:
                val = float(arg)
            except ValueError:
                try:
                    val = int(arg)
                except ValueError:
                    val = arg
            if last in arg_kv:
                if type(arg_kv[last]) == list:
                    arg_kv[last].append(val)
                else:
                    arg_kv[last] = [arg_kv[last], val]
            else:
                arg_kv[last] = val
            last = None
        elif arg[:2] == '--':
            if last:
                arg_kv[last] = True

            last = arg[2:]

    if last:
        arg_kv[last] = True

    if cmd == 'hla_align_ref':
        if 'ref' not in arg_kv or 'hla' not in arg_kv:
            help(cgmhc.hla_align_ref)
        cgmhc.hla_align_ref(**arg_kv)

    elif cmd == 'hla_flanking_fasta':
        if 'ref' not in arg_kv or 'hla' not in arg_kv or 'sam' not in arg_kv:
            help(cgmhc.hla_flanking_fasta)
        cgmhc.hla_flanking_fasta(**arg_kv)

    elif cmd == 'hla_to_bed':
        if 'hla' not in arg_kv or 'sam' not in arg_kv:
            help(cgmhc.hla_to_bed)
        cgmhc.hla_to_bed(**arg_kv)

    elif cmd == 'mask_ref':
        if 'bed' not in arg_kv or 'ref' not in arg_kv:
            help(cgmhc.mask_ref)
        cgmhc.mask_ref(**arg_kv)

    elif cmd == 'align_to_hla':
        if not 'fq1' in arg_kv or 'ref' not in arg_kv:
            help(cgmhc.align_to_hla)
        cgmhc.align_to_hla(**arg_kv)

    elif cmd == 'extract_reads':
        if not 'bam' in arg_kv or 'bed' not in arg_kv or 'out' not in arg_kv:
            help(cgmhc.extract_reads)
        cgmhc.extract_reads(**arg_kv)

    elif cmd == 'score_pairs':
        if not 'sam' in arg_kv or 'hla' not in arg_kv:
            help(cgmhc.score_pairs)
        cgmhc.score_pairs(**arg_kv)

    elif cmd == 'predict':
        if not 'scores' in arg_kv or 'hla' not in arg_kv or 'motifs' not in arg_kv:
            help(cgmhc.predict)
        cgmhc.predict(**arg_kv)

    else:
        usage()
