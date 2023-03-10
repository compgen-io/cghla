#!/usr/bin/env python3 
import cghla
import sys


def usage():
    sys.stderr.write("""Usage: cghla cmd {args}

Valid commands (in order):
    hla_align_ref         Align HLA reference alleles FASTA to a genome reference (bwa)
    hla_flanking_fasta    Add flanking reference sequence to HLA alleles
    hla_to_bed            Calculate BED regions covering HLA alleles
    mask_ref              Mask out the HLA regions from a reference FASTA
    extract_reads         Extract sample reads from an aligned BAM file
    align_to_hla          Align extracted reads to HLA (flanking) FASTA file
    score_pairs           Generate scores for HLA allele pairs
    predict               Find the most likely allele pairs
    similarity            Calculate the cosine similarity between motifs

""")


def help(func, arg_kv=None):
    sys.stderr.write('Usage: cghla %s {args}\n' % func.__name__) 
    sys.stderr.write('%s\n' % func.__doc__)
    if arg_kv:
        sys.stderr.write('\n%s\n' % arg_kv)
    sys.exit(1)


def main(argv):
    cmd = None
    arg_kv = {}
    argv2 = []

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
        else:
            argv2.append(arg)

    if last:
        arg_kv[last] = True

    if cmd == 'hla_align_ref':
        if 'ref' not in arg_kv or 'hla' not in arg_kv:
            help(cghla.hla_align_ref)
        cghla.hla_align_ref(**arg_kv)

    elif cmd == 'hla_flanking_fasta':
        if 'ref' not in arg_kv or 'hla' not in arg_kv or 'sam' not in arg_kv:
            help(cghla.hla_flanking_fasta)
        cghla.hla_flanking_fasta(**arg_kv)

    elif cmd == 'hla_to_bed':
        if 'hla' not in arg_kv or 'sam' not in arg_kv:
            help(cghla.hla_to_bed)
        cghla.hla_to_bed(**arg_kv)

    elif cmd == 'mask_ref':
        if 'bed' not in arg_kv or 'ref' not in arg_kv:
            help(cghla.mask_ref)
        cghla.mask_ref(**arg_kv)

    elif cmd == 'align_to_hla':
        if not 'fq1' in arg_kv or 'ref' not in arg_kv:
            help(cghla.align_to_hla)
        cghla.align_to_hla(**arg_kv)

    elif cmd == 'extract_reads':
        if not 'bam' in arg_kv or 'bed' not in arg_kv or 'out' not in arg_kv:
            help(cghla.extract_reads)
        cghla.extract_reads(**arg_kv)

    elif cmd == 'score_pairs':
        if not 'sam' in arg_kv or 'hla' not in arg_kv:
            help(cghla.score_pairs)
        cghla.score_pairs(**arg_kv)

    elif cmd == 'predict':
        if not 'scores' in arg_kv or 'hla' not in arg_kv or 'motifs' not in arg_kv:
            help(cghla.predict)
        cghla.predict(**arg_kv)

    elif cmd == 'similarity':
        if 'motifs' not in arg_kv:
            help(cghla.similarity, arg_kv)

        if 'allele' in arg_kv and type(arg_kv) != 'list':
            help(cghla.similarity, arg_kv)

        cghla.similarity(**arg_kv)

    else:
        sys.stderr.write("Unknown command: %s\n" % (argv))
        usage()


if __name__ == '__main__':
    main(sys.argv[1:])
