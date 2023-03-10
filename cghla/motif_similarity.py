#
# Given a list of MHC alleles and their binding weight matrices, calculate pair-wise cosine similarities between alleles.
#
#

import sys
import gzip

def read_matrix(fname, length=9, cutoff_fraction=0.01, org=None, valid=None):
    counts = {}
    header = None
    aa_list = None
    f = gzip.open(fname, 'rt')
    for line in f:
        cols = line.strip('\n').split(',')
        if not header:
            header = {}
            aa_list = []
            for i,v in enumerate(cols):
                header[v] = i
                if len(v) == 1:
                    aa_list.append(v)
            continue

        allele1 = cols[header['allele']]
        length1 = int(cols[header['length']])
        cutoff_fraction1 = float(cols[header['cutoff_fraction']])

        if org and allele1.split('-')[0] != org:
            continue

        if valid and allele1 not in valid:
            continue

        if length and length1 != length:
            continue

        if cutoff_fraction and cutoff_fraction1 != cutoff_fraction:
            continue

        gene = allele1.split('*')[0]
        if allele1 not in counts:
            counts[allele1] = []

        aa_freq = []
        for aa in aa_list:
            aa_freq.append(float(cols[header[aa]]))

        counts[allele1].append(aa_freq)


    return counts


def main(freq_fname, length=9, cutoff_fraction=0.01, valid=None):
    counts = read_matrix(freq_fname, length, cutoff_fraction, valid=valid)
    calc_similarity(counts)


def calc_similarity(counts):
    alleles = []
    for k in counts:
        gene = k.split('*')[0]
        alleles.append(k)

    sys.stdout.write('gene1\tgene2\tallele1\tallele2\tsimilarity\n')

    for i, allele1 in enumerate(alleles):
        gene1 = allele1.split('*')[0]
        for allele2 in alleles[i+1:]:
            gene2 = allele2.split('*')[0]
            list1 = []
            list2 = []
            for pos_freq in counts[allele1]:
                for aa_freq in pos_freq:
                    list1.append(aa_freq)
            for pos_freq in counts[allele2]:
                for aa_freq in pos_freq:
                    list2.append(aa_freq)

            sim = cosine_similarity(list1, list2)

            sys.stdout.write('%s\t%s\t%s\t%s\t%s\n' % (gene1, gene2, allele1, allele2, sim))


def cosine_similarity(list1, list2):
    AB_acc = 0
    for a, b in zip(list1, list2):
        AB_acc += (a * b)

    a2_acc = 0
    for a in list1:
        a2_acc += (a * a)

    b2_acc = 0
    for b in list2:
        b2_acc += (b * b)

    a_sqrt = a2_acc ** 0.5
    b_sqrt = b2_acc ** 0.5

    return (AB_acc / (a_sqrt * b_sqrt))