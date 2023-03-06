import sys
import cgmhc


def score_sam(sam, hla, min_as = 0):
    alleles = cgmhc.HLAalleles(hla)

    read_scores = {}
    logger = cgmhc.Logger()

    sys.stderr.write('Reading SAM file\n')
    with cgmhc.open_file(sam, 'rt') as f:
        for line in f:
            if not line.strip() or line[0] == '@':
                continue

            cols = line.strip('\n').split('\t')
            qname = cols[0]
            flag = int(cols[1])
            allele = cols[2]

            unmapped = (flag & 0x4 == 0x4)
            rev = (flag & 0x10 == 0x10)
            nextrev = (flag & 0x20 == 0x20)
            first = (flag & 0x40 == 0x40)

            if unmapped:
                continue

            if not allele in alleles.allele_accn:
                sys.stderr.write("Unknown allele: %s\n" % allele)
                sys.exit(1)

            logger.write(qname)

            as_i = -1

            for col in cols[11:]:
                if col[:4] == 'AS:i':
                    as_i = int(col[5:])

            if as_i < min_as:
                # set a hard threshold
                continue

            if not qname in read_scores:
                read_scores[qname] = {}

            if not allele in read_scores[qname]:
                read_scores[qname][allele] = {}
            
            if first:
                read_scores[qname][allele]['R1'] = (as_i, '-' if rev else '+')
            else:
                read_scores[qname][allele]['R2'] = (as_i, '-' if rev else '+')

    logger.write('Done\n')

    sys.stderr.write('Finding best scores for reads\n')

    allele_scores = {}
    for allele in alleles.allele_accn:
        allele_scores[allele] = set()

    # find best R1/R2 scores
    for qname in read_scores:
        logger.write('%s, matches: %s' % (qname, len(read_scores[qname])))
        best_R1 = 0
        best_R2 = 0
        for allele in read_scores[qname]:
            if 'R1' in read_scores[qname][allele]:
                if read_scores[qname][allele]['R1'][0] > best_R1:
                    best_R1 = read_scores[qname][allele]['R1'][0]
            if 'R2' in read_scores[qname][allele]:
                if read_scores[qname][allele]['R2'][0] > best_R2:
                    best_R2 = read_scores[qname][allele]['R2'][0]

        for allele in read_scores[qname]:
            if 'R1' not in read_scores[qname][allele] or 'R2' not in read_scores[qname][allele]:
                # R1 and R2 must both align
                continue

            if read_scores[qname][allele]['R1'][1] == read_scores[qname][allele]['R2'][1]:
                # must be in reverse orientation
                continue

            if read_scores[qname][allele]['R1'][0] == best_R1 and read_scores[qname][allele]['R2'][0] == best_R2:
                # both sides must be the best match for R1/R2
                allele_scores[allele].add(qname)

    # for gene in alleles.genes:
    #     for allele in alleles.genes[gene]:
    #         logger.write(allele)
    #         if allele in read_scores[qname]:
    #             if 'R1' not in read_scores[qname][allele] or 'R2' not in read_scores[qname][allele]:
    #                 # R1 and R2 must both align
    #                 continue

    #             if read_scores[qname][allele]['R1'][1] == read_scores[qname][allele]['R2'][1]:
    #                 # must be in reverse orientation
    #                 continue

    #             if read_scores[qname][allele]['R1'][0] == best_R1 and read_scores[qname][allele]['R2'][0] == best_R2:
    #                 # both sides must be the best match for R1/R2
    #                 allele_scores[allele].add(qname)

    logger.write("Done\n")
    sys.stdout.write('gene\tallele1\tallele2\tcount1\tcount2\tboth\ttotal\n')

    for gene in alleles.genes:
        sys.stderr.write("Calculating pairs for HLA-%s\n" % gene)
        # this is done to give a sequential order to the alleles (all A*01:01 are next to each other)
        allele_list = alleles.sorted_alleles(gene=gene)

        best_total = 0
        best_alleles = None

        for i, allele1 in enumerate(allele_list):
            for allele2 in allele_list[i:]:
                one = 0
                two = 0
                both = 0

                for q1 in allele_scores[allele1]:
                    if allele2 in allele_scores and q1 in allele_scores[allele2]:
                        both += 1
                    else:
                        one += 1
                if allele2 in allele_scores:
                    for q2 in allele_scores[allele2]:
                        if allele1 in allele_scores and not q2 in allele_scores[allele1]:
                            two += 1

                logger.write("%s/%s (%s, %s, %s)" % (allele1, allele2, one, two, both))

                total = one + two + both

                if total > best_total:
                    best_total = total
                    best_alleles = '%s/%s (%s/%s) => (%s, %s, %s)' % (allele1, allele2, alleles.allele_accn[allele1], alleles.allele_accn[allele2], one, two, both)

                sys.stdout.write('%s\n' % '\t'.join([str(x) for x in [gene ,allele1, allele2, one, two, both, total]]))
        logger.write("    best: %s %s\n" % (best_alleles, best_total))



def predict(scores, hla, motifs, thres=0.95):
    alleles = cgmhc.HLAalleles(hla)

    # for an allele, what is the motif
    motifs_allele = {}
    # for a motif, what are the alleles (first allele is the canonical)
    rev_motifs = {}

    with cgmhc.open_file(motifs, 'rt') as f:
        for line in f:
            cols = line.strip('\n').split(',')
            if cols[0] == 'allele':
                continue

            allele = cols[0][4:] # remove "HLA-"
            motif = cols[1]

            motifs_allele[allele] = motif
            if not motif in rev_motifs:
                rev_motifs[motif] = []
            rev_motifs[motif].append(allele)

    best_scores = {} # keyed by gene
    best_gene_scores = {} # keyed by gene

    with cgmhc.open_file(scores, 'rt') as f:
        for line in f:
            cols = line.strip('\n').split('\t')

            if cols[0] == 'gene':
                continue
                
            gene = cols[0]
            allele1 = cols[1]
            allele2 = cols[2]
            one = int(cols[3])
            two = int(cols[4])
            both = int(cols[5])
            total = one + two + both

            if not gene in best_scores:
                best_scores[gene] = []
                best_gene_scores[gene] = 0

            if total < best_gene_scores[gene] * thres:
                continue

            allele1_4digit = cgmhc.HLAalleles.to_4digit(alleles.allele_accn[allele1])
            allele2_4digit = cgmhc.HLAalleles.to_4digit(alleles.allele_accn[allele2])

            if not allele1_4digit in motifs_allele or not allele2_4digit in motifs_allele:
                # we will only return potential alleles where we have binding motifs
                continue

            best_scores[gene].append((total, both, one, two, allele1, allele2))
            best_scores[gene] = sorted(best_scores[gene], reverse=True)
            best_gene_scores[gene] = best_scores[gene][0][0]
            tmp = []

            for tup in best_scores[gene]:
                if tup[0] >= best_gene_scores[gene] * thres:
                    tmp.append(tup)
            best_scores[gene] = tmp

    sys.stdout.write('gene\tallele\tallele_4digit\tpartner\tpartner_4digit\tcount\tshared\tpartner_count\tclass\tpseudo_motif\tcanonical_motif_allele\n')
    for gene in best_gene_scores:
        primary = {}
        secondary = {}

        for total, both, one, two, allele1, allele2 in best_scores[gene]:
            allele1_4digit = cgmhc.HLAalleles.to_4digit(alleles.allele_accn[allele1])
            allele2_4digit = cgmhc.HLAalleles.to_4digit(alleles.allele_accn[allele2])

            allele1_motif = motifs_allele[allele1_4digit]
            allele2_motif = motifs_allele[allele2_4digit]

            if total == best_gene_scores[gene]:
                if not allele1_motif in primary:
                    primary[allele1_motif] = (allele1, allele1_4digit, allele2, allele2_4digit, one, both, two)
                if not allele2_motif in primary:
                    primary[allele2_motif] = (allele2, allele2_4digit, allele1, allele1_4digit, two, both, one)
            else:
                if not allele1_motif in primary and not allele1_motif in secondary:
                    secondary[allele1_motif] = (allele1, allele1_4digit, allele2, allele2_4digit, one, both, two)
                if not allele2_motif in primary and not allele2_motif in secondary:
                    secondary[allele2_motif] = (allele2, allele2_4digit, allele1, allele1_4digit, two, both, one)

        for allele_motif in primary:
            allele, allele4, partner, partner4, me, both, other = primary[allele_motif]
            sys.stdout.write('HLA-%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tPRIMARY\t%s\t%s\n' % (gene, allele, allele4, partner, partner4, me, both, other, allele_motif, rev_motifs[allele_motif][0]))
        for allele_motif in secondary:
            allele, allele4, partner, partner4, me, both, other = secondary[allele_motif]
            sys.stdout.write('HLA-%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tSECONDARY\t%s\t%s\n' % (gene, allele, allele4, partner, partner4, me, both, other, allele_motif, rev_motifs[allele_motif][0]))

#            sys.stdout.write('%s\n' % '\t'.join([str(x) for x in [gene, allele1, allele2, one, two, both, total, 'PRIMARY' if total == best_gene_scores[gene] else 'SECONDARY']]))
