import sys
import gzip
import subprocess
import tempfile

import cghla
import cghla.samfile


class ReadWriter(object):
    def __init__(self, out, read_limit=1000000):
        self.r1_out = gzip.open('%s_R1.fastq.gz' % out, 'wt')
        self.r2_out = gzip.open('%s_R2.fastq.gz' % out, 'wt')    

        self.read1 = {}
        self.read2 = {}

        self.written = set()

        self.read_limit = read_limit

        self._r1_tmpfiles = []
        self._r2_tmpfiles = []


    def close(self):
        buf1 = {}
        buf2 = {}

        for i, f1 in enumerate(self._r1_tmpfiles):
            f1.seek(0)
            name, seq, qual = self._read_tmp_rec(f1)
            buf1[name] = (i, seq, qual)
            
        for i, f2 in enumerate(self._r2_tmpfiles):
            f2.seek(0)
            name, seq, qual = self._read_tmp_rec(f2)
            buf2[name] = (i, seq, qual)

        while buf1 and buf2:
            found = False
            for read in buf1:
                if read in buf2:
                    found = True

                    self.r1_out.write('@%s\n%s\n+\n%s\n' % (name, buf1[name][1], buf1[name][2]))
                    self.r2_out.write('@%s\n%s\n+\n%s\n' % (name, buf2[name][1], buf2[name][2]))

                    read1 = self._read_tmp_rec(self._r1_tmpfiles[buf1[name][0]])
                    if read1:
                        buf1[read1[0]] = (buf1[name][0], read1[1], read1[2])

                    read2 = self._read_tmp_rec(self._r2_tmpfiles[buf2[name][0]])
                    if read2:
                        buf2[read2[0]] = (buf2[name][0], read2[1], read2[2])

                    del buf1[name]
                    del buf2[name]

            if not found:
                # we didn't find a match, so remove the first (sorted) element from buf1 or buf2
                # this should have matched, so if it didn't, we can assume it is unpaired.

                name1 = sorted(buf1)[0]
                name2 = sorted(buf2)[0]

                # sys.stderr.write('name1: %s => %s\n' % (name1, buf1[name1]))
                # sys.stderr.write('name2: %s => %s\n' % (name2, buf2[name2]))

                if name1 < name2:
                    read1 = self._read_tmp_rec(self._r1_tmpfiles[buf1[name1][0]])
                    if read1:
                        buf1[read1[0]] = (buf1[name1][0], read1[1], read1[2])
                    del(buf1[name1])
                else:
                    read2 = self._read_tmp_rec(self._r2_tmpfiles[buf2[name2][0]])
                    if read2:
                        buf2[read2[0]] = (buf2[name2][0], read2[1], read2[2])
                    del(buf2[name2])


        for f1 in self._r1_tmpfiles:
            f1.close()

        for f2 in self._r2_tmpfiles:
            f2.close()

        self.r1_out.close()
        self.r2_out.close()


    def _read_tmp_rec(self, fobj):
        name = fobj.readline().strip()
        seq = fobj.readline().strip()
        qual = fobj.readline().strip()

        if not name:
            return None

        return (name, seq, qual)


    def trim(self):
        if len(self.read1) > self.read_limit or len(self.read2) > self.read_limit:
            tmpfile1 = tempfile.NamedTemporaryFile('w+t')
            tmpfile2 = tempfile.NamedTemporaryFile('w+t')
            sys.stderr.write("Reads already written: %s, dumping reads in memory to tmpfiles: %s %s\n" % (len(self.written), tmpfile1.name, tmpfile2.name))

            for name in sorted(self.read1):
                tmpfile1.write('%s\n%s\n%s\n' % (name, self.read1[name][0], self.read1[name][1]))

            for name in sorted(self.read2):
                tmpfile2.write('%s\n%s\n%s\n' % (name, self.read2[name][0], self.read2[name][1]))

            self.read1 = {}
            self.read2 = {}

            self._r1_tmpfiles.append(tmpfile1)
            self._r2_tmpfiles.append(tmpfile2)

    def has_pair(self, qname):
        if qname in self.written:
            return False
        if qname in self.read1:
            return True
        if qname in self.read2:
            return True

        return False


    def add_read1(self, name, seq, qual):
        if name in self.written:
            return

        self.read1[name] = (seq, qual)
        if name in self.read2:
            self.write(name)
        else:
            self.trim()


    def add_read2(self, name, seq, qual):
        if name in self.written:
            return

        self.read2[name] = (seq, qual)
        if name in self.read1:
            self.write(name)
        else:
            self.trim()


    def write(self, name):
        self.r1_out.write('@%s\n%s\n+\n%s\n' % (name, self.read1[name][0], self.read1[name][1]))
        self.r2_out.write('@%s\n%s\n+\n%s\n' % (name, self.read2[name][0], self.read2[name][1]))

        self.written.add(name)
        del self.read1[name]
        del self.read2[name]


def extract_reads(bam, bed, out, bam2=None, unmapped=None):
    writer = ReadWriter(out)

    sys.stderr.write("Extracting reads from: %s\n" % bam)
    args = ['ngsutilsj', 'bam-extract', '--paired', '--sam', '--bed', bed]
    if bam2:
        sys.stderr.write("     extra reads from: %s\n" %  bam2)
        args.append('--bam2')
        args.append(bam2)

    args.append(bam)
    args.append('-')
    proc = subprocess.Popen(args, bufsize=1, text=True, stdout=subprocess.PIPE, stderr=None)

    for line in iter(proc.stdout.readline,''):
        if not line.strip() or line[0] == '@':
            continue
        cols = line.strip('\n').split('\t')
        qname = cols[0]
        flags = cghla.samfile.parse_flags(int(cols[1]))
        seq = cols[9]
        qual = cols[10]

        if flags.supplementary:
            # remove supplementary reads (they have truncated seqs)
            continue

        if flags.first:
            if flags.revcomp:
                writer.add_read1(qname, cghla.revcomp(seq), qual[::-1])
            else:
                writer.add_read1(qname, seq, qual)
        else:
            if flags.revcomp:
                writer.add_read2(qname, cghla.revcomp(seq), qual[::-1])
            else:
                writer.add_read2(qname, seq, qual)

    proc.stdout.close()
    proc.wait()
    sys.stderr.write("Done\n")

    # t = 0

    if unmapped:
        #
        # If we have a read that is unmapped, or otherwise failed filters but is withing a BED region,
        # we will pull those reads out too.
        #

        bed_regions = []
        with open(bed, 'rt') as f:
            for line in f:
                if not line.strip() or line[0] == '#':
                    continue

                cols = line.strip('\n').split('\t')
                bed_regions.append((cols[0], int(cols[1]), int(cols[2])))

        sys.stderr.write("Extracting remaining unmapped/failed reads from: %s\n" % unmapped)
        proc = subprocess.Popen(['samtools', 'view', unmapped], bufsize=1, text=True, stdout=subprocess.PIPE, stderr=None)

        in_region = 0
        is_unmapped = 0
        has_pair = 0

        for line in iter(proc.stdout.readline,''):
            # t += 1
            # if t > 10000000:
            #     sys.stderr.write("Skipping...\n")
            #     break

            cols = line.strip('\n').split('\t')
            qname = cols[0]
            flags = cghla.samfile.parse_flags(int(cols[1]))
            chrom = cols[2]
            pos = int(cols[3])
            next_chrom = cols[6]
            next_pos = int(cols[7])
            seq = cols[9]
            qual = cols[10]

            if next_chrom == '=':
                next_chrom = chrom

            if flags.supplementary:
                # remove supplementary reads (they have truncated seqs)
                continue

            good = False

            for bchrom, bstart, bend in bed_regions:
                # keep a read if it falls w/in the HLA BED regions
                if chrom == bchrom and pos > bstart and pos < bend:
                    # bed is zero based and SAM is one based, so > and <, not >=
                    good = True
                    in_region += 1
                elif next_chrom == bchrom and next_pos > bstart and next_pos < bend:
                    # also save a read if the pair matches
                    good = True
                    in_region += 1

            if flags.unmapped or flags.next_unmapped:
                # also keep a read if it (or the pair) is unmapped)
                good = True
                is_unmapped += 1

            if writer.has_pair(qname):
                # also keep a read if we have the R1/R2 pair from the original file (rare)
                good = True
                has_pair += 1

            if not good:
                continue

            if flags.first:
                if flags.revcomp:
                    writer.add_read1(qname, cghla.revcomp(seq), qual[::-1])
                else:
                    writer.add_read1(qname, seq, qual)
            else:
                if flags.revcomp:
                    writer.add_read2(qname, cghla.revcomp(seq), qual[::-1])
                else:
                    writer.add_read2(qname, seq, qual)

        proc.stdout.close()
        proc.wait()
        sys.stderr.write("Done\n")
        sys.stderr.write("Unmapped read criteria -- in bed region: %s, unmapped: %s, had pair %s\n" % (in_region, is_unmapped, has_pair))

    writer.close()

