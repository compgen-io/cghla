import sys
import gzip
import subprocess

import cgmhc
import cgmhc.samfile



class ReadWriter(object):
    def __init__(self, out):
        self.r1_out = gzip.open('%s_R1.fastq.gz' % out, 'wt')
        self.r2_out = gzip.open('%s_R2.fastq.gz' % out, 'wt')    

        self.read1 = {}
        self.read2 = {}

        self.written = set()


    def close(self):
        self.r1_out.close()
        self.r2_out.close()


    def add_read1(self, name, seq, qual):
        if name in self.written:
            return

        self.read1[name] = (seq, qual)
        if name in self.read2:
            self.write(name)


    def add_read2(self, name, seq, qual):
        if name in self.written:
            return

        self.read2[name] = (seq, qual)
        if name in self.read1:
            self.write(name)


    def write(self, name):
        self.r1_out.write('@%s\n%s\n+\n%s\n' % (name, self.read1[name][0], self.read1[name][1]))
        self.r2_out.write('@%s\n%s\n+\n%s\n' % (name, self.read2[name][0], self.read2[name][1]))

        self.written.add(name)
        del self.read1[name]
        del self.read2[name]


def extract_reads(bam, bed, out, unmapped=None):
    writer = ReadWriter(out)

    sys.stderr.write("Extracting reads from: %s\n" % bam)
    proc = subprocess.Popen(['ngsutilsj', 'bam-extract', '--paired', '--sam', '--bed', bed, bam, '-'], bufsize=1, text=True, stdout=subprocess.PIPE, stderr=None)

    for line in iter(proc.stdout.readline,''):
        if not line.strip() or line[0] == '@':
            continue
        cols = line.strip('\n').split('\t')
        qname = cols[0]
        flags = cgmhc.samfile.parse_flags(int(cols[1]))
        seq = cols[9]
        qual = cols[10]

        if flags.supplementary:
            # remove supplementary reads (they have truncated seqs)
            continue

        if flags.first:
            if flags.revcomp:
                writer.add_read1(qname, cgmhc.revcomp(seq), qual[::-1])
            else:
                writer.add_read1(qname, seq, qual)
        else:
            if flags.revcomp:
                writer.add_read2(qname, cgmhc.revcomp(seq), qual[::-1])
            else:
                writer.add_read2(qname, seq, qual)

    proc.stdout.close()
    proc.wait()
    sys.stderr.write("Done\n")

    if unmapped:
        sys.stderr.write("Extracting unmapped reads from: %s\n" % unmapped)
        proc = subprocess.Popen(['samtools', 'view', unmapped], bufsize=1, text=True, stdout=subprocess.PIPE, stderr=None)

        for line in iter(proc.stdout.readline,''):
            cols = line.strip('\n').split('\t')
            qname = cols[0]
            flags = cgmhc.samfile.parse_flags(int(cols[1]))
            seq = cols[9]
            qual = cols[10]

            if flags.supplementary:
                # remove supplementary reads (they have truncated seqs)
                continue

            if flags.first:
                if flags.revcomp:
                    writer.add_read1(qname, cgmhc.revcomp(seq), qual[::-1])
                else:
                    writer.add_read1(qname, seq, qual)
            else:
                if flags.revcomp:
                    writer.add_read2(qname, cgmhc.revcomp(seq), qual[::-1])
                else:
                    writer.add_read2(qname, seq, qual)

        proc.stdout.close()
        proc.wait()
        sys.stderr.write("Done\n")

    writer.close()

