from collections import namedtuple

SamReadFlags = namedtuple("SamReadFlags", "paired properly_aligned unmapped next_unmapped revcomp next_revcomp first last secondary qcfail pcrdup supplementary")

def parse_flags(flags: int) -> SamReadFlags:
    paired = flags & 0x1 > 0
    proper = flags & 0x2 > 0
    unmapped = flags & 0x4 > 0
    next_unmapped = flags & 0x8 > 0
    
    revcomp = flags & 0x10 > 0
    next_revcomp = flags & 0x20 > 0
    first = flags & 0x40 > 0
    last = flags & 0x80 > 0
    
    secondary = flags & 0x100 > 0
    qcfail = flags & 0x200 > 0
    pcrdup = flags & 0x400 > 0
    supplementary = flags & 0x800 > 0

    return SamReadFlags(paired, proper, unmapped, next_unmapped, revcomp, next_revcomp, first, last, secondary, qcfail, pcrdup, supplementary)


def cigar_to_reflen(cigar: str):
    refpos = 0
    for size, op in parse_cigar_op(cigar):
        if op in ['M', '=', 'X', 'D', 'N']:
            refpos += size
    
    return refpos

def cigar_to_readlen(cigar: str):
    readpos = 0
    for size, op in parse_cigar_op(cigar):
        if op in ['M', '=', 'X', 'I']:
            readpos += size
    
    return readpos

def parse_cigar_op(cigar: str) -> tuple[int, str]:
    b=''
    while cigar:
        if cigar[0] in '0123456789':
            b += cigar[0]
            cigar = cigar[1:]
        else:
            size = int(b)
            b = ''
            op = cigar[0]
            cigar = cigar[1:]
            yield (size, op)
