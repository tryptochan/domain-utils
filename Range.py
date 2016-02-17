"""This module contains Range classes
that parse multiple formats of range definition strings.
and facilitate the extraction of residues from
PDB or FASTA.
"""
import re
from . import ParseError

verbose = 0
debug = 0


class Range(object):
    """Range class is designed to provide conveninet interfaces
    to parse domain definitions of different fomats.

    Domain boundaries or the "range" can be divided into
    one or more contiguous fragments in sequence.

    Each contiguous fragment is modelled by ContiguousRange objects.

    Note that Range is a abstract class provides ideas about interfaces
    provided to the specific type of dervied classes of Range
    like SequenceRange or StructureRange.

    """

    def __init__(self):
        self.contiguous_ranges = []
        self.parsing_function_list = []

    def __getitem__(self, i):
        return self.contiguous_ranges[i]

    def __iter__(self):
        return self.contiguous_ranges.__iter__()

    def __str__(self):
        return ','.join([str(fragment) for fragment in self.contiguous_ranges])

    def __len__(self):
        #return len(self.contiguous_ranges)
        raise NotImplementedError

    def get_breaking_positions(self):
        """This function will return a list of breaking points,
        i.e. from range "1-100,101-220,221-300,310-400"
        the returned list will be
        [101,221].

        Note that neither 300 nor 310 will not be marked as breaking
        points becase there are insertions between 300 to 310.

        """
        breaking_points = []
        for r, l in zip(self[:-1], self[1:]):
            if verbose:
                print 'checking breaking points'
                print r, l
                print r.end, l.start

            if r.end + 1 == l.start:
                breaking_points.append(l.start)
                if verbose:
                    print l.start, "added!"

        return breaking_points

    def parse(self, string_of_certain_range_format):
        """Try to parse the given range input format
        using the list of parsing functions defined in parsing_function_list.

        Note that the parsing function list has priority in parsing,
        if the earlier one is successful in parsing, the later functions will
        not be even tried.

        """
        for parse_fn in self.parsing_function_list:
            try:
                parse_fn(string_of_certain_range_format)
            except RangeParseError:
                #if error happens check the next function
                continue
            else:
                #if no error happens
                #accept the result.
                break
        else:
            #visiting this block means that the
            #parsing function list is empty or
            #all functions in the funtion list
            #failed to parse the given input range format.

            raise RangeParseError('All tried parsing methods failed; %s'
                                  % str(self.parsing_function_list))


class ContiguousRange(object):
    """Pretty much like slice object defined in standard library.

    Abstract class.

    """

    def __str__(self):
        raise NotImplementedError

    def __len__(self):
        raise NotImplementedError

    def parse(self, start_stop_string):
        raise NotImplementedError

    def get_start(self):
        raise NotImplementedError

    def get_end(self):
        raise NotImplementedError


class SequenceContiguousRange(ContiguousRange):
    """ContiguousRange for Sequences.

    Start and end class variable can be used for the start and end positions
    for the contiguous region. And has the optional chain ID to use with
    multiple chain sequence domain ranges.

    """
    def __init__(self, start=None, end=None, chain=None):
        if start is not None:
            self.start = int(start)
        else:
            self.start = None
        if end is not None:
            self.end = int(end)
        else:
            self.end = None
        self.chain_id = chain

    def __str__(self):
        if not self.chain_id:
            return '%s-%s' % (self.start, self.end)
        else:
            return '%s:%s-%s' % (self.chain_id, self.start, self.end)

    def parse(self, range_string):
        if not range_string:
            return

        if ':' in range_string:  # still assumes use : as separator
            c = range_string.index(':')
            self.chain_id = range_string[c - 1]
            range_string = range_string[c + 1:]
        try:
            i, j = range_string.split('-')
            self.start, self.end = int(i), int(j)
        except:
            raise RangeParseError('Error in parsing range string')

    def get_start(self):
        return self.start

    def get_end(self):
        return self.end

    def get_overlap(self, contiguous_range):
        """Get overlap bewteen SequenceContiguousRanges."""
        if (self.chain_id is not None and contiguous_range.chain_id is not
                None and self.chain_id != contiguous_range.chain_id):
            return None

        i = self.get_start()
        j = self.get_end()

        n = contiguous_range.get_start()
        m = contiguous_range.get_end()

        x = max(i, n)  # max of starts
        y = min(j, m)  # min of ends

        if x > y:
            return None  # no overlap found!

        new_contig = SequenceContiguousRange()
        new_contig.start = x
        new_contig.end = y
        new_contig.chain_id = self.chain_id

        return new_contig

    def __len__(self):
        if self.get_end() is None or self.get_start() is None:
            return 0
        else:
            return self.get_end() - self.get_start() + 1


class SequenceRange(Range):
    """Range class that deal with Sequence."""
    def __init__(self, range_str=None):
        Range.__init__(self)

        ############################################
        #the function should be added into the list
        #for parsing!!
        #further parsing functions should be added very simply.
        self.parsing_function_list.append(self.parse_scoplike_sequence_format)
        if range_str:
            self.parse(range_str)

    def __len__(self):
        return reduce(lambda x, y: x + y,
                      [len(i) for i in self.contiguous_ranges], 0)

    def add_contiguous_range(self, contig):
        self.contiguous_ranges.append(contig)

    def get_start(self):
        return min([contig.get_start() for contig in self.contiguous_ranges])

    def get_end(self):
        return max([contig.get_end() for contig in self.contiguous_ranges])

    def is_sorted(self):
        """Returns True if the fragement is in sorted order in the sequence.
        Otherwise, returns False.

        Note that this function does not check if the ranges are overlapping or not.

        """
        for contig1, contig2 in zip(self[:-1], self[1:]):
            if contig1.chain_id != contig2.chain_id:
                continue
            if contig1.get_start() < contig2.get_start() and contig1.get_end() < contig2.get_end():
                pass
            else:
                return False
        else:
            return True


    def get_overlap(self, seqrange):
        """Get the intersection between two SequenceRange."""
        # TODO: extend to mutiple ranges
        ol = []
        for contig1 in self.contiguous_ranges:
            for contig2 in seqrange.contiguous_ranges:
                cross = contig1.get_overlap(contig2)
                if cross is not None:
                    ol.append(cross)
        if ol:
            new_r = SequenceRange()
            for c in ol:
                new_r.add_contiguous_range(c)
            return new_r
        else:
            return None

    def get_union(self, *args):
        """Get union of SequenceRanges."""
        contigs = []
        for c1 in self.contiguous_ranges:
            contigs.append(c1)
        for seq_range in args:
            for c2 in seq_range.contiguous_ranges:
                contigs.append(c2)

        def _union2contigs(ctg1, ctg2):
            if ctg1.chain_id != ctg2.chain_id:
                return [ctg1, ctg2]
            if min(ctg1.end, ctg2.end) < max(ctg1.start, ctg2.start) - 1:
                return [ctg1, ctg2]
            contig = SequenceContiguousRange()
            contig.chain_id = ctg1.chain_id
            contig.end = max(ctg1.end, ctg2.end)
            contig.start = min(ctg1.start, ctg2.start)
            return [contig]

        union = SequenceRange()
        try:
            merged = [contigs[0]]
        except IndexError:
            return union
        for ctg in contigs[1:]:
            for m in merged:
                tmp = _union2contigs(ctg, m)
                if len(tmp) == 1:
                    merged.remove(m)
                    merged.extend(tmp)
                    break
            if len(tmp) == 2:  # not merged
                merged.append(ctg)
        for ctg in sorted(merged, key=lambda i: i.start):
            union.add_contiguous_range(ctg)
        return union

    def parse_scoplike_sequence_format(self, range_string):
        """Parses SCOP like sequence range format parser.

        Example:
        "1-100,112-200"

        The above example means the sequence boundary is defined from
        position 1 to 100 and then start from 112 to 200.
        Each contiguous fragment is divided by ','.

        It supports for ranges like '344,355-380', also with chain IDs like
        'A:1-48, B:56-100'.

        """
        contiguous_divider = ','
        chain_divider = ':'
        for contiguous_string in range_string.split(contiguous_divider):
            contig = SequenceContiguousRange()
            contiguous_string = contiguous_string.strip()
            if chain_divider in contiguous_string:
                c = contiguous_string.index(chain_divider)
                contig.chain_id = contiguous_string[c - 1]
                contiguous_string = contiguous_string[c + 1:]
            try:
                if contiguous_string.isdigit():
                    contig.start = contig.end = int(contiguous_string)
                else:
                    start, end = contiguous_string.split('-')
                    contig.start = int(start)
                    contig.end = int(end)

                self.add_contiguous_range(contig)
            except:
                raise RangeParseError(
                    "Error in parsing simple scoplike format")


class PDBContiguousRange(ContiguousRange):
    """ContiguousRange for Sequences.

    Start and end class variable can be used for the start and end positions
    for the contiguous region.

    One difference from SequenceRange is length cannot be determined due to icode.

    """

    def __init__(self):
        self.chain_id = None

        self.start_resnum = None  # stores char
        self.start_icode = ''

        self.end_resnum = None
        self.end_icode = ''

    def __str__(self):
        s = ''
        if self.chain_id != ' ' and self.chain_id:
            s = '%s:' % self.chain_id

        if (self.start_resnum is None and self.end_resnum is None):
            return s
        else:
            if self.start_icode == ' ':
                start_i = ''
            else:
                start_i = self.start_icode
            if self.end_icode == ' ':
                end_i = ''
            else:
                end_i = self.end_icode

            return '%s%s%s-%s%s' % (s, self.start_resnum,
                                    start_i, self.end_resnum, end_i)

    def get_start(self):
        """Returns start residue id.

        if start_residue_id was not defined, returns None

        """

        if self.start_resnum is not None:
            if self.start_icode:
                return (self.chain_id, int(self.start_resnum), self.start_icode)
            else:
                start_icode = ' '
                return (self.chain_id, int(self.start_resnum), start_icode)
        else:
                return None

    def get_end(self):
        """Returns end residue id. if end was not defined, returns None."""

        if self.end_resnum is not None:
            if self.end_icode:
                return (self.chain_id, int(self.end_resnum), self.end_icode)
            else:
                end_icode = ' '
                return (self.chain_id, int(self.end_resnum), end_icode)
        else:
            return None


class PDBRange(Range):
    """Range class that deals with structure Range."""
    def __init__(self, range_str=None):
        Range.__init__(self)
        #the function should be added into the list for parsing!!
        #further parsing functions should be added very simply.
        self.parsing_function_list.append(self.parse_scoprange_format)
        if range_str:
            self.parse(range_str)

    def get_unique_chain_ids(self):
        uniq_ids = []
        for contig in self:
            if not contig.chain_id in uniq_ids:
                uniq_ids.append(contig.chain_id)
        return uniq_ids

    def add_contiguous_range(self, contig):
        self.contiguous_ranges.append(contig)

    def parse_scoprange_format(self, range_string):
        """String format parser for the scop range."""
        #full range starting with chain id
        scoprange_re1 = r'^(.):(-?[0-9]+)([a-zA-Z]?)-(-?[0-9]+)([a-zA-Z]?)$'
        #full range without chain id
        scoprange_re2 = r'^(-?[0-9]+)([a-zA-Z]?)-(-?[0-9]+)([a-zA-Z]?)$'
        #full chain with chainid
        scoprange_re3 = r'^(.):$'
        scoprange_re4 = r'^-$'
        #with single digit range
        scoprange_re5 = r'^(.):(-?[0-9]+)([a-zA-Z]?)$'
        scoprange_re6 = r'^(-?[0-9]+)([a-zA-Z]?)$'

        #first split into fragments using comma, ,
        fragments = range_string.strip().split(',')

        for s in fragments:
            contig = PDBContiguousRange()
            #matching template #1
            match = re.match(scoprange_re1, s)
            if match:
                (contig.chain_id, contig.start_resnum, contig.start_icode,
                    contig.end_resnum, contig.end_icode) = match.groups()
                self.add_contiguous_range(contig)
                continue

            #matching template #2
            match = re.match(scoprange_re2, s)
            if match:
                #no chain id need to be set!
                (contig.start_resnum, contig.start_icode,
                    contig.end_resnum, contig.end_icode) = match.groups()
                self.add_contiguous_range(contig)
                continue

            match = re.match(scoprange_re3, s)
            if match:
                contig.chain_id = match.group(1)
                #no start and end residue information is needed!
                self.add_contiguous_range(contig)
                continue

            match = re.match(scoprange_re4, s)
            if match:
                if len(fragments) == 1:
                    #nothing to be done!
                    pass
                else:
                    raise RangeParseError("SCOP range re4 has problem!")

                self.add_contiguous_range(contig)
                continue

            match = re.match(scoprange_re5, s)
            if match:
                (contig.chain_id, contig.start_resnum,
                 contig.start_icode) = match.groups()
                contig.end_resnum = contig.start_resnum
                contig.end_icode = contig.start_icode
                self.add_contiguous_range(contig)
                continue

            match = re.match(scoprange_re6, s)
            if match:
                (contig.start_resnum, contig.start_icode) = match.groups()
                contig.end_resnum = contig.start_resnum
                contig.end_icode = contig.start_icode
                self.add_contiguous_range(contig)
                continue

            raise RangeParseError(
                "SCOP range parsing has problem!", range_string)


def trim_seq_range(seq_range, length=5):
    """Returns a shallow copy of SequenceRange without small contigs.

    It is useful for structual alignments which usually have some fragments.
    It's better to remove them before converting aligned postions back to PDB
    range or anything.

    """
    if not isinstance(seq_range, SequenceRange):
        return None
    new_range = SequenceRange()
    for contig in seq_range:
        if len(contig) > length:
            new_range.add_contiguous_range(contig)
    return new_range

def aln2range(query, hit, query_start=1, hit_start=1, gap_tolerance=20, case_sensitive=False):
    """Convert a alignment with gaps to sequence ranges.

    Consider paired upper cases as aligned positions.
    Returned sequence range is 1 based.

    case_sensitive: if lower case letters are considered aligned.
    """
    assert len(query) == len(hit)
    qseqr = SequenceRange()
    hseqr = SequenceRange()
    qcontig_start = None
    hcontig_start = None
    qlet_cnt = 0
    hlet_cnt = 0
    qlast_pos = 0
    hlast_pos = 0
    for aaq, aah in zip(query, hit):
        if aaq.isalpha():
            qlet_cnt += 1
        if aah.isalpha():
            hlet_cnt += 1
        if (aaq.isalpha() and aah.isalpha() and (case_sensitive and aaq.isupper() and
                aah.isupper() or not case_sensitive)):
            if qlet_cnt - qlast_pos - 1 > gap_tolerance:
                contig = SequenceContiguousRange(query_start+qcontig_start-1, query_start+qlast_pos-1)
                qseqr.add_contiguous_range(contig)
                qcontig_start = None
            if hlet_cnt - hlast_pos - 1 > gap_tolerance:
                contig = SequenceContiguousRange(hit_start+hcontig_start-1, hit_start+hlast_pos-1)
                hseqr.add_contiguous_range(contig)
                hcontig_start = None
            if not qcontig_start:
                qcontig_start = qlet_cnt
            if not hcontig_start:
                hcontig_start = hlet_cnt
            qlast_pos = qlet_cnt
            hlast_pos = hlet_cnt

    if qcontig_start:
        contig = SequenceContiguousRange(query_start+qcontig_start-1, query_start+qlast_pos-1)
        qseqr.add_contiguous_range(contig)
    if hcontig_start:
        contig = SequenceContiguousRange(hit_start+hcontig_start-1, hit_start+hlast_pos-1)
        hseqr.add_contiguous_range(contig)
    return (qseqr, hseqr)

def map_aligned_range(seq_range, query, hit, query_start=1, hit_start=1, case_sensitive=False, **kwargs):
    """Map corresponding range from query to hit, either represented by aligned
    sequences or equivalent indexes.

    query, hit: string or list of indexes
    start indexes are ignored when indexes are given directly
    case_sensitive: if lower case letters are considered aligned.
    """
    assert len(query) == len(hit)
    q_res = set()
    for contig in seq_range:
        q_res = q_res.union(range(contig.start-query_start+1, contig.end-query_start+2))

    h_res = []
    if isinstance(query, basestring) and isinstance(hit, basestring):
        qlet_cnt = 0
        hlet_cnt = 0
        for aaq, aah in zip(query, hit):
            if aaq.isalpha():
                qlet_cnt += 1
            if aah.isalpha():
                hlet_cnt += 1
            if (aaq.isalpha() and aah.isalpha() and (case_sensitive and aaq.isupper() and
                    aah.isupper() or not case_sensitive)):
                if qlet_cnt in q_res:
                    h_res.append(hlet_cnt+hit_start-1)
    elif isinstance(query, list) and isinstance(hit, list):
        for iq, ih in zip(query, hit):
            if iq in q_res:
                h_res.append(ih)
    else:
        raise TypeError('Query or hit should be string of sequences or list of indexes')
    return index2range(h_res, **kwargs)

def index2range(index, gap_tolerance=20, **kwargs):
    """Convert aligned index of equivalent positions to sequence range."""
    if not index:
        return None
    seq_range = SequenceRange()
    start = int(index[0])
    end = start
    for i in index:
        i = int(i)
        if i - end - 1 <= gap_tolerance:
            end = i
        else:
            contig = SequenceContiguousRange(start, end)
            seq_range.add_contiguous_range(contig)
            start = i
            end = i
    contig = SequenceContiguousRange(start, end)
    seq_range.add_contiguous_range(contig)
    return seq_range


class RangeParseError(ParseError):
    """This error means that all the available parsers in the
    Range class failed to correctly parse the given range format.
    """
    pass


class SequenceRangeError(ParseError):
    """This error means that the SequenceRange is not right."""
    pass
