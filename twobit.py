"""
This following code is a modified version of twobitreader (which is under Perl Artistic License 2.0).
As per license restrictions, the code below indicates what has been modified in relation to the
standard version (retrieved from https://bitbucket.org/thesylex/twobitreader on the 16 May 2012).
No warranty is provided, express or implied

Modifications to package:
- removed download.py and __main__ because they were not used and __main__ had errors.
- removed command-line interface because the BED file functionality is implemented more extensively elsewhere
"""
"""
twobitreader
Licensed under Perl Artistic License 2.0
No warranty is provided, express or implied
"""
from array import array
from bisect import bisect_right
from errno import ENOENT, EACCES
from os import R_OK, access

try:
    from os import strerror
except ImportError:
    strerror = lambda x: 'strerror not supported'
from os.path import exists, getsize
import logging
import textwrap
import sys

if sys.version_info > (3,):
    izip = zip
    xrange = range
    _CHAR_CODE = 'u'
    iteritems = dict.items
else:
    from itertools import izip

    _CHAR_CODE = 'c'
    iteritems = dict.iteritems


def safe_tostring(ary):
    """
    convert arrays to strings in a Python 2.x / 3.x safe way
    """
    if sys.version_info > (3,):
        return ary.tounicode().encode("ascii").decode()
    else:
        return ary.tostring()


def true_long_type():
    """
    OS X uses an 8-byte long, so make sure L (long) is the right size
    and switch to I (int) if needed
    """
    for type_ in ['L', 'I']:
        test_array = array(type_, [0])
        long_size = test_array.itemsize
        if long_size == 4:
            return type_
    raise ImportError("Couldn't determine a valid 4-byte long type to use \
as equivalent to LONG")


LONG = true_long_type()


def byte_to_bases(x):
    """convert one byte to the four bases it encodes"""
    c = (x >> 4) & 0xf
    f = x & 0xf
    cc = (c >> 2) & 0x3
    cf = c & 0x3
    fc = (f >> 2) & 0x3
    ff = f & 0x3
    return [bits_to_base(X) for X in (cc, cf, fc, ff)]


def bits_to_base(x):
    """convert integer representation of two bits to correct base"""
    if x is 0:
        return 'T'
    elif x is 1:
        return 'C'
    elif x is 2:
        return 'A'
    elif x is 3:
        return 'G'
    else:
        raise ValueError('Only integers 0-3 are valid inputs')


def base_to_bin(x):
    """
    provided for user convenience
    convert a nucleotide to its bit representation
    """
    if x == 'T':
        return '00'
    elif x == 'C':
        return '01'
    elif x == 'A':
        return '10'
    elif x == 'G':
        return '11'
    else:
        raise ValueError('Only characters \'ATGC\' are valid inputs')


def create_byte_table():
    """create BYTE_TABLE"""
    d = {}
    for x in xrange(2 ** 8):
        d[x] = byte_to_bases(x)
    return d


def split16(x):
    """
    split a 16-bit number into integer representation
    of its course and fine parts in binary representation
    """
    c = (x >> 8) & 0xff
    f = x & 0xff
    return c, f


def create_twobyte_table():
    """create TWOBYTE_TABLE"""
    d = {}
    for x in xrange(2 ** 16):
        c, f = split16(x)
        d[x] = list(byte_to_bases(c)) + list(byte_to_bases(f))
    return d


BYTE_TABLE = create_byte_table()
TWOBYTE_TABLE = create_twobyte_table()


def longs_to_char_array(longs, first_base_offset, last_base_offset, array_size,
                        more_bytes=None):
    """
    takes in an array of longs (4 bytes) and converts them to bases in
    a char array
    you must also provide the offset in the first and last block
    (note these offsets are pythonic. last_offset is not included)
    and the desired array_size
    If you have less than a long worth of bases at the end, you can provide
    them as a string with more_bytes=

    NOTE: last_base_offset is inside more_bytes not the last long, if more_bytes
          is not None
    returns the correct subset of the array based on provided offsets
    """
    if array_size == 0:
        return array(_CHAR_CODE)
    elif array_size < 0:
        raise ValueError('array_size must be at least 0')

    if not first_base_offset in range(16):
        raise ValueError('first_base_offset must be in range(16)')
    if not last_base_offset in range(1, 17):
        raise ValueError('last_base_offset must be in range(1, 17)')

    longs_len = len(longs)
    if more_bytes is None:
        shorts_length = 0
    else:
        shorts_length = len(more_bytes)
    if array_size > longs_len * 16 + 4 * shorts_length:
        raise ValueError('array_size exceeds maximum possible for input')

    dna = array(_CHAR_CODE, 'N' * (longs_len * 16 + 4 * shorts_length))
    # translate from 32-bit blocks to bytes
    # this method ensures correct endianess (byteswap as neeed)
    i = 0
    if longs_len > 0:
        bytes_ = array('B')
        bytes_.fromstring(longs.tostring())
        # first block
        first_block = ''.join([''.join(BYTE_TABLE[bytes_[x]]) for x in range(4)])
        i = 16 - first_base_offset
        if array_size < i:
            i = array_size
        dna[0:i] = array(_CHAR_CODE, first_block[first_base_offset:first_base_offset + i])
    if longs_len > 1:
        # middle blocks (implicitly skipped if they don't exist)
        for byte in bytes_[4:-4]:
            dna[i:i + 4] = array(_CHAR_CODE, BYTE_TABLE[byte])
            i += 4
        # last block
        last_block = array(_CHAR_CODE, ''.join([''.join(BYTE_TABLE[bytes_[x]])
                                                for x in range(-4, 0)]))
        if more_bytes is None:
            dna[i:i + last_base_offset] = last_block[0:last_base_offset]
        else:  # if there are more bytes, we need the whole last block
            dna[i:i + 16] = last_block[0:16]
        i += 16
    if more_bytes is not None:
        bytes_ = array('B')
        bytes_.fromstring(more_bytes)
        j = i
        for byte in bytes_:
            j = i + 4
            if j > array_size:
                dnabytes = array(_CHAR_CODE, BYTE_TABLE[byte])[0:(array_size - i)]
                dna[i:array_size] = dnabytes
                break
            dna[i:i + last_base_offset] = array(_CHAR_CODE, BYTE_TABLE[byte])
            i += 4
    return dna[0:array_size]


class TwoBitFile(dict):
    """
python-level reader for .2bit files (i.e., from UCSC genome browser)
(note: no writing support)
TwoBitFile inherits from dict
You may access sequences by name, e.g.
>>> genome = TwoBitFile('hg18.2bit')
>>> chr20 = genome['chr20']
Sequences are returned as TwoBitSequence objects
You may access intervals by slicing or using str() to dump the entire entry
e.g.
>>> chr20[100100:100120]
'ttttcctctaagataatttttgccttaaatactattttgttcaatactaagaagtaagataacttccttttgttggta
tttgcatgttaagtttttttcc'
>>> whole_chr20 = str(chr20)
Fair warning: dumping the entire chromosome requires a lot of memory
See TwoBitSequence for more info
    """

    def __init__(self, foo):
        super(TwoBitFile, self).__init__()
        if not exists(foo):
            raise IOError(ENOENT, strerror(ENOENT), foo)
        if not access(foo, R_OK):
            raise IOError(EACCES, strerror(EACCES), foo)
        self._filename = foo
        self._file_size = getsize(foo)
        self._file_handle = open(foo, 'rb')
        self._load_header()
        self._load_index()
        for name, offset in iteritems(self._offset_dict):
            self[name] = TwoBitSequence(self._file_handle, offset,
                                        self._file_size,
                                        self._byteswapped)
        return

    def __reduce__(self):  # enables pickling
        return (TwoBitFile, (self._filename,))

    def _load_header(self):
        file_handle = self._file_handle
        header = array(LONG)
        header.fromfile(file_handle, 4)
        # check signature -- must be 0x1A412743
        # if not, swap bytes
        byteswapped = False
        (signature, version, sequence_count, reserved) = header
        if not signature == 0x1A412743:
            byteswapped = True
            header.byteswap()
            (signature2, version, sequence_count, reserved) = header
            if not signature2 == 0x1A412743:
                raise TwoBitFileError('Signature in header should be ' +
                                      '0x1A412743, instead found 0x%X' %
                                      signature)
        if not version == 0:
            raise TwoBitFileError('File version in header should be 0.')
        if not reserved == 0:
            raise TwoBitFileError('Reserved field in header should be 0.')
        self._byteswapped = byteswapped
        self._sequence_count = sequence_count

    def _load_index(self):
        file_handle = self._file_handle
        byteswapped = self._byteswapped
        remaining = self._sequence_count
        sequence_offsets = []
        file_handle.seek(16)
        while remaining > 0:
            name_size = array('B')
            name_size.fromfile(file_handle, 1)
            if byteswapped:
                name_size.byteswap()
            # name = array(_CHAR_CODE)
            name = array('B')
            name.fromfile(file_handle, name_size[0])
            name = "".join([chr(X) for X in name])

            if byteswapped:
                name.byteswap()
            offset = array(LONG)
            offset.fromfile(file_handle, 1)
            if byteswapped:
                offset.byteswap()
            sequence_offsets.append((name, offset[0]))
            remaining -= 1
        self._sequence_offsets = sequence_offsets
        self._offset_dict = dict(sequence_offsets)

    def sequence_sizes(self):
        """returns a dictionary with the sizes of each sequence"""
        d = {}
        file_handle = self._file_handle
        byteswapped = self._byteswapped
        for name, offset in iteritems(self._offset_dict):
            file_handle.seek(offset)
            dna_size = array(LONG)
            dna_size.fromfile(file_handle, 1)
            if byteswapped:
                dna_size.byteswap()
            d[name] = dna_size[0]
        return d


class TwoBitSequence(object):
    """
A TwoBitSequence object refers to an entry in a TwoBitFile
You may access intervals by slicing or using str() to dump the entire entry
e.g.
>>> genome = TwoBitFile('hg18.2bit')
>>> chr20 = genome['chr20']
>>> chr20[100100:100200] # slicing returns a string
'ttttcctctaagataatttttgccttaaatactattttgttcaatactaagaagtaagataacttccttttgttggta
tttgcatgttaagtttttttcc'
>>> whole_chr20 = str(chr20) # get whole chr as string
Fair warning: dumping the entire chromosome requires a lot of memory
Note that we follow python/UCSC conventions:
Coordinates are 0-based, end-open
(Note: The UCSC web-based genome browser uses 1-based closed coordinates)
If you attempt to access a slice past the end of the sequence,
it will be truncated at the end.
Your computer probably doesn't have enough memory to load a whole genome
but if you want to string-ize your TwoBitFile, here's a recipe:
x = TwoBitFile('my.2bit')
d = x.dict()
for k,v in d.items(): d[k] = str(v)
    """

    def __init__(self, file_handle, offset, file_size, byteswapped=False):
        self._file_size = file_size
        self._file_handle = file_handle
        self._original_offset = offset
        self._byteswapped = byteswapped
        file_handle.seek(offset)
        header = array(LONG)
        header.fromfile(file_handle, 2)
        if byteswapped:
            header.byteswap()
        dna_size, n_block_count = header
        self._dna_size = dna_size  # number of characters, 2 bits each
        self._n_bytes = (dna_size + 3) / 4  # number of bytes
        # number of 32-bit fragments
        self._packed_dna_size = (dna_size + 15) / 16
        n_block_starts = array(LONG)
        n_block_sizes = array(LONG)
        n_block_starts.fromfile(file_handle, n_block_count)
        if byteswapped:
            n_block_starts.byteswap()
        n_block_sizes.fromfile(file_handle, n_block_count)
        if byteswapped:
            n_block_sizes.byteswap()
        self._n_block_starts = n_block_starts
        self._n_block_sizes = n_block_sizes
        mask_rawc = array(LONG)
        mask_rawc.fromfile(file_handle, 1)
        if byteswapped:
            mask_rawc.byteswap()
        mask_block_count = mask_rawc[0]
        mask_block_starts = array(LONG)
        mask_block_starts.fromfile(file_handle, mask_block_count)
        if byteswapped:
            mask_block_starts.byteswap()
        mask_block_sizes = array(LONG)
        mask_block_sizes.fromfile(file_handle, mask_block_count)
        if byteswapped:
            mask_block_sizes.byteswap()
        self._mask_block_starts = mask_block_starts
        self._mask_block_sizes = mask_block_sizes
        file_handle.read(4)
        self._offset = file_handle.tell()

    def __len__(self):
        return self._dna_size

    def __getitem__(self, slice_or_key):
        """
        return a sub-sequence, given a slice object
        """
        step = None
        if isinstance(slice_or_key, slice):
            step = slice_or_key.step
            if step is not None:
                raise ValueError("Slicing by step not currently supported")
            return self.get_slice(min_=slice_or_key.start, max_=slice_or_key.stop)

        elif isinstance(slice_or_key, int):
            max_ = slice_or_key + 1
            if max_ == 0:
                max_ = None
            return self.get_slice(min_=slice_or_key, max_=max_)

    def get_slice(self, min_, max_=None):
        """
        get_slice returns only a sub-sequence
        """
        # handle negative coordinates
        dna_size = self._dna_size
        if min_ is None:  # for slicing e.g. [:]
            min_ = 0
        if max_ is not None and max_ < 0:
            if max_ < -dna_size:
                raise IndexError('index out of range')
            max_ = dna_size + max_
        if min_ is not None and min_ < 0:
            if min_ < -dna_size:
                raise IndexError('index out of range')
            min_ = dna_size + min_
        # make sure there's a proper range
        if max_ is not None and min_ > max_:
            return ''
        if max_ == 0 or max_ == min_:
            return ''

        # load all the data
        if max_ is None or max_ > dna_size:
            max_ = dna_size

        file_handle = self._file_handle
        byteswapped = self._byteswapped
        n_block_starts = self._n_block_starts
        n_block_sizes = self._n_block_sizes
        mask_block_starts = self._mask_block_starts
        mask_block_sizes = self._mask_block_sizes
        offset = self._offset
        packed_dna_size = self._packed_dna_size
        # n_bytes = self._n_bytes

        # region_size is how many bases the region is
        if max_ is None:
            region_size = dna_size - min_
        else:
            region_size = max_ - min_

        # start_block, end_block are the first/last 32-bit blocks we need
        # blocks start at 0
        start_block = min_ // 16
        # jump directly to desired file location
        local_offset = offset + (start_block * 4)
        end_block = (max_ - 1 + 16) // 16
        # don't read past seq end

        file_handle.seek(local_offset)

        # note we won't actually read the last base
        # this is a python slice first_base_offset:16*blocks+last_base_offset
        first_base_offset = min_ % 16
        last_base_offset = max_ % 16
        if last_base_offset == 0:
            last_base_offset = 16
        # +1 we still need to read end_block maybe
        blocks_to_read = end_block - start_block
        if (blocks_to_read + start_block) > packed_dna_size:
            blocks_to_read = packed_dna_size - start_block

        fourbyte_dna = array(LONG)
        # remainder_seq = None
        if (blocks_to_read * 4 + local_offset) > self._file_size:
            fourbyte_dna.fromfile(file_handle, blocks_to_read - 1)
            morebytes = file_handle.read()  # read the remaining characters
        # if byteswapped:
        #                morebytes = ''.join(reversed(morebytes))
        else:
            fourbyte_dna.fromfile(file_handle, blocks_to_read)
            morebytes = None
        if byteswapped:
            fourbyte_dna.byteswap()
        str_as_array = longs_to_char_array(fourbyte_dna, first_base_offset,
                                           last_base_offset, region_size,
                                           more_bytes=morebytes)
        for start, size in izip(n_block_starts, n_block_sizes):
            end = start + size
            if end <= min_:
                continue
            if start > max_:
                break
            if start < min_:
                start = min_
            if end > max_:
                end = max_
            start -= min_
            end -= min_
            # this should actually be decoded, 00=N, 01=n
            str_as_array[start:end] = array(_CHAR_CODE, 'N' * (end - start))
        lower = str.lower
        first_masked_region = max(0,
                                  bisect_right(mask_block_starts, min_) - 1)
        last_masked_region = min(len(mask_block_starts),
                                 1 + bisect_right(mask_block_starts, max_,
                                                  lo=first_masked_region))
        for start, size in izip(mask_block_starts[first_masked_region:
        last_masked_region],
                                mask_block_sizes[first_masked_region:
                                last_masked_region]):
            end = start + size
            if end <= min_:
                continue
            if start > max_:
                break
            if start < min_:
                start = min_
            if end > max_:
                end = max_
            start -= min_
            end -= min_
            str_as_array[start:end] = array(_CHAR_CODE,
                                            lower(safe_tostring(str_as_array[start:end])))
        if not len(str_as_array) == max_ - min_:
            raise RuntimeError("Sequence was the wrong size")
        return safe_tostring(str_as_array)

    def __str__(self):
        """
        returns the entire chromosome
        """
        return self.get_slice(0, None)


class TwoBitFileError(Exception):
    """
    Base exception for TwoBit module
    """

    def __init__(self, msg):
        errtext = 'Invalid 2-bit file. ' + msg
        return super(TwoBitFileError, self).__init__(errtext)
