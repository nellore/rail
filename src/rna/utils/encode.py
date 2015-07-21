"""
encode.py
Part of Rail-RNA

Contains functions for encoding and decoding integers and base sequences
as shorter ASCII chars.
"""
import string

_alphabet = string.digits + string.ascii_lowercase + string.ascii_uppercase + \
            '-_'
_seq_translation_table = string.maketrans('ATCGN', '12345')
_reverse_translation_table = string.maketrans('12345', 'ATCGN')

def encode(value, base=64):
    """ Encodes an integer in base (len alphabet) to shorten it.

        This function is used to assign short names to reads by record.
        Based on http://stackoverflow.com/questions/561486/
        how-to-convert-an-integer-to-the-shortest-url-safe-string-in-python/
        561534#561534 . Does not reverse encode string because it's 
        unnecessary.

        value: integer

        Return value: base (len alphabet) encoded value
    """
    assert value >= 0
    s = []
    while True:
        value, r = divmod(value, base)
        s.append(_alphabet[r])
        if value == 0: break
    return ''.join(reversed(s))

def encode_sequence(seq):
    """ Encodes a sequence of ATCGNs in base 36.

        seq: sequences of ATCGNs

        Return value: encoded string
    """
    return encode(int(seq.translate(_seq_translation_table), 6), base=36)

def decode_sequence(seq):
    """ Decodes a sequence encoded by encode_sequence().

        seq: sequence to be decoded

        Return value: decoded string
    """
    return encode(int(seq, 36), 6).translate(_reverse_translation_table)