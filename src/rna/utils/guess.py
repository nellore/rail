#!/usr/bin/env python
"""
guess.py
Part of Rail-RNA

Tools for inferring the properties of samples.
"""
import itertools
import math

_RANGES = {
    'Sanger': (33, 93),
    'Phred64': (64, 104)
}

_uniques = {}
for version in _RANGES:
    _uniques[version] = set(range(_RANGES[version][0], _RANGES[version][1]))
    for compared_version in _RANGES:
        if compared_version != version:
            _uniques[version] -= set(
                range(
                    _RANGES[compared_version][0], _RANGES[compared_version][1]
                )
            )
for version in _uniques:
    _uniques[version] = set([chr(char) for char in _uniques[version]])

def inferred_phred_format(fastq_stream, at_once=500):
    """ Studies a selection of reads from a sample to determine Phred format.

        fastq_stream: where to read input fastq lines or None if format is
            provided
        at_once: number of quality records to read in at once before checking
            for distinguishing characters

        Return value: one of {Sanger, Solexa, Illumina-1.3, Illumina-1.5};
            assumes Sanger if no distinguishing characters are found
    """
    chars = set()
    for i, qual in enumerate(itertools.islice(fastq_stream, None, None, 3)):
        if not (i % at_once):
            for char in chars:
                for version in _uniques:
                    if char in _uniques[version]:
                        return version
            chars = set()
        chars.add(qual)
    for char in qual:
        for version in _uniques:
            if char in _uniques[version]:
                return version
    # Still nada? Assume Sanger
    return 'Sanger'

def phred_converter(fastq_stream=None, phred_format=None, at_once=500):
    """ Provides a function that converts a quality string to Sanger format

        Inspired by https://github.com/brentp/bio-playground/blob/master/
        reads-utils/guess-encoding.py . Returns a function

        fastq_stream: where to read input fastq lines or None if format is
            provided
        phred_format: Phred format from _RANGES or None if fastq_stream is
            provided
        at_once: number of quality records to read in at once before checking
            for distinguishing characters

        Return value: one of (Sanger, Solexa, Illumina-1.3, Illumina-1.5)
    """
    assert fastq_stream is not None or platform is not None, ('Either a '
        'fastq stream must be provided to infer platform or a platform '
        'must be provided directly.')
    assert phred_format in _RANGES, ('Platform must be Sanger, Solexa, '
        'Illumina-1.3, Illumina-1.5')
    if phred_format is None:
        phred_format = inferred_phred_format(fastq_stream, at_once)
    # if phred_format == 'Solexa':
    #    def final_converter(qual):
    #        return ''.join([
    #                            chr(round(
    #                        10*math.log(1+10**((ord(char)-64)/10.0),10)
    #                )+33) for char in qual]
    #            )
    if phred_format == 'Sanger':
        def final_converter(qual):
            return qual
    else:
        # It's Phred64
        def final_converter(qual):
            return ''.join([chr(ord(char) - 31) for char in qual])
    return final_converter