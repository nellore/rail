#!/usr/bin/env python
"""
guess.py
Part of Rail-RNA

Tools for inferring the properties of samples.
"""
import itertools

_RANGES = {
    'Sanger': (33, 93),
    'Solexa': (59, 104),
    'Illumina-1.3': (64, 104),
    'Illumina-1.5': (67, 104)
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

def phred_format(fastq_stream):
    """ Studies a selection of reads from a sample to determine Phred format.

        Inspired by https://github.com/brentp/bio-playground/blob/master/
        reads-utils/guess-encoding.py

        fastq_stream: where to read input fastq lines

        Return value: one of (Sanger, Solexa, Illumina-1.3, Illumina-1.5)
    """
    for qual in itertools.islice(fastq_stream, 0, None, 3):
        for char in qual:
            for version in _uniques:
                if char in _uniques[version]:
                    return version