"""
check.py

Sanity check a MANIFEST file; basically, checking that all the URLs in the
file exist and we haven't left out an unpaired/paired-end counterpart.

I recommend you run once without --check-urls and, if that passes, again with
--check-urls.
"""

import os
import sys
import string
import time
import subprocess
import argparse

parser = argparse.ArgumentParser(description=\
    'Scan manifest file, summarize it, and check for obvious mistakes.')
parser.add_argument(\
    '--manifest', type=str, required=True, help='Manifest file')
parser.add_argument(\
    '--check-urls', action='store_const', const=True, default=False,
    help='Check whether files at URLs exist')

args = parser.parse_args()

bytes = 0

def exists(url, attempts=5, timeout=5):
    global bytes
    while attempts > 0:
        pipe = subprocess.Popen('curl -I %s 2>/dev/null' % url, shell=True, stdout=subprocess.PIPE)
        for ln in pipe.stdout:
            if ln.startswith('Content-Length'):
                toks = ln.split()
                bytes += int(toks[-1])
        if pipe.wait() == 0:
            return True
        attempts -= 1
        time.sleep(timeout)
        timeout *= 2
        print >> sys.stderr, 'WARNING: URL "%s" yielded non-zero return from curl' % url
    return False

def urlToSrr(url):
    toks = string.split(url, '/')
    return toks[-2]

def nameToGroup(name):
    return name[:name.index('-')]

srrs = set()
urls = set()
namesP, namesU = set(), set()
groups = set()
groups = set()

with open(args.manifest) as fh:
    print >> sys.stderr, 'Processing manifest file "%s"' % args.manifest
    for ln in fh:
        ln = ln.rstrip()
        if len(ln) == 0:
            continue
        if ln[0] == '#':
            continue
        toks = string.split(ln, '\t')
        if len(toks) == 3:
            url, md5, name = toks
            assert url not in urls, url
            urls.add(url)
            srrs.add(urlToSrr(url))
            assert name.count('-') == 2
            assert name not in namesU, name
            namesU.add(name)
            groups.add(nameToGroup(name))
            if args.check_urls:
                if not exists(url):
                    raise RuntimeError('No such unpaired URL: "%s"' % url)
        else:
            assert len(toks) == 5, toks
            url_1, md5_1, url_2, md5_2, name = toks
            assert url_1 not in urls, url_1
            urls.add(url_1)
            srrs.add(urlToSrr(url_1))
            assert url_2 not in urls, url_2
            urls.add(url_2)
            srrs.add(urlToSrr(url_2))
            assert name.count('-') == 2, name
            assert name not in namesP, name
            namesP.add(name)
            groups.add(nameToGroup(name))
            if args.check_urls:
                if not exists(url_1):
                    raise RuntimeError('No such mate-1 URL: "%s"' % url_1)
                if not exists(url_2):
                    raise RuntimeError('No such mate-2 URL: "%s"' % url_2)
        
        if args.check_urls:
            print >> sys.stderr, 'GB so far: %0.3f' % (bytes / float(1024 * 1024 * 1024))

print >> sys.stderr, "Saw %d URLs" % len(urls)
print >> sys.stderr, "Saw %d SRRs" % len(srrs)
print >> sys.stderr, "Saw %d groups (%d names unp, %d names paired)" % (len(groups), len(namesU), len(namesP))
