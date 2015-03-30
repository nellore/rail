import sys

with open('dmel1.fastq', 'w') as f1:
    with open('dmel2.fastq', 'w') as f2:
        for i, line in enumerate(sys.stdin):
            if i % 8 == 0:
                cfh = f1
            elif (i + 4) % 8 == 0:
                cfh = f2
            cfh.write(line)
