"""
log_analyzer.py
Part of Rail-RNA

Takes a job flow ID log directory and tests whether counters
(# map input records and # reduce output records) are consistent. Automatically
takes into account multiple outputs as long as they're in different
directories.
"""

base_path = os.path.abspath(
                    os.path.dirname(os.path.dirname(os.path.dirname(
                        os.path.realpath(__file__)))
                    )
                )
utils_path = os.path.join(base_path, 'rna', 'utils')
site.addsitedir(utils_path)
site.addsitedir(base_path)

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description=__doc__, 
                    formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_args('-l', '--log', type=str, required=True,
            help=('Job flow ID subdirectory of a log directory; '
                  'example: s3://my-bucket/log/j-17GPQ51QPR5YV'
        )