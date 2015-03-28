"""
harvest_spliced_read_performance.py

Puts all spliced-read performance results in the same directory, labeling
files appropriately.
"""
import os
import shutil

def files_in_dir(path):
    """ Returns full paths of all files in dir and its subdirectories.

        WARNING: recursive function; could choke in Python, but works for this
            purpose.

        path: a directory

        Return value: list of files
    """
    all_files = []
    for a_file in os.listdir(path):
        a_file = os.path.join(path, a_file)
        if os.path.isfile(a_file):
            all_files.append(os.path.abspath(a_file))
        elif os.path.isdir(a_file):
            all_files.extend(files_in_dir(os.path.abspath(a_file)))
    return all_files

if __name__ == '__main__':
    import argparse
    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__, 
                formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--input', '-i', type=str, required=True,
            help='Where to find performance files. Can be the directory '
                 'containing sample names or a sample name directory itself'
        )
    parser.add_argument('--output', '-o', type=str, required=True,
            help='Where to copy performance files'
        )
    args = parser.parse_args()
    try:
        os.makedirs(args.output)
    except OSError:
        # Already exists
        pass
    performance_files = [a_file for a_file 
                            in files_in_dir(os.path.abspath(args.input))
                            if 'perform' in a_file]
    args.output = os.path.abspath(args.output)
    root_size = len(args.output)
    for a_file in performance_files:
        new_name = os.path.join(args.output,
                                    a_file.replace('/', '.'))[root_size+1:]
        shutil.copy(a_file, new_name)