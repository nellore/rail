#!/usr/bin/env python
"""
tempdel.py
Part of Rail-RNA

For deleting temporary directories on exit from Python script.
"""
import shutil

def remove_temporary_directories(temp_dir_paths):
    """ Deletes temporary directory.

        temp_dir_paths: iterable of paths of temporary directories

        No return value.
    """
    for temp_dir_path in temp_dir_paths:
        shutil.rmtree(temp_dir_path)
