#!/usr/bin/env python2.7
import sys
import os
import unittest
import emr_simulator as esim

def setUpModule():
    pass

def tearDownModule():
    pass

class Test_presorted_tasks(unittest.TestCase):
    '''
    tests the original presorted_tasks method
    and it's refactored sub-methods
    '''
    
    def setUp(self):
        self.test_presort_input_file = 'sort_2_liner.tsv'
        with open(self.test_presort_input_file,"w") as f:
            f.write("1\thttp://verve.webfactional.com/dm3_example_1_left.fastq\t0\thttp://verve.webfactional.com/dm3_example_1_right.fastq\t0\tdm3_example-1-1\n")
            f.write("0\thttp://verve.webfactional.com/dm3_example_1_left.fastq\t0\thttp://verve.webfactional.com/dm3_example_1_right.fastq\t0\tdm3_example-1-2\n")
        self.test_presort_output_dir = './'
        pass

    def test_presorted_tasks_toplevel(self):
        '''
        basic test to see if presorted_tasks(...) returns and error or not
        '''
        input_files = ["%s" % self.test_presort_input_file]
        process_id = 0
        sort_options = '-k1,1'
        output_dir = "%s" % self.test_presort_output_dir
        key_fields = 1
        separator = ' '
        partition_options = '-k1,1'
        task_count = 4
        memcap = 307200
        gzip = False
        gzip_level = 3
        scratch = '-'
        direct_write = False
        sort = 'sort'
        mod_partition = False
        max_attempts = 0
        error = esim.presorted_tasks(input_files,process_id,sort_options,output_dir,key_fields,separator,partition_options,task_count,memcap,gzip=gzip,gzip_level=gzip_level,scratch=scratch,direct_write=direct_write,sort=sort,mod_partition=mod_partition,max_attempts=max_attempts)
        self.assertEqual(error,None)
        
        gzip = True
        error = esim.presorted_tasks(input_files,process_id,sort_options,output_dir,key_fields,separator,partition_options,task_count,memcap,gzip=gzip,gzip_level=gzip_level,scratch=scratch,direct_write=direct_write,sort=sort,mod_partition=mod_partition,max_attempts=max_attempts)
        self.assertEqual(error,None)
        
if __name__ == '__main__':
    unittest.main()
