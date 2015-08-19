#!/usr/bin/env python2.7

import sys
import os
SCRIPT_PATH = os.path.dirname(os.path.realpath(os.path.abspath(__file__)))
sys.path.append(os.path.dirname(SCRIPT_PATH))
import unittest
import socket
import tempfile
from dooplicity import emr_simulator

class TestNewParallelCode(unittest.TestCase):

    def setUp(self):
       (temp_file_handle_low_level, 
            self.temp_file_path) = tempfile.mkstemp(
                          prefix="test_emr_simulator_",text=True)
       temp_file_handle = os.fdopen(temp_file_handle_low_level,"w")
       temp_file_handle.write("test123")
       temp_file_handle.close()
       self.iface = emr_simulator.create_iface(None, None) 
       hosts = ["test.localhost", "test.stringray"]
       self.file_tracker = emr_simulator.FileTracker(hosts, self.iface)

    def tearDown(self):
        os.unlink(self.temp_file_path)
        new_temp_file = "%s_new" % self.temp_file_path
        if os.path.exists(new_temp_file):
            os.unlink(new_temp_file)

    def test_FileTracker_add_file_mapping(self):
       self.file_tracker.add_file_mapping("test.txt", "test.localhost", "/tmp") 

    def test_copy_file_locally(self):
       error = emr_simulator.copy_file(self.temp_file_path, 
                                       "%s_new" % self.temp_file_path)
       self.assertEqual(error, None)

#for the copy_file_remotely test to work a public key 
#must be in the authorized_keys file for this host for this user
    def test_copy_file_remotely(self):
       localhost = socket.gethostname()
       error = emr_simulator.copy_file(self.temp_file_path, 
                                       "%s_new" % self.temp_file_path,
                                       remote_host=localhost)
       self.assertEqual(error, None)

if __name__ == '__main__':
    unittest.main()
