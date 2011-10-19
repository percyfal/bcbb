"""This directory is setup with configurations to run the main functional test.

It exercises a samplebased analysis pipeline on a smaller subset of data, as implemented at SciLife.
"""
import os
import sys
import subprocess
import unittest
import shutil
import contextlib
import glob
import yaml

@contextlib.contextmanager
def make_workdir():
    dirname = os.path.join(os.path.dirname(__file__), "projects", "j_doe_00_01", "intermediate", "nobackup", "110106_FC70BUKAAXX")
    if os.path.exists(dirname):
        if os.path.islink(dirname):
            os.remove(dirname)
        else:
            shutil.rmtree(dirname)
    os.makedirs(dirname)
    orig_dir = os.getcwd()
    try:
        os.chdir(dirname)
        yield
    finally:
        os.chdir(orig_dir)

def init_flowcell_dir():
    dirname = os.path.join(os.path.dirname(__file__), "110106_FC70BUKAAXX")
    fcdir = os.path.join(os.path.dirname(__file__), "test_automated_output")
    if os.path.exists(dirname):
        if os.path.islink(dirname):
            os.remove(dirname)
        else:
            shutil.rmtree(dirname)
    os.symlink(fcdir, dirname)
    

class SampleBasedAnalysisTest(unittest.TestCase):
    """Setup a sample based scilife analysis
    """
    def setUp(self):
        self.file_dir = os.path.join(os.path.dirname(__file__))
        self.fc_dir = os.path.join(self.file_dir, "110106_FC70BUKAAXX")
        self.proj_dir = os.path.join(self.file_dir, "projects", "j_doe_00_01")
        ##self.fcdir = os.path.join(os.path.dirname(__file__), "test_automated_output")
        self.run_info = os.path.join(self.fc_dir, "run_info.yaml")
        self.archive_base_dir  = os.path.join(self.file_dir)
        self.analysis_base_dir = os.path.join(self.file_dir)

        # Remove fcdir if exists and setup new link
        init_flowcell_dir()
        if not os.path.exists(self.proj_dir):
            os.makedirs(self.proj_dir)
        if not os.path.exists(os.path.join(self.file_dir, "test_automated_output", "run_info.yaml")):
            os.symlink(os.path.join(self.file_dir, "data", "automated", "run_info-project.yaml"), os.path.join(self.file_dir, "test_automated_output", "run_info.yaml"))
        if not os.path.exists(os.path.join(self.file_dir, "test_automated_output", "tool-data")):
            os.symlink(os.path.join(self.file_dir, "data", "automated", "tool-data"), os.path.join(self.file_dir, "test_automated_output", "tool-data"))

        # Post_process.yaml
        with open(os.path.join(self.file_dir, "data", "automated", "post_process.yaml"), "r") as fh:
            post_process = yaml.load(fh)
        post_process["analysis"]["store_dir"] = os.path.join(self.archive_base_dir)
        post_process["analysis"]["base_dir"] = os.path.join(self.analysis_base_dir)
        post_process["algorithm"]["snpcall"] = "true"
        post_process["algorithm"]["dbsnp"] = os.path.join("data", "genomes", "hg19", "variation", "dbsnp_132.vcf")
        with open(os.path.join(self.fc_dir, "post_process.yaml"), "w") as fh:
            yaml.dump(post_process, stream=fh)

    def _deliver_data(self):
        print "Delivering data"
        cl = ["project_analysis_setup.py",
              os.path.join(self.fc_dir, "post_process.yaml"),
              self.fc_dir, self.proj_dir,
              self.run_info,
              "--data_prefix=intermediate/nobackup",
              "--project_desc=%s" % "J.Doe_00_01",
              "--flowcell_alias=20000101A_hiseq2000"]
        subprocess.check_call(cl)
        print "Finished delivering data..."

    def test_targeted_resequencing_pipeline(self):
        """Test a sample based targeted resequencing pipeline"""
        with make_workdir():
            print "Going to deliver data"
            self._deliver_data()
            cl = ["project_exome_pipeline.py",
                  os.path.join(self.fc_dir, "post_process.yaml"),
                  os.path.join(self.proj_dir, "intermediate", "nobackup", "110106_FC70BUKAAXX"),
                  os.path.join(self.proj_dir, "intermediate", "nobackup", "20000101A_hiseq2000", "project_run_info.yaml")]

            subprocess.check_call(cl)
