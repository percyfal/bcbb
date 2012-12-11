#!/usr/bin/env python
"""Perform an automated analysis on a sequencing run using Galaxy information.

Given a directory of solexa output, this retrieves details about the sequencing
run from the Galaxy description, and uses this to perform an initial alignment
and analysis.

Usage:
    automated_initial_analysis.py <YAML config file> <flow cell dir>
                                  [<YAML run information>]

The optional <YAML run information> file specifies details about the flowcell
lanes, instead of retrieving it from Galaxy. An example configuration file is
located in 'config/run_info.yaml'

Workflow:
    - Retrieve details on a run.
    - Align fastq files to reference genome.
    - Perform secondary analyses like SNP calling.
    - Generate summary report.
"""
import os
import sys
from optparse import OptionParser
import datetime

import yaml
import logbook

from bcbio.galaxy.api import GalaxyApiAccess
from bcbio.solexa.flowcell import get_fastq_dir
from bcbio import utils
from bcbio.log import logger2 as logger
from bcbio.log import setup_logging
from bcbio.log import version
from bcbio.log import create_log_handler
from bcbio.log import RecordProgress
from bcbio.distributed.messaging import parallel_runner
from bcbio.pipeline.run_info import get_run_info
from bcbio.pipeline.demultiplex import add_multiplex_across_lanes
from bcbio.pipeline.merge import organize_samples
from bcbio.pipeline.qcsummary import write_metrics, write_project_summary
from bcbio.variation.realign import parallel_realign_sample
from bcbio.variation.genotype import parallel_variantcall
from bcbio.pipeline.config_loader import load_config


def main(config_file, fc_dir, run_info_yaml=None):
    config = load_config(config_file)
    work_dir = os.getcwd()
    if config.get("log_dir", None) is None:
        config["log_dir"] = os.path.join(work_dir, "log")

    def insert_command(record):
        record.extra["command"] = sys.argv
        record.extra["version"] = version.get_pipeline_version()

    setup_logging(config)
    handler = create_log_handler(config)
    with handler, \
         logbook.Processor(insert_command):

        run_main(config, config_file, fc_dir, work_dir, run_info_yaml)


def run_main(config, config_file, fc_dir, work_dir, run_info_yaml):
    _record_sw_versions(config, os.path.join(work_dir, "bcbb_software_versions.txt"))
    prog = RecordProgress(work_dir)
    to_compress = set()
    prog.progress("analysis_start")

    align_dir = os.path.join(work_dir, "alignments")
    run_module = "bcbio.distributed"
    fc_name, fc_date, run_info = get_run_info(fc_dir, config, run_info_yaml)
    fastq_dir, galaxy_dir, config_dir = _get_full_paths(get_fastq_dir(fc_dir),
                                                        config, config_file)

    config_file = os.path.join(config_dir, os.path.basename(config_file))
    dirs = {"fastq": fastq_dir, "galaxy": galaxy_dir, "align": align_dir,
            "work": work_dir, "flowcell": fc_dir, "config": config_dir}

    run_parallel = parallel_runner(run_module, dirs, config, config_file)
    run_items = add_multiplex_across_lanes(run_info["details"], dirs["fastq"], fc_name)

    lanes = ((info, fc_name, fc_date, dirs, config) for info in run_items)
    lane_items = run_parallel("process_lane", lanes)
    [to_compress.add(f) for f in lane_items[0][0:2]]
    prog.dummy()
    prog.progress("process_lane")

    # Remove spiked in controls, contaminants etc.
    lane_items = run_parallel("remove_contaminants", lane_items)
    [to_compress.add(f) for f in lane_items[0][0:2]]
    prog.dummy()
    prog.progress("remove_contaminants")
    align_items = run_parallel("process_alignment", lane_items)
    [to_compress.add(f) for f in align_items[0]['fastq']]
    prog.dummy()
    prog.progress("process_alignment")

    # process samples, potentially multiplexed across multiple lanes
    samples = organize_samples(align_items, dirs, config_file)
    samples = run_parallel("merge_sample", samples)
    to_compress.add(samples[0][0]['fastq1'])
    to_compress.add(samples[0][0]['fastq2'])
    prog.dummy()
    prog.progress("merge_sample")
    samples = run_parallel("mark_duplicates_sample", samples)
    to_compress.add(samples[0][0]['fastq1'])
    to_compress.add(samples[0][0]['fastq2'])
    prog.dummy()
    prog.progress("mark_duplicates_sample")
    run_parallel("screen_sample_contaminants", samples)
    prog.dummy()
    prog.progress("screen_sample_contaminants")
    samples = run_parallel("recalibrate_sample", samples)
    prog.dummy()
    prog.progress("recalibrate_sample")
    samples = parallel_realign_sample(samples, run_parallel)
    prog.dummy()
    prog.progress("realign_sample")
    samples = parallel_variantcall(samples, run_parallel)
    prog.dummy()
    prog.progress("variantcall")
    samples = run_parallel("detect_sv", samples)
    prog.dummy()
    prog.progress("detect_sv")
    samples = run_parallel("process_sample", samples)
    prog.dummy()
    prog.progress("process_sample")
    samples = run_parallel("generate_bigwig", samples, {"programs": ["ucsc_bigwig"]})
    prog.dummy()
    prog.progress("generate_bigwig")
    write_project_summary(samples)
    write_metrics(run_info, fc_name, fc_date, dirs)
    prog.dummy()
    prog.progress("write_metrics")

    # Compress all files in to_compress
    if config['algorithm'].get('compress_files', True):
        (before, after) = utils.compress_files(to_compress)
        logger.info("Space used by the files before compressing (in bytes): " \
                     + str(before))
        logger.info("Space used by the files after compressing (in bytes): " \
                     + str(after))
        logger.info("Saved space (in bytes): " + str(before - after))


# Utility functions

def _record_sw_versions(config, sw_version_file):
    """Get the versions of software used in the pipeline and output to
       log and text file in working directory
    """
    sw_versions = version.get_versions(config)
    sw_versions['bcbb'] = version._get_git_commit()

    logger.info("bcbb pipeline is running with software versions: %s" % sw_versions)

    with open(sw_version_file, 'w') as fh:
        fh.write("%s\n" % datetime.datetime.now().isoformat())
        for sw, ver in sw_versions.items():
            fh.write("%s\t%s\n" % (sw, ver))


def _get_full_paths(fastq_dir, config, config_file):
    """Retrieve full paths for directories in the case of relative locations.
    """
    fastq_dir = utils.add_full_path(fastq_dir)
    config_dir = utils.add_full_path(os.path.dirname(config_file))
    galaxy_config_file = utils.add_full_path(config["galaxy_config"], config_dir)

    return fastq_dir, os.path.dirname(galaxy_config_file), config_dir


def _get_run_info(fc_name, fc_date, config, run_info_yaml):
    """Retrieve run information from a passed YAML file or the Galaxy API.
    """
    if run_info_yaml and os.path.exists(run_info_yaml):
        logger.info("Found YAML samplesheet, using %s instead of Galaxy API" % run_info_yaml)
        with open(run_info_yaml) as in_handle:
            run_details = yaml.load(in_handle)

        return dict(details=run_details, run_id="")

    else:
        logger.info("Fetching run details from Galaxy instance")
        galaxy_api = GalaxyApiAccess(config['galaxy_url'], config['galaxy_api_key'])

        return galaxy_api.run_details(fc_name, fc_date)

if __name__ == "__main__":
    parser = OptionParser()
    (options, args) = parser.parse_args()
    if len(args) < 2:
        print "Incorrect arguments"
        print __doc__
        sys.exit()

    kwargs = dict()
    main(*args, **kwargs)
