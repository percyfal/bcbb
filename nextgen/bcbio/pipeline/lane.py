"""Top level driver functionality for processing a sequencing lane.
"""
import os
import copy
import csv

from bcbio.pipeline import log
from bcbio.pipeline.fastq import get_fastq_files
from bcbio.pipeline.demultiplex import split_by_barcode
from bcbio.pipeline.alignment import align_to_sort_bam

def process_lane(info, fc_name, fc_date, dirs, config):
    """Prepare lanes, potentially splitting based on barcodes.
    """
    config = _update_config_w_custom(config, info)

    sample_name = info.get("description", "")
    if (config["algorithm"].get("include_short_name", True) and
            info.get("name", "")):
        sample_name = "%s---%s" % (info.get("name", ""), sample_name)
    genome_build = info.get("genome_build", None)
    multiplex = info.get("multiplex", None)

    log.info("Processing sample: %s; lane %s; reference genome %s; " \
             "researcher %s; analysis method %s" %
             (sample_name, info["lane"], genome_build,
              info.get("researcher", ""), info.get("analysis", "")))
    if multiplex:
        log.debug("Sample %s is multiplexed as: %s" % (sample_name, multiplex))

    full_fastq1, full_fastq2 = get_fastq_files(dirs["fastq"], info['lane'], fc_name)
    lane_name = "%s_%s_%s" % (info['lane'], fc_date, fc_name)
    lane_items = []
    for mname, msample, fastq1, fastq2 in split_by_barcode(full_fastq1,
            full_fastq2, multiplex, lane_name, dirs, config):
        mlane_name = "%s_%s" % (lane_name, mname) if mname else lane_name
        if msample is None:
            msample = "%s---%s" % (sample_name, mname)
        lane_items.append((fastq1, fastq2, genome_build, mlane_name, msample,
                           dirs, config))
    
    # Append the demultiplexing results for this lane to the report file
    if multiplex:
        metrics_file = os.path.join(dirs["work"], "%s_barcode" % lane_name, "%s_bc.metrics" % lane_name)
        dmplx_report_file = os.path.join(dirs["work"], "demultiplexed_read_counts.txt")
        if os.path.exists(metrics_file):
            
            # Lookup the name and sequence of the barcode index
            barcodes = dict()
            for m in multiplex:
                barcodes[m.get('barcode_id','')] = [m.get('name','.'),m.get('sequence','.')]
            
            # Parse the demultiplexed barcode counts and store them and the metadata in a list
            dmplx = []
            with open(metrics_file,"rb") as mfr:
                csvr = csv.reader(mfr,dialect='excel-tab')
                for row in csvr:
                    d = [fc_date,fc_name,info['lane'],info.get('description','Lane '+str(info['lane'])),row[0]]
                    d.extend(barcodes.get(row[0],['.','.']))
                    d.append(row[1])
                    dmplx.append(d)
            
            # Append the results to the report file 
            with open(dmplx_report_file,"ab") as mfw:
                csvw = csv.writer(mfw,dialect='excel-tab')
                csvw.writerows(dmplx)    
        
    return lane_items

def process_alignment(fastq1, fastq2, genome_build, lane_name, sample, dirs, config):
    """Do an alignment of fastq files, preparing a sorted BAM output file.
    """
    aligner = config["algorithm"].get("aligner", None)
    if os.path.exists(fastq1) and aligner:
        log.info("Aligning lane %s with %s aligner" % (lane_name, aligner))
        align_to_sort_bam(fastq1, fastq2, genome_build, aligner,
                          lane_name, sample, dirs, config)

def _update_config_w_custom(config, lane_info):
    """Update the configuration for this lane if a custom analysis is specified.
    """
    config = copy.deepcopy(config)
    analysis_type = lane_info.get("analysis", "")
    custom = config["custom_algorithms"].get(analysis_type, None)
    if custom:
        for key, val in custom.iteritems():
            config["algorithm"][key] = val
    return config

