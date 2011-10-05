"""
Mako templates for rst output
"""

import os
import sys
import glob
import json
from texttable import *

def teqc_config(indir, samples):
    """Simple teqc configuration if no run_info.yaml file"""
    jdata = {}
    fc = os.path.basename(indir)
    for s in samples:
        jdata[fc] = {}
        infiles = glob.glob(os.path.join(indir, "*.json"))
        for f in infiles:
            fp = open(f)
            jd = json.load(fp)
            fp.close()
            jdata[fc][os.path.basename(f)] = jd
    tdata = {}
    png = {'chrom-barplot': [os.path.relpath(x) for x in glob.glob(os.path.join(indir, "*-chrom-barplot.png"))],
           'coverage-hist': [os.path.relpath(x) for x in glob.glob(os.path.join(indir, "*-coverage-hist.png"))],
           'coverage-targetlength-plot-avgCoverage': [os.path.relpath(x) for x in glob.glob(os.path.join(indir, "*-coverage-targetlength-plot-avgCoverage.png"))],
           'coverage-targetlength-plot-nReads': [os.path.relpath(x) for x in glob.glob(os.path.join(indir, "*-coverage-targetlength-plot-nReads.png"))],
           'coverage-uniformity': [os.path.relpath(x) for x in glob.glob(os.path.join(indir, "*-coverage-uniformity.png"))],
           'duplicates-barplot': [os.path.relpath(x) for x in glob.glob(os.path.join(indir, "*-duplicates-barplot.png"))],
           'insert-size-hist': [os.path.relpath(x) for x in glob.glob(os.path.join(indir, "*-insert-size-hist.png"))]
           }
    tdata[fc] = png
    return [jdata, tdata]
        
    
def program_info(proj_conf):
    d = proj_conf['program']
    tab = Texttable()
    tab.set_cols_align(["l", "l"])
    tab.set_cols_align(["b", "b"])
    tab.header(["Program", "Value"])
    for k in d.keys():
        tab.add_row([k, d[k]])
    return tab.draw()

def duplication_metrics(d):
    if d==None:
        return
    tab = Texttable()
    add = False
    first = True
    for fc in d.keys():
        md = d[fc]
        for lab in md.keys():
            data = md[lab]
            for row in data:
                if row[0] == "LIBRARY":
                    add = True
                    if first:
                        tab.add_row(row)
                    first = False
                elif row[0] == "BIN":
                    add = False
                tab.add_row(row)
    return tab.draw()
# Has 19 columns, need better way of representing it
def align_metrics(d):
    if d==None:
        return
    tab = Texttable()
    first = True
    for fc in d.keys():
        md = d[fc]
        for lab in md.keys():
            data = md[lab]
            for row in data:
                if row[0] == "CATEGORY": 
                    if first:
                        tab.header(row)
                        tw = []
                        for td in row:
                            tw.append(len(td))
                        tab.set_cols_width(tw)
                        first = False
                else:
                    tab.add_row(row)
    return tab.draw()


##################################################
## Illumina raw data
##################################################
def image(fp, width):
    res = ".. figure:: %s\n    :width: %s\n\n" % (fp, width)
    return res

##################################################
## Output for TEQC
##################################################
def teqc_json(d):
    if d is None:
        return
    res = []
    for fc in d.keys():
        js = d[fc]
        for s in js.keys():
            tab = Texttable()
            tab.add_rows([["key", "value"],
                          ['enrichment', d[fc][s]['enrichment']],
                          ['max theoretical enrichment', "%.1f" % ( 1.0/float(d[fc][s]['target']['fraction']))],
                          ['target fraction (%)', d[fc][s]['target']['fraction'] * 100],
                          ['target width (Mb)', "%.2f" % (int(d[fc][s]['target']['width']) / int(1000000))],
                          ['mean coverage', d[fc][s]['coverage']['avg']],
                          ['coverage sd', d[fc][s]['coverage']['sd']]
                          ])
                         
            res.append("%s\n^^^^^^^^^^^^^^^^^\n" % (s))
            res.append(tab.draw())

            tab = Texttable()
            tab.add_rows([["flanking region", "capture specificity"],
                          ["0", d[fc][s]['capture_specificity']['flank_0']],
                          ["50", d[fc][s]['capture_specificity']['flank_50']],
                          ["100", d[fc][s]['capture_specificity']['flank_100']]])
            res.append(tab.draw())
            
            tab = Texttable()
            ck = d[fc][s]['coverage']['k']
            tab.add_rows([sorted(ck.keys(),key=int),
                          [ck[x] for x in sorted(ck.keys(), key=int)]])
            res.append(tab.draw())

            tab = Texttable()
            tab

    return "\n\n".join(res)

def teqc_graphics(d, which="chrom-barplot", width="65%", as.table=True, columns=2):
    res = []
    if d==None:
        return
    for fc in d.keys():
        flowcell = d[fc]
        png = flowcell[which]
        if len(png) == 0:
            next
        else:
            for grf in png:
                res.append(".. figure:: %s\n    :width: %s\n\n" % (grf, width))
    return "\n".join(res)

    
