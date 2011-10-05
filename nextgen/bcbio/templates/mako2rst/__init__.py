"""
Mako templates for rst output
"""

import os
import sys
import glob
import json
from texttable import *

##################################################
## Misc inputs
########################################
def insert_images(png, columns=2):
    tab = Texttable()
    rows = []
    for i in range(0, len(png), columns):
        rows.append([".. image:: %s" % x for x in png[i:(i+columns)]])
    maxl = 0
    for row in rows:
        for col in row:
            if len(col) > maxl:
                maxl = len(col)
    tab.add_rows(rows, header=False)
    tab.set_cols_width([maxl for i in range(0,columns)])
    return tab.draw()

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
        

def teqc_json(d, by_sample=True):
    if d is None:
        return
    res = []
    if by_sample:
        enrichment_table = Texttable()
        flanking_table = Texttable()
        coverage_table = Texttable()
        for fc in d.keys():
            js = d[fc]
            samples = js.keys()
            enrichment_table.add_rows([["key"] + samples,
                                       ["enrichment"] + [d[fc][x]['enrichment'] for x in samples],
                                       ["max theoretical enrichment"] + ["%.1f" % (1.0/float(d[fc][x]['target']['fraction'])) for x in samples],
                                       ["target fraction (%)"] + [(d[fc][x]['target']['fraction'] * 100) for x in samples],
                                       ["target width (Mb)"] + [ ("%.2f" % (int(d[fc][x]['target']['width']) / int(100000))) for x in samples],
                                       ["mean coverage"] + [d[fc][x]['coverage']['avg'] for x in samples],
                                       ["mean coverage"] + [d[fc][x]['coverage']['sd'] for x in samples]
                                       ])
            res.append(enrichment_table.draw())

            flanking_table.add_rows([["flanking region \ capture specificity"] + samples,
                                     ["0"] + [d[fc][s]['capture_specificity']['flank_0'] for s in samples],
                                     ["50"] + [d[fc][s]['capture_specificity']['flank_50'] for s in samples],
                                     ["100"] + [d[fc][s]['capture_specificity']['flank_100'] for s in samples],
                                     ])
            res.append("\n")
            res.append(flanking_table.draw())
            ck = d[fc][samples[0]]['coverage']['k']
            def get_coverage(s):
                tmp = [s] + [d[fc][s]['coverage']['k'][x] for x in sorted(ck.keys(), key=int)]
                return tmp

            coverage_table.add_rows([["sample"] + [x for x in sorted(ck.keys(), key=int)]])
            coverage_table.add_rows([get_coverage(s) for s in samples ], header=False)

            res.append("\n")
            res.append(coverage_table.draw())


    else:
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
                
    return "\n\n".join(res)

def teqc_graphics(d, which="chrom-barplot", width="65%", as_table=True, columns=2):
    res = []
    tab = Texttable()
    if d==None:
        return
    for fc in d.keys():
        flowcell = d[fc]
        png = flowcell[which]
        if len(png) == 0:
            next
        else:
            if as_table:
                rows = []
                for i in range(0, len(png), columns):
                    rows.append([".. image:: %s" % x for x in png[i:(i+columns)]])
                maxl = 0
                for row in rows:
                    for col in row:
                        if len(col) > maxl:
                            maxl = len(col)
                tab.add_rows(rows, header=False)
                tab.set_cols_width([maxl for i in range(0,columns)])
                res = tab.draw()
            else:
                for grf in png:
                    res.append(".. figure:: %s\n    :width: %s\n\n" % (grf, width))
                res = "\n".join(res)
    return res

    
