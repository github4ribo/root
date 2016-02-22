#! /usr/bin/python
# -*- coding: utf-8 -*-
###################################################
# Program: exome_annalysis_prepare.py
# Function: preparation for exom analysis
# Author: Ouyang Guojun <guojun.ouyang@sirna.cn>
# Version: V1.0
# Date: Fri Jun 27 09:51:48 2014
###################################################
import sys
import os
import glob
import time
BIN = os.path.abspath(os.path.dirname(sys.argv[0]))

anodb     =   BIN + "/humandb/"
chrbed    =   BIN + "/locexon/hg19.bed"
#exom database
exombed   = BIN + "/humandb/Trueseq/TruSeq_exome_targeted_regions.exome.hg19.bed"
nearbed   = BIN + "/humandb/Trueseq/TruSeq_exome_targeted_regions.near.hg19.bed"
targetbed = BIN + "/humandb/Trueseq/TruSeq_exome_target_flanking.merged.bed"

Reseq_exe = BIN + "/resequencing3.py"
#Ref_fa    = BIN +  "/humandb/hg19.fa"
Ref_fa    = BIN +  "/GATK_resource_bundle/ucsc.hg19.fasta"
sh_head   = "#! /bin/bash\n"

def main(Args):
    if not Args:
        print ("Usage: python %s <rawdata_dir> <analysis_dir> <sample1..N>" % sys.argv[0])
        sys.exit(0)

    in_dir = os.path.abspath(Args[0])
    analysis_dir = os.path.abspath(Args[1])
    samples = Args[2:]

    #-----------------------make shell----------------------------------------#
    for sample in samples:
        sample_analysis_dir = analysis_dir + "/analysis/" + sample
        if os.path.exists(in_dir + "/" + 'cleandata'):
            if not os.path.exists(sample_analysis_dir):
                os.makedirs(sample_analysis_dir)
            try:
                raw_data_dir = in_dir + "/" + 'cleandata/' + sample
                os.system("ln -fs %s/01.DataCleaning/%s.read*_Clean.fastq.gz %s" % (raw_data_dir, sample, sample_analysis_dir))
                reseq_sh = open(sample_analysis_dir + "/" + sample + ".sh", "w")
                reseq_sh.write(sh_head)
                reseq_sh.write("#! /bin/sh\n#BSUB -n 1\n#BSUB -oo %s.%%J.o\n#BSUB -e %s.%%J.e\n" % (sample, sample))
                reseq_sh.write("%s -i %s -o %s -r %s -M True -E True -R True -C True --qcCheck False\
 --insertsize 306  --mismatch 4 --gapopen 5 --gapextention 20 \
 -A True --anodb %s --species hg19 --chrbed %s --readlen 100 \
 --exombed %s --nearbed %s --targetbed %s --basename %s"\
 %(Reseq_exe, sample_analysis_dir, sample_analysis_dir, Ref_fa, anodb, chrbed, exombed, nearbed, targetbed, sample))
                reseq_sh.close()
                os.system("sh %s&" % (sample_analysis_dir + "/" + sample + ".sh"))
                time.sleep(2)
            except:
                print >> sys.stderr, "sample %s is not in %s!" % (sample, raw_data_dir)
                sys.exit(0)
        else:
            if not os.path.exists(sample_analysis_dir):
                os.makedirs(sample_analysis_dir)
            try:
                os.system("ln -fs %s/%s_R*.fastq.gz %s/" % (in_dir, sample, sample_analysis_dir))
                reseq_sh = open(sample_analysis_dir + "/" + sample + ".sh", "w")
                reseq_sh.write(sh_head)
                reseq_sh.write("#! /bin/sh\n#BSUB -n 1\n#BSUB -oo %s.%%J.o\n#BSUB -e %s.%%J.e\n" % (sample, sample))
                reseq_sh.write("%s -i %s -o %s -r %s --qcCheck True\
 -M True -E True -R True -C True\
 --insertsize 306  --mismatch 4 --gapopen 5 --gapextention 20 \
 -A True --anodb %s --species hg19 --chrbed %s --readlen 100 \
 --exombed %s --nearbed %s --targetbed %s --basename %s"\
 %(Reseq_exe, sample_analysis_dir, sample_analysis_dir, Ref_fa, anodb, chrbed, exombed, nearbed, targetbed, sample))
                reseq_sh.close()
                os.system("sh %s&" % (sample_analysis_dir + "/" + sample + ".sh"))
                time.sleep(2)
            except:
                print >> sys.stderr, "sample %s is not in %s!" % (sample, in_dir)
                sys.exit(0)

###################################################
USAGE = """
python %s <rawdata_dir> <analysis_dir> <sample1..N>
    sample1..N: replaced by sample names in use

Version: V1.0
  Author: Ouyang Guojun <guojun.ouyang@sirna.cn>
  Date: Fri Jun 27 09:51:48 2014
        """ % sys.argv[0]
###################################################
if __name__ == '__main__':
    if len(sys.argv) > 3:
        main(sys.argv[1:])
    else:
        print USAGE
        sys.exit(0)
