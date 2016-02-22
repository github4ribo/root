#! /usr/bin/python
########################################################################
#Program: resequencing3.py
#Fuction: Main program of Exomes sequencing annalysis 
#Author: Ouyang Guojun
#Version: 1.1.0.20150918_Alpha
#Date: Fri Sep 18 16:50:44 2015
########################################################################
import re
import os
import sys
import glob
import time
import commands
import optparse

BIN = os.path.dirname(sys.argv[0])
BIN = os.path.abspath(BIN)

GATK            = BIN + '/GenomeAnalysisTK.jar'
VARSCAN         = BIN + '/VarScan.v2.3.3.jar'
ANNOVAR         = BIN + '/annotate_variation.pl'
FILTANNO        = BIN + '/variants_reduction.pl'
sumanno         = BIN + '/table_annovar.pl'
convert2annovar = BIN + '/convert2annovar.pl'
samtools        = BIN + "/samtools-0.1.18/samtools"
picardDir       = BIN + "/picard-tools-1.119/"
GATK_resource   = BIN + "/GATK_resource_bundle/"
circos          = BIN + "/circos/circos-0.63-4/bin/circos"

#databases
humandb         = BIN + '/humandb'
ratdb           = BIN + '/rn5'
mousedb         = BIN + '/mm10'
pigdb           = BIN + '/susScr3'
exon            = BIN + '/humandb/Trueseq/TruSeq_exome_targeted_regions.exome.hg19.bed'
exonplus200     = BIN + '/humandb/Trueseq/TruSeq_exome_targeted_regions.near.hg19.bed'
exon_flank      = BIN + '/humandb/Trueseq/TruSeq_exome_target_flanking.merged.bed'

sortCover       = BIN + '/sortcover.py'
vcfutils        = BIN + '/vcfutils.pl'
ReadMap         = BIN + '/ReadMap.R'
maploc          = BIN + '/maploc.py'
editPlot        = BIN + '/editplot.py'
coverDepth      = BIN + '/cover_depth.py'
svg_head        = BIN + "/svg.head"
Samtools        = BIN + "/samtools-0.1.18/samtools"
DataClean       = BIN + "/Qualitycheck/quality_check3.py"
coverexome      = BIN + "/coverexome.py"
easydraw        = BIN + "/easydraw"
sam2bam         = BIN + "/sam2bam.o"
variantcalling  = BIN + "/variantcalling.o"

AllJobids = list()
#-----------------general shells-----------------------------#
def makeShell(shell_name, shell_options, output_dir):
	"""make shell file from shell name and shell script"""
	outfile = open(output_dir+ '/' + shell_name, 'w')
	outfile.write('%s\n' % shell_options)
	outfile.close()

def getFile(path, regex):
	"""get target file path"""
        return glob.glob('%s/%s' %(path, regex))

def mkDir(path):
	"""make dir"""
	if not os.path.exists(path):
		os.mkdir(path)

def subjobs(jobs, depend_cond):
	"""sub jobs with dependence"""
	jobids = list()
	if depend_cond:
	    depend_cond = "&&".join(["done(%s)" % dep for dep in depend_cond])
	for job in jobs:
		if depend_cond:
			jobid = commands.getoutput("bsub -w \"%s\" < %s" % (depend_cond, job))
		else:
			jobid = commands.getoutput("bsub < %s" % job)
		jobids.append(re.split("<|>", jobid)[1])
	        AllJobids.append(re.split("<|>", jobid)[1])
		time.sleep(1)
	return jobids

def JobStat(jobids):
	"""detect jobs runing status"""
	jobstats = list()
	for jobid in jobids:
		jobstat = commands.getoutput("bjobs -l %s" % jobid)
		if re.search("Done successfully", jobstat, re.I):
			jobstats.append("done")
		elif re.search("exit", jobstat, re.I):
			os.system("bkill %s" % " ".join(AllJobids))
			sys.exit(0)
		else:
			jobstats.append("running")
	return jobstats

#------------------Mapping-------------------------#

def mapping(fq1, fq2, ref, mismatch, gap_open, gap_extention, insize, out_dir, basename):
	"""make alignment module"""
	out_aln = out_dir + '/02.Alignment'
	mkDir(out_aln)
	mapshell = ''
	fq1_name = re.split("\.fq|\.fastq", os.path.basename(fq1))[0]
	if fq2:fq2_name = re.split("\.fq|\.fastq", os.path.basename(fq2))[0]
	sai_1 = out_aln + '/' + fq1_name + ".sai"
	if fq2:sai_2 = out_aln + '/' + fq2_name + ".sai"
	sam = "%s/%s.sam" % (out_aln, basename)
	shell_aln_1 = 'bwa aln  -n %s -o %s -e %s -m 100000 -t 4 -q 4 -f %s %s %s' % (mismatch, gap_open, gap_extention, sai_1, ref, fq1)
	if fq2:
		shell_aln_2 = 'bwa aln  -n %s -o %s -e %s -m 100000 -t 4 -q 4 -f %s %s %s' % (mismatch, gap_open, gap_extention, sai_2, ref, fq2)
		shell_map = 'bwa sampe -r "@RG\\tID:%s\\tLB:%s\\tPL:illumina\\tSM:%s\\tPI:%s" -f %s %s %s %s %s %s\n' %(
				basename, basename, basename, insize, sam, ref, sai_1, sai_2,fq1, fq2)
	else:
		shell_map = 'bwa samse -r "@RG\\tID:%s\\tPL:illumina\\tSM:%s" -f %s %s %s %s\n' %(basename, basename, sam, ref, sai_1, fq1)
	
	shell_samtools = "%s %s %s %s %s %s %s/ucsc.hg19.fasta %s %s %s %s" % (sam2bam, samtools, picardDir, GATK, GATK_resource, sam, GATK_resource, exon, exonplus200, out_aln, basename)

	#------call Variants------------#
	out_var = out_dir  + '/03.Variants'
	mkDir(out_var)
	shell_cns = "%s %s %s/%s.uniq.rmdup.sort.bam %s/ucsc.hg19.fasta %s %s" % (variantcalling, GATK, out_aln, basename, GATK_resource, out_var, basename)

	if fq2:
		return shell_aln_1, shell_aln_2, shell_map, shell_samtools, shell_cns
	else:
		return shell_aln_1, shell_map, shell_samtools, shell_cns


def annovar(snpVcf, indelVcf, annodb, species, out_dir, basename):
	"""Variants annotation module"""
	out_anno = out_dir + '/04.Annotation'
	mkDir(out_anno)
	annovarInput = "%s/%s.annovar" % (out_anno, basename)
	out_anno = "%s/%s" % (out_anno, basename)
	covert2ano = '%s %s --includeinfo --allallele -format vcf4 >%s\n' %(
			convert2annovar, snpVcf, annovarInput)
	covert2ano += '%s %s --includeinfo --allallele -format vcf4 >>%s\n' %(
			convert2annovar, indelVcf, annovarInput)
	if species == 'hg19':
		#anovar = 'perl %s %s %s --buildver hg19 -protocol refGene,phastConsElements46way,genomicSuperDups,esp6500si_all,1000g2012apr_all,snp137,avsift,ljb_all -operation g,r,r,f,f,f,f,f -nastring NA --outfile %s\n' %(sumanno, annovarInput, annodb, out_anno)
		anovar = 'perl %s %s %s --buildver hg19 -protocol refGene,phastConsElements46way,genomicSuperDups,esp6500si_all,1000g2012apr_asn,snp138,avsift,ljb_all -operation g,r,r,f,f,f,f,f -nastring NA --outfile %s\n' %(sumanno, annovarInput, annodb, out_anno)
	#	anovar = 'perl %s %s %s --buildver hg19 -protocol refGene -operation g -nastring NA --outfile %s\n' %(sumanno, annovarInput, annodb, out_anno)
	elif species == 'rn5':
		anovar = 'perl %s %s %s --buildver rn5 -protocol refGene -operation g -nastring NA --outfile %s\n' %(sumanno, annovarInput, annodb, out_anno)
	elif species == 'mm10':
		anovar = 'perl %s %s %s --buildver mm10 -protocol refGene,snp138 -operation g,f -nastring NA --outfile %s\n'%(sumanno, annovarInput, annodb, out_anno)
	elif species == 'susScr3':
		anovar = 'perl %s %s %s --buildver susScr3 -protocol refGene -operation g -nastring NA --outfile %s\n' %(sumanno, annovarInput, annodb, out_anno)
	else:
		reftxt = getFile(annodb, '*refGene.txt')
		redfa = getFile(annodb, '*refGeneMrna.fa')
		refgff = getFile(annodb, '*gff')
		refffn = getFile(annodb, '*ffn')
		if len(reftxt) == len(redfa) and len(reftxt) == 1:
			anovar = 'perl %s %s %s --buildver %s -protocol refGene -operation g -nastring NA --outfile %s\n' %(sumanno, annovarInput, annodb, species, out_anno)
		elif len(refgff) > 0 and len(refgff) == len(reffn):
			gff2refGene(annodb, species)
			anovar = 'perl %s %s %s --buildver %s -protocol refGene -operation g -nastring NA --outfile %s\n' %(sumanno, annovarInput, annodb, species, out_anno)
		else:
			print >>sys.stderr, 'Annotation files:\n\t1.*.refGene.txt and *.refGeneMrna.fa files from UCSC;\n\t2.*.gff3 file and *.ffn\n'
			sys.exit(0)
	anoshell = covert2ano + anovar
	return anoshell

def readMap(rmdupbam, chrbed, windowsize, species, annodb, out_dir, basename):
	out_loc = out_dir + '/05.Readloc'
	mkDir(out_loc)
	locshell = ''
	rmdup  = os.path.basename(rmdupbam)
	readbed = out_loc + '/' +  rmdup.replace('.bam', '.bed')
	editbed = out_loc + '/' +  rmdup.replace('.bam', '.edit.bed')
	bam2bed = 'bedtools bamtobed -i %s >%s\n' %(rmdupbam, readbed)

	#-----chromosome location-----#
	chrname = os.path.basename(chrbed)
	winbed = out_loc + '/' + chrname.replace('.bed', '.win.bed')
	cover_win = out_loc + '/cover.' + rmdup.replace('.bam', '.win.txt')
	sort_cover = out_loc + '/sorted.cover.' + rmdup.replace('.bam', '.win.txt')
	shell_makewindow = 'bedtools makewindows -b %s -w %s >%s\n' %(chrbed, windowsize, winbed)
	os.system(shell_makewindow)

	# reads distribution plot
	edit_winbed = out_loc + '/edit_%s' %os.path.basename(winbed)
	#shell_cover = 'edtools coverage -abam %s -b %s >%s\n' %(rmdupbam, winbed, cover_win)
	shell_cover  = 'bedtools coverage -abam %s -b %s |cut -f 1-4 > %s\n' %(rmdupbam, winbed, cover_win)
	shell_plot  = '%s --cyto=%s/hg19_cytoBand.txt --depth2=%s --basename=%s/%s.readsmap\n' % (easydraw, humandb, cover_win, out_loc, basename)
	shell_readmap = bam2bed + shell_makewindow + shell_cover + shell_plot
	locshell += shell_readmap
	#-----Fuction location------#
	bedsub = "awk '{print $1\"\\t\"$2\"\\t\"$3\"\\t0\\t0\\t\"$6}' %s >%s\n" %(readbed, editbed)
	readano = '%s %s %s -buildver %s\n' %(ANNOVAR, editbed, annodb, species)
	varfuction = editbed + '.variant_function'
	dat = out_loc + '/dat.txt'
	staread = "awk '{print $1}' %s |sort |uniq -c > %s\n" %(varfuction, dat)
	editfile = out_loc + '/plot.dat'
	editplot = 'python %s %s %s\n' %(editPlot, dat, editfile)
	staloc =  bedsub + readano + staread + editplot
	locshell += staloc

	bamname = os.path.basename(rmdupbam)
	bamname = bamname.split('.')[1]
	Rstrip = open(out_loc + '/plotpie.R', 'w')
	script = 'library("plotrix")\n'
	script += "pdf('%s/readsmap_%s.pdf')\n"  %(out_loc, bamname)
	script += 'dat=read.table("%s")\n' %editfile
	script += 'ratio=sprintf("%.2f",100*dat[,2]/sum(dat[,2]))\n'
	script += 'ratio=paste(ratio,"%",sep="")\n'
	script += 'label=paste(dat[,1],ratio,sep="\\n")\n'
	script += 'pie3D(dat[,2],col=rainbow(10),main="", border="black", labels=label,font=2,labelcex=1,,explode=0.1,radius=0.95)\n'
	Rstrip.write('%s\n' %script)
	Rstrip.close()
	Rstrip = out_loc + '/plotpie.R'
	Rplot = 'Rscript %s\n' %Rstrip
	locshell += Rplot
	pdf2png = "/usr/bin/mogrify -density 250 -format png %s/readsmap_%s.pdf\n" % (out_loc, bamname) #added by Ouyang Guojun [2014-07-01]
	locshell += pdf2png #added by Ouyang Guojun [2014-07-01]
	return locshell

def Circos(rmdupbam, chrbed, filename, exonfunc, out_dir):
	out_dir = out_dir + '/06.Circos'
	mkDir(out_dir)
	shell_circos = ''
	chrname=os.path.basename(chrbed)
	rmdup = os.path.basename(rmdupbam)

	# make window
	winbed = out_dir + '/' + chrname.replace('.bed', '.win2000000.bed')
	shell = 'bedtools makewindows -b %s -w 2000000 >%s\n' %(chrbed, winbed)

	# edit windows bed
	edit_winbed = out_dir + '/edit_%s' %os.path.basename(winbed)
	shell += 'awk \'{print $1"\\t"$2-1"\\t"$3-2}\' %s >%s\n' %(winbed, edit_winbed)

	# coverage
	readmap = out_dir + '/cover.' + rmdup.replace('.bam', '.win.txt')
	shell += 'bedtools coverage -abam %s -b %s >%s\n' %(rmdupbam, edit_winbed, readmap)
	coverbed = out_dir + '/coverage.txt'
	shell += 'awk \'{sub("chr", "hs")}$4!=0{print $1" "$2" "$3" "$4}\' %s >%s\n' %(readmap, coverbed)
	confile = out_dir + '/%s.conf' % filename

	shell += "awk -F \"\\t\" 'split($NF, arr, \":\")split(arr[1], arr2, \"/\"){if (arr2[1]==arr2[2]) {print $4\"\\t\"$5\"\\t\"$6}}' %s > %s/homozygous.bed\n" % (exonfunc, out_dir)
	shell += "awk -F \"\\t\" 'split($NF, arr, \":\")split(arr[1], arr2, \"/\"){if (arr2[1]!=arr2[2]) {print $4\"\\t\"$5\"\\t\"$6}}' %s > %s/heterozygous.bed\n" % (exonfunc, out_dir)
	shell += "awk -F \"\\t\" '$2~/synonymous/ {print $4\"\\t\"$5\"\\t\"$6}' %s > %s/synonymous.bed\n" % (exonfunc, out_dir)
	shell += "awk -F \"\\t\" '$2~/nonsynonymous/ {print $4\"\\t\"$5\"\\t\"$6}' %s > %s/nonsynonymous.bed\n" % (exonfunc, out_dir)
	shell += "awk -F \"\\t\" '$2~/insertion/ {print $4\"\\t\"$5\"\\t\"$6}' %s > %s/insert.bed\n" % (exonfunc, out_dir)
	shell += "awk -F \"\\t\" '$2~/deletion/ {print $4\"\\t\"$5\"\\t\"$6}' %s > %s/delete.bed\n" % (exonfunc, out_dir)
	shell += "awk -F \"\\t\" '$2~/stopgain/ {print $4\"\\t\"$5\"\\t\"$6}' %s > %s/stopgain.bed\n" % (exonfunc, out_dir)
	shell += "awk -F \"\\t\" '$2~/stoploss/ {print $4\"\\t\"$5\"\\t\"$6}' %s > %s/stoploss.bed\n" % (exonfunc, out_dir)
	
	# beds coverage
	for bed in [out_dir + '/homozygous.bed', out_dir + '/heterozygous.bed', out_dir + '/insert.bed', out_dir + '/delete.bed', out_dir + '/nonsynonymous.bed', out_dir + '/synonymous.bed', out_dir + '/stoploss.bed', out_dir + '/stopgain.bed']:
		bedout = bed.replace('bed', 'txt')
		shell += 'bedtools coverage -a %s -b %s |awk  \'{sub("chr", "hs")}$4!=0{print $1" "$2" "$3" "$4}\' - >%s\n' %(bed, edit_winbed, bedout)
	conf = """
<colors>
<<include /RiboBio/home/bi.ouyangguojun/pipeline/exome_pipeline/circos/circos-0.63-4/etc/colors.conf>>
</colors>

<fonts>
<<include /RiboBio/home/bi.ouyangguojun/pipeline/exome_pipeline/circos/circos-0.63-4/etc/fonts.conf>>
</fonts>

<<include /RiboBio/home/bi.ouyangguojun/pipeline/exome_pipeline/ideogram.conf>>
<<include /RiboBio/home/bi.ouyangguojun/pipeline/exome_pipeline/ticks.conf>>

karyotype = /RiboBio/home/bi.ouyangguojun/pipeline/exome_pipeline/circos/circos-0.63-4/data/karyotype/karyotype.human.txt

<image>
dir = %s
file  = %s.variants.svg
radius         = 1500p
background     = white
angle_offset   = -62
24bit = yes
</image>

<<include /RiboBio/home/bi.ouyangguojun/pipeline/exome_pipeline/circos/circos-0.63-4/etc/colors_fonts_patterns.conf>>

<<include /RiboBio/home/bi.ouyangguojun/pipeline/exome_pipeline/circos/circos-0.63-4/etc/housekeeping.conf>>

chromosomes_units           = 5000000

<plots>

<plot>
type = histogram
file    = %s/coverage.txt
r1   = 0.97r
r0   = 0.85r
stroke_type = bin
thickness   = 2
color = black
orientation = out
fill_color = lgrey
extend_bin = no
</plot>

<plot>
type             = scatter
file             = %s/insert.txt
r0               = 0.97r
r1               = 0.95r
fill_color  = green
scale_log_base   = 2
</plot>

<plot>
type             = scatter
file             = %s/delete.txt
r0               = 0.95r
r1               = 0.93r
fill_color  = vdgreen
scale_log_base   = 2
</plot>

<plot>
type = histogram
stroke_type = bin
thickness   = 2
color = black
file    = %s/homozygous.txt
r1   = 0.85r
r0   = 0.78r
orientation = out
fill_color = orange
extend_bin = no
</plot>

<plot>
type = histogram
stroke_type = bin
thickness   = 2
color = black
file    = %s/heterozygous.txt
r1   = 0.78r
r0   = 0.71r
orientation = in
fill_color = yellow
extend_bin = no
</plot>

<plot>
type = scatter
file    = %s/synonymous.txt
r1   = 0.7r
r0   = 0.55r
thickness = 2
color = red
extend_bin = no
</plot>

<plot>
type = scatter
file    = %s/nonsynonymous.txt
r1   = 0.54r
r0   = 0.39r
orientation = out
color = purple
extend_bin = no
</plot>

<plot>
type = scatter
file    = %s/stopgain.txt
r1   = 0.65r
r0   = 0.65r
orientation = out
fill_color = blue
extend_bin = no
</plot>

<plot>
type = scatter
file    = %s/stoploss.txt
r1   = 0.68r
r0   = 0.68r
orientation = out
fill_color = black
extend_bin = no
</plot>
</plots>
	""" %(out_dir, filename, out_dir, out_dir, out_dir, out_dir, out_dir, out_dir, out_dir, out_dir, out_dir)
	outconf = open(out_dir + '/%s.conf' %filename, 'w')
	outconf.write('%s\n' %conf)
	outconf.close()
	
	shell += '%s -conf %s\n' % (circos, confile)
	return shell


def gff2refGene(annodb, species):
	out_1 = open(annodb + '/' + species + '_refGene.txt', 'w')
	out_2 = open(annodb + '/' + species + '_refGeneMrna.fa', 'w')
	out_3 = open(annodb + '/check.txt', 'w')
	gffs = getFile(annodb, '*gff')
	gfflist = sorted(gffs)
	ffns = getFile(annodb, '*ffn')
	ffnlist = sorted(ffns)
	num = 0

	if len(gfflist) == 0 or len(ffnlist) == 0:
		print >>sys.stderr, 'Annotation files:\n\t1.speices_refGene.txt and species_refGeneMrna.fa files from UCSC;\n\t2.*.gff3 file and *.ffn\n'
		sys.exit(0)
	elif len(gfflist) != len(ffnlist):
		print >>sys.stderr, 'The gff files and ffn files must be equivalent!!!'
		sys.exit(0)

	for i in range(len(gfflist)):
		gname = {}
		gff = open(gfflist[i], 'r')
		ffn = open(ffnlist[i], 'r')
		name = ''
		cdsnum = 0
		cdssta = ''
		cdsend = ''
		exonframs = ''
		cdsstaloc = 0
		cdsendloc = 0
		check = 0

		for line in gff:
			line = line.rstrip()
			if len(line.split()) > 2:
				if line.split()[2] == 'region':
					if re.search('chromosome=', line):
						chrid = line.split('chromosome=')[1].split(';')[0]
						if not re.match('chr', chrid):
							chrid = 'chr%s' %chrid
					else:
						chrid = line.split()[0]

				elif line.split()[2] == 'gene':
					if num == 0:
						geneid = line.split('GeneID:')[1].split(';')[0]
						name2 = line.split('Name=')[1].split(';')[0]
						strand = line.split()[6]
					else:
						if (not cdsstaloc == 0) and (not cdsendloc == 0):
							if check == 0:
								txsta = cdsstaloc
								txend = cdsendloc
							gene = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tunk\tunk\t%s' %(
									num, name, chrid, strand, txsta-1, txend, cdsstaloc-1, cdsendloc,cdsnum, cdssta-1, cdsend, geneid, name2, exonframs)
							if not len(gene.split('\t')) == 16:
								continue
							group = gene.split('\t')
							for i in gene.split('\t'):
								if len(i) == 0:
									continue
							out_1.write('%s\n' %gene)
						cdssta = ''
						cdsend = ''
						cdsnum = 0
						cdsstaloc = 0
						cdsendloc = 0
						geneid = line.split('GeneID:')[1].split(';')[0]
						name2 = line.split('Name=')[1].split(';')[0]
						strand = line.split()[6]
						exonframs = ''
					num += 1
				elif line.split()[2] == 'mRNA':
					txsta = line.split()[3]
					txend = line.split()[4]
					check = 1
				elif line.split()[2] == 'CDS':
					cdssta += '%s,' % line.split()[3]
					cdsend += '%s,' % line.split()[4]
					cdsnum += 1
					if cdsstaloc == 0:
						cdsstaloc = int(line.split()[3])
						cdsendloc = int(line.split()[4])
					else:
						if int(line.split()[4]) > cdsendloc:
							cdsendloc = int(line.split()[4])
					exonframs += '%s,' %line.split()[7]
					name = line.split('Name=')[1].split(';')[0]
					genestrand = line.split()[6]
					if genestrand == '+':
						mark = '%s-%s' %(cdsstaloc, cdsendloc)
					else:
						mark = '%s-%s' %(cdsendloc, cdsstaloc)
					if not mark in gname:
						gname[mark] = name
				else:
					continue
		if (not cdsstaloc == 0) and (not cdsendloc == 0):
			out_1.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tunk\tunk\t%s\n' %(
				num, name, chrid, strand, txsta-1, txend, cdsstaloc-1, cdsendloc-1, cdsnum, cdssta-1, cdsend, geneid, name2, exonframs))
		gff.close()
		
		isprint = 0
		seq = ''
		for line in ffn:
			line = line.rstrip()
			if re.		tch('>', line):
				locid = line.split()[0].split(':')[1]
				pattern = re.compile('\d+-\d+')
				if pattern.search(line):
					geneid = pattern.search(line).group()
				else:
					print line
					sys.exit(0)
				isprint = 0
				if len(seq) > 0:
					out_2.write('%s\n' %seq)
				if geneid in gname:
					isprint = 1
					out_2.write('>%s\n' % gname[geneid])
				else:
					out_3.write('%s\n' % line)
				seq = ''
			else:
				if isprint == 1:
					seq += line
		out_2.write('%s\n' %seq)
		ffn.close()

	out_1.close()
	out_2.close()
	out_3.close()

def Stavar(outAnnovar, out_dir, basename):
	"""statistics variants"""
	annofile = "%s/%s.hg19_multianno.txt" % (outAnnovar, basename)
	varfun = "%s/%s.refGene.variant_function" % (outAnnovar, basename)

	#-------------------statistic variant function------------------------#
	var_indel = "%s/%s.indel.txt" % (out_dir, basename)
	var_snp = "%s/%s.snp.txt" % (out_dir, basename)

	os.system('head -1 %s >%s' %(annofile, var_indel))
	os.system('head -1 %s >%s' %(annofile, var_snp))
	shell_split_variant = "awk -F \"\\t\" '$1!~/Chr/&&($4~/-/||$5~/-/||length($4)!=length($5))' %s >>%s\n" %(annofile, var_indel)
	shell_split_variant += "awk -F \"\\t\" '$1!~/Chr/&&$4!~/-/&&$5!~/-/&&length($4)==length($5)' %s >>%s\n" %(annofile, var_snp)
	sta_var_indel = out_dir + '/staloc.' + os.path.basename(var_indel)
	sta_var_snp = out_dir + '/staloc.' + os.path.basename(var_snp)
	shell_sta_variant = "awk '$0!~/Chr/{if($6~/ncRNA/){split($6,array, \"_\");print array[1]}else{print $6}}' %s|sort|uniq -c - > %s\n" %(var_indel, sta_var_indel) 
	shell_sta_variant += "awk '$0!~/Chr/{if($6~/ncRNA/){split($6,array, \"_\");print array[1]}else{print $6}}' %s|sort|uniq -c - > %s\n" %(var_snp, sta_var_snp)
	#------------------statistic exome variant function-------------------#
	exon_var_indel = out_dir + '/stafun.' + os.path.basename(var_indel)
	exon_var_snp = out_dir + '/stafun.' + os.path.basename(var_snp)
	shell_sta_exon = "awk -F \"\\t\" '$0!~/Chr/&&$8!~/NA|unknown/{print $8}' %s|sort|uniq -c  - > %s\n" %(var_indel, exon_var_indel)
	shell_sta_exon += "awk -F \"\\t\" '$0!~/Chr/&&$8!~/NA|unknown/{print $8}' %s|sort|uniq -c - > %s\n" %(var_snp, exon_var_snp)

	sta_exon_syno_nosyno = out_dir + '/sta.syno.nosyno.' + os.path.basename(varfun) + '.txt' 
	shell_sta_syno_nosyno = "awk  -F \"\\t\" '$1!~/Chr/&&$6!~/INDEL/ {print $NF}' %s|awk -F \":\" '{print $1}'|sort|uniq -c - >%s\n" %(varfun,sta_exon_syno_nosyno)
	shell_stavar = shell_split_variant + "\n" + shell_sta_variant + "\n" + shell_sta_exon +"\n" + shell_sta_syno_nosyno + "\n"
	return shell_stavar

def editano(annodir, basename):
	"""Edits snp.txt and indel.txt; gets snp.xls and indel.xls"""
	anno = {}
	anofile = '%s/%s.refGene.variant_function' % (annodir, basename)
	snpfile = "%s/%s.snp.txt" % (annodir, basename)
	indfile = "%s/%s.indel.txt" % (annodir, basename)

	in_1 = open(anofile, 'r')
	in_2 = open(snpfile, 'r')
	in_3 = open(indfile, 'r')
	out_1 = open("%s/%s.snp.xls" % (annodir, basename), 'w')
	out_2 = open("%s/%s.indel.xls" % (annodir, basename), 'w')

	for line in in_1:
		line = line.rstrip()
		line = line.replace(' SNV', '')
		chrid = line.split('\t')[2]
		loc = line.split('\t')[3]
		ref = line.split('\t')[5]
		var = line.split('\t')[6]
		depth = line.split('DP=')[1].split(';')[0]
		heter = line.split('\t')[-1].split(':')[0]
		if heter == '1/1':
			heter = 'Hom'
		elif heter == '0/1' or heter == '1/2':
			heter = 'Het'
		else:
			print heter
			sys.exit(0)
		mark = '%s\t%s\t%s\t%s' %(chrid, loc, ref, var)
		if mark not in anno:
			anno[mark] = [depth, heter]
		else:
			print mark
			sys.exit(0)
	in_1.close()

	for line in in_2:
		line = line.rstrip()
		line = line.replace(' SNV', '')
		if re.search('Func.refGene', line):
			out_1.write('%s\tDepth\tZygosity\n' %line)
		else:
			chrid = line.split('\t')[0]
			locid = line.split('\t')[1]
			ref = line.split('\t')[3]
			var = line.split('\t')[4]
			mark_2 = '%s\t%s\t%s\t%s' %(chrid, locid, ref, var)
			if mark_2 in anno:
				out_1.write('%s\t%s\t%s\n' %(line, anno[mark_2][0], anno[mark_2][1]))
			else:
				print line
				sys.exit(0)
	in_2.close()
	out_1.close()

	for line in in_3:
		line = line.rstrip()
		if re.search('Func.refGene', line):
			out_2.write('%s\tDepth\tZygosity\n' %('\t'.join(line.split('\t')[:14]))) 
		else:
			chrid = line.split('\t')[0]
			locid = line.split('\t')[1]
			ref = line.split('\t')[3]
			var = line.split('\t')[4]
			mark_2 = '%s\t%s\t%s\t%s' %(chrid, locid, ref, var)
			edit_line = '\t'.join(line.split('\t')[:14])
			if mark_2 in anno:
				out_2.write('%s\t%s\t%s\n' %(edit_line, anno[mark_2][0], anno[mark_2][1]))
			else:
				print line
				sys.exit(0)
	in_3.close()
	out_2.close()

#################################################################################
def main(Args):
	op = optparse.OptionParser()
	op.add_option('-i', '--indir', help='fastq files path')
	op.add_option('-s', '--label', help='sample name')
	op.add_option('', '--qcCheck', help='quality check[True]', default="True")
	op.add_option('', '--qcArgs', help='quality_check args.[-t 2 -l AGATCGGAAGAGC -r AGATCGGAAGAGC]', default="-t 2 -l AGATCGGAAGAGC -r AGATCGGAAGAGC")
	op.add_option('-o', '--outdir', help='output path.[./]', default="./")
	op.add_option('-r', '--ref', help='fasta to use as reference')
	op.add_option('-M', '--map', help='read mapping [true]', default = 'true')
	op.add_option('','--insertsize', help = 'insert size [400]', default = '400')
	op.add_option('', '--mismatch', help = 'max #diff (int) or missing prob under 0.02 err rate (float) [0.04]', default = 0.04)
	op.add_option('', '--gapopen', help = 'maximum number or fraction of gap opens [1]', default = '1')
	op.add_option('', '--gapextention', help = 'maximum number of gap extensions, -1 for disabling long gaps [-1]', default = -1)
	op.add_option('-E', '--exome', help = 'coverage of exome and near 200bp region [false]', default ='false')
	op.add_option('', '--readlen', help = 'read length [100]', default = 100)
	op.add_option('','--exombed', help = 'exome bed file')
	op.add_option('','--nearbed', help = 'flanking region bed file')
	op.add_option('','--targetbed', help = 'exome and flanking region bed file')
	op.add_option('-A', '--annotation', help = 'variants annotaion [false]', default = 'false')
	op.add_option('', '--anodb', help = 'needed by -A and -R options, annotation database')
	op.add_option('', '--species', help = 'needed by -A and -R options, hg19,mm10,rn5, other species')
	op.add_option('-R', '--readloc', help = 'read mapping plot [false]', default = 'false')
	op.add_option('', '--chrbed', help = 'needed by  "-R" option, bed file for chromosome.')
	op.add_option('', '--winsize', help = 'needed with "-R" option, plot windowsize [500,000]', default = 500000)
	op.add_option('-C', '--circos', help = 'variants plot [false], just for hg19', default = 'false')
	op.add_option('-b', '--basename', help = 'basename.[false]', default = 'false')

	if not Args or Args[0] == '--help' or Args[0] == '-h':
		print op.print_help()
		sys.exit(0)
	opts, args = op.parse_args()
	out_dir = opts.outdir
	mkDir(out_dir)
	sample = opts.label
	fqdir = opts.indir
        basename = opts.basename
	fqlist = getFile(fqdir, '*_R*.fastq.gz')
	if not fqlist:
		fqlist = getFile(fqdir, '*read*_Clean.fastq.gz')
	if not fqlist or len(fqlist)>2:
		print >>sys.stderr, 'There must be one or two paired fastq files for qcheck'
		sys.exit(0)
	fq1 = sorted(fqlist)[0]
	fq2 = "" if len(fqlist)==1 else sorted(fqlist)[1]
        if re.search("false", basename, re.I):
            basename = re.split("_R1|\.read1_", os.path.basename(fq1), re.I)[0]

	jobids = list()
	shellHead = "#!/bin/bash\n#BSUB -n 5\n#BSUB -R \"span[ptile=5]\"\n#BSUB -o %s/01.qcCheck.%%J.o\n#BSUB -e %s/01.qcCheck.%%J.e\n" % (out_dir, out_dir)
	if re.search('true', opts.qcCheck, re.I):
		if fq2:
			qcShell = shellHead + "python %s -1 %s -2 %s %s -o %s" % (DataClean, fq1, fq2, opts.qcArgs, out_dir)
		else:
			qcShell = shellHead + "python %s -1 %s %s -o %s" % (DataClean, fq1, opts.qcArgs, out_dir)
		makeShell("01.qcCheck.sh", qcShell, out_dir)
		jobids = subjobs([out_dir+"/01.qcCheck.sh"], [])
		while True:
		    jobstats = JobStat(jobids)
		    if jobstats.count("running"):
		        time.sleep(2)
		    else:
			break
		try:
			fq1 = getFile(out_dir + "/01.DataCleaning/", "*.read1_Clean.fastq.gz")[0]
			if fq2: fq2 = getFile(out_dir + "/01.DataCleaning/", "*.read2_Clean.fastq.gz")[0]
		except:
			print >> sys.stderr, "*.read*_Clean.fastq.gz not exists!!!"

	if re.search('true', opts.map, re.I):
		if not opts.ref:
			print '-r option is needed!!!\n'
			sys.exit(0)
		ref = opts.ref
		mis = opts.mismatch
		gapo = opts.gapopen
		gape = opts.gapextention
		insize = opts.insertsize
		depend_cond = jobids
		if fq2:
			(shell_aln_1, shell_aln_2, shell_map, shell_samtools, shell_cns) = mapping(fq1, fq2, ref, mis, gapo, gape, insize, out_dir, basename)
		else:
			(shell_aln_1, shell_map, shell_samtools, shell_cns) = mapping(fq1, fq2, ref, mis, gapo, gape, insize, out_dir, basename)
		shellHead = "#!/bin/bash\n#BSUB -n 4\n#BSUB -R \"span[ptile=4]\"\n#BSUB -o %s/02.sai1.sh.%%J.o\n#BSUB -e %s/02.sai1.sh.%%J.e\n" % (out_dir, out_dir)
		makeShell("02.sai1.sh", shellHead+shell_aln_1, out_dir)
		shellHead = "#!/bin/bash\n#BSUB -n 4\n#BSUB -R \"span[ptile=4]\"\n#BSUB -o %s/02.sai2.sh.%%J.o\n#BSUB -e %s/02.sai2.sh.%%J.e\n" % (out_dir, out_dir)
		jobids = subjobs([out_dir+"/02.sai1.sh"], depend_cond)
		if fq2:
			makeShell("02.sai2.sh", shellHead+shell_aln_2, out_dir)
			jobids.extend(subjobs([out_dir+"/02.sai2.sh"], depend_cond))
		shellHead = "#!/bin/bash\n#BSUB -n 8\n#BSUB -R \"span[ptile=8]\"\n#BSUB -o %s/02.mapping.sh.%%J.o\n#BSUB -e %s/02.mapping.sh.%%J.e\n" % (out_dir, out_dir)
		makeShell("02.mapping.sh", shellHead+shell_map, out_dir)
		shellHead = "#!/bin/bash\n#BSUB -n 4\n#BSUB -o %s/02.sam2bam.sh.%%J.o\n#BSUB -e %s/02.sam2bam.sh.%%J.e\n" % (out_dir, out_dir)
		makeShell("02.sam2bam.sh", shellHead+shell_samtools, out_dir)
		shellHead = "#!/bin/bash\n#BSUB -n 4\n#BSUB -o %s/02.cns.sh.%%J.o\n#BSUB -e %s/02.cns.sh.%%J.e\n" % (out_dir, out_dir)
		makeShell("02.cns.sh", shellHead+shell_cns, out_dir)
		jobids = subjobs([out_dir+"/02.mapping.sh"], jobids)
		jobids = subjobs([out_dir+"/02.sam2bam.sh"], jobids)
		jobids = subjobs([out_dir+"/02.cns.sh"], jobids)
		while True:
		    jobstats =  JobStat(jobids)
		    if jobstats.count("running"):
		        time.sleep(2)
		    else:
			break
	if re.search('true', opts.exome, re.I):
		exome_bed = opts.exombed
		near_bed = opts.nearbed
		target_bed = opts.targetbed
		if not target_bed or (not near_bed):
			print 'The exomebed file and nearbed file are needed for exome option!!!'
			sys.exit(0)
		out_cover = out_dir + '/07.CoverTarget'
		out_map = out_dir + '/02.Alignment'
		mkDir(out_cover)
		bamfile = "%s/%s.sort.bam" % (out_map, basename)
		seq_len=opts.readlen
		cover_exon = out_cover + '/cover.exon.' + os.path.basename(bamfile).replace('bam', 'txt')
		cover_near = out_cover + '/cover.near200.' + os.path.basename(bamfile).replace('bam', 'txt')
		cover_exon_plus_near = out_cover + '/cover.target.' + os.path.basename(bamfile).replace('bam', 'txt')
		shell_cover_exome = 'bedtools coverage -abam %s -b %s >%s\n' %(bamfile, exome_bed, cover_exon)
		shell_cover_near = 'bedtools coverage -abam %s -b %s >%s\n' %(bamfile, near_bed, cover_near)
		shell_cover_target = 'bedtools coverage -abam %s -b %s >%s\n' %(bamfile,target_bed, cover_exon_plus_near)
		shellHead = "#!/bin/bash\n#BSUB -n 1\n#BSUB -o %s/07.CoverTarget.exome.sh.%%J.o\n#BSUB -e %s/07.CoverTarget.exome.sh.%%J.e\n" % (out_dir, out_dir)
		makeShell("07.CoverTarget.exome.sh", shellHead+shell_cover_exome, out_dir)
		shellHead = "#!/bin/bash\n#BSUB -n 1\n#BSUB -o %s/07.CoverTarget.near.sh.%%J.o\n#BSUB -e %s/07.CoverTarget.near.sh.%%J.e\n" % (out_dir, out_dir)
		makeShell("07.CoverTarget.near.sh", shellHead+shell_cover_near, out_dir)
		shellHead = "#!/bin/bash\n#BSUB -n 1\n#BSUB -o %s/07.CoverTarget.target.sh.%%J.o\n#BSUB -e %s/07.CoverTarget.target.sh.%%J.e\n" % (out_dir, out_dir)
		makeShell("07.CoverTarget.target.sh", shellHead+shell_cover_target, out_dir)
		jobids = subjobs([out_dir+"/07.CoverTarget.exome.sh", out_dir+"/07.CoverTarget.near.sh", out_dir+"/07.CoverTarget.target.sh"], jobids)
		shell = "python %s %s %s %s %s %s %s" % (coverexome, bamfile, seq_len, exome_bed, near_bed, target_bed,out_cover)
		shellHead = "#!/bin/bash\n#BSUB -n 1\n#BSUB -o %s/07.CoverTarget.coverexome.sh.%%J.o\n#BSUB -e %s/07.CoverTarget.coverexome.sh.%%J.e\n" % (out_dir, out_dir)
		makeShell("07.CoverTarget.coverexome.sh", shellHead+shell, out_dir)
		jobids = subjobs([out_dir+"/07.CoverTarget.coverexome.sh"], jobids)
		while True:
		    jobstats =  JobStat(jobids)
		    if jobstats.count("running"):
		        time.sleep(2)
		    else:
			break

	if re.search('true', opts.readloc, re.I):
		out_map = out_dir + '/02.Alignment'
		chr_bed = opts.chrbed
		window = opts.winsize
		spe = opts.species
		annotadb = opts.anodb
		if (not opts.chrbed) or (not opts.species) or (not opts.anodb):
			 print >>sys.stderr, 'if "readloc" is "true", "chrbed","species" and "anodb" options are needed!!!'
			 sys.exit(0)
		rmdupbam = "%s/%s.rmdup.sort.bam" % (out_map, basename)
		if rmdupbam:
			shell = readMap(rmdupbam, chr_bed, window, spe, annotadb, out_dir, basename)
			shellHead = "#!/bin/bash\n#BSUB -n 1\n#BSUB -o %s/05.Readloc.sh.%%J.o\n#BSUB -e %s/05.Readloc.sh.%%J.e\n" % (out_dir, out_dir)
			makeShell("05.Readloc.sh", shellHead+shell, out_dir)
		else:
			print >>sys.stderr, ('uniq.rmdup*bam not exists!!! %s\n' %commands.getoutput('date'))
			sys.exit(0)
		jobids = subjobs([out_dir+"/05.Readloc.sh"], jobids)
		while True:
		    jobstats =  JobStat(jobids)
		    if jobstats.count("running"):
		        time.sleep(2)
		    else:
			break

	if re.search('true', opts.annotation, re.I):
		out_vcf = out_dir + '/03.Variants' 
		mkDir(out_vcf)
		snpVcf = "%s/%s.final.snp.vcf" % (out_vcf, basename)
		indelVcf = "%s/%s.final.indel.vcf" % (out_vcf, basename)
		annodb = opts.anodb
		spe = opts.species
		if (not opts.anodb) or (not opts.species):
			print 'anodb and species options are needed!!'
			sys.exit(0)
		shellHead = "#!/bin/bash\n#BSUB -o %s/04.Annotation.sh.%%J.o\n#BSUB -e %s/04.Annotation.sh.%%J.e\n" % (out_dir, out_dir)
		makeShell("04.Annotation.sh", shellHead + annovar(snpVcf, indelVcf, annodb, spe, out_dir, basename), out_dir)
		jobids = subjobs([out_dir+"/04.Annotation.sh"], jobids)
		while True:
		    jobstats =  JobStat(jobids)
		    if jobstats.count("running"):
		        time.sleep(2)
		    else:
			break
		out_anovar = out_dir + '/04.Annotation'
		shell_stavar = Stavar(out_anovar, out_anovar, basename)
		shellHead = "#!/bin/bash\n#BSUB -n 1\n#BSUB -o %s/04.Stavar.sh.%%J.o\n#BSUB -e %s/04.Stavar.sh.%%J.e\n" % (out_dir, out_dir)
		makeShell("04.Stavar.sh", shellHead+shell_stavar, out_dir)
		jobids = subjobs([out_dir+"/04.Stavar.sh"], jobids)
		while True:
		    jobstats =  JobStat(jobids)
		    if jobstats.count("running"):
		        time.sleep(2)
		    else:
			break
		editano(out_anovar, basename)

	if re.search('true', opts.circos, re.I):
		outano = out_dir + '/04.Annotation'
		out_map = out_dir + '/02.Alignment'
		chr_bed = opts.chrbed
		if getFile(out_map, '*uniq.rmdup*bam') and chr_bed and getFile(outano, '*exonic_variant_function'):
			rmdupbam =  getFile(out_map, '*uniq.rmdup*bam')[0]
			filename = os.path.basename(rmdupbam).split('.')[0]
			if not getFile(outano, '*exonic_variant_function'):
				print 'their is no *.exonic_variant_function file'
				sys.exit(0)
			exonfunc = getFile(outano, '*exonic_variant_function')[0]
			shell_circos = Circos(rmdupbam, chr_bed, filename, exonfunc, out_dir)
			shellHead = "#!/bin/bash\n#BSUB -n 1\n#BSUB -o %s/06.Circos.sh.%%J.o\n#BSUB -e %s/06.Circos.sh.%%J.e\nsource /etc/profile.d/apps.sh\n" % (out_dir, out_dir)
			makeShell("06.Circos.sh", shellHead+shell_circos, out_dir)
			jobids = subjobs([out_dir+"/06.Circos.sh"], jobids)
			while True:
		    		jobstats =  JobStat(jobids)
		    		if jobstats.count("running"):
		        		time.sleep(2)
		    		else:
					break
			svg_file = "%s/%s.variants.svg" % (out_dir+"/06.Circos", filename)
			FH = open(svg_file)
			svg_raw_lines = FH.readlines()
			FH.close()
			head = open(svg_head).readlines()
			out_svg = open(svg_file, "w")
			out_svg.write("".join(head))
			out_svg.write("".join(svg_raw_lines[3:]))
			out_svg.close()
			os.system("mogrify -format png %s\n" % svg_file)
		else:
			print >> sys.stderr, ('circos plot location error!!! %s\n' %commands.getoutput('date'))
			

if __name__ == '__main__':
	try:
		main(sys.argv[1:])
	except KeyboardInterrupt:
		os.system("bkill %s" % " ".join(AllJobids))
		sys.exit(0)
