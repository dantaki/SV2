#!/usr/bin/env python2
'''
Copyright <2017> <Danny Antaki>

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
'''
__version__='1.4.3'
from sv2_backend import make_dir,rand_id,report_time,slash_check
from Bam import bam_init
from sv2Config import Config
from FeatureExtraction import extract_feats
from Output import sv2_train_output
from Ped import ped_init
from Preprocess import preprocess
from Snv import snv_init
from SV import sv_init
from argparse import RawTextHelpFormatter
from glob import glob
from time import time
import shutil,sys,os,argparse
splash='                       ____\n  _____________   ___ |___ \\\n /   _____/\   \ /   // ___/\n \_____  \  \   Y   //_____)\n /        \  \     /\n/_________/   \___/\n'
__useage___="""
Support Vector Structural Variation Genotyper Train: generate a training set
Version {}    Author: Danny Antaki <dantaki at ucsd dot edu>    github.com/dantaki/SV2

  sv2 -i <bam ...> -v <vcf ...> -b <bed ...> -snv <snv vcf ...> -p <ped ...> [OPTIONS]

input arguments: github.com/dantaki/SV2/wiki/Options#input-arguments

  Note: input arguments can take multiple files, separated by space <-i BAM1 BAM2 ...>

  -i, -bam    ...     bam file(s)
  -v, -vcf    ...     vcf files(s) of SVs
  -b, -bed    ...     bed files(s) of SVs
  -snv        ...     snv vcf files(s), must be bgzipped and tabixed
  -p, -ped    ...     ped files(s)

genotype arguments: github.com/dantaki/SV2/wiki/Options#genotype-arguments

  -g, -genome         reference genome build <hg19, hg38, mm10> [default: hg19]
  -pcrfree            GC content normalization for pcr free sequences
  -M                  bwa mem -M compatibility, split-reads flagged as secondary instead of supplementary

  -pre                preprocessing output directory, skips preprocessing
  -feats              feature extraction output directory, skips feature extraction

optional arguments:
 
  -L, -log            log file for standard error messages [default: STDOUT]
  -T, -tmp-dir        directory for temporary files [default: working directory]
  -s, -seed           random seed for preprocessing shuffling [default: 42]
  -o, -out            output prefix [default: sv2_training_features]
  -O, -odir           output path, location for sv2 output directories [default: working directory]
  -h, -help           show this message and exit
""".format(__version__)

def main():
	init_time = int(time())
 	parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter,usage=splash.replace("       ","",1)+__useage___,add_help=False)
	inArgs,genoArgs,optArgs = parser.add_argument_group('input arguments'),parser.add_argument_group('genotype arguments'),parser.add_argument_group('optional arguments')
	inArgs.add_argument('-i','-bam',type=str,default=None,nargs='*')
	inArgs.add_argument('-b','-bed',type=str,default=None,nargs='*')
	inArgs.add_argument('-v','-vcf',type=str,default=None,nargs='*')
	inArgs.add_argument('-snv',type=str,default=None,nargs='*')
	inArgs.add_argument('-p','-ped',type=str,default=None,nargs='*')
	genoArgs.add_argument('-g','-genome',required=False,default='hg19',type=str)
	genoArgs.add_argument('-pcrfree',required=False,default=False,action="store_true")
	genoArgs.add_argument('-M',default=False,required=False,action="store_true")
	genoArgs.add_argument('-pre',required=False,default=None)
	genoArgs.add_argument('-feats',required=False,default=None)
	optArgs.add_argument('-L','-log',default=None,required=False)
	optArgs.add_argument('-T','-tmp-dir',default=os.getcwd()+'/sv2_tmp_'+rand_id(),required=False)
	optArgs.add_argument('-s','-seed',required=False,default=42,type=int)
	optArgs.add_argument('-o','-out',required=False,default="sv2_training_features",type=str)
	optArgs.add_argument('-O','-odir',required=False,default=os.getcwd(),type=str)
	optArgs.add_argument('-h','-help',required=False,action="store_true",default=False)
	args = parser.parse_args()
	bams,bed,vcf,snv,ped = args.i,args.b,args.v,args.snv,args.p
	gen,pcrfree,legacy_m,predir,featsdir= args.g,args.pcrfree,args.M,args.pre,args.feats
	logfh, tmp_dir, seed, ofh, odir = args.L,args.T,args.s,args.o,args.O
	_help = args.h
	if (_help==True or len(sys.argv)==1):
		print splash+__useage___
		sys.exit(0)
	if logfh!=None:
		lfh = open(logfh,'w')
		sys.stderr=lfh
	preprocess_files,feats_files={},{}
	gens = ['hg19','hg38','mm10']
	olog = logfh
	if olog == None: olog = 'STDOUT'
	print 'sv2 version:{}    report bugs to <dantaki at ucsd dot edu>       error messages located in {}'.format(__version__,olog)
	Confs=Config()
	if bams==None and predir==None and featsdir==None:
		print 'FATAL ERROR: No BAM file specified <-i, -bam  FILE ...>'
		sys.stderr.write('FATAL ERROR: No BAM file specified <-i, -bam  FILE ...>\n')
		sys.exit(1)
	if snv==None and predir==None and featsdir==None:
		print 'FATAL ERROR: No SNV VCF file specified <-snv  FILE ...>'
		sys.stderr.write('FATAL ERROR: No SNV VCF file specified <-snv  FILE ...>\n')
		sys.exit(1)
	if ped==None:
		print 'FATAL ERROR: No PED file specified <-p, -ped  FILE ...>'
		sys.stderr.write('FATAL ERROR: No PED file specified <-p, -ped  FILE ...>\n')
		sys.exit(1)
	if bed==None and vcf==None:
		print 'FATAL ERROR: No SVs provided <-b, -bed  BED ...> <-v,-vcf  VCF ...>'
		sys.stderr.write('FATAL ERROR: No SVs provided <-b, -bed  BED ...> <-v,-vcf  VCF ...>\n')
		sys.exit(1)
	if gen not in gens:
		print 'FATAL ERROR -g must be hg19 or hg38. NOT {}'.format(gen)
		sys.stderr.write('FATAL ERROR -g must be hg19 or hg38. NOT {}\n'.format(gen))
		sys.exit(1)
	Peds=ped_init(ped)
	if bams!=None: Bams=bam_init(bams,Peds,snv_init(snv),gen)
	SV = sv_init(bed,vcf,gen)
	ofh = ofh.replace('.vcf','').replace('.out','').replace('.txt','')
	make_dir(tmp_dir)
	tmp_dir=slash_check(tmp_dir)
	if not odir.endswith('/'): odir = odir+'/'
	make_dir(odir)
	"""
	PREPROCESSING
	"""
	if predir == None:
		outdir = odir+'sv2_preprocessing/'
		make_dir(outdir)
		for bam in Bams:
			preofh = outdir+bam.id+'_sv2_preprocessing.txt'
			preprocess_files[bam.id]=preofh
			preprocess(bam,preofh,seed,gen,tmp_dir)
	else:
		predir=slash_check(predir)
		for fh in glob(predir+'*sv2_preprocessing.txt'):
			f = open(fh)
			if sum(1 for l in open(fh)) <= 1: continue
			else:
				preids=[]
				for l in f:
					if l.startswith('#'):continue
					preids.append(l.rstrip('\n').split('\t').pop(0))
			f.close()
			for iid in set(preids):
				if iid in Peds.ids : preprocess_files[iid]=fh
	report_time(init_time,'PREPROCESSING COMPLETE')
	""""
	FEATURE EXTRACTION
	"""
	if featsdir == None:
		outdir = odir+'sv2_features/'
		make_dir(outdir)
		for bam in Bams:
			if preprocess_files.get(bam.id) == None:
				sys.stderr.write('WARNING: BAM sample id {} not found in preprocessing files. Skipping ...\n'.format(bam.id))
				continue
			prefh = preprocess_files[bam.id]
			featfh = outdir+bam.id+'_sv2_features.txt'
			feats_files[bam.id]=featfh
			extract_feats(bam,SV.raw,prefh,featfh,gen,pcrfree,legacy_m,Confs,tmp_dir)
	else:
		featsdir=slash_check(featsdir)
		for fh in glob(featsdir+'*sv2_features.txt'):
			f = open(fh)
			if sum(1 for l in open(fh)) <= 1: continue
			else:
				featsid=[]
				for l in f:
					if l.startswith('#'):continue
					featsid.append(l.rstrip('\n').split('\t').pop(5))
			f.close()
			for iid in set(featsid):
				if iid in Peds.ids : feats_files[iid]=fh
	feats=[]
	train_dir = odir+'sv2_training_features/'
	make_dir(train_dir)
	for iid in feats_files:
		with open(feats_files[iid]) as f:
			for l in f: feats.append(tuple(l.rstrip('\n').split('\t')))
	sv2_train_output(feats,Peds,gen,train_dir+ofh)
	shutil.rmtree(tmp_dir)
	lfh.close()
	report_time(init_time,'FEATURE EXTRACTION COMPLETE')