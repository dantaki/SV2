'''
Copyright <2017> <Danny Antaki>

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
'''
__version__='1.3.2'
import sys,os,argparse
from core import writeConfig,checkConfig,check_in,Bed,check_sv,errFH,reportTime,preprocess,extract_feats,outputTrain
from argparse import RawTextHelpFormatter
from multiprocessing import Pool
from glob import glob
from time import time
def main():
	init_time = int(time())
	splash='                       ____\n  _____________   ___ |___ \\\n /   _____/\   \ /   // ___/\n \_____  \  \   Y   //_____)\n /        \  \     /\n/_________/   \___/\nSupport Vector Structural Variation Genotyper Train: Generate a training set\nVersion {}    Author: Danny Antaki <dantaki@ucsd.edu>    github.com/dantaki/SV2\n\n'.format(__version__)
	parser = argparse.ArgumentParser(description=splash,formatter_class=RawTextHelpFormatter)
	inArgs,genoArgs = parser.add_argument_group('input arguments'),parser.add_argument_group('genotype arguments')
	inArgs.add_argument('-i','-in', help='Tab delimited input [ ID, BAM-PATH, VCF-PATH, M/F ]',default=None,type=str)
	inArgs.add_argument('-b','-bed',help='BED file(s) of SVs',type=str,default=None,nargs='*')
	inArgs.add_argument('-v','-vcf',help='VCF file(s) of SVs',type=str,default=None,nargs='*')
	genoArgs.add_argument('-c','-cpu', help='Parallelize sample-wise. 1 per cpu',required=False,default=1,type=int)
	genoArgs.add_argument('-g','-genome',  help='Reference genome build [ hg19, hg38 ]',required=False,default='hg19',type=str)
	genoArgs.add_argument('-pcrfree',  help='GC content normalization for PCR free libraries',required=False,default=False,action="store_true")
	genoArgs.add_argument('-M', help='bwa mem -M compatibility. Split-reads flagged as secondary instead of supplementary',default=False,required=False,action="store_true")
	genoArgs.add_argument('-s','-seed', help='Preprocessing: integer seed for genome shuffling',required=False,default=42,type=int)
	genoArgs.add_argument('-pre', help='Preprocessing output directory',required=False,default=None)
	genoArgs.add_argument('-feats', help='Feature extraction output directory',required=False,default=None)
	genoArgs.add_argument('-o','-out', help='Output prefix',required=False,default="sv2",type=str)
	args = parser.parse_args()
	infh,bed,vcf = args.i,args.b,args.v
	cores,gen, pcrfree,legacyM,seed,ofh = args.c,args.g,args.pcrfree,args.M,args.s,args.o
	predir,featsdir = args.pre,args.feats
	preprocess_files={}
	feats_files={}
	gens = ['hg19','hg38']
	if infh==None: 
		sys.stderr.write('ERROR sample information file <-i | -in> not defined.\n')
		sys.exit(1)
	if gen not in gens: 
		sys.stderr.write('ERROR -g must be either hg19 or hg38. NOT {}\n'.format(gen))
		sys.exit(1)
	checkConfig()
	raw,sv= [],[]
	if bed!=None:
		for x in bed: raw,sv = check_sv(Bed(x,False),gen,raw,sv)
	if vcf!=None:
		for x in vcf: raw,sv = check_sv(Bed(x,True),gen,raw,sv)
	bam_dict,vcf_dict,gender_dict=check_in(infh)
	"""
	PREPROCESSING
	"""
	if predir == None:
		if cores > 1 :
			pool = Pool(processes=cores)
	       		for bam_id in bam_dict:
				preofh = bam_id+'_sv2_preprocessing.txt'
	      			preprocess_files[bam_id]=os.getcwd()+'/sv2_preprocessing/'+preofh
				pool.apply_async(preprocess, args=(bam_id,bam_dict[bam_id],vcf_dict[bam_id],preofh,gen,seed) )
			pool.close()
			pool.join()
		else:
			for bam_id in bam_dict:
				preofh = bam_id+'_sv2_preprocessing.txt'
				preprocess_files[bam_id]=os.getcwd()+'/sv2_preprocessing/'+preofh
				preprocess(bam_id,bam_dict[bam_id],vcf_dict[bam_id],preofh,gen,seed)	
	else: 
		if not predir.endswith('/'): predir = predir+'/'
		if not os.path.isdir(predir): errFH(predir)
		for fh in glob(predir+'*sv2_preprocessing.txt'):
			f = open(fh)
			if sum(1 for l in open(fh)) <= 1: continue
			else:
				preids=[]
				head=f.next().rstrip('\n')
				for l in f: preids.append(l.rstrip('\n').split('\t').pop(0))
			f.close()
			for iid in set(preids): 
				if gender_dict.get(iid) != None: preprocess_files[iid]=fh
	reportTime(init_time,'PREPROCESSING COMPLETE')
	""""
	FEATURE EXTRACTION
	"""
	if featsdir == None:
		if cores > 1 :
			pool = Pool(processes=cores)
      			for bam_id in bam_dict:
				if preprocess_files.get(bam_id) == None:
					print 'WARNING: preprocessing iid {} does not match inlist iid'.format(bam_id)
					continue
				prefh = preprocess_files[bam_id]
				gtofh = bam_id+'_sv2_features.txt'
				feats_files[bam_id]=os.getcwd()+'/sv2_features/'+gtofh
				pool.apply_async(extract_feats, args=(bam_id,bam_dict[bam_id],vcf_dict[bam_id],sv,prefh,gender_dict[bam_id],gtofh,gen,pcrfree,legacyM) )
			pool.close()
			pool.join()
		else:
			for bam_id in bam_dict:
				if preprocess_files.get(bam_id) == None:
					print 'WARNING: preprocessing iid {} does not match inlist iid'.format(bam_id)
					continue
				prefh = preprocess_files[bam_id]
				gtofh = bam_id+'_sv2_features.txt'
				feats_files[bam_id]=os.getcwd()+'/sv2_features/'+gtofh
				extract_feats(bam_id,bam_dict[bam_id],vcf_dict[bam_id],sv,prefh,gender_dict[bam_id],gtofh,gen,pcrfree,legacyM)
	else: 
		if not featsdir.endswith('/'): featsdir=featsdir+'/'
		if not os.path.isdir(featsdir): errFH(featsdir)
		for fh in glob(featsdir+'*sv2_features.txt'):
			f = open(fh)
			if sum(1 for l in open(fh)) <= 1: continue
			else:
				featsid=[]
				head=f.next().rstrip('\n')
				for l in f: featsid.append(l.rstrip('\n').split('\t').pop(5))
			f.close()
			for iid in set(featsid): 
				if gender_dict.get(iid) != None: feats_files[iid]=fh
	feats=[]
	for iid in feats_files:
		with open(feats_files[iid]) as f:
			for l in f: feats.append(tuple(l.rstrip('\n').split('\t')))
	outputTrain(feats,gender_dict,gen,ofh)
	reportTime(init_time,'FEATURE EXTRACTION COMPLETE')