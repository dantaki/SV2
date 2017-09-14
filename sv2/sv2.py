__version__='1.3.2'
import sys,os,argparse
from .core import loadClf,writeConfig,checkConfig,check_in,Bed,check_sv,errFH,reportTime,preprocess,extract_feats,genotype,mergeSV,annotate
from argparse import RawTextHelpFormatter
from multiprocessing import Pool
from glob import glob
from time import time
def main():
	init_time = int(time())
	splash='                       ____\n  _____________   ___ |___ \\\n /   _____/\   \ /   // ___/\n \_____  \  \   Y   //_____)\n /        \  \     /\n/_________/   \___/\nSupport Vector Structural Variation Genotyper\nVersion 1.3.2    Author: Danny Antaki <dantaki@ucsd.edu>    github.com/dantaki/SV2\n\n'
	parser = argparse.ArgumentParser(description=splash,formatter_class=RawTextHelpFormatter)
	inArgs, genoArgs,clfArgs,configArgs= parser.add_argument_group('input arguments'),parser.add_argument_group('genotype arguments'), parser.add_argument_group('classifier arguments'),parser.add_argument_group('config arguments')
	inArgs.add_argument('-i','-in', help='Tab or space delimited input [ID, BAM-PATH, VCF-PATH, M/F]',default=None,type=str)
	inArgs.add_argument('-b','-bed',help='BED file(s)',type=str,default=None,nargs='*')
	inArgs.add_argument('-v','-vcf',help='VCF file(s)',type=str,default=None,nargs='*')
	genoArgs.add_argument('-c','-cpu', help='Parallelize sample-wise. 1 per cpu',required=False,default=1,type=int)
	genoArgs.add_argument('-g','-genome',  help='Reference genome build [hg19, hg38]',required=False,default='hg19',type=str)
	genoArgs.add_argument('-pcrfree',  help='GC content normalization for PCR free chemistries',required=False,default=False,action="store_true")
	genoArgs.add_argument('-M', help='bwa mem -M compatibility. Split-reads flagged as secondary instead of supplementary',default=False,required=False,action="store_true")
	genoArgs.add_argument('-s','-seed', help='Preprocessing integer seed for genome shuffling',required=False,default=42,type=int)
	genoArgs.add_argument('-o','-out', help='Output',required=False,default="sv2_genotypes.vcf",type=str)
	genoArgs.add_argument('-merge',help='Merge SV after genotyping',required=False,default=False,action="store_true")
	genoArgs.add_argument('-min-ovr',help='Minimum reciprocal overlap for merging SVs [0.8]',required=False,default=None,type=float)
	genoArgs.add_argument('-pre', help='Preprocessing output directory',required=False,default=None)
	genoArgs.add_argument('-feats', help='Feature extraction output directory',required=False,default=None)
	clfArgs.add_argument('-load-clf',help='Add custom classifiers. -load-clf <clf.JSON>', required=False, default=None,type=str)
	clfArgs.add_argument('-clf',help='Specify classifiers for genotyping [default]',required=False,default='default',type=str)
	configArgs.add_argument('-hg19',default=None,help='hg19 FASTA',required=False)
	configArgs.add_argument('-hg38',default=None,help='hg38 FASTA',required=False)	
	args = parser.parse_args()
	infh,bed,vcf = args.i,args.b,args.v
	cores,gen, pcrfree,legacyM,ofh,seed = args.c,args.g,args.pcrfree,args.M,args.o,args.s
	mergeFlag,minOvr = args.merge,args.min_ovr
	predir,featsdir = args.pre,args.feats
	clfLoad, CLF = args.load_clf,args.clf
	conf_hg19,conf_hg38=args.hg19,args.hg38
	preprocess_files={}
	feats_files={}
	gens = ['hg19','hg38']
	if clfLoad!=None: 
		loadClf(clfLoad)
		sys.exit(0)
	if conf_hg19 != None or conf_hg38 != None: 
		writeConfig(conf_hg19,conf_hg38)
		sys.exit(0)
	if infh==None: 
		sys.stderr.write('ERROR sample information file <-i | -in> not defined.\n')
		sys.exit(1)
	if gen not in gens: 
		sys.stderr.write('ERROR -g must be hg19 or hg38. NOT {}\n'.format(gen))
		sys.exit(1)
	checkConfig()
	bam_dict,vcf_dict,gender_dict=check_in(infh)
	raw,sv= [],[]
	if bed!=None:
		for x in bed: raw,sv = check_sv(Bed(x,False),gen,raw,sv)
	if vcf!=None:
		for x in vcf: raw,sv = check_sv(Bed(x,True),gen,raw,sv)
	ofh = ofh.replace('.txt','.vcf').replace('.out','.vcf')
	if not ofh.endswith('.vcf'): ofh=ofh+'.vcf'
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
					print 'WARNING: preprocessing iid {} does not match sample information iid'.format(bam_id)
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
					print 'WARNING: preprocessing iid {} does not match sample information iid'.format(bam_id)
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
	reportTime(init_time,'FEATURE EXTRACTION COMPLETE')	
	"""
	GENOTYPING
	"""
	if mergeFlag==True and minOvr==None: minOvr=0.8
	if minOvr!=None: mergeFlag=True
	feats=[]
	for iid in feats_files:
		with open(feats_files[iid]) as f:
			for l in f: feats.append(tuple(l.rstrip('\n').split('\t')))
	genos,REF,NON,GQ,HEMI,STD_FILT,DNM_FILT = genotype(CLF,raw,feats,gender_dict,gen,ofh.replace('.vcf','.txt'))
	if mergeFlag==True: raw = mergeSV(raw,NON,minOvr)
	annotate(raw,genos,gen,REF,NON,GQ,ofh,gender_dict,HEMI,STD_FILT,DNM_FILT)
	reportTime(init_time,'GENOTYPING COMPLETE')	
