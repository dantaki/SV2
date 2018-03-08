#!/usr/bin/env python2
from datetime import timedelta
from pybedtools import BedTool
from time import time
import json,os,random,string,sys
def accepted_chrom(flag,gen):
	chroms,prefix=[],''
	if flag==True: prefix='chr'
	if 'hg' in gen:
		for i in range(22): chroms.append(prefix+str(i+1))
	else:
		for i in range(19): chroms.append(prefix+str(i+1))
	chroms.append(prefix+'X')
	chroms.append(prefix+'Y')
	return chroms
def check_ids(ped,feats):
	missing_feats, missing_ped = list(set(ped)-set(feats)),list(set(feats)-set(ped))
	if len(missing_feats)>0:
		for x in missing_feats: sys.stderr.write('WARNING: sample {} in ped file is missing features. Skipping {}\n'.format(x,x))
	if len(missing_ped)>0:
		for x in missing_ped: sys.stderr.write('WARNING: sample {} has features but is missing in the ped file. Skipping {}\n'.format(x,x))
	return sorted(list(set(feats).intersection(ped)),key=str.lower)
def check_PAR(sv,gen):
	from sv2Config import Config
	in_par = False
	if len(BedTool(sv,from_string=True).intersect(BedTool('{}par_{}.bed'.format(Config().resource_path(),gen)),f=0.5)) > 0: in_par=True
	return in_par
def dump_json(o,d):
	ofh=open(o,'w')
	json.dump(d,ofh,indent=4)
	ofh.close()
def errFH(fh):
	if not os.path.isfile(fh):
		print 'FATAL ERROR {} does not exist'.format(fh)
		sys.stderr.write('FATAL ERROR {} does not exist\n'.format(fh))
		sys.exit(1)
def format_chrom(c):
	if 'chr' not in c: c='chr'+str(c)
	c = c.replace('chr23','chrX').replace('chr24','chrY')
	return c
def get_path():
	try:
		root = __file__
		if os.path.islink(root): root = os.path.realpath(root)
		return os.path.dirname(os.path.abspath(root))
	except:
		print 'FATAL ERROR {} NOT FOUND\n'.format(__file__)
		sys.stderr.write('FATAL ERROR {} NOT FOUND\n'.format(__file__))
		sys.exit(1)
def make_dir(mydir):
	if not os.path.exists(mydir): os.makedirs(mydir)
def mask_bed(fh,gen):
	from sv2Config import Config
	return BedTool(fh).subtract('{}{}_excluded.bed.gz'.format(Config().resource_path(),gen))
def match_chrom_key(call,Bam,gen):
	chrom_key = format_chrom(str(call[0]))
	if Bam.sex == 1 and (chrom_key=='chrX' or chrom_key=='chrY') and check_PAR('{} {} {}'.format(chrom_key,call[1],call[2]),gen)== True: chrom_key='GENOME'
	return chrom_key
def match_chrom_prefix(c,flag):
	if flag==False and c.startswith('chr'): c=c.replace('chr','')
	elif flag==True and not c.startswith('chr'): c='chr'+c
	return c
def query_yes_no(question, default="yes"):
	"""https://code.activestate.com/recipes/577058/"""
	valid = {"yes": True, "y": True, "ye": True,"no": False, "n": False}
	if default is None: prompt = " [y/n] "
	elif default == "yes": prompt = " [Y/n] "
	elif default == "no": prompt = " [y/N] "
	else: raise ValueError("invalid default answer: '%s'" % default)
	while True:
		sys.stdout.write(question + prompt)
		choice = raw_input().lower()
		if default is not None and choice == '': return valid[default]
		elif choice in valid: return valid[choice]
		else: sys.stdout.write("Please respond with 'yes' or 'no' (or 'y' or 'n')\n")
def rand_id(size=7, chars=string.ascii_uppercase + string.digits):
	return ''.join(random.choice(chars) for _ in range(size))
def reciprocal_overlap(x):
	(s1,s2,e1,e2) = x
	sz1 = e1-s1
	sz2 = e2-s2
	smax=max(s1,s2)
	emin=min(e1,e2)
	ovr = emin-smax
	return min((ovr/sz1),(ovr/sz2))
def report_time(init_time,s):
	f=[]
	s = ' '+s+' <time elapsed: {}>'.format(timedelta(seconds=int(time())-init_time))
	for _ in range(len(s)+1): f.append('*')
	print ''.join(f)+'\n'+s+'\n'+''.join(f)
def slash_check(mydir):
	if not mydir.endswith('/'): mydir = mydir+'/'
	if not os.path.isdir(mydir):
		print 'FATAL ERROR: {} not found.\n'.format(mydir)
		sys.stderr.write('FATAL ERROR: {} not found.\n'.format(mydir))
		sys.exit(1)
	return mydir
def tokenize(l):
	r=l.rstrip().split('\t')
	if len(r)==1: r=l.rstrip().split(' ')
	if len(r)==1:
		sys.stderr.write('WARNING: {} not tokenized correctly. Check if {} is Tab or Space delimited\n'.format(l,l))
		return 0
	else: return r