#!/usr/bin/env python2
from sv2_backend import dump_json,errFH,get_path,query_yes_no
try:from configparser import ConfigParser
except ImportError:from ConfigParser import ConfigParser
import json,os,sys,wget,zipfile
__memsize__='238M'
__resource_url__ = 'https://github.com/dantaki/SV2/releases/download/sv2v1.4.0/sv2_resources.zip'
class Config():
	def __init__(self):
		self.fh=get_path()+'/config/sv2.ini'
		self.json=get_path()+'/config/sv2_clf.json'
		self.fasta_chr_flag=False
		self.clfs=None
		self.gtclf={} # classifiers for genotyping
		self.gcfh=get_path()+'/resources/GC_content_reference.txt'
		self.gc_dict={}
		self.resource=None # directory to resource files
		if not os.path.isfile(self.fh):
			print 'FATAL ERROR: {} config file is missing! Please reinstall sv2: [$ pip uninstall sv2 -y && pip install sv2-VERSION.tar.gz]'.format(self.fh)
			sys.stderr.write('FATAL ERROR: {} config file is missing! Please reinstall sv2: [$ pip uninstall sv2 -y && pip install sv2-VERSION.tar.gz]\n'.format(self.fh))
			sys.exit(1)
		if not os.path.isfile(self.json):
			print 'FATAL ERROR: {} json file is missing! Please reinstall sv2: [$ pip uninstall sv2 -y && pip install sv2-VERSION.tar.gz]'.format(self.json)
			sys.stderr.write('FATAL ERROR: {} json file is missing! Please reinstall sv2: [$ pip uninstall sv2 -y && pip install sv2-VERSION.tar.gz]\n'.format(self.json)) 		
			sys.exit(1)
		if not os.path.isfile(self.gcfh):
			print 'FATAL ERROR: {} file is missing! Please reinstall sv2: [$ pip uninstall sv2 -y && pip install sv2-VERSION.tar.gz]'.format(self.gcfh)
			sys.stderr.write('FATAL ERROR: {} file is missing! Please reinstall sv2: [$ pip uninstall sv2 -y && pip install sv2-VERSION.tar.gz]\n'.format(self.gcfh))
			sys.exit(1)
		with open(self.json) as data_file: self.clfs = json.load(data_file)
	def check_config(self):
		conf=ConfigParser()
		conf.read(self.fh)
		if conf.get('RESOURCE_DIR','sv2_resource')==None or not os.path.isdir(conf.get('RESOURCE_DIR','sv2_resource')):
			print 'FATAL ERROR: sv2 resource files not found.\nDownload the resource files with [ $ sv2 -download ] github.com/dantaki/SV2/wiki/Installation#download-resources'
			sys.stderr.write('FATAL ERROR: sv2 resource files not found.\nDownload the resource files with [ $ sv2 -download ] github.com/dantaki/SV2/wiki/Installation#download-resources\n')
			sys.exit(1)
		if not os.path.isfile(conf.get('FASTA_PATHS','hg19')) and not os.path.isfile(conf.get('FASTA_PATHS','hg38')) and not os.path.isfile(conf.get('FASTA_PATHS','mm10')):
			print 'FATAL ERROR: No reference fasta specified.\nConfigure sv2 before genotyping [ $ sv2 -hg19 FASTA -hg38 FASTA -mm10 FASTA ] github.com/dantaki/SV2/wiki/Installation#configure'
			sys.stderr.write('FATAL ERROR: No reference fasta specified.\nConfigure sv2 before genotyping [ $ sv2 -hg19 FASTA -hg38 FASTA -mm10 FASTA ] github.com/dantaki/SV2/wiki/Installation#configure\n')
			sys.exit(1)
	def download_resource(self,):
		install_loc = get_path()+'/resources/'
		change_loc = query_yes_no('Install sv2 resource files to the default location? < {} >'.format(install_loc))
		if change_loc == False:
			sys.stdout.write('Enter sv2 resource file destination:\n')
			install_loc = raw_input()
			if not install_loc.endswith('/'): install_loc=install_loc+'/'
			if not os.path.isdir(install_loc):
				sys.stderr.write('FATAL ERROR: {} not a directory\n'.format(install_loc))
				sys.stdout.write('FATAL ERROR: {} not a directory\n'.format(install_loc))
				sys.exit(1)
		install_ok = query_yes_no('sv2 resource files will require {} of free space in {}. Proceed with installation?'.format(__memsize__,install_loc))	
		if install_ok==False: sys.exit(1)
		else:
			zip_file = install_loc+'sv2_resources.zip'
			wget.download(__resource_url__,zip_file)
			zip_fh= zipfile.ZipFile(zip_file)
			if zip_fh.testzip() is not None:
				sys.stdout.write('FATAL ERROR: {} corrupted. Try <sv2 -download> again\n'.format(zip_file))
				sys.stderr.write('FATAL ERROR: {} corrupted. Try <sv2 -download> again\n'.format(zip_file))
				sys.exit(1)
			zip_fh.extractall(install_loc)
			zip_fh.close()
			os.remove(zip_file)
			self.resource = install_loc
			self.write_config(None,None,None)
		sys.stdout.write('INSTALL COMPLETE\nRun < sv2 -hg19 FASTA [-hg38,-mm10] > to configure sv2\n')
	def fasta(self,gen=None):
		conf = ConfigParser()
		conf.read(self.fh)
		fasta_path=conf.get('FASTA_PATHS',gen)
		if str(fasta_path)=='None':
			print 'FATAL ERROR: {} fasta for reference {} not found. Please configure SV2\n[ $ sv2 -hg19 FASTA -hg38 FASTA -mm10 FASTA ] github.com/dantaki/SV2/wiki/Installation#configure'.format(fasta_path,gen)
			sys.stderr.write('FATAL ERROR: {} fasta for reference {} not found. Please configure SV2\n[ $ sv2 -hg19 FASTA -hg38 FASTA -mm10 FASTA ] github.com/dantaki/SV2/wiki/Installation#configure\n'.format(fasta_path,gen))
			sys.exit(1)
		with open(fasta_path,'r') as f:
			l = f.next()
			if not l.startswith('>'):
				print 'FATAL ERROR: {} headers must contain ">"'.format(fasta_path)
				sys.stderr.write('FATAL ERROR: {} headers must contain ">"\n'.format(fasta_path))
				sys.exit(1)
			if l.startswith('>chr'): self.fasta_chr_flag=True
		return fasta_path
	def gc_norm(self,pcrfree):
		with open(self.gcfh,'r') as f:
			for l in f:
				tag,GC,norm = l.rstrip('\n').split('\t')
				key = 'DOC'
				if 'READCOUNT' in tag: key ='RC'
				if pcrfree == True:
					if 'PCRFREE' in tag: self.gc_dict[(key,int(GC))]=float(norm)
				else:
					if 'PCRFREE' not in tag: self.gc_dict[(key,int(GC))]=float(norm)
	def get_clf(self,classifier_name=None):
		names=[x for x in self.clfs['default']]
		if self.clfs.get(classifier_name)==None:
			sys.stderr.write('WARNING: {} classifier name not found. Using default classifiers. Please check {}\n'.format(classifier_name,self.json))
			classifier_name='default'
		if self.clfs.get(classifier_name)!=None:
			for x in names:
				if self.clfs[classifier_name].get(x)==None:
					sys.stderr.write('WARNING {} {} classifier not found. Using default {} classifier. Please check {}\n'.format(classifier_name,x,x,self.json))
					self.gtclf[x]=self.clfs['default'][x]
				else:
					if not os.path.isfile(self.clfs[classifier_name][x]):
						sys.stderr.write('WARNING {} {} pickle file not found. Using default {} classifier. Please check {} and reload the classifier json file <sv2 -load-clf <clf.json>\n'.format(classifier_name,x,x,self.json))
					else: self.gtclf[x]=self.clfs[classifier_name][x]
		else:
			print 'FATAL ERROR: {} classifier not found. Please check {}. If the default classifiers are not found, reinstall sv2 <pip uninstall sv2 -y && pip install sv2-VERSION-.tar.gz'.format(classifier_name,self.json)
			sys.stderr.write('FATAL ERROR: {} classifier not found. Please check {}. If the default classifiers are not found, reinstall sv2 <pip uninstall sv2 -y && pip install sv2-VERSION-.tar.gz\n'.format(classifier_name,self.json))
			sys.exit(1)
	def load_clf(self,jsonfh=None):
		import shutil
		errFH(jsonfh)
		clfs={}
		with open(jsonfh) as f: clfs = json.load(f)
		realpaths=[]
		for name in clfs:
			print 'loading {} classifier ...'.format(name)
			clf_dir = get_path()+'/resources/training_sets/'+name+'/'
			if not os.path.exists(clf_dir): os.makedirs(clf_dir)		
			for x in clfs[name]:
				clffh = str(clfs[name][x])
				if not os.path.isfile(clffh): sys.stderr.write('WARNING: {} does not exist. Please check the paths in {}\n'.format(clffh,jsonfh))
				else:
					clfname = clffh.split('/').pop()
					newpath = clf_dir+clfname
					clfreplace=True
					if os.path.isfile(newpath): clfreplace= query_yes_no('WARNING: {} exists... Replace?'.format(newpath),'no')
					if clfreplace==True:
						if self.clfs.get(name)==None: self.clfs[name]={}
						self.clfs[name][x]=newpath
						realpaths.append(newpath)
						shutil.copyfile(clffh, newpath)
		for x in realpaths: print 'installed classifier {}'.format(x)
		print 'appending to {}... DO NOT ALTER {}'.format(self.json,self.json)
		dump_json(self.json,self.clfs)
	def resource_path(self):
		conf = ConfigParser()
		conf.read(self.fh)
		resource_loc = conf.get('RESOURCE_DIR','sv2_resource')
		if not resource_loc.endswith('/'): resource_loc=resource_loc+'/'
		if not os.path.isdir(resource_loc):
			sys.stdout('FATAL ERROR: {} not found. Install resource files < sv2 -download >\n'.format(resource_loc))
			sys.stderr('FATAL ERROR: {} not found. Install resource files < sv2 -download >\n'.format(resource_loc))
			sys.exit(1)
		return resource_loc
	def write_config(self,hg19=None,hg38=None,mm10=None):
		conf=ConfigParser()
		for x in self.clfs['default']:
			clf_fh = str(self.clfs['default'][x])
			if not os.path.isfile(clf_fh):
				realpath = get_path()+'/resources/training_sets/'+clf_fh
				if os.path.isfile(realpath): self.clfs['default'][x]=realpath
				else:
					print 'FATAL ERROR: {} pickle file not found. If this file is missing, reinstall sv2: pip uninstall sv2 -y && pip install sv2-VERSION.tar.gz'.format(realpath)
					sys.stderr.write('FATAL ERROR: {} pickle file not found. If this file is missing, reinstall sv2: pip uninstall sv2 -y && pip install sv2-VERSION.tar.gz\n'.format(realpath))
					sys.exit(1)
		dump_json(self.json,self.clfs)
		if not os.path.isfile(self.fh):
			conf_fh = open(self.fh,'w')
			conf.add_section('FASTA_PATHS')
			conf.set('FASTA_PATHS','hg19',hg19)
			conf.set('FASTA_PATHS','hg38',hg38)
			conf.set('FASTA_PATHS','mm10',mm10)
			conf.set('RESOURCE_DIR','sv2_resource',self.resource)
			conf.write(conf_fh)
			conf_fh.close()
		else:
			conf.read(self.fh)
			if hg19 != None:
				errFH(hg19)
				conf.set('FASTA_PATHS','hg19',hg19)
			if hg38 != None:
				errFH(hg38)
				conf.set('FASTA_PATHS','hg38',hg38)
			if mm10 != None:
				errFH(mm10)
				conf.set('FASTA_PATHS','mm10',mm10)
			if self.resource != None:
				conf.set('RESOURCE_DIR','sv2_resource',self.resource)
			with open(self.fh,'w') as conf_fh: conf.write(conf_fh)