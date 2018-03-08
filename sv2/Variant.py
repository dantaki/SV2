#!/usr/bin/env python
from sv2_backend import format_chrom
import numpy as np
class Variant:
	def __init__(self,data=None):
		self.locus=(format_chrom(data[0]),int(data[1]),int(data[2]),str(data[3])) # (chrom,start,end,type)
		self.ref={} # [classifier_name]=[reference genotype likelihoods]
		self.med_ref='NA'
		self.alt={} # [classifier_name]=[alternate genotype likelihoods]
		self.med_alt='NA'
		self.format={} # [locus+id]='GT:CN:PE:SR:SC:SP:AR:HT:SQ:GL'
		self.gt={} # gt[locus+id]=genotype
		self.gq={} # gq[locus+id]=(ref,het,hom) likelihoods
		self.breakpoint_feats=0 # breakpoint features: discordant paired-ends + split-reads
		self.clf=[] # classifiers
		self.standard_filter='PASS'
		self.denovo_filter='PASS'
	def add_denovo_filters(self):
		yea,nay = 0,0
		for clfname in self.alt:
			non,feat,svlen =  self.alt[clfname],self.breakpoint_feats,self.locus[2]-self.locus[1]
			passflag=0
			ref=20
			if self.ref.get(clfname)!=None: ref=self.ref[clfname]
			if 'DEL' in clfname and 'male_sex_chrom' not in clfname:
				if feat != 0 and non >=12 and ref >=12: passflag=1
				elif feat==0 and 1000<=svlen<3000 and non >=24 and ref >=20: passflag=1
				elif feat==0 and 3000<=svlen<5000 and non >=20 and ref >=20: passflag=1
				elif feat==0 and svlen>= 5000 and non >=20 and ref >=18: passflag=1	
			elif 'DUP_breakpoint' in clfname:
				if feat != 0 and non >=11 and ref >=11: passflag=1
				if feat==0 and 3000<=svlen<100000 and non >=10 and ref >=15: passflag=1
				if feat==0 and  svlen>=100000 and non >=10 and ref >=13: passflag=1
			elif 'DUP_snv' in clfname and non >=10 and ref >=10 and 'DUP_breakpoint' in self.clf: passflag=1
			elif 'DUP_snv' in clfname and 'DUP_breakpoint' not in self.clf:
				if  svlen<5000 and ref >=18 and non >=18: passflag=1
				elif 5000<= svlen<100000 and ref >=15 and non >=10: passflag=1
				elif svlen>=100000 and ref >=13 and non >=10: passflag=1
			elif 'DEL_male_sex_chrom' in clfname and non>=8 and ref>=8: passflag=1
			elif 'DUP_male_sex_chrom' in clfname and  svlen>=5000 and ref>=10 and non>=10: passflag=1
			if passflag==0: nay+=1
			else: yea+=1
		if nay>yea: self.denovo_filter='FAIL'
	def add_standard_filters(self):
		yea,nay=0,0
		for clfname in self.alt:
			lik,feat,svlen = self.alt[clfname],self.breakpoint_feats,self.locus[2]-self.locus[1]
			passflag=1
			if 'DEL' in clfname and 'male_sex_chrom' not in clfname:
				if feat!=0 and lik<8: passflag=0
				elif feat==0:
					if 'DEL_lt1kb' in clfname and lik <18: passflag=0
					if 'DEL_gt1kb' in clfname:
						if  svlen < 3000: passflag=0
						elif 3000<= svlen<5000 and lik <20: passflag=0
						elif  svlen >=5000 and lik <18: passflag=0
			elif 'DUP_breakpoint' in clfname:
				if feat != 0:
					if  svlen < 1000 and lik < 12: passflag=0
					if  svlen >= 1000 and lik < 10: passflag=0
				elif feat==0:
					if  svlen < 3000: passflag=0
					elif  svlen >=3000 and lik < 12: passflag=0
			elif 'DUP_snv' in clfname:
				if 'DUP_breakpoint' in self.clf:
					if lik < 8: passflag=0
					if  svlen<3000 and lik < 18: passflag=0
					if  svlen>=3000 and lik < 13: passflag=0
			elif 'DEL_male_sex_chrom' in clfname:
				if lik<8 and feat!=0: passflag=0
				if lik<10 and feat==0 and  svlen<=1000: passflag=0
				if  svlen<1000 and feat==0: passflag=0
			elif 'DUP_male_sex_chrom' in clfname:
				if  svlen<5000: passflag=0
				elif  svlen>=5000 and lik<10: passflag=0
				elif feat==0: passflag=0
			if passflag==1: yea+=1
			else: nay+=1
		if nay>yea: self.standard_filter='FAIL'
	def genotype(self,data=None,label=None,copy_number=None,lik=None,biallele=True):
		self.breakpoint_feats+=float(data[7])+float(data[8])
		if biallele==True:
			gq=','.join((format(float(data[14]),'.2f'),format(float(data[15]),'.2f'),format(float(data[16]),'.2f')))
		else:
			gq=','.join((format(float(data[14]),'.2f'),format(float(data[15]),'.2f')))
		_format = '{}:{}:{}'.format(label,copy_number,':'.join(map(str,(format(float(data[7]),'.2f',), format(float(data[8]),'.2f',), format(float(data[9]),'.2f',), int(data[10]), format(float(data[11]),'.2f'), int(np.rint(float(data[12]))),int(np.rint(lik)),gq))))		
		self.gt[self.locus+(str(data[5]),)]=_format
	def median_likelihood(self):
		if len(self.ref)>0:
			self.med_ref=[]
			for clfname in self.ref:
				self.med_ref+=self.ref[clfname]
				self.ref[clfname]=int(np.rint(np.median(self.ref[clfname])))
			self.med_ref=int(np.rint(np.median(self.med_ref)))
		if len(self.alt)>0:
			self.med_alt=[]
			for clfname in self.alt:
				self.med_alt+=self.alt[clfname]
				self.alt[clfname]=int(np.rint(np.median(self.alt[clfname])))
			self.med_alt=int(np.rint(np.median(self.med_alt)))