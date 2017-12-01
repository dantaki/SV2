#!/usr/bin/env python2
import pandas as pd
import numpy as np
try: import cPickle as pickle
except ImportError: import pickle
def biallelic_del_svm(df,gtclf):
	gt_df,gt_mat,lt_df,lt_mat = prep_biallelic_dels(df)
	final=pd.DataFrame()
	if len(gt_mat) > 0:
		clf=None
		with open(gtclf['delgt'],'rb') as f: clf=pickle.load(f)
		final=final.append(biallelic_del_tabulate(gt_df,clf.predict_proba(gt_mat),'DEL_gt1kb'))
	if len(lt_mat) > 0:
		clf=None
		with open(gtclf['dellt'],'rb') as f: clf=pickle.load(f)
		final=final.append(biallelic_del_tabulate(lt_df,clf.predict_proba(lt_mat),'DEL_lt1kb'))
	return final
def biallelic_dup_svm(df,gtclf):
	snv_df,brk_df = prep_biallelic_dups(df)
	final=pd.DataFrame()
	if len(snv_df) >0:
		clf=None
		snv_mat = tabulate_snv(snv_df)
		with open(gtclf['duphar'],'rb') as f: clf=pickle.load(f)
		final=final.append(biallelic_dup_tabulate(snv_df,clf.predict_proba(snv_mat),'DUP_har'))
	if len(brk_df) >0:
		clf=None
		brk_mat= tabulate_3feats(brk_df)
		with open(gtclf['dupbrk'],'rb') as f: clf=pickle.load(f)
		final=final.append(biallelic_dup_tabulate(brk_df,clf.predict_proba(brk_mat),'DUP_breakpoint'))
	return final
def male_sex_chrom_svm(df,gtclf,clfname):
	df[['covr_GC','dpe','sr']]=df[['covr_GC','dpe','sr']].astype(float)
	df['covr']=df['covr_GC']
	df = df[df['covr'] < 5.0]
	mat =  tabulate_3feats(df)
	final=pd.DataFrame()
	if len(mat)>0:
		clf=None
		with open(gtclf[clfname],'rb') as f: clf=pickle.load(f)
		if clfname=='delmsc':
			final=final.append(male_sex_chrom_del_tabulate(df,clf.predict_proba(mat),'DEL_male_sex_chrom'))
		else:
			final=final.append(male_sex_chrom_dup_tabulate(df,clf.predict_proba(mat),'DUP_male_sex_chrom'))
	return final
def training_features(df,case,out): # case 0:dels; 1:dups; 2:del_msc; 3:dup_msc
	pd.options.mode.chained_assignment = None
	if case==0:
		gt_df,gt_mat,lt_df,lt_mat = prep_biallelic_dels(df)
		if len(gt_df)>0:
			ofh = ofh_train('{}_deletion_gt1kb.txt'.format(out))
			for x in train_reformat(gt_df).values: ofh.write('\t'.join(map(str,x))+'\tNA\tdeletion_gt1kb\n')
			ofh.close()
		if len(lt_df)>0:
			ofh = ofh_train('{}_deletion_lt1kb.txt'.format(out))
			for x in train_reformat(lt_df).values: ofh.write('\t'.join(map(str,x))+'\tNA\tdeletion_lt1kb\n')
	if case==1:
		snv_df,brk_df = prep_biallelic_dups(df)
		if len(brk_df)>0:
			ofh = ofh_train('{}_duplication_breakpoint.txt'.format(out))
			for x in train_reformat(brk_df).values: ofh.write('\t'.join(map(str,x))+'\tNA\tduplication_breakpoint\n')
			ofh.close()
		if len(snv_df)>0:
			ofh= ofh_train('{}_duplication_har.txt'.format(out))
			for x in train_reformat(snv_df).values: ofh.write('\t'.join(map(str,x))+'\tNA\tduplication_har\n')
			ofh.close()
	if case==2 or case==3:
		df[['covr_GC','dpe','sr']]=df[['covr_GC','dpe','sr']].astype(float)
		df['covr']=df['covr_GC']
		df = df[df['covr'] < 5.0]
		out_suffix = 'deletion_male_sex_chrom'
		if case==3: out_suffix='duplication_male_sex_chrom'
		if len(df)>0:
			ofh = ofh_train('{}_{}.txt'.format(out,out_suffix))
			for x in train_reformat(df).values: ofh.write('\t'.join(map(str,x))+'\tNA\t{}\n'.format(out_suffix))
			ofh.close()
cdef data_table(df,lik,cn,biallelic=True):
	df['copy_number']=cn
	if biallelic==True:
		df['likHOM'],df['likHET'],df['likREF']=lik[:,0],lik[:,1],lik[:,2]
	else:
		df['likHOM'],df['likHET'],df['likREF']=lik[:,0],lik[:,0],lik[:,1]
	df['NONREF'] = 1 - df['likREF']
	pref = 1.0 - df['likREF']
	palt = 1.0 - df['NONREF']
	pref[pref==1]=1-1e-12
	pref[pref==0]=0+1e-12
	palt[palt==1]=1-1e-12
	palt[palt==0]=0+1e-12
	df['PHRED_REF'] = -10.0 * np.log10(pref)
	df['PHRED_NONREF'] = -10.0 * np.log10(palt)
	df= df[['chr','start','end','type','size','id','covr','dpe','sr','SNP_coverage','SNPs','HET_ratio','HET_SNPs','copy_number','likREF','likHET','likHOM','NONREF','PHRED_REF','PHRED_NONREF','classif']]	
	return df
cdef biallelic_del_tabulate(df,lik,clfname):
	df['classif']=clfname
	cn=[]
	for x in lik:
		if x[0] > x[1] and x[0] > x[2]: cn.append(0)
		elif x[1] > x[0] and x[1] > x[2]: cn.append(1)
		elif x[2] > x[0] and x[2] > x[1]: cn.append(2)
		else: cn.append(2)
	return data_table(df,lik,cn)
cdef biallelic_dup_tabulate(df,lik,clfname):
	df['classif']=clfname
	cn=[]
	for x in lik:
		if x[0] > x[1] and x[0] > x[2]: cn.append(4)
		elif x[1] > x[0] and x[1] > x[2]: cn.append(3)
		elif x[2] > x[0] and x[2] > x[1]: cn.append(2)
		else: cn.append(2)	
	return data_table(df,lik,cn)
cdef male_sex_chrom_del_tabulate(df,lik,clfname):
	df['classif']=clfname
	cn=[]
	for x in lik:
		if x[0] > x[1]: cn.append(0)
		else: cn.append(1)
	return data_table(df,lik,cn,False)
cdef male_sex_chrom_dup_tabulate(df,lik,clfname):
	df['classif']=clfname
	cn=[]
	for x in lik:
		if x[0] > x[1]: cn.append(2)
		else: cn.append(1)
	return data_table(df,lik,cn,False)
cdef ofh_train(o):
	ofh = open(o,'w')
	ofh.write('\t'.join(('chr','start','end','type','size','id','covr','dpe_rat','sr_rat','SNP_coverage','SNPs','HET_ratio','HET_SNPs','copy_number','classifier'))+'\n')
	return ofh
cdef prep_biallelic_dels(df):
	df[['covr_GC','dpe','sr']]=df[['covr_GC','dpe','sr']].astype(float)
	df[['size']]=df[['size']].astype(int)
	df = df[df['covr_GC'] < 5.0]
	df['covr']=df['covr_GC']
	gt_df = df[df['size'] > 1000]
	lt_df = df[df['size'] <= 1000]
	gt_mat = tabulate_3feats(gt_df)
	lt_mat = tabulate_3feats(lt_df)
	return gt_df,gt_mat,lt_df,lt_mat
cdef prep_biallelic_dups(df):
	snv_df,brk_df= pd.DataFrame(),pd.DataFrame()
	df[['SNPs','HET_SNPs']]=df[['SNPs','HET_SNPs']].astype(int)
	df[['covr_GC','dpe','sr']]=df[['covr_GC','dpe','sr']].astype(float)
	mean_covr=[]
	for ind,ele in df.iterrows():
		if ele['SNPs'] > 0: mean_covr.append(np.mean((ele['covr_GC'],float(ele['SNP_coverage']))))
		else: mean_covr.append(ele['covr_GC'])
	df['covr']=mean_covr
	df=df[df['covr'] < 5.0]
	snvs=[]
	brks=[]
	for ind,ele in df.iterrows():
		if ele['dpe']==0 and ele['sr']==0 and ele['HET_SNPs'] > 0: snvs.append(ind)
		else: brks.append(ind)
	snv_df=df.loc[snvs]
	snv_df[['HET_ratio']]=snv_df[['HET_ratio']].astype(float)
	brk_df=df.loc[brks]
	return snv_df,brk_df
cdef tabulate_3feats(df): return df[['covr','dpe','sr']].as_matrix()
cdef tabulate_snv(df): return df[['covr','HET_ratio']].as_matrix()
cdef train_reformat(df): return df[['chr','start','end','type','size','id','covr','dpe','sr','SNP_coverage','SNPs','HET_ratio','HET_SNPs']]	