#!/usr/bin/env python
from sv2_backend import format_chrom
from pybedtools import BedTool
import itertools
def check_merge(sv,minOvr):
	to_merge=[]
	for x in sv.intersect(sv,f=minOvr,F=minOvr,wao=True):
		x=tuple(x)
		if skip_sv(x)==False:  to_merge.append(x)
	return  to_merge
def first_merge(sv,minOvr):
	to_merge=[]
	for x in sv.intersect(sv,f=minOvr,F=minOvr,wao=True):
		x=tuple(x)
		k1,k2 = tokenize_sv(x),tokenize_sv(x,False)
		if k1 != k2:
			 to_merge.append(k1)
			 to_merge.append(k2)
	return [x for x in sv if tokenize_sv(tuple(x)) not in  to_merge]
def merging(sv,non):
	merged,trash=[],[]
	for x in sv:
		k1,k2 = tokenize_sv(x),tokenize_sv(x,False)
		non1,non2=0,0
		if non.get(k1)!=None: non1=non[k1]
		if non.get(k2)!=None: non2=non[k2]
		if non1=='NA':non1=0
		if non2=='NA':non2=0
		if non1 > non2:
			if k1 not in trash: merged.append(k1)
			trash.append(k2)
		elif non2 > non1:
			if k2 not in trash: merged.append(k2)
			trash.append(k1)
		else:
			if k1 not in trash: merged.append(k1)
			trash.append(k2)
	return BedTool(list(set(merged)))
def merge_sv(raw,SVs,minOvr):
	non,result,merged={},{},[]
	for locus in SVs: non[locus]=SVs[locus].med_alt
	dels = BedTool([(format_chrom(x[0]),x[1],x[2],x[3]) for x in raw if 'DEL' in x[3]])
	dups = BedTool([(format_chrom(x[0]),x[1],x[2],x[3]) for x in raw if 'DUP' in x[3]])
	dels_nomerge = first_merge(dels,minOvr)
	while len(check_merge(dels,minOvr))>0:
		dels = merging(check_merge(dels,minOvr),non)
	dups_nomerge = first_merge(dups,minOvr)
	while len(check_merge(dups,minOvr))>0:
		dups = merging(check_merge(dups,minOvr),non)
	for sv in list(itertools.chain(dels,dels_nomerge,dups,dups_nomerge)): 
		sv=tokenize_sv(tuple(sv))
		result[sv]=1
	for sv in raw:
		chrom,start,end,svtype = sv
                _sv = (format_chrom(chrom),start,end,svtype)
                if result.get(_sv)!=None: merged.append(sv) 
	return merged
def skip_sv(x):
	k1,k2 = tokenize_sv(x),tokenize_sv(x,False)
	if k1==k2: return True
	else: return False
def tokenize_sv(x,first_entry=True):
	if first_entry==True: return (str(x[0]),int(x[1]),int(x[2]),str(x[3]))
	else: return (str(x[4]),int(x[5]),int(x[6]),str(x[7]))
