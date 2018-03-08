#!/usr/bin/env python
from sv2_backend import errFH,tokenize
import sys
class Ped():
	def __init__(self,peds=None):
		self.sex={}
		self.males={}
		self.ids=[]
		for fh in peds:
			errFH(fh)
			with open(fh,'r') as f:
				for l in f:
					r = tokenize(l)
					if r!=0:
						if len(r)<5: sys.stderr.write('WARNING {} does not contain 5 elements. PED files are formatted as:\nFamily ID  Individual ID  Paternal ID  Maternal ID  Sex(1=male;2=female;)\n'.format(l))
						else:
							sample_id,sex = str(r[1]),int(r[4])
							if sex != 1 and sex != 2: sys.stderr.write('WARNING {} is not an accepted sex entry. Accepted entries: 1=male; 2=female. Error found here:{}\n'.format(sex,l))
							else:
								self.sex[sample_id]=sex
								self.ids.append(sample_id)
								if sex==1: self.males[sample_id]=1
		self.ids=sorted(list(set(self.ids)),key=str.lower)
def ped_init(peds=None): return Ped(peds)