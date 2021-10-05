import os
import sys
from io import BytesIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from collections import defaultdict
from Bio import SeqIO
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq, UnknownSeq
from StringIO import StringIO
from Bio import Alphabet
from Bio.Align import MultipleSeqAlignment
from BCBio.GFF import GFFExaminer
from BCBio import GFF
import pprint
from Bio import SeqIO
alll=dict()
all_blocks=dict()
fo = open("/homes/biertank/halima/Downloads/Python-3.4.1/halima_saker_project/foo.txt", "a")
limit_info = dict(
        gff_id = ["1765255"],
        )
#for rec in GFF.parse(in_handle, limit_info=limit_info):
#
files=('xar','xbk','xcd','xcw','xdp','xei', 'xfb','xfu','xgn'  'xhg''xas' , 'xbl' , 'xce' , 'xcx' , 'xdq'  ,'xej' , 'xfc' , 'xfv'  ,'xgo'  ,'xhh','xaa','xat','xbm','xcf' ,'xcy' ,
'xdr' ,'xek','xfd' ,'xfw' ,'xgp','xhi','xab','xau' , 'xbn' , 'xcg' , 'xcz' , 'xds',  'xel' , 'xfe' , 'xfx' , 'xgq ' ,'xhj','xac','xav','xbo',  'xch'  ,'xda' , 'xdt' , 'xem'  ,'xff' ,
'xfy','xgr' , 'xhk','xad','xaw' , 'xbp' , 'xci' , 'xdb','xdu','xen','xfg' , 'xfz','xgs' ,'xhl','xae','xax' , 'xbq'  ,'xcj' , 'xdc',  'xdv' , 'xeo' , 'xfh',  'xga' , 'xgt' , 'xhm',
'xaf','xay' , 'xbr',  'xck' , 'xdd',  'xdw' , 'xep', 'xfi' , 'xgb' , 'xgu' , 'xhn','xag','xaz' , 'xbs' , 'xcl' , 'xde' , 'xdx' , 'xeq' , 'xfj' , 'xgc' , 'xgv' , 'xho','xah','xba' ,
'xbt' , 'xcm' , 'xdf' , 'xdy' , 'xer' , 'xfk' , 'xgd' , 'xgw' , 'xhp','xai','xbb' , 'xbu',  'xcn' , 'xdg' , 'xdz' , 'xes',  'xfl'  ,'xge',  'xgx','xaj','xbc' , 'xbv' , 'xco' , 'xdh',
'xea' , 'xet' , 'xfm' , 'xgf' , 'xgy','xak','xbd' , 'xbw' , 'xcp',  'xdi' , 'xeb',  'xeu',  'xfn', 'xgg'  ,'xgz','xal','xbe' , 'xbx' , 'xcq' , 'xdj' , 'xec' ,'xev' , 'xfo'  ,'xgh' ,
'xha','xam','xbf' , 'xby' , 'xcr' , 'xdk' , 'xed'  ,'xew' , 'xfp' ,'xgi' ,'xhb','xan','xbg' ,'xbz' , 'xcs'  ,'xdl' , 'xee'  ,'xex' , 'xfq'  ,'xgj' , 'xhc','xao','xbh' , 'xca' , 'xct',
'xdm',  'xef' , 'xey', 'xfr' , 'xgk' , 'xhd','xap','xbi' , 'xcb' , 'xcu' , 'xdn' , 'xeg',  'xez'  ,'xfs' , 'xgl',  'xhe','xaq','xbj',  'xcc' , 'xcv' , 'xdo' , 'xeh'  ,'xfa' , 'xft' , 
'xgm' , 'xhf')
for d in files:
	in_file = "/homes/biertank/halima/Downloads/Python-3.4.1/halima_saker_project/Untitled Folder 2/"+d
	#in_handle = open(in_file)	
	with open(in_file) as in_handle:
		for rec in GFF.parse(in_handle):
			t=0
			for record in rec.features:
				elem_metadatas = list()
				keyss=('clst_id', 'SubjectScore','SubjectOrganism')
		#print record.type
				for key in sorted(record.qualifiers.iterkeys()):
					if key in keyss:
						if key=="SubjectOrganism":	
							load_profile = open('/homes/biertank/halima/Downloads/Python-3.4.1/halima_saker_project/speciess.gff')
				
							read_it = load_profile.read()
							myLine =list()
							myscore=list()
							purse=dict()
							for val in record.qualifiers[key]:
								i=1
								for line in read_it.splitlines():
									if line == val:
										myLine.append(i)
									#print (val, "  -------------------------  ",i)
										break
									i=i+1	
							elem_metadatas.append(str(key) + '=' + str(myLine))	
						else:				
							elem_metadatas.append(str(key) + '=' + str(record.qualifiers[key]))
				load_profile.close()
				myscore=record.qualifiers["SubjectScore"]
				for i in range(0,len(myLine)):						
					purse[myLine[i]]=str(myscore[i])	
				alll[rec.id,record.qualifiers["clst_id"][0]]=purse		
		
				for i in range(0,len(myLine)):						
					t=t+float(myscore[i])	
		
		
		
		
				min=0
				max=1
				u=0
				z=1
				lst=sorted(purse.iterkeys())
				blocks=dict()
				list_elmt=list()
				for u in range(0,len(lst)-1):
					#print lst[u]
					#print min, "--",  max
					if lst[max]-lst[max-1]<=2:		
					#print u		
				
						if max-min==1:
							#print (u,"-----------------------")
							list_elmt.append(lst[min])
							list_elmt.append(lst[max])
						else:
							list_elmt.append(lst[max])
						max=max+1						
					else:
				
						if max-min>1:
							blocks[z]=list_elmt
							list_elmt=list()
							z=z+1
						min=max
						max=max+1			
				all_blocks[rec.id]=blocks
			print(all_blocks)
			print d,"\n"	
	#print alll
			fo.write(str(all_blocks))
			all_blocks=dict()
		in_handle.close()
fo.close()
