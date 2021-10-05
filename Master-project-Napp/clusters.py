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
import fileinput
from Bio import SeqIO

def species(val):
	clus_id=list()
	#fo = open("/homes/biertank/halima/Downloads/halima_saker_project_/"+val+"xxxx", "a")
	for line in fileinput.input(['/homes/biertank/halima/Downloads/halima_saker_project_/All_Napp_Pred_Profils.gff']):
		values = line.split("\t")
		clus_id.append(values[3])
	#print clus_id
	mynewlist = list(set(clus_id))	
	print mynewlist
spec="test.gff"
	
species(spec)  
