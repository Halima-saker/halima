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
	fo = open("/homes/biertank/halima/halima__/rec_spe/"+val+".gff", "a")
	for line in fileinput.input(['/homes/biertank/halima/halima__/All_Napp_Pred_Profils.gff']):
		print " halima"
		values = line.split("\t")
		if values[1]==val:
			fo.write(str(line))

load_profile = open('/homes/biertank/halima/halima__/speciess.gff')				
read_it = load_profile.read()		

for lines in read_it.splitlines():
	species(str(lines))
	print "halima"	

