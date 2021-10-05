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

def species():
	with open('/homes/biertank/halima/Downloads/halima_saker_project_/Escherichia_coli_BL21_DE3-Clusters.gff') as in_handle:
		for rec in GFF.parse(in_handle):
			for record in rec.features:
				#print record.qualifiers["cluster"]
				if record.qualifiers["cluster"]==["7"]:
					print rec

	
species()  
