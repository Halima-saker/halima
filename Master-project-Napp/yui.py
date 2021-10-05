from Bio import Entrez
from Bio import SeqIO
Entrez.email = 'A.N.Other@example.com'
handle = Entrez.efetch(db="nuccore",
                    id="CP001618",
                    rettype="genbank",
                    seq_start=70731,
                    seq_stop=70780)
whole_sequence = SeqIO.read(handle, "genbank")
print (whole_sequence.seq)
