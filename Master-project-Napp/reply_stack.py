from Bio import Entrez
from Bio import SeqIO

Entrez.email = "sample@example.org"

handle = Entrez.efetch(db="nuccore",
                       id="CP001665",
                       rettype="gb",
                       retmode="text",
                       seq_start=6373, 
                       seq_stop=6422
                       )

whole_sequence = SeqIO.read(handle, "genbank")

print whole_sequence.seq
