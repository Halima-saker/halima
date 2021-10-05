from Bio import Entrez, SeqIO
from Bio.Blast import NCBIXML
from Bio.Blast import NCBIWWW

#help(NCBIWWW.qblast)



#result_handle = NCBIWWW.qblast("blastn","genbank", "CACTTATTTAGTTAGCTTGCAACCCTGGATTTTTGTTTACTGGAGAGGCC",entrez_query='Beutenbergia_cavernae_DSM_12333[org]')
result_handle = NCBIWWW.qblast(
    "blastn",
    "nr",
    "CACTTATTTAGTTAGCTTGCAACCCTGGATTTTTGTTTACTGGAGAGGCC",
    megablast=False,
    expect=1000,
    word_size=11,
   # query_to=49,
    nucl_reward=1,
    nucl_penalty=-1,
    gapcosts="5 2",
    entrez_query='Beutenbergia cavernae DSM 12333 [Organism]',
    alignments=10)
#records = NCBIXML.parse(result_handle)
blast_records = NCBIXML.parse(result_handle)
#blast_record = next(blast_records)
#blast_record = NCBIXML.read(result_handle)
#blast_record = next(records)
#blast_record = next(blast_records)
#print blast_records.gi_code()
for blast_record in blast_records:
	for alignment in blast_record.alignments:
		for hsp in alignment.hsps:
			print('****Alignment****')
			print('sequence:', alignment.title)
			print('length:', alignment.length)
			print('e value:', hsp.expect)
			print(hsp.query + '...')
			print(hsp.match + '...')
			print(hsp.sbjct + '...')
			print hsp.sbjct_end
			print hsp.sbjct_start
			print hsp.score


print dir(hsp)
#print dmel_qblast.query
#result_handle = NCBIWWW.qblast("blastn", "nt", "8332116")


#print result_handle.readlines()
#print dir(result_handle)

