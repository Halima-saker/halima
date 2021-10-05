from Bio import Entrez, SeqIO
import urllib
import webbrowser
import fileinput
#webbrowser.open('http://napp.u-psud.fr/SqlGff.php?specie=355&SpecieName=Pyrococcus_horikoshii')
import os.path
fo = open("G:/master_2/2eme semestre_project/halima__/fasta_file/testt.fas", "w")
in_file = open("C:/Users/User/Desktop/seq.txt")
in_handle = in_file.read()                  
for lines in in_handle.splitlines():
    values = lines.split("\t")
    print ("sequence__________::::")
    for line in fileinput.input(['C:/Users/User/Desktop/Desktop/org.gff']):
        valuees = line.split("\t")
        if valuees[0]==values[3]:#and float(values[2])>0.5:                             
            chrom=valuees[2][:-4]
            start_end=values[1].split("..")
            #if start_end[0]>start_end[1]:
            startt=int(start_end[0])-30
            endd=int(start_end[1])+30
            print (chrom)
            print (startt)
            print (endd)
                            # else:
                                # startt=start_end[0]
                                # endd=start_end[1]
            handle = Entrez.efetch(db="nuccore",
							id=chrom,
							rettype="gb",
							retmode="text",
							seq_start=startt, 
							seq_stop=endd
							)
            whole_sequence = SeqIO.read(handle, "genbank")
            print (whole_sequence.seq)
            fo.write(">"+chrom+"\n"+str(whole_sequence.seq)+"\n")
                                 
fo.close()
