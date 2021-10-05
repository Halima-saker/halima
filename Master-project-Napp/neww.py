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
#from StringIO import StringIO
from Bio import Alphabet
from Bio.Align import MultipleSeqAlignment
from BCBio.GFF import GFFExaminer
from BCBio import GFF
import pprint
import os.path
from Bio import SeqIO
from Bio import Entrez
import urllib
import webbrowser
import fileinput
#import requests
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO
#import Bio.Clustalw
from Bio.Alphabet import IUPAC
from sys import *
#import urllib2
import urllib.request
from Bio.Align import MultipleSeqAlignment
from Bio import Entrez
Entrez.email = "sample@example.org"

def get_seq(d):
    fo = open("G:/master_2/2eme semestre_project/halima__/infos/"+d, "w")
    c=0
    erra=dict()
    id_spe=0
    alll=dict()
    all_blocks=dict()
    #in_file = "/homes/biertank/halima/Downloads/halima_saker_project_/database/"+d
    in_file = "G:/master_2/2eme semestre_project/halima__/"+d
    examiner = GFFExaminer()
    in_handle = open(in_file)       
    for rec in GFF.parse(in_handle):
        t=0
        #print c
        for record in rec.features:
            elem_metadatas = list()
            keyss=('clst_id', 'SubjectScore','SubjectOrganism')
        #print record.type
            for key in sorted(record.qualifiers.keys()):
                if key in keyss:
                    if key=="SubjectOrganism":  
                        load_profile = open('G:/master_2/2eme semestre_project/halima__/speciess.gff')
                        len_org=len(record.qualifiers["SubjectOrganism"])               
                        read_it = load_profile.read()
                        myLine =list()
                        myscore=list()
                        purse=dict()
                        for val in record.qualifiers[key]:
                            i=1
                            for line in read_it.splitlines():
                                if line == record.__dict__['qualifiers']['source'][0]:
                                    id_spe=i
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
                if float(myscore[i])>0.5:
                    purse[myLine[i]]=str(myscore[i])    
            alll[rec.id,record.qualifiers["clst_id"][0]]=purse              
            for i in range(0,len(myLine)):                      
                t=t+float(myscore[i])
            print (purse)    
            min=0
            max=1
            u=0
            z=1
            lst=sorted(purse.keys())
            blocks=dict()
            list_elmt=list()
            for u in range(0,len(lst)-1):
                #print lst[u]
                #print min, "--",  max
                if lst[max]-lst[max-1]<=4:      
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
            #pprint.pprint(purse)
        c=c+1
        
        
        
        
        
        
        
        
        
        
        #print (c,"\n",rec.id,"\n",d,"\n",id_spe)
        #print blocks
        bboolean="true"
        for w in range (1,len(blocks)+1):       
            if id_spe in blocks[w]:
                ind_block_spe=w
                len_block_spe=len(blocks[w])
                break
            else:
                w=0
        pp=0
        #print rec.id
        #print w
        #print len_block_spe,"_____",len_org
        tot=0
        inter_block=list()
        inter_org=list()        
        if w==0:
            bboolean="false"
            #print id_spe,"_____________________________________________________________"
        else:
            if float(len_block_spe)/float(len_org) >= 0.5:
                bboolean="true"
            #   print id_spe,"_____________________________________________________________"
            else:
                for f in range(1,len(blocks)+1):
                    if f != w:
                        for j in blocks[f]:
                            #print float(purse[j])
                            tot=tot+float(purse[j])
                        
                        if tot/len(blocks[f])>0.45:
                            #print tot/len(blocks[f])
                            #print id_spe,"_____________________________________________________________"
                            #print tot/len(blocks[f])
                            inter_block.append(str(f))
                            #print blocks[f]
                        #if sum(x for x in (myscore[]))/len(blocks[pp]) < 0.45
                        tot=0
                for ff in myLine:
                    
                    lstt=sorted(purse.keys())
                    #print lstt
                    #print len(lstt)
                    #print len(myLine)
                    if ff not in lstt:
                        if purse[ff]>0.45 :
                            
                            inter_org.append(str(ff))
                #print inter_block,"___",inter_org      
                if len(inter_block)==0 and len(inter_org)==0:
                    bboolean="true"
                else:
                    #print id_spe,"_____________________________________________________________"
                    bboolean="false"            
                        
                    
                            
                
                
                    
        #print ind_block_spe    
        #print rec.id,"____",len_org, "__",len_block_spe, "__",id_spe, "__",blocks,"___",sorted(purse.iterkeys()),"____",len(blocks)
        #print alll
        #fo.write(str(all_blocks))
        #fo.write(rec.id)
        #fo.write("")
        if bboolean=="false":
            erra[rec.id]=record.__dict__['qualifiers']['source'][0]    
    #   print bboolean, "____",w
    #print dir(rec)
    #print record.__dict__['qualifiers']['source']
        
    print (erra)   
    if erra:
            for idd in erra:
                org_pos=dict()
                fol = open("G:/master_2/2eme semestre_project/halima__/fasta_file/"+idd+".fas", "w")
                i=0
                print (erra[idd])
                for line in fileinput.input(['C:/Users/User/Desktop/Desktop/profiles.gff']):
                    values = line.split("\t")
                    if values[0]==idd:# in values:
                        org_pos[values[3]]=values[1]
                        #fo.write(idd+">\n"+values[1])
                        
                        print (line)
                        #print ("___________"+str(i)+"____________")
                        #break
        
                print (org_pos)   
                for za in org_pos:#range(0, len(org_pos)):
                    print ("sequence__________::::")
                    for line in fileinput.input(['C:/Users/User/Desktop/Desktop/org.gff']):
                         valuees = line.split("\t")
                         if valuees[0]==za:                             
                             chrom=valuees[3][:-4]
                             start_end=org_pos[za].split("..")
                             #if start_end[0]>start_end[1]:
                             startt=start_end[0]
                             endd=start_end[1]
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
                             fol.write(">"+org_pos[za]+"\n"+whole_sequence.seq+"\n")
                fol.close()
                

    in_handle.close()
    fo.close()


        
#for rec in GFF.parse(in_handle, limit_info=limit_info):
# 
#files=('xaa', 'xab','xac', 'xad', 'xae', 'xaf', 'xag', 'xah', 'xai', 'xaj', 'xak', 'xal', 'xam', 'xan', 'xao', 'xap', 'xaq', 'xar', 'xas', 'xat')

k=os.listdir("G:/master_2/2eme semestre_project/halima__/rec_spe/")
files=("kfkk.gff",)
#print (k)
for d in files:
    get_seq(d)
    

