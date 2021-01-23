#Python script to pick probes from the user provided consensus region of a DNA seq
#Inatall primer3, Biopython, re libraries
#input filename should be test.fasta
#to add decremental from 5' end
import primer3
from Bio.Seq import Seq
from Bio.SeqUtils import GC
from Bio import SeqIO
import re
import collections
import os
import sys

ipname = str(sys.argv)
wkdir = os.getcwd()
#IPfile1 = wkdir + '/' + str(sys.argv[1])
IPfile1 = str(sys.argv[1])
print ("Input File:",IPfile1)

output = open("output_4.txt", "a")


min_len = 15
value = 0
prefix_nt = "T"
poly_nt = 15

for seq_record in SeqIO.parse(IPfile1, "fasta"):
    probe_bind = str(seq_record.seq.reverse_complement())
    probe_non_nt = re.search(r"[^ATGC]", probe_bind)
    if probe_non_nt:
        print('probe is not a nucleotide sequence')
        non_nt = probe_confirm.group()
        print("The non nt is " + non_nt)
    else:
        probe_list=[]
        long_seq_list=[]
        probe_bind = str(seq_record.seq.reverse_complement())
        probe_rc = probe_bind
        probe_fc = probe_bind
        probe_bind_len = len(probe_rc)
        probe_bind_len2 = len(probe_rc)
        probe_list.append(probe_bind)
        if(probe_bind_len <= 60):
           idx_sl = 0
           idx_rc = 0
           idx_fc = 0
           for i in range(0, len(probe_bind)-min_len, 1):
               probe_list.append(probe_bind[i:i+20])
               idx_sl = idx_sl + 1
           while probe_bind_len >= min_len:
               probe_rc = probe_rc[:-1]
               probe_list.append(probe_rc)
               probe_bind_len = len(probe_rc)
               idx_rc = idx_rc + 1
           while probe_bind_len2 >= min_len:
               probe_fc = probe_fc[1:]
               probe_list.append(probe_fc)
               probe_bind_len2 = len(probe_fc)
               idx_fc = idx_fc + 1
        else:
           long_seq_list.append(seq_record.id)
    
    idx = 1
    short = 0
    for seq_probe in probe_list:    
        probe = (prefix_nt * poly_nt) + str(seq_probe)
        probe_length = len(probe)
        #Tm_ori = primer3.calcTm(probe)
        #Tm_ori2 = ("%.2f" % Tm_ori)
        #GC_percent_ori = GC(probe)
        #GC_per_2 = ("%.2f" % GC_percent_ori)
        #print(seq_record.id + "_" + str(idx) +  "\t" + probe + "\t" + str(GC_per_2) + "\t" +  str(Tm_ori2) + "\t", end = '')
        #if(probe_length <=40):
        #    short = short+1
       # else:
        #    print("\n")
        if(probe_length <= 60):
            Tm = primer3.calcTm(probe)
            Tm2 = ("%.2f" % Tm)
            GC_percent = GC(probe)
            GC_per = ("%.2f" % GC_percent)
            Homodimer = primer3.calcHomodimer(probe)
            Hairpin  =  primer3.calcHairpin(probe)
            print(seq_record.id + "_" + str(idx) +  "\t" + probe + "\t" + "\t" + str(GC_per) + "\t" + str(Tm2) + "\t" + str(Hairpin.structure_found) + "\t" + str(Homodimer.structure_found))
        idx = idx+1




output.close()
