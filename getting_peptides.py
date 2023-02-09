# -*- coding: utf-8 -*-

import glob 
import gzip
from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-input_folders", type=str, nargs='+', required=True)
parser.add_argument("-cont_path", type=str, required=False, default = None)
args = parser.parse_args()


# finding the names of the folders with the mzID files and their common path
folders = []
windows_like = 0
for i in args.input_folders:
    tmp = i.split("/")
    if tmp[0] == i:
        tmp = i.split("\\")
        windows_like = 1
    if len(tmp) == 1:
        folders.append(i)
        path = ""
    else:
        folders.append(tmp[-1])
        if windows_like == 0:
            path = i.replace("/"+tmp[-1], "")
        else:
            path = i.replace("\\"+tmp[-1], "")

print("folders:", folders)
print("path:", path)
if type(args.cont_path) == str:
    print("path to contaminants file:", args.cont_path)


# Amino acids accepted by MultiPep
vocab = ['A',
    'R',
    'N',
    'D',
    'C',
    'Q',
    'E',
    'G',
    'H',
    'I',
    'L',
    'K',
    'M',
    'F',
    'P',
    'S',
    'T',
    'W',
    'Y',
    'V',]




# function to read gz files
def read_gz_files(inp):
    f=gzip.open(inp,'rb')
    file_content=f.read()
    f.close()
    file_content = file_content.decode("utf-8")
    file_content = file_content.split("\n")

    return file_content


# 'master function' to retrieve peptides from raw mZid files
def get_relevant_peps(folders):

    dct = {}
    for i in range(len(folders)):
        print(folders[i])
        dct[folders[i]] = {}
        tmp = glob.glob1(path + folders[i]+"/", "*mzid*")
        for ii in tmp:
            print("\t"+ii)
            if ii[-3:] == ".gz":
                lines = read_gz_files(path + folders[i]+"/"+ii)
                dct[folders[i]][ii] = read_peaks_studio_mascot_mzid(lines)
            else:
                with open(path + folders[i] + "/" + ii, "r") as f:
                    lines = f.read()
                    lines = lines.split("\n")
                    dct[folders[i]][ii] = read_peaks_studio_mascot_mzid(lines)
        print("\n")
    return dct




# function that gets the protein name
def find_prot(lin):
    if "HUMAN" in lin:
        tmp = lin.split("_")
        indx = tmp.index("HUMAN")
        tmp = "_".join(tmp[indx - 1:indx + 1])
    else:
        tmp = lin.split("_")[-1]
    return tmp


# The function runs through the mZID files and extracts peptides that
# pass the given threshold
def read_peaks_studio_mascot_mzid(inp, condition = False):

    if condition != False:
        conds = []

    thrs = ""
    spec = ""

    seqsp = {}
    for j,line in enumerate(inp):
        line = line.strip()
        # retrieve peptide ID and peptide sequence
        if '<Peptide id=' in line:
            idd = line.split('"')[-2]
            if idd in seqsp:
                pass
            else:
                seqsp[idd] = ''
            # retrieve associated peptide
            line_tmp = inp[j + 1].strip()
            seqsp[idd] += line_tmp.split("<PeptideSequence>")[-1].split("</PeptideSequence>")[0]


        # retrieving decoy status
        if '<PeptideEvidence ' in line:
            q = line.split(" ")
            idd = [i for i in q if 'peptide_ref=' in i][0].split('"')[-2]
            decoy = [i for i in q if 'isDecoy=' in i][0]
            if decoy not in seqsp[idd]:
                seqsp[idd] += "," + decoy


        # retrieve condition data    
        if condition != False:
            if "<SpectraData" in line:
                #tmp = [i+"_"+line.split('"')[-2] for i in condition if i in line]
                tmp = [line.split('"')[-2] for i in condition if i in line]
                if len(tmp) > 0:
                    conds.append( tmp[0] )


        if "<SpectrumIdentificationItem" in line:
            thrs = ""
            spec = ""
            if "passThreshold" in line and "peptide_ref" in line:
                if condition != False:
                    line_tmp = inp[j-1].strip()
                    if "<SpectrumIdentificationResult" in line_tmp:
                        spec = [i for i in line_tmp.split(" ") if "spectraData_ref=" in i][0].split('"')[1]
                        if spec in conds:
                            q = line.split(" ")
                            idd = [i for i in q if 'peptide_ref=' in i][0].split('"')[-2]
                            thrs = [i for i in q if 'passThreshold=' in i][0].strip(">")
                            #seqsp[idd] += "-" + thrs

                else:
                    q = line.split(" ")
                    idd = [i for i in q if 'peptide_ref=' in i][0].split('"')[-2]
                    thrs = [i for i in q if 'passThreshold=' in i][0].strip(">")


        if '<PeptideEvidenceRef peptideEvidence_ref=' in line:
            if condition != False:
                if spec in conds:
                    tmp = line.split('"')[-2]
                    tmp = find_prot(tmp)
                    if len(tmp) >= 6:
                        thrs += "," + tmp

            else:
                tmp = line.split('"')[-2]
                tmp = find_prot(tmp)
                if len(tmp) >= 6:
                    thrs += "," + tmp

        if '<PeptideEvidenceRef peptideEvidence_ref=' in inp[j-1] and '<cvParam accession' in line:
            if thrs not in seqsp[idd] and thrs != "":
                seqsp[idd] += "," + thrs



    if condition != False:
        return seqsp, conds

    if condition == False:
        return seqsp



# function that reads the "contamination.fasta" file and extracs sequences
def getting_cont_seqs(path):
    
    fasta_sequences = SeqIO.parse(open(path),'fasta')
    fs = []
    for fasta in fasta_sequences:
        fs.append(str(fasta.seq))
    return "_".join(fs)
    


# function that ensures that only peptides that pass the trheshold are used
def extract_info(lst, kk, k):
    prots = set()
    key = 0
    for i in lst:
        if i == 'passThreshold="true"':
            key = 1
        elif i == 'passThreshold="false"':
            key = 0
        elif key == 1:
            prots.add(i) 

    prots = list(prots)
    prots.sort()
    return ",".join([";".join(prots), kk, k])


# add relevant information to extracted peptide: protein, folder, file name, etc.
def tissue_group_peptides(dct, cont):
    
    dct_new = {}
    not_in_vocab = 0
    long_short = 0
    decoys = 0
    not_pass_thr = 0
    noise = 0
    for k,v in dct.items():
        for kk,vv in dct[k].items():
            for kkk,vvv in dct[k][kk].items():
                
                sq = vvv.split(",")[0]
                    
                if len(set(sq).difference(vocab)) == 0:
                    if len(sq) >= 2 and len(sq) <= 200:
                        
                        if 'isDecoy="true"' not in vvv and 'isDecoy="false"' in vvv:
                            if 'passThreshold="true"' in vvv:
                                tmp = vvv.split(",")
                                sq = tmp[0]
                                if sq not in cont and sq[::-1] not in cont: 
                                    if sq not in dct_new:
                                        dct_new[sq] = [extract_info(tmp[2:], kk, k)]
                                    else:
                                        #try:
                                        dct_new[sq].append(extract_info(tmp[2:], kk, k))
                                        #except:
                                        #    print(dct_new[sq], kk, k)
                                        #    assert 1 == 2
                                else:
                                    noise += 1
                            else:
                                not_pass_thr += 1
                        else:
                            decoys += 1 
                    else:
                        long_short += 1
                else:
                    not_in_vocab += 1
                    
    return dct_new, not_in_vocab, long_short,decoys,not_pass_thr,noise



if __name__ == "__main__":

    print("Extracting all relevant information from the raw files...")
    d = get_relevant_peps(folders)

    cont = ""
    if type(args.cont_path) == str:
        print("reading the contamination fasta file...")
        cont = getting_cont_seqs(args.cont_path)

    print("dividing the peptides and associated info based on tissue type...")
    dct_new, not_in_vocab, long_short,decoys,not_pass_thr,noise = tissue_group_peptides(d, cont)

    
    n = ["not_in_vocab", "long_short","decoys","not_pass_thr","noise"]
    n1  = [not_in_vocab, long_short,decoys,not_pass_thr,noise]
    print("removed peptides:")
    for i in range(len(n)):
        print("\t",n[i], n1[i])
    
    
    print("extracting peptides based on tissue..")
    cl = {}
    for k, v in dct_new.items():
        for i in v:
            folder_name = i.split(",")[-1]
            if folder_name not in cl:
                cl[folder_name] = {k}
            else:
                cl[folder_name].add(k)
    
    print("sizes of classes")
    for k,v in cl.items():
        print("\t", k, len(v))
    
    print("Finding tissue-unique peptides...")
    cl_unique = {}
    cl_intersection = {}
    tmp_inter = set()
    cl_keys = list(cl.keys())
    for k, v in cl.items():
        tmp = set()
        others = [i for i in cl_keys if i != k]
        for ii in others:
            tmp.update(cl[ii].copy())
            
        cl_unique[k] = v.difference(tmp)
            
        if len(tmp_inter) == 0:
            tmp_inter = v.copy()
        else:
            tmp_inter = tmp_inter.intersection(v)
            
    cl_intersection["__".join(sorted(cl_keys))] = tmp_inter
        
    
    print("sizes of unique-classes")
    for k,v in cl_unique.items():
        print("\t",k, len(v))
    
    
    print("sizes of intersection-classes")
    for k,v in cl_intersection.items():
        print("\t", k, len(v))
    
    
    print("writing the unique-classes to TXT file...")
    for k, v in cl_unique.items():
        with open(k + ".txt", "w") as f:
            f.write("\n".join(list(v)))
    
    
    print("writing the intersection-classes to TXT file...")
    for k, v in cl_intersection.items():
        with open(k + ".txt", "w") as f:
            f.write("\n".join(list(v)))

