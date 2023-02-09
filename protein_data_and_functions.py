# -*- coding: utf-8 -*-
import requests
import os
import pandas as pd
import pickle
import numpy as np
import re
import gzip
from Bio import SeqIO
import time
from multi_key_dict import multi_key_dict






def get_uniprot_data(path, rev=True, iso=False, organism="HUMAN"):

    name = "uniprot_"

    if not os.path.exists(path):
        os.makedirs(path)


    if iso == False:
        if rev == True:
            print("getting reviewed proteins...")
            name += "rev_true.csv"

            if os.path.exists(path+name):
                return pd.read_csv(path+name, sep=None, engine="python")
            else:
                url = "https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Creviewed%2Cid%2Cprotein_name%2Cgene_names%2Corganism_name%2Ccc_tissue_specificity%2Cft_topo_dom%2Cft_intramem%2Ccc_subcellular_location%2Cft_transmem%2Cft_propep%2Cft_signal%2Cft_transit%2Ccc_ptm%2Cft_init_met%2Cft_lipid%2Cft_mod_res%2Cft_peptide%2Cft_carbohyd%2Cft_disulfid%2Cft_crosslnk%2Cft_chain%2Cgo_p%2Cgo_c%2Cgo_f%2Cgo%2Cgo_id%2Csequence%2Ccc_alternative_products%2Cft_var_seq&format=tsv&query=%28%28organism_id%3A9606%29%20AND%20%28reviewed%3Atrue%29%29%20AND%20%28reviewed%3Atrue%29"

                with open(path+name, "w") as f:
                    f.write(requests.get(url).text)

                return pd.read_csv(path+name, sep=None, engine="python")

        else:
            print("getting un-reviewed proteins...")
            name += "rev_false.csv"

            if os.path.exists(path+name):
                return pd.read_csv(path+name, sep=None, engine="python")

            else:
                url = "https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Creviewed%2Cid%2Cprotein_name%2Cgene_names%2Corganism_name%2Ccc_tissue_specificity%2Cft_topo_dom%2Cft_intramem%2Ccc_subcellular_location%2Cft_transmem%2Cft_propep%2Cft_signal%2Cft_transit%2Ccc_ptm%2Cft_init_met%2Cft_lipid%2Cft_mod_res%2Cft_peptide%2Cft_carbohyd%2Cft_disulfid%2Cft_crosslnk%2Cft_chain%2Cgo_p%2Cgo_c%2Cgo_f%2Cgo%2Cgo_id%2Csequence%2Ccc_alternative_products%2Cft_var_seq&format=tsv&query=%28%28organism_id%3A9606%29%20AND%20%28reviewed%3Afalse%29%29%20AND%20%28reviewed%3Afalse%29"

                with open(path+name, "w") as f:
                    f.write(requests.get(url).text)

                return pd.read_csv(path+name, sep=None, engine="python")


    else:
        print("getting isoforms...")
        name += "isoform_all_organisms.fasta"

        if os.path.exists(path+name):
            fasta_sequences = SeqIO.parse(open(path+name),'fasta')
            fs = {}
            for fasta in fasta_sequences:
                if organism in fasta.id.split("|")[-1]: 
                    fs[fasta.id.split("|")[1]] = str(fasta.seq)
            return fs
        else:
            url = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot_varsplic.fasta.gz"
            iso = requests.get(url)
            iso = gzip.decompress(iso.content).decode("utf-8")
            with open(path+name, "w") as f:
                f.write(iso)

            fasta_sequences = SeqIO.parse(open(path+name),'fasta')
            fs = {}
            for fasta in fasta_sequences:
                if organism in fasta.id.split("|")[-1]:
                    fs[fasta.id.split("|")[1]] = str(fasta.seq)
            return fs



def prep_info_for_mapper(dct):
    long_sq = "_"
    lookup = multi_key_dict()

    for k,v in dct.items():
        lensq = len(long_sq)
        lookup[tuple(range(lensq, lensq + len(v) + 1))] =  [k,lensq]
        long_sq += v + "_"

    return long_sq, lookup
    

def mapper(peptide_list, longsq, lookup, out, not_found, 
           seen_before=False, add=0):

    #out = {}
    #not_found = []

    for k in peptide_list:
        if seen_before:
            if k in out:
                continue
        ln_k = len(k)
        mapz = [m.start() for m in re.finditer(k, longsq)]
        found = []
        
        for i in mapz:
            tmp = lookup[i]
            start = i-tmp[1]
            stop = start + ln_k - 1
            found.append(tmp[0] + "|{}_{}".format(start+add,stop+add))
        
        if len(found) == 0:
            not_found.append(k)
        
        elif len(found) > 0:
            if k in out:
                diff = list(set(found).difference(out[k]))
                if len(diff) > 0:
                    out[k] += diff
            else:
                out[k] = found
        
        else:
            not_found.append(k)

    # udating 'not_found' list
    not_found = set(not_found)
    not_found = not_found.difference( set(out) )
    return out, list(not_found)


def global_mapping(peptides, path, prot_rev_dct, #iso, prot_unrev_dct,
                   add=1, seen_rev=False, seen_iso=False,
                   seen_unrev=False):
    
    names = ["long_sq_lookup_rev.pkl", 
             "long_sq_lookup_iso.pkl",
             "long_sq_lookup_unrev.pkl",
             "mapper_dct.pkl"]
    
    
    if not os.path.exists(path):
        print("creating the path: '{}'".format(path))
        os.makedirs(path)
        
    if not os.path.exists(path+names[-1]):
        print("Creating the mapper dictionary and list...\n")
        pep_map, pep_not_found = {}, []
        
    if os.path.exists(path+names[-1]):
        print("Retrieving the mapper dictionary...")
        with open(path+names[-1], "rb") as f:
            pep_map, pep_not_found= pickle.load(f)
            seen_rev=True 
            seen_iso=True
            seen_unrev=True
    


    print("Handling the reviewed proteins...")
    if os.path.exists(path+names[0]):
        print("\tRetrieving the rev_lookup dictionary...")
        with open(path+names[0], "rb") as f:
            long_sq, lookup = pickle.load(f)    

    if not os.path.exists(path+names[0]):
        print("\tCreating the rev_lookup dictionary...")
        long_sq, lookup = prep_info_for_mapper(prot_rev_dct)
        with open(path+names[0], "wb") as f:
            pickle.dump([long_sq, lookup], f)

    print("\nMapping to reviewed proteins...")
    pep_map, pep_not_found = mapper(peptides, long_sq, lookup, 
                                            pep_map, 
                                            pep_not_found, 
                                            seen_before=seen_rev,
                                            add=add)
    
    
    #print("\nHandling the isoform proteins...")
    #if os.path.exists(path+names[1]):
    #    print("\tRetrieving the iso_lookup dictionary...")
    #    with open(path+names[1], "rb") as f:
    #        long_sq, lookup = pickle.load(f)

    #if not os.path.exists(path+names[1]):
    #    print("\tCreating the iso_lookup dictionary...")
    #    long_sq, lookup = prep_info_for_mapper(iso)
    #    with open(path+names[1], "wb") as f:
    #        pickle.dump([long_sq, lookup], f)
    
    #print("\nMapping to isoforms...")
    #pep_map, pep_not_found = mapper(peptides, long_sq, lookup, 
    #                                    pep_map, 
    #                                    pep_not_found, 
    #                                    seen_before=seen_iso,
    #                                    add=add)
    
    ##if len(pep_not_found) > 0:
    #print("\nHandling the unreviewed proteins...")
    #if os.path.exists(path+names[2]):
    #    print("\tRetrieving the unrev_lookup dictionary...")
    #    with open(path+names[2], "rb") as f:
    #            long_sq, lookup = pickle.load(f)
            
    #if not os.path.exists(path+names[2]):
    #    print("\tCreating the unrev_lookup dictionary...")
    #    long_sq, lookup = prep_info_for_mapper(prot_unrev_dct)
    #    with open(path+names[2], "wb") as f:
    #            pickle.dump([long_sq, lookup], f)
        
    #print("\nMapping to unreviewed proteins...")
    #pep_map, pep_not_found = mapper(peptides, long_sq, lookup, 
    #                                pep_map, 
    #                                pep_not_found, 
    #                                seen_before=seen_unrev,
    #                                add=add)


    print("\nSaving the mapped data...")
    with open(path+names[-1], "wb") as f:
        pickle.dump([pep_map, pep_not_found], f)

    return pep_map, pep_not_found
    
    #todo: make a global_mapping_seen_before that can 
          #be used when a run has been made once
          # do this with the for the unrev proteins

def compile_dicts_from_uniprot_data(path, name, lst1, lst2):

    if not os.path.exists(path):
        os.makedirs(path)

    if os.path.exists(path+name):
        with open(path+name, "rb") as f:
            new_dct = pickle.load(f)
        return new_dct

    else:
        tmp = {}
        for j,i in enumerate(lst1):
            tmp[i] = lst2[j]

        with open(path+name, "wb") as f:
            pickle.dump(tmp, f)

        return tmp





"""
if __name__ == "__main__":
    prot_rev = get_uniprot_data(protein_data, rev=True, iso=False)
    #prot_unrev = get_uniprot_data(protein_data, rev=False, iso=False)
    #iso = get_uniprot_data(protein_data, rev=True, iso=True)
    prot_rev_dct = compile_dicts_from_uniprot_data(protein_data_pickle, 
                                                   "prot_rev_dct.pkl", 
                                                   list(prot_rev["Entry"]), 
                                                   list(prot_rev["Sequence"]))

    pep_map, pep_not_found = global_mapping(sqs, protein_data_pickle, prot_rev_dct)
"""