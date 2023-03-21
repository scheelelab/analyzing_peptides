# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import seaborn as sns
sns.set_theme(style="darkgrid", font_scale=5)
from matplotlib import pyplot as plt
import pickle
import argparse
from protein_data_and_functions import get_uniprot_data, prep_info_for_mapper, mapper, global_mapping, compile_dicts_from_uniprot_data

parser = argparse.ArgumentParser()
parser.add_argument("-input_unique", type=str, nargs='+', required=True)
parser.add_argument("-input_inter", type=str, nargs='+', required=True)
parser.add_argument("-multipep_seqs", type=str, required=False, default=None)
parser.add_argument("-toxic", type=str, nargs='+', required=False, default=['hemolytic', 'toxic'], choices=['hemolytic', 'toxic', "None"])
args = parser.parse_args()


# total_classes_seq32.pkl file from the MultiPep repository: https://github.com/scheelelab/MultiPep
if type(args.multipep_seqs) == str:
    with open(args.multipep_seqs,"rb") as f:
        mlpep_tr_data = pickle.load(f)
else:
    mlpep_tr_data = {}

protein_data = "protein_data/"
protein_data_pickle = "protein_data/pickle/"

mulp_files = args.input_unique
mulp_files_inter = args.input_inter


thrs = [0.5, 0.9]

# A function that removes peptides already in the MultiPep training data
def remove_train_data(df, train_data):
    pred_already = []
    for i in list(df.index):
        if i in train_data:
            pred_already.append(0)
        else:
            pred_already.append(1)
            
    return np.array(pred_already) == 1

# A function that removes peptides not classified as any class
def find_relevant_peptides(out_dfs, th=0.5):

    sq_dct = {}
    for k,v in out_dfs.items():
        tmp = v.values
        tmp = np.sum(tmp > th,axis=-1) > 0
        sq_dct[k] = list(v[tmp].index)
    return sq_dct

# getting MultiPep predicted peptides based on thresholds ('thrs')
def retrieve_data(mulp_files, thrs):
    out_dct = {}
    out_dct_prcnt = {}
    out_dfs = {}
    out_no_peps = {}
    for i in mulp_files:
        name = i.replace(".csv","")
        # removing peptides in MultiPep training data
        preds = pd.read_csv(i, index_col=0)
        bool_ = remove_train_data(preds, mlpep_tr_data)
        preds = preds[bool_]
        out_dfs[name] = preds
        
        # Getting column names (bioactivity classes)
        if i == mulp_files[0]:
            cols = list(preds.columns)
        
        obv = preds.values
        out_dct[name] = {}
        out_dct_prcnt[name] = {}
        out_no_peps[name] = {}
        for ii in thrs:
            
            obv_tmp = obv > ii
            no_peps = np.sum(np.sum(obv_tmp, axis=-1) > 0)
            sc = np.sum(obv_tmp, axis=0)
            out_dct[name][ii] = sc
            out_dct_prcnt[name][ii] = (sc / np.sum(sc)) * 100
            out_no_peps[name][ii] = no_peps

    return out_dct, out_dct_prcnt, out_dfs, cols, out_no_peps


def covert_to_df(dct, cols, names, remove=0):
    df = pd.DataFrame(dct)
    clss = cols.copy()
    clss[-8] = 'cytokines/growthf.'
    clss[-1] = 'dipeptidyl pept. inh.'
    df["Bioactivity classes"] = clss
    bool_ = np.sum(df[names].values,axis=-1) > remove
    return df, bool_
    

# A very hard-coded plotting function
def plot_predictions(dct, cols, names, thrs, dct_pno, remove=0, extra_name = False, color = False):
    
    fig, axes = plt.subplots(1, 2, figsize=(40, 20), sharey=True)
    for j,i in enumerate(thrs):
        

        tmp = {}
        for ii in names:
            tmp[ii] = dct[ii][i]

        if type(extra_name) == list:
            legend_names = ["{} ({})".format(extra_name[ii], dct_pno[names[ii]][i]) for ii in range(len(names))]
        else:
            legend_names = ["{} ({})".format(names[ii].split("_")[0], dct_pno[names[ii]][i]) for ii in range(len(names))]

        if i == 0.5:
            df, b = covert_to_df(tmp, cols, names, remove)
        else:
            df, _ = covert_to_df(tmp, cols, names, remove)

        df1 = pd.melt(df[b], id_vars="Bioactivity classes", var_name="Tissue", value_name="Percent")

        if type(color) == str:
            f1 = sns.barplot(data=df1, ax=axes[j], x="Bioactivity classes", y="Percent", hue="Tissue",
                     palette=[color] )
        else:
            f1 = sns.barplot(data=df1, ax=axes[j], x="Bioactivity classes", y="Percent", hue="Tissue",
                        palette="dark")

        f1.set_xticklabels(
            f1.get_xticklabels(), 
            rotation=45, 
            horizontalalignment='right',
            )
            

        if i != 0.5:
            cols_tmp = list(np.array(cols)[b])
            _tmp = _[b]
            to_red = [cols_tmp.index(clas) for clas in np.array(cols_tmp)[_tmp == False]]
            for clas in to_red:
                axes[j].get_xticklabels()[clas].set_color("red")
            axes[j].set(ylabel=None)

        axes[j].set_title("Prediction scores > {}".format(i), fontsize="x-large")

        leg = axes[j].get_legend()
        for t, l in zip(leg.texts, legend_names):
            t.set_text(l)
        leg.set_title("")
    fig.tight_layout()
    plt.subplots_adjust(wspace=0.1)

    plt.savefig("_".join(names)+".tif", bbox_inches='tight', dpi=350)



# remove peptide classified as toxic or hemolytic:
#tox = ['hemolytic', 'toxic']
tox = args.toxic

# the classes that should be analyzed further:
bio = ['neuropeptide', 'peptidehormone']

# UniProt info that should be extracted:
relevant = ['Entry', 'Gene Names', 'Tissue specificity', 'Subcellular location [CC]', 
            'Propeptide','Signal peptide', 'Transit peptide', 'Peptide', 
            'Gene Ontology (biological process)']

# names of columns with protein position information:
st_st = ["Start", "Stop"]




def make_small_table(out_dfs, lst, toxic, relevant, pep_map, th = 0.5, one_match=True):

    """
    1. Predicted toxic peptides are removed. 
    2. Looking at predicted 'neuropeptide' and 'peptidehormone' ppetides
    3. Finding protein information of peptides that mach to s single protein (if 'one_match=True')
    4. generating a table with relevant information
    """    


    out_lst = {}
    for k,v in out_dfs.items():

        # removing toxic peptides
        if "None" in toxic:
            tmp_df = v[lst]
        else:
            tmp_tox = np.sum(v[toxic].values <= 0.5, axis=-1) == len(toxic)
            tmp_df = v[tmp_tox][lst]

        # getting sorted peptides with prediction above 'th'
        tmp_val = tmp_df.values
        tmp_indx = np.array(tmp_df.index)
        tmp_mean_arg = np.argsort( np.mean(tmp_val,axis=-1) )[::-1]
        tmp_indx = tmp_indx[tmp_mean_arg]
        tmp_df = np.round(tmp_df.loc[tmp_indx],2)
        tmp_preds = np.sum(tmp_df.values > th,axis=-1) > 0
        tmp_df = tmp_df[tmp_preds]
        print(k, tmp_df.shape)
        use = np.ones(len(tmp_df))
        # getting protein information
        to_df = []
        for j,i in enumerate(list(tmp_df.index)):
            # Only using the first matched protein, attaching all matches in the final df
            if i in pep_map:    
                prot = pep_map[i]
            else:
                use[j] = 0
                continue
            if one_match == True:
                if len(prot) > 1:
                    use[j] = 0
                    continue
                else:
                    all_matches = [";".join([ii.split("|")[0] for ii in prot])]
                    prot, start_stop = prot[0].split("|")
                    start, stop = start_stop.split("_")
                    to_df.append(list(prot_rev[prot_rev["Entry"] == prot][relevant].values[0]) + [start, stop] + all_matches)
            else:
                all_matches = [";".join([ii.split("|")[0] for ii in prot])]
                prot, start_stop = prot[0].split("|")
                start, stop = start_stop.split("_")
                to_df.append(list(prot_rev[prot_rev["Entry"] == prot][relevant].values[0]) + [start, stop] + all_matches)
        df = pd.DataFrame(to_df, columns=relevant + st_st + ["all_matches"], index=tmp_df.index[use == 1] )
        print(k, df.shape)
        out_lst[k] = pd.concat([tmp_df[use == 1], df], axis=1)
    return out_lst
    

# The function finds the top X predicted peptides. the 'bio' list contain the relavant peptides
# In addition, information from the UniProt database is handled
def get_top_x(df, clmns, X, st_st):

    df1 = df[clmns].loc[df.index[:X]].copy()
    dfst = df[st_st].loc[df.index[:X]].copy()
    
    peps = list(df1["Peptide"])
    new_peps = []
    for i in peps:
        if type(i) == float:
            new_peps.append("Not known")
        else:
            new_peps.append( ";".join([ii.replace("PEPTIDE","") for ii in i.split(";") if "PEPTIDE" in ii]) )
    df1["Peptide"] = new_peps
    
    subc = list(df1['Subcellular location [CC]'])
    new_sub = []
    for i in subc:
        if type(i) == float:
            new_sub.append("Not known")
        else:
            tmp = i.split("SUBCELLULAR LOCATION: ")
            new_sub.append(";".join([ii.split("{")[0] for ii in tmp if len(ii) > 0]))
    df1['Subcellular location [CC]'] = new_sub
    
    ss = []
    for i in zip(dfst["Start"], dfst["Stop"]):
        ss.append("..".join(i))
    df1["Start..stop"] = ss
    
    renm = {"neuropeptide":"NE", "peptidehormone":"PH"}
    df1 = df1.rename(columns=renm)
    
    reorder = ['Entry',
     "Start..stop",
     'NE',
     'PH',
     'Gene Names',
     'Subcellular location [CC]',
     'Peptide']
    
    df1 = df1[reorder]
    #sig = list(df1['Signal peptide'])
    #new_sig = []
    #for i in sig:
    #    if type(i) == float:
    #        new_sig.append("Not known")
    #    else:
    #        new_sig.append(";".join([ii for ii in i.split(";") if "SIGNAL" in ii]))
    #df1['Signal peptide'] = new_sig

    return df1
    



if __name__ == "__main__":
    
    print("Retieving data for unique brain and plasma peptides...")
    out_dct, out_dct_prcnt, out_dfs, cols, out_no_peps = retrieve_data(mulp_files, thrs)
    names = list(out_dct_prcnt.keys())
    plot_predictions(out_dct_prcnt, cols, names, thrs, out_no_peps, remove=0)

    print("Retieving data for peptides in both the plasma and brain datasets...")
    out_dct_bp, out_dct_prcnt_bp, out_dfs_bp, cols_bp, out_no_peps_bp = retrieve_data(mulp_files_inter, thrs)
    names_bp = list(out_dct_prcnt_bp.keys())
    plot_predictions(out_dct_prcnt_bp, cols_bp, names_bp, thrs, out_no_peps_bp, remove=0, extra_name=["brain_plasma"], color="green" )
    
    print("Finding peptides with a prediction score > 0.5...")
    sq_dct = find_relevant_peptides(out_dfs, th=0.5)
    
    sq_dct_bp = find_relevant_peptides(out_dfs_bp, th=0.5)
    
    print("Readying UniProt protein information...")
    prot_rev = get_uniprot_data(protein_data, rev=True, iso=False)
    #prot_unrev = get_uniprot_data(protein_data, rev=False, iso=False)
    #iso = get_uniprot_data(protein_data, rev=True, iso=True)
    prot_rev_dct = compile_dicts_from_uniprot_data(protein_data_pickle, 
                                                   "prot_rev_dct.pkl", 
                                                   list(prot_rev["Entry"]), 
                                                   list(prot_rev["Sequence"]))

    print("Starting the peptide-protein mapping...")
    mapped_peps = {}
    for k,v in sq_dct.items():
        pep_map, pep_not_found = global_mapping(v, protein_data_pickle, prot_rev_dct)
        #print(k, "found {}, not found {}".format(len(pep_map), len(pep_not_found)))
        #mapped_peps[k] = [pep_map, pep_not_found]
        
    mapped_peps_bp = {}
    for k,v in sq_dct_bp.items():
        pep_map, pep_not_found = global_mapping(v, protein_data_pickle, prot_rev_dct)
        #print(k, "found {}, not found {}".format(len(pep_map), len(pep_not_found)))
        #mapped_peps[k] = [pep_map, pep_not_found]
    
    
    print("Retrieving the relevant information and filtering away protein mapped to more than 1 protein...")
    dct_dfs = make_small_table(out_dfs, bio, tox, relevant, pep_map)
    dct_dfs_bp = make_small_table(out_dfs_bp, bio, tox, relevant, pep_map)
    
    
    relevant_top_x = ["Entry"] + bio + ['Gene Names', 'Subcellular location [CC]', 
                'Peptide']
    
    X = 5
    print("Generating the top {} tables".format(X))
    for k,v in dct_dfs.items():
        df1 = get_top_x(v, relevant_top_x, X, st_st)
        df1.to_excel("top_{}_{}.xlsx".format(X, k))
    
    for k,v in dct_dfs_bp.items():
        df1 = get_top_x(v, relevant_top_x, X, st_st)
        df1.to_excel("top_{}_{}.xlsx".format(X, k))
    
    
