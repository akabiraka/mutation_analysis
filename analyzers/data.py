import sys
sys.path.append("../mutation_analysis")

import pandas as pd


def print_class_distribution(df):
    n_mutation_rows = df.shape[0]
    print("#-mutations: {}".format(n_mutation_rows))

    n_pdb_ids = len(df["pdb_id"].dropna().unique().tolist())
    n_uniprot_ids = len(df["uniprot_id"].dropna().unique().tolist())
    n_pubmed_ids = len(df["pubmed_id"].dropna().unique().tolist())
    n_proteins = len(df["protein"].dropna().unique().tolist())
    n_inv_pdb_ids = len(df["inverse_pdb_id"].dropna().unique().tolist())
    
    print("n_pdb_ids: {}, n_uniprot_ids: {}, n_pubmed_ids:{}, n_proteins:{}, n_inv_pdb_ids:{}".format(n_pdb_ids, n_uniprot_ids, n_pubmed_ids, n_proteins, n_inv_pdb_ids))
    
    n_destabilizing = df[df["ddg"]<-1.0].shape[0]
    n_destabilizing_hard = df[df["ddg"]<=-5.0].shape[0]
    print("n_destabilizing: {}, n_destabilizing_hard: {}".format(n_destabilizing, n_destabilizing_hard))

    n_neutral = df[(df["ddg"]<=1.0) & (df["ddg"]>=-1.0)].shape[0]
    print("n_neutral: {}".format(n_neutral))

    n_stabilizing = df[df["ddg"]>1.0].shape[0]
    n_stabilizing_hard = df[df["ddg"]>=5.0].shape[0]
    print("n_stabilizing: {}, n_stabilizing-hard: {}".format(n_stabilizing, n_stabilizing_hard))
    
    
    print("#-mutations: {}".format(n_destabilizing+n_neutral+n_stabilizing))
    

print_class_distribution(pd.read_csv("data/clean_1/ProThermDB.csv"))
# print_class_distribution(pd.read_csv("data/clean_1/PON_TStab.csv"))
# print_class_distribution(pd.read_csv("data/clean_1/I_Mutant_2_seq.csv"))
# print_class_distribution(pd.read_csv("data/clean_1/I_Mutant_2_structure.csv"))
# print_class_distribution(pd.read_csv("data/clean_1/PoPMuSiC_2.csv"))
# print_class_distribution(pd.read_csv("data/clean_1/Ssym.csv"))
# print_class_distribution(pd.read_csv("data/clean_1/Saraboji_S1396.csv"))
# print_class_distribution(pd.read_csv("data/clean_1/Saraboji_S2204.csv"))
# print_class_distribution(pd.read_csv("data/clean_1/iPTREE_STAB.csv"))
# print_class_distribution(pd.read_csv("data/clean_1/AUTOMUTE_S1925.csv"))
# print_class_distribution(pd.read_csv("data/clean_1/AUTOMUTE_S1962.csv"))
# print_class_distribution(pd.read_csv("data/clean_1/Broom.csv"))
# print_class_distribution(pd.read_csv("data/clean_1/ThermoMutDB_single.csv"))


