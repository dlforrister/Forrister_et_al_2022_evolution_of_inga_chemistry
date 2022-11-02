import pandas as pd
import numpy as np
from sklearn.neighbors.kde import KernelDensity
from scipy.signal import argrelextrema
from bisect import bisect_left
from datetime import date
import os

#set working directory
os.chdir("K://Lab_Map//git_repositories//Evolution_Of_Inga_Chemistry-master//data")

# define names of your two columns / methods -- for example, "C18" and "amide"
method_1 = "C18"
method_2 = "amide"

# set paths to compound feature tables for both methods
mgf_path_1 = ".//C18_one_per_comp.mgf"
mgf_path_2 = ".//polar_one_per_comp.mgf"

# set paths to compound_feature tables for both methods
cft_path_1 = ".//C18//c18_compound_feature_table.csv"
cft_path_2 = ".//polar//polar_compound_feature_table.csv"

# load mgf for first method
with open(mgf_path_1) as temp:
    mgf_1 = []
    for line in temp:
        mgf_1.append(line.rstrip('\n'))

# load compound_feature table for first method
with open(cft_path_1) as temp:
    cft_1 = {}
    for line in temp:
        comp_num, feat_num, TIC, rel_abund, mz, rt = line.rstrip('\n').split(',')
        if comp_num == "compound_number":
            continue
        if comp_num in cft_1:
            cft_1[comp_num]["mz"] += [mz]
            cft_1[comp_num]["TIC"] += [TIC]
        else:
            cft_1[comp_num] = {
                "mz": [mz],
                "TIC": [TIC]
            }

# convert to dictionary containing both ms features (from compound_feature table) and MSMS features (from .mgf)
compound = {}
for comp in cft_1.keys():
    current_compound = method_1 + "_" + comp
    compound[current_compound] = {
        "ms_features": cft_1[comp]['mz'],
        "ms_feat_abund": cft_1[comp]['TIC'],
        "msms_features":[],
        "msms_feat_abund":[],
        "parent_ion":[]
    }
    for i,line in enumerate(mgf_1):
        if line.split('=')[0] == "SCANS" and line.split('=')[1] == comp:
            compound[current_compound]["parent_ion"] = float(mgf_1[i-1].split('=')[1])
            ion = i + 4
            while True:
                if mgf_1[ion] == "END IONS":
                    break
                else:
                    peak, abund = mgf_1[ion].split('\t')
                    compound[current_compound]["msms_features"] += [str(round(float(peak), 5))]
                    compound[current_compound]["msms_feat_abund"] += [float(abund)]
                    ion += 1
        
# load mgf for second method
with open(mgf_path_2) as temp:
    mgf_2 = []
    for line in temp:
        mgf_2.append(line.rstrip('\n'))

# load compound feature table for second method
with open(cft_path_2) as temp:
    cft_2 = {}
    for line in temp:
        comp_num, feat_num, TIC, rel_abund, mz, rt = line.rstrip('\n').split(',')
        if comp_num == "compound_number":
            continue
        if comp_num in cft_2:
            cft_2[comp_num]["mz"] += [mz]
            cft_2[comp_num]["TIC"] += [TIC]
        else:
            cft_2[comp_num] = {
                "mz": [mz],
                "TIC": [TIC]
            }

# add to same compound dictionary created to store compounds from first method
for comp in cft_2.keys():
    current_compound = method_2 + "_" + comp
    compound[current_compound] = {
        "ms_features": cft_2[comp]['mz'],
        "ms_feat_abund": cft_2[comp]['TIC'],
        "msms_features":[],
        "msms_feat_abund":[],
        "parent_ion":[]
    }
    for i,line in enumerate(mgf_2):
        if line.split('=')[0] == "SCANS" and line.split('=')[1] == comp:
            compound[current_compound]["parent_ion"] = float(mgf_2[i-1].split('=')[1])
            ion = i + 4
            while True:
                if mgf_2[ion] == "END IONS":
                    break
                else:
                    peak, abund = mgf_2[ion].split('\t')
                    compound[current_compound]["msms_features"] += [str(round(float(peak), 5))]
                    compound[current_compound]["msms_feat_abund"] += [float(abund)]
                    ion += 1

# define cosine score function
def cosine_score(compound_1, compound_2, abundance="abundance", features="features"):
    """Return cosine score of two vectors.
    Arguments:
    compound_1 -- dictionary entry containing named lists with feature and
     abundance values of components of first vector to be compared
    compound_2 -- dictionary entry containing named lists with feature and
     abundance values of components of second vector to be compared
    abundance -- name of list containing abundance values
    features -- name of list containing feature names
    """
    tics_1 = []
    tics_2 = []
    if type(compound_1[abundance][0]) is not float:
        compound_1[abundance] = [float(i) for i in compound_1[abundance]]
    if type(compound_2[abundance][0]) is not float:    
        compound_2[abundance] = [float(i) for i in compound_2[abundance]]
    if type(compound_1[features][0]) is not float:
        compound_1[features] = [float(i) for i in compound_1[features]]
    if type(compound_2[features][0]) is not float:
        compound_2[features] = [float(i) for i in compound_2[features]]
    for feature in compound_1[features]:
        if feature in compound_2[features]:
            tics_1 += [compound_1[abundance][compound_1[features].index(feature)]]
            tics_2 += [compound_2[abundance][compound_2[features].index(feature)]]
    if len(tics_1) > 0:
        cosine_score = sum([tics_1[x] * tics_2[x] for x in range(len(tics_1))]) / \
         ((sum([x ** 2 for x in compound_1[abundance]]) ** (0.5)) * \
         (sum([x ** 2 for x in compound_2[abundance]]) ** (0.5)))
    else:
        cosine_score = 0

    return cosine_score

# calculate cosine score between all pairs of compound from each method. If cosine score for both 
# MS features and MSMS features is above the set cutoff, they are considered to be the same.
compound_matches = {}
for comp1 in compound.keys():
    if method_1 in comp1:
        continue
    for comp2 in compound.keys():
        if method_2 in comp2:
            continue
        # if cosine score between MS peaks > 0.2
        if cosine_score(compound[comp1], compound[comp2], "ms_feat_abund", "ms_features") >= 0.2:
            #if no MSMS data present for EITHER compound being compared
            if len(compound[comp1]['msms_feat_abund']) < 1 or len(compound[comp2]['msms_feat_abund']) < 1:
                #if no MSMS data present - comp 1
                if len(compound[comp1]['msms_feat_abund']) < 1:
                    #match only things with ms cos score > 0.5
                    if cosine_score(compound[comp1], compound[comp2], "ms_feat_abund", "ms_features") >= 0.5:
                        if comp1 in compound_matches:
                            compound_matches[comp1] += [("matched",comp2,cosine_score(compound[comp1], compound[comp2], \
                                "ms_feat_abund", "ms_features"),len(compound[comp1]["ms_features"]),len(compound[comp2]["ms_features"]),\
                            len(set(compound[comp1]["ms_features"]).intersection(compound[comp2]["ms_features"])),comp1 + " no MSMS")]
                        else:
                            compound_matches[comp1] = [("matched",comp2,cosine_score(compound[comp1], compound[comp2], \
                                "ms_feat_abund", "ms_features"),len(compound[comp1]["ms_features"]),len(compound[comp2]["ms_features"]),\
                            len(set(compound[comp1]["ms_features"]).intersection(compound[comp2]["ms_features"])),comp1 + " no MSMS")]
                    #keep comparisons with ms cos score < 0.5, but label as mismatch
                    else:
                        if comp1 in compound_matches:
                            compound_matches[comp1] += [("not matched",comp2,cosine_score(compound[comp1], compound[comp2], \
                            "ms_feat_abund", "ms_features"),len(compound[comp1]["ms_features"]),len(compound[comp2]["ms_features"]),\
                            len(set(compound[comp1]["ms_features"]).intersection(compound[comp2]["ms_features"])),comp1 + " no MSMS")]
                        else:
                            compound_matches[comp1] = [("not matched",comp2,cosine_score(compound[comp1], compound[comp2], \
                            "ms_feat_abund", "ms_features"),len(compound[comp1]["ms_features"]),len(compound[comp2]["ms_features"]),\
                            len(set(compound[comp1]["ms_features"]).intersection(compound[comp2]["ms_features"])),comp1 + " no MSMS")]
                #if no MSMS data present - comp 2 
                if len(compound[comp2]['msms_feat_abund']) < 1:
                     #match only things with ms cos score > 0.5
                    if cosine_score(compound[comp1], compound[comp2], "ms_feat_abund", "ms_features") >= 0.5:
                        if comp1 in compound_matches:
                            compound_matches[comp1] += [("matched",comp2,cosine_score(compound[comp1], compound[comp2], \
                                "ms_feat_abund", "ms_features"),len(compound[comp1]["ms_features"]),len(compound[comp2]["ms_features"]),\
                            len(set(compound[comp1]["ms_features"]).intersection(compound[comp2]["ms_features"])),comp2 + " no MSMS")]
                        else:
                            compound_matches[comp1] = [("matched",comp2,cosine_score(compound[comp1], compound[comp2], \
                                "ms_feat_abund", "ms_features"),len(compound[comp1]["ms_features"]),len(compound[comp2]["ms_features"]),\
                            len(set(compound[comp1]["ms_features"]).intersection(compound[comp2]["ms_features"])),comp2 + " no MSMS")]
                    #keep comparisons with ms cos score < 0.5, but label as mismatch
                    else:
                        if comp1 in compound_matches:
                            compound_matches[comp1] += [("not matched",comp2,cosine_score(compound[comp1], compound[comp2], \
                            "ms_feat_abund", "ms_features"),len(compound[comp1]["ms_features"]),len(compound[comp2]["ms_features"]),\
                            len(set(compound[comp1]["ms_features"]).intersection(compound[comp2]["ms_features"])),comp2 + " no MSMS")]
                        else:
                            compound_matches[comp1] = [("not matched",comp2,cosine_score(compound[comp1], compound[comp2], \
                            "ms_feat_abund", "ms_features"),len(compound[comp1]["ms_features"]),len(compound[comp2]["ms_features"]),\
                            len(set(compound[comp1]["ms_features"]).intersection(compound[comp2]["ms_features"])),comp2+ " no MSMS")]
                #if no MSMS data present for BOTH compounds being compared
                if len(compound[comp1]['msms_feat_abund']) < 1 and len(compound[comp2]['msms_feat_abund']) < 1:
                    #match only things with ms cos score > 0.5
                    if cosine_score(compound[comp1], compound[comp2], "ms_feat_abund", "ms_features") >= 0.5:
                        if comp1 in compound_matches:
                            compound_matches[comp1] += [("matched",comp2,cosine_score(compound[comp1], compound[comp2], \
                            "ms_feat_abund", "ms_features"),len(compound[comp1]["ms_features"]),len(compound[comp2]["ms_features"]),\
                            len(set(compound[comp1]["ms_features"]).intersection(compound[comp2]["ms_features"])),"BOTH no MSMS")]
                        else:
                            compound_matches[comp1] = [("matched",comp2,cosine_score(compound[comp1], compound[comp2], \
                            "ms_feat_abund", "ms_features"),len(compound[comp1]["ms_features"]),len(compound[comp2]["ms_features"]),\
                            len(set(compound[comp1]["ms_features"]).intersection(compound[comp2]["ms_features"])),"BOTH no MSMS")]
                    #keep comparisons with ms cos score < 0.5, but label as mismatch
                    else:
                        if comp1 in compound_matches:
                            compound_matches[comp1] += [("not matched",comp2,cosine_score(compound[comp1], compound[comp2], \
                            "ms_feat_abund", "ms_features"),len(compound[comp1]["ms_features"]),len(compound[comp2]["ms_features"]),\
                            len(set(compound[comp1]["ms_features"]).intersection(compound[comp2]["ms_features"])),"BOTH no MSMS")]
                        else:
                            compound_matches[comp1] = [("not matched",comp2,cosine_score(compound[comp1], compound[comp2], \
                            "ms_feat_abund", "ms_features"),len(compound[comp1]["ms_features"]),len(compound[comp2]["ms_features"]),\
                            len(set(compound[comp1]["ms_features"]).intersection(compound[comp2]["ms_features"])),"BOTH no MSMS")]
            #if MSMS data present for both compounds
            else:
                #match compounds with MS cosine > 0.2 and MSMS cosine greater than 0.25
                if cosine_score(compound[comp1], compound[comp2], "msms_feat_abund", "msms_features") >= 0.25:
                    #label them as matches in dictionary
                    if comp1 in compound_matches:
                        compound_matches[comp1] += [("matched",comp2,cosine_score(compound[comp1], compound[comp2], \
                            "ms_feat_abund", "ms_features"),len(compound[comp1]["ms_features"]),len(compound[comp2]["ms_features"]),\
                            len(set(compound[comp1]["ms_features"]).intersection(compound[comp2]["ms_features"])),\
                            cosine_score(compound[comp1], compound[comp2],"msms_feat_abund", "msms_features"))]
                    else:
                        compound_matches[comp1] = [("matched",comp2,cosine_score(compound[comp1], compound[comp2], \
                            "ms_feat_abund", "ms_features"),len(compound[comp1]["ms_features"]),len(compound[comp2]["ms_features"]),\
                            len(set(compound[comp1]["ms_features"]).intersection(compound[comp2]["ms_features"])),\
                            cosine_score(compound[comp1], compound[comp2],"msms_feat_abund", "msms_features"))]
                # keep comparisons < 0.25 msms too
                else:
                    if comp1 in compound_matches:
                        compound_matches[comp1] += [("not matched",comp2,cosine_score(compound[comp1], compound[comp2], \
                            "ms_feat_abund", "ms_features"),len(compound[comp1]["ms_features"]),len(compound[comp2]["ms_features"]),\
                            len(set(compound[comp1]["ms_features"]).intersection(compound[comp2]["ms_features"])),\
                            cosine_score(compound[comp1], compound[comp2],"msms_feat_abund", "msms_features"))]
                    else:
                        compound_matches[comp1] = [("not matched",comp2,cosine_score(compound[comp1], compound[comp2], \
                            "ms_feat_abund", "ms_features"),len(compound[comp1]["ms_features"]),len(compound[comp2]["ms_features"]),\
                            len(set(compound[comp1]["ms_features"]).intersection(compound[comp2]["ms_features"])),\
                            cosine_score(compound[comp1], compound[comp2],"msms_feat_abund", "msms_features"))]  
        # if ms cos score < 0.2, still keep comparison
        else:
            #if no MSMS data present for EITHER compound being compared
            if len(compound[comp1]['msms_feat_abund']) < 1 or len(compound[comp2]['msms_feat_abund']) < 1:
                #if no MSMS data present - comp 1
                if len(compound[comp1]['msms_feat_abund']) < 1:
                    if comp1 in compound_matches:
                        compound_matches[comp1] += [("not matched",comp2,cosine_score(compound[comp1], compound[comp2], \
                        "ms_feat_abund", "ms_features"),len(compound[comp1]["ms_features"]),len(compound[comp2]["ms_features"]),\
                        len(set(compound[comp1]["ms_features"]).intersection(compound[comp2]["ms_features"])),comp1 + " no MSMS")]
                    else:
                        compound_matches[comp1] = [("not matched",comp2,cosine_score(compound[comp1], compound[comp2], \
                        "ms_feat_abund", "ms_features"),len(compound[comp1]["ms_features"]),len(compound[comp2]["ms_features"]),\
                        len(set(compound[comp1]["ms_features"]).intersection(compound[comp2]["ms_features"])),comp1 + " no MSMS")]
                #if no MSMS data present - comp 2 
                if len(compound[comp2]['msms_feat_abund']) < 1:
                    if comp1 in compound_matches:
                        compound_matches[comp1] += [("not matched",comp2,cosine_score(compound[comp1], compound[comp2], \
                            "ms_feat_abund", "ms_features"),len(compound[comp1]["ms_features"]),len(compound[comp2]["ms_features"]),\
                            len(set(compound[comp1]["ms_features"]).intersection(compound[comp2]["ms_features"])),comp2+ " no MSMS")]
                    else:
                        compound_matches[comp1] = [("not matched",comp2,cosine_score(compound[comp1], compound[comp2], \
                            "ms_feat_abund", "ms_features"),len(compound[comp1]["ms_features"]),len(compound[comp2]["ms_features"]),\
                            len(set(compound[comp1]["ms_features"]).intersection(compound[comp2]["ms_features"])),comp2+ " no MSMS")]
                # if MSMS data missing for BOTH compounds being compared
                if len(compound[comp1]['msms_feat_abund']) < 1 and len(compound[comp2]['msms_feat_abund']) < 1:
                    if comp1 in compound_matches:
                        compound_matches[comp1] += [("not matched",comp2,cosine_score(compound[comp1], compound[comp2], \
                        "ms_feat_abund", "ms_features"),len(compound[comp1]["ms_features"]),len(compound[comp2]["ms_features"]),\
                        len(set(compound[comp1]["ms_features"]).intersection(compound[comp2]["ms_features"])),"BOTH no MSMS")]
                    else:
                        compound_matches[comp1] = [("not matched",comp2,cosine_score(compound[comp1], compound[comp2], \
                        "ms_feat_abund", "ms_features"),len(compound[comp1]["ms_features"]),len(compound[comp2]["ms_features"]),\
                        len(set(compound[comp1]["ms_features"]).intersection(compound[comp2]["ms_features"])),"BOTH no MSMS")]
            #if MSMS data present for both compounds
            else:
                #provide MSMS cosine in mismatch dictionary
                if comp1 in compound_matches:
                    compound_matches[comp1] += [("not matched",comp2,cosine_score(compound[comp1], compound[comp2], \
                        "ms_feat_abund", "ms_features"),len(compound[comp1]["ms_features"]),len(compound[comp2]["ms_features"]),\
                        len(set(compound[comp1]["ms_features"]).intersection(compound[comp2]["ms_features"])),cosine_score(compound[comp1], \
                            compound[comp2], "msms_feat_abund", "msms_features"))]
                else:
                    compound_matches[comp1] = [("not matched",comp2,cosine_score(compound[comp1], compound[comp2], \
                        "ms_feat_abund", "ms_features"),len(compound[comp1]["ms_features"]),len(compound[comp2]["ms_features"]),\
                        len(set(compound[comp1]["ms_features"]).intersection(compound[comp2]["ms_features"])),cosine_score(compound[comp1], \
                            compound[comp2],"msms_feat_abund", "msms_features"))]

# write file with match and mismatch information as CSV output
with open(".//matches_across_columns.csv",'w') as f2:
            f2.write("compound_1, compound_2,match,ms_cosine,comp1_feat,comp2_feat,feat_overlap,msms_cosine\n")
            for comp1 in compound_matches.keys():
                for x in range(len(compound_matches[comp1])):
                    if type(compound_matches[comp1][x][6]) is int or type(compound_matches[comp1][x][6]) is float:
                        f2.write("%s, %s, %s, %f, %d, %d, %d, %f \n" % (comp1,compound_matches[comp1][x][1], compound_matches[comp1][x][0], \
                            compound_matches[comp1][x][2], compound_matches[comp1][x][3],compound_matches[comp1][x][4],\
                                compound_matches[comp1][x][5],compound_matches[comp1][x][6]))
                    else:
                        f2.write("%s, %s, %s, %f, %d, %d, %d, %s \n" % (comp1,compound_matches[comp1][x][1], compound_matches[comp1][x][0], \
                            compound_matches[comp1][x][2], compound_matches[comp1][x][3],compound_matches[comp1][x][4],\
                                compound_matches[comp1][x][5],compound_matches[comp1][x][6]))
            
