import requests
import numpy as np
import pandas as pd
import sys
from sklearn.neighbors.kde import KernelDensity
from scipy.signal import argrelextrema
from bisect import bisect_left
import math
from operator import itemgetter
import base64
from datetime import date
import csv

class Compound:

    def __init__(self):
        """
        Constructor - Description of the instance variables
        self.polarity         :: ion polarity (positive or negative)
        self.column           :: type of column used
        self.method           :: UPLC method
        self.primary_mz       :: primary mass to charge value
        self.ms_features      :: Pandas dataframe containing on each row
                                    a (MZ, RT, abundance, sample, mz_cluster_id, 
                                    mz_cluster_mass, rt_cluster_id, rt_cluster_time) entry
        self.msms_features    :: Pandas dataframe containing on each row
                                    a (mass, abundance, msms_feature_id) entry
        self.compound_number  :: number assigned to each unique compound once 
                                    once matched across samples
        self.pcid_idx         :: If compound is component of Compounds object, the index 
                                    of the pcid in Compounds.pcids list
        """
        self.polarity = 0
        self.column = ""
        self.method = ""
        self.primary_mz = -1
        self.ms_features = None
        self.msms_features = None
        self.compound_number = -1
        self.pcid_idx = -1


def cosine_score(compound_1, compound_2, abundance="rel_abund", features="features",min_cos_score=0.5):
    """Return cosine score of two vectors.

    Arguments:
    compound_1 -- Pandas dataframe containing columns with feature and
     abundance values of components of first vector to be compared
    compound_2 -- Pandas dataframe containing columns with feature and
     abundance values of components of first vector to be compared
    abundance -- name of column containing abundance values
    features -- name of column containing feature names
    """

    if len(list(set(compound_1[features]).intersection(compound_2[features]))) == 0:
        return 0
    else:
        cos_temp = compound_1[[abundance]].join(compound_2[[abundance]], how="outer", lsuffix="_1", rsuffix="_2")
        cos_temp = cos_temp.fillna(0)
        
        cosine_score = sum([x[abundance + "_1"] * x[abundance + "_2"] for i, x in cos_temp.iterrows()]) / ((sum([x ** 2 for x in cos_temp[abundance + "_1"]]) ** (0.5)) * (sum([x ** 2 for x in cos_temp[abundance + "_2"]]) ** (0.5)))
        return cosine_score


class Compounds(dict):

    def __init__(self, description=""):
        """
        Constructor - Description of the instance variables
          self.description        :: desription of this group of compounds
          self.total_ms_features  :: total number of ms features imported from frtd table
          self.pcids              :: number of unique pcids imported from ftrd table
          self.all_ms_features    :: pandas dataframe with MZ, RT, pcid index, mz_cluster_id, cluster_mz, 
                                     rt_cluster_id, and cluster_rt for each feature in each sample
          self.compounds          :: dictionary where keys are compound numbers and values are pandas
                                     dataframes containing feature_id, MZ, RT, and TIC for each feature
                                     associated with that compound.
          self.samples            :: dictionary where keys are sample names and values are pandas dataframes
                                     containing sample name, compound number, and TIC for each compound present
                                     in that sample.
        """        
        self.description = description
        self.total_ms_features = 0
        self.pcids = []
        self.all_ms_features = pd.DataFrame(data=None, columns=["MZ", "RT", "pcid_idx", "mz_cluster_id", "cluster_mz",
                            "rt_cluster_id", "cluster_rt"], dtype=float)
        self.compounds = {}
        self.compounds_2 = {}
        self.samples = {}

    def add_compounds(self, ftrd):
        """
        Method which adds compounds from raw feature table to Compounds object
        as elements of dictionary.
        For each unique PC_ID in ftrd table, this function adds:
          a. Compound object to dictionary, with PC_ID as compound name
          b. Pandas dataframe as self[pcid].ms_features with RT, MZ, and TIC of each ms feature
        :param ftrd:REQUIRED  :: pandas dataframe containing the columns RT, PC_ID, TIC, MZ
        """
        if set(["PC_ID", "MZ", "RT", "TIC", "sample"]) <= set(ftrd.columns):
            self.total_ms_features += len(ftrd)
            self.pcids += list(set(ftrd["PC_ID"]))
            for pcid in set(ftrd["PC_ID"]):
                self[pcid] = Compound()
                self[pcid].ms_features = ftrd.ix[ftrd["PC_ID"]==pcid, ["MZ","RT","TIC", "sample"]]
                self[pcid].pcid_idx = self.pcids.index(pcid)
        else:
            print( '   ERROR ftrd must contain columns named')
            print( '   PC_ID, MZ, RT, TIC, sample')
            return
    
    def group_features(self, mz_bandwidth=0.004, rt_bandwidth=0.3):
        """
        Method which clusters ms features across samples first by mass, then by retention time. 
        Uses kernel density estimation to find local maxima and minima. Clusters are defined by 
        everything falling between adjacent minima. The maximum for each cluster is used as the
        value for mass or retention time of that cluster.
        For each ms feature, adds to the ms_features numpy array:
          a. mz cluster index as 5th column
          b. mz value for that cluster as 6th column
          c. rt cluster index as 7th column
          d. rt value for that cluster as 8th column
        :param mz_bandwidth: OPTIONAL, default 0.002
        :param rt_bandwidth: OPTIONAL, default 0.2
        """
        # make array of unique mz_rt values for each PC_ID
        if len(self.pcids) == 0:
            print( "No compounds associated with this Compounds object")
            sys.exit()
        else:
            print( "Collecting ms feature data for " + str(len(self.pcids)) + " PC_IDs")

        for comp in self:
            self.all_ms_features = self.all_ms_features.append(self[comp].ms_features[["MZ", "RT"]].drop_duplicates())
            self.all_ms_features.loc[[str(x)=="nan" for x in self.all_ms_features["pcid_idx"]], "pcid_idx"] = self[comp].pcid_idx

        # kernel density estimation for mz
        print( "Estimating mz kernel density with a bandwidth of " + str(mz_bandwidth))
        kde = KernelDensity(kernel='gaussian', bandwidth=mz_bandwidth).fit(self.all_ms_features["MZ"].values.reshape(-1,1))
        s = np.linspace(85, 2001, num=400000)
        e = kde.score_samples(s.reshape(-1,1))
        mi, ma = argrelextrema(e, np.less)[0], argrelextrema(e, np.greater)[0]
        self.all_ms_features["mz_cluster_id"] = [bisect_left(s[mi], x) for x in self.all_ms_features["MZ"]]
        self.all_ms_features["cluster_mz"] = [s[ma[x]] for x in self.all_ms_features["mz_cluster_id"]]
        print( "mz values grouped into " + str(len(set(self.all_ms_features["mz_cluster_id"]))) + " clusters")

        # use kde to group within every mz cluster
        print( "Estimating rt kernel density with a bandwidth of " + str(rt_bandwidth))
        rt_min, rt_max = float(min(self.all_ms_features["RT"])) - 0.5, float(max(self.all_ms_features["RT"])) + 0.5
        n_breaks = round((rt_max - rt_min)*20)
        s = np.linspace(rt_min, rt_max, num=n_breaks)
        for mzclust in list(set(self.all_ms_features["mz_cluster_id"])):
            if mzclust % 100 == 0:
                print( "    processing mz cluster " + str(mzclust))
            temp1 = self.all_ms_features[self.all_ms_features["mz_cluster_id"] == mzclust].copy()
            kde = KernelDensity(kernel='gaussian', bandwidth=rt_bandwidth).fit(temp1["RT"].values.reshape(-1,1))
            e = kde.score_samples(s.reshape(-1,1))
            mi, ma = argrelextrema(e, np.less)[0], argrelextrema(e, np.greater)[0]
            temp1["rt_cluster_id"] = [bisect_left(s[mi], x) for x in temp1["RT"]]
            temp1["cluster_rt"] = [s[ma[x]] for x in temp1["rt_cluster_id"]]
            self.all_ms_features[self.all_ms_features["mz_cluster_id"] == mzclust] = temp1
        self.all_ms_features["mz_cluster_id"] = [str(int(x)) for x in self.all_ms_features["mz_cluster_id"]]
        self.all_ms_features["rt_cluster_id"] = [str(int(x)) for x in self.all_ms_features["rt_cluster_id"]]
        self.all_ms_features["feature_id"] = [str(x["mz_cluster_id"]) + '_' + str(x["rt_cluster_id"]) for i, x \
        in self.all_ms_features.iterrows()]

        # insert MZ_RT ids into ms feature arrays for each compound 
        for comp_idx in set(self.all_ms_features["pcid_idx"]):
            to_insert = self.all_ms_features[self.all_ms_features["pcid_idx"] == comp_idx]
            self[self.pcids[int(comp_idx)]].ms_features = pd.merge(self[self.pcids[int(comp_idx)]].ms_features, \
            to_insert, how="left", left_on=["MZ", "RT"], right_on=["MZ","RT"])

    def feature_table_to_csv(self, filepath=None):
        if filepath is None:
            print("Please specify filepath for output file")
            return
        with open(filepath, 'w') as file1:
            file1.write("feature_id,cluster_rt,cluster_mz\n")
            for i, row in self.all_ms_features[["feature_id", "cluster_rt", "cluster_mz"]].drop_duplicates().iterrows():
                file1.write("%s,%f,%f\n" % (row["feature_id"], row["cluster_rt"], row["cluster_mz"]))
  
    def group_compounds(self,feature_table,column_type,min_cos_score=0.5):
        """
        Groups pcids based on shared features. 
        Arguments:
        pc_id_table -- dictionary of pcgroups and associated features
        built with build_pc_id function
        min_cos_score -- pcgroups with cosine score above this cutoff
        are grouped as the same compound (default = 0.5, best for our data)
        """
        sorted_pc_ids = sorted(self.pcids, key=lambda k: self[k].ms_features.shape[0],reverse=True)
        feature_table = feature_table.set_index('feature_id')
        compound_rt = {}
        cosine_comp = dict.fromkeys(sorted_pc_ids)

        if len(self.compounds) == 0:
            current_compound_number = 0
        else:
            current_compound_number = max(self.compounds)

        # compare PCIDs to each other   
        for i, pcid in enumerate(sorted_pc_ids):
            if i % 10 == 0:
                print("   Grouping pcid " + str(i) + " of " + str(len(sorted_pc_ids)))
            matched = False
            pcid_temp = pd.pivot_table(self[pcid].ms_features, values=["TIC"], index="feature_id", aggfunc=[sum])
            pcid_temp = pcid_temp["sum"][["TIC"]]
            pcid_temp["feature_id"] = pcid_temp.index.values
            pcid_temp = pd.merge(pcid_temp,feature_table,how="left", left_index = True, right_index = True)
            pcid_temp = pcid_temp.rename(columns={"cluster_mz":"MZ","cluster_rt":"RT"})
            # only compare PCID features if they constitute the top 25% of TIC
            pcid_temp_top_25 = pcid_temp.sort_values(ascending=False, by = 'TIC')
            max_abun = max(pcid_temp_top_25["TIC"])
            pcid_temp_top_25["rel_abund"] = (pcid_temp_top_25["TIC"]/max_abun)*100
            pcid_temp_top_25 =  pcid_temp_top_25[(pcid_temp_top_25['rel_abund'] > 1) & (pcid_temp_top_25["TIC"] > 500)]
            pcid_rt_med = round(np.median(pcid_temp_top_25["RT"]),2)
        
            if current_compound_number == 0:
                current_compound_number += 1
                self.compounds[current_compound_number] = pcid_temp_top_25
                self[pcid].compound_number = current_compound_number
                compound_rt[current_compound_number] = pcid_rt_med
                continue
            
            # compare PCIDs to compounds within a 1 minute window of their median retention time
            comps_to_compare = {key:value for key, value in compound_rt.items() if value >= pcid_rt_med-0.5 and value <= pcid_rt_med+0.5}
            if len(comps_to_compare) == 0:
                current_compound_number += 1
                self.compounds[current_compound_number] = pcid_temp_top_25
                self[pcid].compound_number = current_compound_number
                compound_rt[current_compound_number] = pcid_rt_med
                continue
    
            # make a list of tuples containing the PCID-compound comparisons and cosine scores of each one
            comparisons = []
            for compound in comps_to_compare:
                comparisons.append(tuple([compound,round(cosine_score(pcid_temp_top_25, self.compounds[compound],abundance="rel_abund", features="feature_id"),3)]))

            # add list of tuples to cosine_comp dictionary
            cosine_comp[pcid] = comparisons
            pcid_comps_comp = cosine_comp[pcid]

            # combine compounds with cosine score above the minimum threshold (0.5), append new PCID data to compound grouping
            # if not matched (cosine score < 0.5), add data to distinct compound entry in dictionary
            if len(pcid_comps_comp) > 0:
                if max(map(itemgetter(1),pcid_comps_comp)) >= min_cos_score:
                    compound = [i[0] for i in sorted(pcid_comps_comp,key = lambda x: (x[1],x[1]), reverse=True)][0]
                    temp = pd.concat([self.compounds[compound], pcid_temp_top_25])
                    temp.index.name = None
                    self.compounds[compound] = pd.pivot_table(temp, values=["MZ", "RT", "TIC"], index="feature_id", aggfunc=[sum, np.mean])
                    self.compounds[compound] = self.compounds[compound]["sum"][["TIC"]].join(self.compounds[compound]["mean"][["MZ","RT"]]).sort_values(ascending=False, by = 'TIC')
                    self.compounds[compound]["feature_id"] = self.compounds[compound].index.values
                    max_abun = max(self.compounds[compound]["TIC"])
                    self.compounds[compound]["rel_abund"] = (self.compounds[compound]["TIC"]/max_abun)*100
                    self[pcid].compound_number = compound
                    compound_rt[compound] = round(np.median(self.compounds[compound]["RT"]),2)
                    matched = True
                if not matched:
                    current_compound_number += 1
                    self.compounds[current_compound_number] = pcid_temp_top_25
                    self[pcid].compound_number = current_compound_number
                    compound_rt[current_compound_number] = pcid_rt_med
        
        # create CSV with cosine comparisons
        with open("./" + column_type + "_cosine_comparison_" + str(date.today()) +".csv",'w') as f:
            for k,v in cosine_comp.items():
                if v is not None:
                    for j in v:
                        f.write("{0}: {1}\n".format(k, ", ".join(str(k) for k in j)))                       
        
    def compound_table_to_csv(self, filepath=None):
        """
        Writes compound table created by group_compounds() to comma separated file.
        Requires filepath (either absolute or relative).
        """
        if filepath is None:
            print("Please specify filepath for output file")
            return
        with open(filepath, "w") as file1:
            file1.write("compound_number,feature_number,TIC,rel_abundance,mz,rt\n")
            for compound in self.compounds:
                tot_tic = sum(self.compounds[compound]["TIC"])
                for i, row in self.compounds[compound].iterrows():
                    file1.write("%d,%s,%f,%f,%f,%f\n" % (compound, row["feature_id"], row["TIC"], 
                    row["TIC"] / tot_tic, row["MZ"], row["RT"]))


    def fill_compounds(self, filled_features_filepath=None):
        """
        Determines which compounds are present in a sample based on filled features from R.
        """
        if filled_features_filepath is None:
            print("Please specify filepath of filled features table.")
            return
        filled_features = pd.read_csv(filled_features_filepath,encoding="utf-8-sig",low_memory=False)
        if not set(["sample_name", "feature_id", "actual_rt", "actual_mz", "TIC","norm_TIC"]) <= set(filled_features.columns):
            print('ERROR: filled features file must contain columns named')
            print('       "feature_id", "TIC", "norm_TIC", "actual_mz", "actual_rt", "sample_name"')
            return

        compound_comparisons = {}
        cc_key = 0
        for sample in list(filled_features.sample_name.drop_duplicates()):
            sample = str(sample)
            if sample in self.samples.keys():
                if len(self.samples[sample]) > 0:
                    print("Compounds for " + sample + " have already been filled.")
                    continue
            print("Filling compounds for " + sample)
            self.samples[sample] = pd.DataFrame(data=None, columns=["sample", "compound", "TIC", "norm_TIC","Cosine_Shared_FT","Nshared","Nshared_Major","Percent_Shared"], dtype=float)
            temp_features = filled_features[filled_features.sample_name == sample].drop_duplicates()
            temp_features["feature_available"] = [True] * temp_features.shape[0]
            

            for compound in self.compounds:
                match = False
                temp_comp = self.compounds[compound]
                                
                # get list of features shared between sample and current compound
                shared_features = list(set(temp_features["feature_id"]).intersection(temp_comp["feature_id"]))            
                if len(shared_features) == 0:
                    continue
    
                # get list of major features of compound (within 2 orders of magnitude or as abundant as most abundant peak)
                top_feature = temp_comp["TIC"].idxmax()
                major_features = temp_comp.loc[temp_comp.TIC >= 0.05 * temp_comp.TIC[top_feature], "feature_id"].tolist()
                # get list of major features present in sample
                shared_major_features = [x for x in shared_features if x in major_features]
                # if sample contains at least half of major features:
                shared_per = round(len(shared_major_features)/len(major_features),3)
                
                # sample compound compares the shared features between a sample and a given compound
                sample_compound = temp_features.loc[temp_features["feature_id"].isin(list(set(temp_features["feature_id"]).intersection(temp_comp["feature_id"])))]
                sample_compound = sample_compound.set_index(sample_compound.feature_id)
                top_feature_tic = float(max(sample_compound["TIC"]))
                sample_compound["rel_abund"] = (sample_compound["TIC"]/top_feature_tic)*100
                temp_comp = temp_comp.set_index("feature_id",drop=False)
                cos_score_all = round(cosine_score(sample_compound, temp_comp, abundance="rel_abund", features="feature_id"),3)

                # classify the sample as having a compound if it shares at least 3 major features or 33% of all features AND cosine score > 0.5
                if len(shared_major_features) >= 3 or shared_per > 0.33:
                    if cos_score_all >= 0.5:
                        match = True
                        matched_features = temp_features[temp_features["feature_id"].isin(shared_features)]
                        self.samples[sample] = self.samples[sample].append(pd.Series([sample, int(compound), int(sum(matched_features["TIC"])), int(sum(matched_features["norm_TIC"])),cos_score_all,len(shared_features),len(shared_major_features),shared_per], index=["sample", "compound", "TIC", "norm_TIC","Cosine_Shared_FT","Nshared","Nshared_Major","Percent_Shared"]), ignore_index=True)
                compound_comparisons[cc_key] = (tuple([sample,int(compound),len(shared_features),len(shared_major_features),shared_per,cos_score_all,match]))
                cc_key += 1
            print("    found {} compounds.".format(str(len(self.samples[sample]))))
            
        # make a CSV file with all the cosine comparisons of compounds in filled samples
        with open("./fill_compounds_all_comparisons_" + str(date.today()) +".csv",'w') as f2:
            f2.write("row.names,sample,compound_number,Nshared,Nshared_Major,Percent_Shared,Cosine_shared,match\n")
            for key,value in compound_comparisons.items():
                f2.write("%s, %s\n" % (key, str(value)[1:-1]))
    
    def filled_compounds_to_csv(self, filepath=None):
        """
        Writes filled compound table created by fill_compounds() to comma separated file.
        Requires filepath (either absolute or relative).
        """
        if filepath is None:
            print("Please specify filepath for output file")
            return
        x=1
        with open(filepath, 'w') as file1:
            file1.write("row.names,sample,compound_number,TIC,norm_TIC,Cosine_Shared,Nshared,Nshared_Major,Percent_Shared\n")
            for sample in self.samples:
                for i, row in (self.samples[sample]).iterrows():
                   s = "%s,%s,%d,%f,%f,%f,%f,%f,%f\n" % (x,sample, int(row["compound"]), row["TIC"], row["norm_TIC"],row["Cosine_Shared_FT"],row["Nshared"], row["Nshared_Major"],row["Percent_Shared"])
                   x+=1
                   file1.write(s)


