import imp
import pickle
import pandas as pd
import base64
import os

# set working directory
os.chdir('K://Lab_Map//git_repositories//Evolution_Of_Inga_Chemistry-master')

# load file with feature + compound grouping and compound filling functions
ct = imp.load_source('ct', './/code//build_compound_table_with_classes.py')

# load feature_table_long output from XCMS run in R (step 1)
# This may be either a csv or may be stored in the database, depending on where you put it.
# Just make sure you load in the file you want.


##### C18 ######
ftrd = pd.read_csv(".//data//C18//xcms_feature_table_long_c18.csv",sep=",")
ftrd.head()

# put compounds from feature table into Compound object
c18comps = ct.Compounds(description="c18 compounds negative mode")
c18comps.add_compounds(ftrd)

# perform clustering on features with default settings (mz_bandwidth=0.004, rt_bandwidth=0.3)
# if you want to use different bandwidths, set them using mz_bandwidth= and rt_bandwidth=
c18comps.group_features()
#knowncomps.merge_fill(filled_features=pd.read_csv(".//data//filled_features_fix_2020_02_17.csv"))

#write pickle before the compounds have been assigned...
pickle.dump(c18comps, open(".//data//C18//c18_compound_object_precomps.pkl", "wb"))

# write feature table to csv
c18comps.feature_table_to_csv(".//data//C18//c18_features_matched_across_samples.csv")

# cluster pcids into compounds by calculating cosine score between features with default settings (min_cos_score=0.5)
# if you want to use a different cosine score cutoff, set using min_cos_score = x
c18comps.group_compounds(min_cos_score=0.5,feature_table=pd.read_csv(".//data//C18//c18_features_matched_across_samples.csv",encoding="utf8"))

# write newly constructed compound table to csv
c18comps.compound_table_to_csv(".//data//C18//c18_compound_feature_table.csv")

# write Compounds object to pickle file:
pickle.dump(c18comps, open(".//data//C18//c18_compound_object_grouped_comps.pkl", "wb"))

# for creating a CSV of entire compound pkl object and PCID (original XCMS compounds) cosine comparisons
c18_df = pd.DataFrame()
for x in range(len(c18comps.pcids)):
    temp_df = c18comps[c18comps.pcids[x]].ms_features
    pcid = c18comps.pcids[x]
    comp_num = c18comps[c18comps.pcids[x]].compound_number
    d = {'pcid': [pcid] * len(temp_df), 'mz':temp_df.MZ, 'rt':temp_df.RT,'TIC':temp_df.TIC,\
        'samp':temp_df['sample'],'cluster_mz':temp_df.cluster_mz,'cluster_rt':temp_df.cluster_rt,\
            'feature_number':temp_df.feature_id,'compound_number':[comp_num]*len(temp_df)}
    comp_df = pd.DataFrame(data=d)
    c18_df = c18_df.append(comp_df)

c18_df.to_csv(".//data//C18//c18_comp_feat_pcid_mzrt_rtwindow.csv")



##### polar ######
ftrd = pd.read_csv(".//data//polar//xcms_feature_table_long_polar.csv",sep=",")
ftrd.head()

# put compounds from feature table into Compound object
polcomps = ct.Compounds(description="polar compounds negative mode")
polcomps.add_compounds(ftrd)

# perform clustering on features with default settings (mz_bandwidth=0.004, rt_bandwidth=0.3)
# if you want to use different bandwidths, set them using mz_bandwidth= and rt_bandwidth=
polcomps.group_features()
#knowncomps.merge_fill(filled_features=pd.read_csv(".//data//filled_features_fix_2020_02_17.csv"))

#write pickle before the compounds have been assigned in case you have to reload
pickle.dump(polcomps, open(".//data//polar//polar_compound_object_precomps.pkl", "wb"))

# write feature table to csv
polcomps.feature_table_to_csv(".//data//polar//polar_features_matched_across_samples.csv")

# cluster pcids into compounds by calculating cosine score between features with default settings (min_cos_score=0.5)
# if you want to use a different cosine score cutoff, set using min_cos_score = x
polcomps.group_compounds(min_cos_score=0.5,feature_table=pd.read_csv(".//data//polar//polar_features_matched_across_samples.csv",encoding="utf8"))

# write newly constructed compound table to csv
polcomps.compound_table_to_csv(".//data//polar//polar_compound_feature_table.csv")

# write Compounds object to pickle file:
pickle.dump(polcomps, open(".//data//polar//polar_compound_object_grouped_comps.pkl", "wb"))

# for creating a CSV of entire compound pkl object and PCID (original XCMS compounds) cosine comparisons
pol_df = pd.DataFrame()
for x in range(len(polcomps.pcids)):
    temp_df = polcomps[polcomps.pcids[x]].ms_features
    pcid = polcomps.pcids[x]
    comp_num = polcomps[polcomps.pcids[x]].compound_number
    d = {'pcid': [pcid] * len(temp_df), 'mz':temp_df.MZ, 'rt':temp_df.RT,'TIC':temp_df.TIC,\
        'samp':temp_df['sample'],'cluster_mz':temp_df.cluster_mz,'cluster_rt':temp_df.cluster_rt,\
            'feature_number':temp_df.feature_id,'compound_number':[comp_num]*len(temp_df)}
    comp_df = pd.DataFrame(data=d)
    pol_df = pol_df.append(comp_df)

pol_df.to_csv(".//data//polar//polar_comp_feat_pcid_mzrt_rtwindow.csv")

