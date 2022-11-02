import imp
import pickle
import pandas as pd
import os
import csv

# set working directory
os.chdir('K://Lab_Map//git_repositories//Evolution_Of_Inga_Chemistry-master')

# load file with feature + compound grouping and compound filling functions
ct = imp.load_source('ct', './/code//build_compound_table_with_classes.py')


##### C18 ######
# load Compounds object created in 2_group_features_compounds.py
file=open(".//data//C18//c18_compound_object_grouped_comps.pkl","rb")
c18comps = pickle.load(file)

c18comps.fill_compounds(filled_features_filepath=".//data//C18//filled_features_c18_allsamples.csv")
c18comps.filled_compounds_to_csv(".//data//C18//filled_compounds_c18.csv")


##### polar ######
# load Compounds object created in 2_group_features_compounds.py
file=open(".//data//polar//polar_compound_object_grouped_comps.pkl","rb")
polcomps = pickle.load(file)

polcomps.fill_compounds(filled_features_filepath=".//data//polar//filled_features_polar_allsamples.csv")
polcomps.filled_compounds_to_csv(".//data//polar//filled_compounds_polar.csv")
    