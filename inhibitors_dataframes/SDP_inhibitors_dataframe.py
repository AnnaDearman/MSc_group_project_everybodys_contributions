#!/usr/bin/env python
# coding: utf-8

# ### BIO727P - Bioinformatics Software Development Group Project (2019/20)

# #### AIM: To retrive information about kinases inhibitors from various databases and compile into two tables

# In[ ]:


# Python version: Python 3.7.4

# import the required packages

get_ipython().system('pip install pandas')

import pandas as pd # import pandas
import re # import regular expression
import urllib.request # import url library module
from datetime import datetime


# __Step 1:__ Load in the list of UniProt accession numbers to put into BindingDB

# In[ ]:


uniprot_IDs=pd.read_csv("uniprot_acc_nums.csv")
uniprot=list(uniprot_IDs.UniProt_ID) # load in uniprot accession numbers and convert to list


# __Step 2:__ Using the BindingDB API system, retrieve the relevant information related to all the ligands that bind to the kinases from the UniProt list

# In[ ]:


# open empty lists to become dataframe columns

UniProtKB=[]
monomerid=[]
bindingDB_ID=[]
failed_uniprot=[]
Ki=[] 
IC50=[]
Kd=[]
EC50=[]

for u in range(len(uniprot)):
    try:
        # open up the bindingDB API URL to retrieve the ligands associated with the human kinases
        binding_url="http://www.bindingdb.org/axis2/services/BDBService/getLigandsByUniprots?uniprot="+uniprot[u]+"&&code=0&response=json"
        binding_webpage=urllib.request.urlopen(binding_url)
        binding_file=binding_webpage.read().decode()
        # using regex, extract a single monomer/ligand entry
        entry=re.findall(r"<affinities>(.*?)</affinities>", binding_file)
        
        for e in range(len(entry)):
            
            try:
                # extract the monomer id and append it to a list - also append the uniprot id to another list
                # this stage is essentially mapping the bindingDB ID to the uniprot id it is associated with
                ID_pattern=re.compile(r"<monomerid>([0-9]+)</monomerid>")
                ID_match=ID_pattern.search(entry[e])
                monomerid.append(ID_match.group(1))
                UniProtKB.append(uniprot[u])
            except:
                monomerid.append("N/A")
            
            try:
                # adding "BDBM" infront of the monomer id converts it into the bindingdb accession number
                ID_pattern=re.compile(r"<monomerid>([0-9]+)</monomerid>")
                ID_match=ID_pattern.search(entry[e])
                bindingDB_ID.append("BDBM"+ID_match.group(1))
            except:
                bindingDB_ID.append("N/A")

            try:
                # extracting the Ki (inhibitor constant) of the ligand
                ki_pattern=re.compile(r"<affinity_type>Ki</affinity_type><affinity>(.*?)</affinity>")
                ki_match=ki_pattern.search(entry[e])
                Ki.append(ki_match.group(1).replace("&gt;", "").replace("&lt;", ""))
            except:
                Ki.append("N/A")

            try:
                # extracting the IC50 of the ligand
                # IC50 is the concentration of an inhibitor where the response (or binding) is reduced by half
                ic50_pattern=re.compile(r"<affinity_type>IC50</affinity_type><affinity>(.*?)</affinity>")
                ic50_match=ic50_pattern.search(entry[e])
                IC50.append(ic50_match.group(1).replace("&gt;", "").replace("&lt;", ""))
            except:
                IC50.append("N/A")
    
            try:
                # extracting the Kd (dissociation constant) of the ligand
                kd_pattern=re.compile(r"<affinity_type>Kd</affinity_type><affinity>(.*?)</affinity>")
                kd_match=kd_pattern.search(entry[e])
                Kd.append(kd_match.group(1).replace("&gt;", "").replace("&lt;", ""))
            except:
                Kd.append("N/A")
    
            try:
                # extracting the EC50 of the ligand
                # EC50 is the concentration of a drug that gives half-maximal response
                ec50_pattern=re.compile(r"<affinity_type>EC50</affinity_type><affinity>(.*?)</affinity>")
                ec50_match=ec50_pattern.search(entry[e])
                EC50.append(ec50_match.group(1).replace("&gt;", "").replace("&lt;", ""))
            except:
                EC50.append("N/A")
    
    except:
        # any uniprot ids that can't be found on bindingdb will be appended to this list
        failed_uniprot.append(uniprot[u])
        continue


# __Step 3:__ Convert the monomer IDs (based on BindingDB IDs) into chEMBL IDs using the UniChem API system and append the resulting IDs to a list

# In[ ]:


# open new lists to append information to

chEMBL_ID=[]
failed_list=[]

counts=0

for b in range(len(monomerid)):
    try:
        # this url takes a query monomerid/bindingdb id and maps it to a chEMBL ID
        unichem_url="https://www.ebi.ac.uk/unichem/rest/src_compound_id/"+monomerid[b]+"/31"
        unichem_webpage=urllib.request.urlopen(unichem_url)
        unichem_file=unichem_webpage.read().decode()
        
        # add a counter to check on the progress of the script
        counts=counts+1
        print (counts)
        
        try:
            # using regex extract the chEMBL ID
            chembl_pattern=re.compile(r"{\"src_id\":\"1\",\"src_compound_id\":\"(.*?)\"}")
            chembl_match=chembl_pattern.search(unichem_file)
            chEMBL_ID.append(chembl_match.group(1))
        except:
            # if bindingdb id can't be mapped to the chembl id then append an N/A
            chEMBL_ID.append("N/A")
        
    except:
        # any errors that occur can be appended to the failed list
        failed_list.append(monomerid[b])
        chEMBL_ID.append("N/A")
        pass


# __Step 4:__ Make a new dataframe and append all the information gathered so far

# In[ ]:


# make a new dataframe and append all the lists that have been created

df1=pd.DataFrame({'UniProtID' : UniProtKB})
df1['BindingDB_ID']=bindingDB_ID
df1['Monomer_ID']=monomerid
df1['chEMBL_ID']=chEMBL_ID
df1['Ki_nM']=Ki
df1['IC50_nM']=IC50
df1['Kd_nM']=Kd
df1['EC50_nM']=EC50


# __Step 5:__ Remove all "N/A"s in the chEMBL ID column as they can't be mapped from BindingDB 

# In[ ]:


chEMBL_inhibs=df1[df1.chEMBL_ID != "N/A"] # remove rows where chEMBL_ID == "N/A", these are not useful for the next phase


# __Step 6:__ Create the basis of the inhibitor_kinases relationships table, i.e. a table that maps UniProt ID to its associated BindingDB ID

# In[ ]:


# removing the monomerid, chembl id, and other information related to the properties of the kinase
# this table acts as a mapping device, it relates uniprot ids to bindingdb ids

inhib_kin=chEMBL_inhibs.drop("Monomer_ID", axis=1)
inhib_kin=inhib_kin.drop("chEMBL_ID", axis=1)
inhib_kin=inhib_kin.drop("Ki_nM", axis=1)
inhib_kin=inhib_kin.drop("IC50_nM", axis=1)
inhib_kin=inhib_kin.drop("Kd_nM", axis=1)
inhib_kin=inhib_kin.drop("EC50_nM", axis=1)


# __Step 7:__ From the original chEMBL_inhibs dataframe, remove the UniProt IDs and monomer IDs to create the basis of the inhibitors table, dropping duplicates of inhibitor IDs as we only want this table to act like a dictionary for inhibitors

# In[ ]:


# create a basis of the inhibitors table by removing uniprot id and monomer id
# for this table, we only require one record of every inhibitor, thus we remove duplicates of chEMBL IDs 

bdb_to_chEMBL=chEMBL_inhibs.drop("UniProt_ID", axis=1)

bdb_to_chEMBL=bdb_to_chEMBL.drop("Monomer_ID", axis=1)

bdb_to_chEMBL=bdb_to_chEMBL.drop_duplicates(subset=['chEMBL_ID'], keep='first')


# __Step 8:__ Create a list of all the chEMBL IDs and insert them into a URL to retrieve the page source code. Using this HTML code, using regular expressions we can parse the page and extract the relevant information for each inhibitor. If a chEMBL entry does not have a name associated with it, then it will be skipped

# In[ ]:


ch_id=list(bdb_to_chEMBL.chEMBL_ID) # extract chEMBL IDs as a list


# In[ ]:


# open lists to append information to

chEMBL=[]
molecule_name=[]
molecule_type=[]
molecular_formula=[]
molecular_weight=[]
synonyms=[]

failed_list=[]

counter = 0

for c in range(len(ch_id)):
    # open the compound id page for each element in the list
    ebi_url="https://www.ebi.ac.uk/chembl/compound_report_card/"+ch_id[c]+"/"
    ebi_webpage=urllib.request.urlopen(ebi_url)
    ebi_file=ebi_webpage.read().decode()
    
    # add a counter to check the progress of the loop
    counter = counter + 1
    print (counter)
    
    try:
        # as some chembl id entries don't have an inhibitor name associated with them in the 'Compound:' field, but do have a 
        # name associated with it in the 'Synonyms:' field, to avoid including inhibitors without an official name, we can 
        # skip those entries as they won't be as useful to us
        
        syn_pattern=re.compile(r"Synonyms: (.*?)\"/>")
        syn_match=syn_pattern.search(ebi_file)
        
        name_pattern=re.compile(r"Compound: (.*?)\"\/")
        name_match=name_pattern.search(ebi_file)
        
        # once the entry has an inhibitor name and synonyms associated with it, we can extract the other information associated
        # with the inhibitor
        info_pattern=re.compile(r"Molecule Type: (.*?), Molecular Formula: (.*?), Molecular Weight: (.*?), Synonyms: (.*?)\"/>")
        info_match=info_pattern.search(ebi_file)
            
        molecule_type.append(info_match.group(1))
        molecular_formula.append(info_match.group(2))
        molecular_weight.append(info_match.group(3))
        synonyms.append(info_match.group(4).replace(";", ", "))
        
        # also append the chembl id to the list so we know which accession number the info is associated with
        chEMBL.append(ch_id[c])
            
        # sometimes the molecule name is just the chEMBL ID so to avoid that, we can append the synonym instead
        if name_match.group(1) == ch_id[c]:
            molecule_name.append(info_match.group(4))
        elif name_match.group(1) != ch_id[c]:
            molecule_name.append(name_match.group(1))
    except:
        # if there's no molecule name associated with the entry, we append it to the failed list
        failed_list.append(ch_id[c])
        pass


# __Step 9:__ Create a new dataframe with the information extracted from chEMBL

# In[ ]:


df2=pd.DataFrame({'chEMBL_ID' : chEMBL})
df2['Molecule_name']=molecule_name
df2['Molecule_type']=molecule_type
df2['Molecular_formula']=molecular_formula
df2['Molecular_weight']=molecular_weight
df2['Synonyms']=synonyms


# __Step 10:__ Download the dataframe and clean it up manually. Some entries will have two names in the 'Molecule_name' column, and some other names are displayed in a corrupted way on the HTML source code, and thus need to be re-entered manually. Also in the 'Molecule_name' column some names with a dash in them, that match the names in the 'Synonyms' column can be removed. Lastly, in 'Synonyms', there is a string 'Trade Names:' which can be removed as well. 

# In[ ]:


df2.to_csv("chembl_inhibitors.csv") # download it and clean it up manually


# __Step 10:__ Re-load the cleaned dataframe back into the script

# In[ ]:


inhibs_df=pd.read_csv("chembl_inhibitors.csv")


# __Step 11:__ Merge the dataframe that contains information from BindingDB with the dataframe that contains information from chEMBL based on the chEMBL ID. This forms the basis of the inhibitors table

# In[ ]:


inhibs=bdb_to_chEMBL.merge(inhibs_df, on="chEMBL_ID")


# __Step 12:__ As some of the chEMBL IDs were skipped due to not having a commercial name given to it, we need to remove the entries in the inhibitor-kinase relationships table that are associated with a BindingDB ID that is no longer found in the inhibitors table

# In[ ]:


binding_list=list(inhibs.BindingDB_ID) # list of bindingdb ids from the inhibitors list


# In[ ]:


inhib_kin_binding_list=list(inhib_kin.BindingDB_ID) # list of bindingdb ids from the inhib-kin list


# In[ ]:


# if the bindingdb id from inhib-kin isn't in the bindingdb ids list from the inhibitors list, then append the index of the
# inhib-kin bindingbd id

unwanted_index=[]

for x in range(len(inhib_kin_binding_list)):
    ik_bdb=inhib_kin_binding_list[x]
    if ik_bdb not in binding_list:
        unwanted_index.append(x)


# In[ ]:


inhib_kin2=inhib_kin.drop(unwanted_index) # drop the bindingdb ids that aren't in the inhibs table


# __Step 13:__ Populate the inhibitor-kinase relationship table by adding the chEMBL ID and inhibitor name for each UniProtID to BindingBDID mapping

# In[ ]:


# extract a list of bindingdb ids that map to the uniprot ids

new_binding=list(inhib_kin2.BindingDB_ID) 


# In[ ]:


# extract the inhibitor information we want to put into the inhib-kin table

b=list(inhibs.BindingDB_ID)
c=list(inhibs.chEMBL_ID)
n=list(inhibs.Molecule_name)


# In[ ]:


# open some new lists to append to the inhib-kin dataframe

c_list=[]
n_list=[]

# this loops acts like dictionary. it takes the bindingdb id of the inhib-kin list (that maps to the uniprot id), and then finds
# the index of it in the inhibitor bindingdb id list, using this index we can acquire the chEMBL ID and molecule name from the
# inhibitors table and append it to the c and n lists to form the columns of the inhib-kin dataframe

for y in range(len(new_binding)):
    iden=new_binding[y]
    i=b.index(iden)
    c_list.append(c[i])
    n_list.append(n[i])


# __Step 14:__ Create unique accession numbers for each entry in the inhibitor-kinase relationships table

# In[ ]:


# doing this by adding numbers 0-end of list to "IK" 

inki=[]

for f in range(len(c_list)):
    inki.append("IK"+str(f))


# __Step 15:__ Append the populated columns to the dataframe

# In[ ]:


# append the information to the inhib-kin dataframe

inhib_kin2['chEMBL_ID']=c_list
inhib_kin2['Molecule_name']=n_list
inhib_kin2['IN_KI']=inki


# __Step 16:__ Create unique accession numbers for each inhibitor in the inhibitor table and append it to the inhibitors dataframe

# In[ ]:


# doing this by adding numbers 0-end of list to "IN" 

i_list=list(inhibs.BindingDB_ID)

inhi=[]

for z in range(len(i_list)):
    inhi.append("IN"+str(z))


# In[ ]:


# add the unique accession numbers to the inhibitors dataframe

inhibs['IN_ID']=inhi


# __Step 17:__ Construct image URLs using the chEMBL API system and append the URL list to the inhibitors dataframe

# In[ ]:


# the URL has the same basis, so we can just sub in the chEMBLID into the template to acquire it

chEMBL_URL=[]

for w in range(len(c)):
    chEMBL_URL.append("https://www.ebi.ac.uk/chembl/api/data/image/"+c[w]+".png?bgColor=white")


# In[ ]:


# append the URL list to the dataframe

inhibs['chEMBL_URL']=chEMBL_URL


# __Step 18:__ Manually add some entries that are known to exist but got excluded due to lack of linking up of resources, and finally export the dataframes as .csv files

# In[ ]:


# make lists and add the relevant information for the inhibitors table

BindingDB_ID=["BDBM50398336", "BDBM17052", "BDBM50427326"]
chEMBL_ID=["CHEMBL2178577", "CHEMBL3916849", "CHEMBL2325697"]
Ki=["", "", ""]
IC50=["0.17", "12", "2400"]
Kd=["", "", ""]
EC50=["", "", ""]
Molecule_name=["AZD5363", "BX912", "AZ20"]
Molecule_type=["Small Molecule", "Small Molecule", "Small Molecule"]
Molecular_formula=["C32H44F2N2O5S2", "C20H23BrN8O", "C21H24N4O3S"]
Molecular_weight=["638.84", "471.36", "412.52"]
Synonyms=["AZD-5363, AZD-5672", "BX-912", "AZ-20"]
IN_ID=["IN1429", "IN1430", "IN1431"]
chEMBL_URL=["https://www.ebi.ac.uk/chembl/api/data/image/CHEMBL2178577.png?bgColor=white", 
            "https://www.ebi.ac.uk/chembl/api/data/image/CHEMBL3916849.png?bgColor=white", 
            "https://www.ebi.ac.uk/chembl/api/data/image/CHEMBL2325697.png?bgColor=white"]


# In[ ]:


# create a new dataframe with the information

extra_i=pd.DataFrame({'BindingDB_ID' : BindingDB_ID})
extra_i['chEMBL_ID']=chEMBL_ID
extra_i['Ki_nM']=Ki
extra_i['IC50_nM']=IC50
extra_i['Kd_nM']=Kd
extra_i['EC50_nM']=EC50
extra_i['Molecule_name']=Molecule_name
extra_i['Molecule_type']=Molecule_type
extra_i['Molecular_formula']=Molecular_formula
extra_i['Molecular_weight']=Molecular_weight
extra_i['Synonyms']=Synonyms
extra_i['IN_ID']=IN_ID
extra_i['chEMBL_URL']=chEMBL_URL


# In[ ]:


# merge the two dataframes together and export as a .csv file

frames=[inhibs, extra_i]
result = pd.concat(frames)
result.to_csv("inhibitors_dataframe.csv", index=False)


# In[ ]:


# make lists and add the relevant information for the inhibitors-kinases relaionships table

UniProtKB=["P31749", "P31751", "Q9Y243", "O14965", "Q13535", "P78527", "Q13315", "P42336", "P42345"]
BindingDB_ID=["BDBM50398336", "BDBM50398336", "BDBM50398336", "BDBM17052", "BDBM50427326", "BDBM50427326", "BDBM50427326", 
              "BDBM50427326", "BDBM50427326"]
chEMBL_ID=["CHEMBL2178577", "CHEMBL2178577", "CHEMBL2178577", "CHEMBL3916849", "CHEMBL2325697", "CHEMBL2325697", "CHEMBL2325697", 
           "CHEMBL2325697", "CHEMBL2325697"]
Molecule_name=["AZD5363", "AZD5363", "AZD5363", "BX912", "AZ20", "AZ20", "AZ20", "AZ20", "AZ20"]
IN_KI=["IK54753", "IK54754", "IK54755", "IK54756", "IK54757", "IK54758", "IK54759", "IK54760", "IK54761"]


# In[ ]:


# create a new dataframe with the information

extra_inhibs=pd.DataFrame({'UniProt_ID' : UniProtKB})
extra_inhibs['BindingDB_ID']=BindingDB_ID
extra_inhibs['chEMBL_ID']=chEMBL_ID
extra_inhibs['Molecule_name']=Molecule_name
extra_inhibs['IN_KI']=IN_KI


# In[ ]:


# merge the two dataframes together and export as a .csv file

frames2=[inhib_kin2, extra_inhibs]
result2 = pd.concat(frames2)
result2.to_csv("inhib_kin_dataframe.csv", index=False)

