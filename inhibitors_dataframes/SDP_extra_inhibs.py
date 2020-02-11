#!/usr/bin/env python
# coding: utf-8

# ### BIO727P - Bioinformatics Software Development Group Project (2019/20)

# #### AIM: To append extra inhibitor information found on http://www.icoa.fr/pkidb/

# In[ ]:


# Python version: Python 3.7.4

# import the required packages

get_ipython().system('pip install chembl_webresource_client')
get_ipython().system('pip install pandas')

import pandas as pd
import csv
from chembl_webresource_client.new_client import new_client
import urllib.request
import re # import regular expression


# __Step 1:__ Load in the chEMBL accession IDs found on http://www.icoa.fr/pkidb/ that weren't present in the original dataframe

# In[ ]:


extra_chembl=pd.read_table("extra_chembls.txt")


# __Step 2:__ Convert the IDs to a list

# In[ ]:


compounds=list(extra_chembl.chEMBL_ID)


# __Step 3:__ Create an empty dictionary and add the chEMBL IDs as the keys

# In[ ]:


compounds2targets = dict()

for x in compounds:
    compounds2targets[x] = set()


# __Step 4:__ Using compound IDs, we can then go and find which proteins the inhibitors act on. These targets found are in terms of as chEMBL IDs

# In[ ]:


# now we have our chEMBL IDs, we can process them in chunks

chunk_size = 50

for i in range(0, len(compounds), chunk_size):
    # we jump from compounds to targets through activities:
    activities = new_client.activity.filter(molecule_chembl_id__in=keys[i:i + chunk_size]).only(
        ['molecule_chembl_id', 'target_chembl_id'])
    # extracting target ChEMBL IDs from activities:
    for act in activities:
        compounds2targets[act['molecule_chembl_id']].add(act['target_chembl_id'])


# __Step 5:__ Convert the chEMBL target IDs into UniProt IDs

# In[ ]:


for key, val in compounds2targets.items():
    # We don't know how many targets are assigned to a given compound so again it's
    # better to process targets in chunks:
    lval = list(val)
    uniprots = set()
    for i in range(0, len(val), chunk_size):
        targets = new_client.target.filter(target_chembl_id__in=lval[i:i + chunk_size]).only(
            ['target_components'])
        uniprots |= set(
            sum([[comp['accession'] for comp in t['target_components']] for t in targets],[]))
    compounds2targets[key] = uniprots


# __Step 6:__ Now we want to map the chEMBL compound IDs to UniProt IDs, i.e. 1-to-1 in the same dataframe, and create a dataframe

# In[ ]:


# extract the chEMBL IDs and uniprot IDs as lists from the dictionary

chemblid=list(compounds2targets.keys())
uniprot=list(compounds2targets.values())


# In[ ]:


# open up new lists to form the dataframe

UniProtKB=[]
chEMBL_ID=[]

# for every uniprot id, append the corresponding chEMBL ID to a list, and the uniprot id, to make a one-on-one mapping table

for y in range(len(chemblid)):
    for j in range(len(uniprot[y])):
        chEMBL_ID.append(chemblid[y])
        UniProtKB.append((list(uniprot[y])[j]))


# In[ ]:


inhib_kin=pd.DataFrame({'UniProt_ID' : UniProtKB})
inhib_kin['chEMBL_ID']=chEMBL_ID


# __Step 7:__ Remove the targets that aren't found in the human kinase list

# In[ ]:


# load in the uniprot accession numbers dataframe and extract the IDs as a list

human_kinases=pd.read_csv("uniprot_acc_nums.csv")
kinase_list=list(human_kinases.UniProt_ID)


# In[ ]:


index=[]

# if the uniprot id from the compound-target mapping isn't in the human kinase uniprot list, then append the index of the 
# compound-mapping uniprot id to a list

for k in range(len(UniProtKB)):
    pos=UniProtKB[k]
    if pos not in kinase_list:
        index.append(k)


# In[ ]:


inhib_kin2=inhib_kin.drop(index) # drop the uniprot ids that aren't in the human kinase list. by doing this, we may remove some 
# of the initial inhibitors because they aren't linked to any of the human kinases in the uniprot list


# __Step 8:__ Convert chEMBL IDs to BindingDB IDs

# In[ ]:


# make a list of individual inhibitors by removing duplicates in the inhib-kin template dataframe

inhib_kin3=inhib_kin2.drop_duplicates(subset=['chEMBL_ID'], keep='first')

inhib_kin3.reset_index() # reset the indexing of the dataframe


# In[ ]:


# make a list of the new chEMBL IDs that are linked to the human kinase inhibitors

new_chembls=list(inhib_kin3.chEMBL_ID)


# In[ ]:


# open up new lists to append new information to

BindingDB=[]
failed_list=[]

for b in range(len(new_chembls)):
    try:
        # this url takes a query chEMBL ID and maps it to the corresponding monomerid/bindingdb id
        unichem_url="https://www.ebi.ac.uk/unichem/rest/src_compound_id/"+new_chembls[b]+"/1"
        unichem_webpage=urllib.request.urlopen(unichem_url)
        unichem_file=unichem_webpage.read().decode()
        try:
            # using regex extract the chEMBL ID
            binding_pattern=re.compile(r"{\"src_id\":\"31\",\"src_compound_id\":\"(.*?)\"}")
            binding_match=binding_pattern.search(unichem_file)
            BindingDB.append(binding_match.group(1))
            
        except:
            # if chEMBL ID can't be mapped to the bindingdb id then append an N/A
            BindingDB.append("N/A")
            
    except:
        # any errors that occur can be appended to the failed list
        failed_list.append(new_chembls[b])
        BindingDB.append("N/A")
        pass


# __Step 9:__ Add the bindingDB information for all entries manually

# In[ ]:


# manually find the bindingDB info for the above monomer IDs

Ki=["2", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "","", "", "", "", "", "", "", "", 
   "", ""]
IC50=["", "", "", "9000", "", "", "20000", "5", "", "", "", "", "27", "4", "", "10000", "", "", "4800", "", "", "", "", "", "", 
     "1.3", "", "5", "3.38", "", "0.8", "", "22.3"]
Kd=["", "", "", "", "", "", "", "", "", "", "", "", "", "2.1", "1200", "", "", "", "", "", "1000", "", "11", "", "", "", "", "", 
   "", "", "", "", ""]
EC50=["", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "",
      "", "", ""]


# __Step 10:__ Iterate through the chEMBL ID list and use HTML source code to extract the relevant information about the inhibitors from chEMBL

# In[ ]:


ch_id=new_chembls

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


# __Step 11:__ Create structure image URLs for the inhibitors from chEMBL

# In[ ]:


chEMBL_URL=[]

for w in range(len(ch_id)):
    chEMBL_URL.append("https://www.ebi.ac.uk/chembl/api/data/image/"+ch_id[w]+".png?bgColor=white")


# __Step 12:__ It is known that all the inhibitors in the extra_chembls.txt file are in chEMBL, but one ID got excluded because it doesn't have any synonyms - this is a minor flaw of the code that could obviously be improved. However, to re-append it to the dataframe, the information was manually gathered up and re-inserted to the lists in the correct position

# In[ ]:


# find the index of the skipped inhibitor

print (ch_id.index("CHEMBL1277072"))

# the position is 12


# In[ ]:


# manually find the relevant information and append it back to the lists at the correct position

chEMBL.insert(12, "CHEMBL1277072")
molecule_name.insert(12, "HENATINIB")
molecule_type.insert(12, "Small molecule")
molecular_formula.insert(12, "C25H29FN4O4")
molecular_weight.insert(12, "468.53")
synonyms.insert(12, "N/A")


# __Step 13:__ Create unique accession numbers for the inhibitors table, starting from the last accession number in the original/main file

# In[ ]:


ID=1431
IN_ID=[]

for x in range(len(chEMBL)):
    ID = ID + 1
    IN_ID.append("IN"+str(ID))


# __Step 14:__ Create a new dataframe of all the information gathered

# In[ ]:


df1=pd.DataFrame({'BindingDB_ID' : BindingDB})
df1['chEMBL_ID']=chEMBL
df1['Ki_nM']=Ki
df1['IC50_nM']=IC50
df1['Kd_nM']=Kd
df1['EC50_nM']=EC50
df1['Molecule_name']=molecule_name
df1['Molecule_type']=molecule_type
df1['Molecular_formula']=molecular_formula
df1['Molecular_weight']=molecular_weight
df1['Synonyms']=synonyms
df1['IN_ID']=IN_ID
df1['chEMBL_URL']=chEMBL_URL


# __Step 15:__ Download the dataframe to check if it needs manual cleaning. The only thing that needed cleaning was in the 'Synonyms' column, whereby it contained a section on 'Trade Names:' for 5 inhibitors, and thus that needed to be removed from the column

# In[ ]:


df1.to_csv("extra_inhibs.csv", index=False)


# __Step 16:__ Load in the extra_inhibitors dataframe and original dataframes created using the 'SDP_inhibitors_dataframe.ipynb' script

# In[ ]:


inhibs_main=pd.read_csv("inhibitors_dataframe.csv")

inhibs_extra=pd.read_csv("extra_inhibs.csv")

inhib_kin_main=pd.read_csv("inhib_kin_dataframe.csv")


# __Step 17:__ Append the original inhibitors dataframe to the extra inhibitors dataframe to create the final inhibitors dataframe, and export is as a final .csv

# In[ ]:


frames=[inhibs_main, inhibs_extra]
result = pd.concat(frames)
result
result.to_csv("final_inhibitors_dataframe.csv", index=False)


# __Step 18:__ Create the extra inhibitors dataframe as outlined below

# In[ ]:


# acquire lists of the uniprot to chEMBL ID mappings

ik_uniprot=list(inhib_kin2.UniProt_ID)
ik_chembl=list(inhib_kin2.chEMBL_ID)

# create the dataframe with the first column being the UniProt ID

inhib_kin_extra=pd.DataFrame({'UniProt_ID' : ik_uniprot})


# In[ ]:


# open new lists to form the basis of the columns in the inhib kin

ik_binding=[]
ik_molecule=[]
ik_id=[]
IK=54761

# this loops acts like dictionary. it takes the chEMBL id of the inhib-kin list (that maps to the uniprot id), and then finds
# the index of it in the inhibitor chEMBL id list, using this index we can acquire the bindingdb ID and molecule name from the
# inhibitors table and append it to the ik_binding and ik_molecule lists to form the columns of the inhib-kin dataframe

# this list also incoporates the accession number creating loop starting from the last accession number in the main dataframe

for e in range(len(ik_chembl)):
    iden=ik_chembl[e]
    i=chEMBL.index(iden)
    ik_binding.append(BindingDB[i])
    ik_molecule.append(molecule_type[i])
    IK=IK + 1
    ik_id.append("IK"+str(IK))


# __Step 19:__ Append the extra inhib_kin dataframe to the original inhib_kin dataframe, and export is as a .csv file

# In[ ]:


inhib_kin_extra['BindingDB_ID']=ik_binding
inhib_kin_extra['chEMBL_ID']=ik_chembl
inhib_kin_extra['Molecule_name']=ik_molecule
inhib_kin_extra['IN_KI']=ik_id


# In[ ]:


frames2=[inhib_kin_main, inhib_kin_extra]
result2 = pd.concat(frames2)
result2.to_csv("final_inhib_kin_dataframe.csv", index=False)

