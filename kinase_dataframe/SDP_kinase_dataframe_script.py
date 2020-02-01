#!/usr/bin/env python
# coding: utf-8

# ### BIO727P - Bioinformatics Software Development Group Project (2019/20)
# 

# #### AIM: To retrive information about human protein kinases from various databases and compile into one table

# In[ ]:


# Python version: Python 3.7.4

# import the required packages

get_ipython().system('pip install pandas')

import pandas as pd # import pandas
import re # import regular expression
import urllib.request # import url library module


# __Step 1:__ Compile a list of human protein kinases

# __Step 1a:__ Use UniProt to gather a list of entry names and accession numbers for the protein kinases listed on the URL specified in the code

# In[ ]:


url="https://www.uniprot.org/docs/pkinfam.txt" # retrieve webpage of human kinases from UniProt
webpage=urllib.request.urlopen(url) # open the URL
myfile=webpage.read().decode() # read the URL contents

hum_kin=re.findall(r"([A-Z0-9]+)_HUMAN", myfile) # create a regular expression to find all human kinases and extract all to a list

uni_id=re.findall(r"_HUMAN\s+\((.*?)\)", myfile)

human_kinases=[] # open an empty list to store UniProt identifiers for all human kinases
uniprot_id=[] # open an empty list to store UniProt accession numbers for all human kinases

for h in range(len(hum_kin)):
    human_kinases.append(hum_kin[h]+"_HUMAN") # append identifiers to new list, re-adding the "_HUMAN" found on all identifiers
    uniprot_id.append(uni_id[h].replace(" ", "")) # append accession numbers to new list


# __Step 1b:__ Append a list of additional kinases found from phosphosite webscraping to existing lists using acquired UniProt entry names and accession numbers 

# In[ ]:


extra_kinases=pd.read_csv("extra_kinases.csv") # open the extra kinase list 
extra_entries_list=list(extra_kinases.Entry_name) # extract extra entry names as a list
extra_accessions_list=list(extra_kinases.uniprot_id) # extract extra entry names as a list
    
for k in range(len(extra_accessions_list)):
    uniprot_id.append(extra_accessions_list[k]) # append extra accession numbers to list
    human_kinases.append(extra_entries_list[k]) # append extra identifiers to list


# __Step 1c:__ Create a dataframe with the UniProt entries and accession numbers - this is the basis for the dataframe

# In[ ]:


df=pd.DataFrame({'Entry_name' : human_kinases}) # create a new table with the header entry name consisting of the human_kinases list
df['UniProt_ID']=uniprot_id # add the uniprot id list as a column to the dataframe


# __Step 1d:__ Create a dataframe with the UniProt accession numbers only - this is so it is possible to convert UniProtKB AC/ID to PDB IDs

# In[ ]:


df1=pd.DataFrame({'UniProt_ID' : uniprot_id})
df1.to_csv("uniprot_acc_nums.csv", index=False)


# __Step 2:__ Extract following kinase information from UniProt:

# __Step 2a:__ Protein names

# In[ ]:


primary_protein_name=[] # open two empty lists for the two different types of names seen
alternative_protein_names=[]

for x in range(len(uniprot_id)):
    protein_url="https://www.uniprot.org/uniprot/?query=id:"+uniprot_id[x]+"&columns=protein%20names&format=tab" # retrieve webpage
    protein_webpage=urllib.request.urlopen(protein_url) # open the URL
    protein_file=protein_webpage.read().decode() # read the URL contents 
        
    protein_list=protein_file.split("\n") # seperate protein name from header
    
    pattern1=re.compile(r"\[(Cleaved into|Includes): (.*?)\]") # regular expression for notes within protein names
    prot_names1=pattern1.sub("", (str(protein_list[1]))) # remove the notes from the string
    
    pattern2=re.compile(r"\(EC(.*?)\)") # regular expression for EC numbers
    prot_names2=pattern2.sub("", prot_names1) # remove EC numbers from the string
    
    # certain names have brackets within their names, and to avoid getting recognised by the regular expression, they need
    # to be replaced by curly brackets
    pattern3=re.compile(r"\((.*?)\)]") 
    prot_names3=pattern3.sub("{acetyl-transferring}]", prot_names2)
    
    pattern4=re.compile(r"\(A\)")
    prot_names4=pattern4.sub("{A}", prot_names3)
    
    pattern5=re.compile(r"\(C\)")
    prot_names5=pattern5.sub("{C}", prot_names4)
    
    pattern6=re.compile(r"\(II\)")
    prot_names6=pattern6.sub("{II}", prot_names5)
    
    pattern7=re.compile(r"\(v-fgr\)")
    prot_names7=pattern7.sub("{v-fgr}", prot_names6)
    
    # exclude all the alternative names that are contained within regular brackets
    pattern8=re.compile(r"\((.*?)\)")
    prot_names8=pattern8.sub("", prot_names7)
    
    prot_names9=prot_names8.replace("{", "(").replace("}", ")").replace("' ", "") # replace punctuation in the string

    primary_protein_name.append(prot_names9) # append primary protein names to list
    
    matches=re.findall(r"\((.*?)\)", prot_names7) # find all alternate protein names contained within brackets 
    alt_prot_names=(str(matches)).replace("[", "").replace("]", "").replace("'", "").replace("{", "(").replace("}", ")").replace('"', "")
        # replace punctuation in the string 
    alternative_protein_names.append(alt_prot_names) # append alternative protein names to list

df['Primary_Protein_Name']=primary_protein_name # add the primary names list as a column to the dataframe
df['Alternative_Protein_Name(s)']=alternative_protein_names # add the alternate names list as a column to the dataframe


# __Step 2b:__ Gene symbols

# In[ ]:


# primary gene symbol

primary_gene_name=[] # open empty list for primary gene names

for p in range(len(uniprot_id)):
    gene_url="https://www.uniprot.org/uniprot/?query=id:"+uniprot_id[p]+"&columns=genes(PREFERRED)&format=tab" # retrieve webpage
    gene_webpage=urllib.request.urlopen(gene_url) # open the URL
    gene_file=gene_webpage.read().decode() # read the URL contents 
    
    gene_list=gene_file.split("\n") # seperate gene name from header
    
    primary_gene_name.append(str(gene_list[1])) # append gene name to the new list

# alternate gene symbol(s)

alternative_gene_names=[] # open new list for alternative gene names

for a in range(len(uniprot_id)):
    gene_alt_url="https://www.uniprot.org/uniprot/?query=id:"+uniprot_id[a]+"&columns=genes(ALTERNATIVE)&format=tab" # retrieve webpage
    gene_alt_webpage=urllib.request.urlopen(gene_alt_url) # open the URL
    gene_alt_file=gene_alt_webpage.read().decode() # read the URL contents 
    
    gene_alt_list=gene_alt_file.split("\n") # seperate gene name from header
    
    gene_alt_names1=(str(gene_alt_list[1])).split(" ") # split by space
    gene_alt_names2=(str(gene_alt_names1)).replace("[", "").replace("]", "").replace("'", "")
    alternative_gene_names.append(gene_alt_names2) # append gene names to the new list

df['Gene_Symbol']=primary_gene_name # add the primary names list as a column to the dataframe
df['Alternative_Gene_Name(s)']=alternative_gene_names # add the alternate names list as a column to the dataframe


# __Step 2c:__ Families

# In[ ]:


family=[] # open new list for family names 

for f in range(len(uniprot_id)):
    fam_url="https://www.uniprot.org/uniprot/?query=id:"+uniprot_id[f]+"&columns=families&format=tab" # retrieve webpage
    fam_webpage=urllib.request.urlopen(fam_url) # open the URL
    fam_file=fam_webpage.read().decode() # read the URL contents 
    
    fam_list=fam_file.split("\n") # seperate family name from header
    
    family.append(fam_list[1]) # append family names to a new list
    
df['Families']=family # add the families names list as a column to the dataframe


# __Step 2d:__ Sequence

# In[ ]:


sequence=[] # open new list for sequence

for s in range(len(uniprot_id)):
    seq_url="https://www.uniprot.org/uniprot/?query=id:"+uniprot_id[s]+"&columns=sequence&format=tab" # retrieve webpage
    seq_webpage=urllib.request.urlopen(seq_url) # open the URL
    seq_file=seq_webpage.read().decode() # read the URL contents 
    
    seq_list=seq_file.split("\n") # seperate sequence name from header
    
    sequence.append(seq_list[1]) # append sequences to a new list

df['AA_Seq']=sequence # add the AA seq list as a column to the dataframe


# __Step 2e:__ Molecular Mass

# In[ ]:


molecular_mass=[] # open new list for molecular mass

for m in range(len(uniprot_id)):
    mass_url="https://www.uniprot.org/uniprot/?query=id:"+uniprot_id[m]+"&columns=mass&format=tab" # retrieve webpage
    mass_webpage=urllib.request.urlopen(mass_url) # open the URL
    mass_file=mass_webpage.read().decode() # read the URL contents 
    
    mass_list=mass_file.split("\n") # seperate sequence name from header
    
    molecular_mass.append(mass_list[1]) # append masses to a new list
    
df['Molecular_Mass_(Da)']=molecular_mass # add the mass list as a column to the dataframe


# __Step 2f:__ Subcellular location

# In[ ]:


sub_cell_locs=pd.read_csv("subcellular_locations.csv") # open the locations list 
sub_cell_locs_list=list(sub_cell_locs.Location) # extract locations as a list

subcellular_location=[] # open new list for subcellular location 

for l in range(len(uniprot_id)):
    cell_url="https://www.uniprot.org/uniprot/?query=id:"+uniprot_id[l]+"&columns=comment(SUBCELLULAR%20LOCATION)&format=tab" # retrieve webpage
    cell_webpage=urllib.request.urlopen(cell_url) # open the URL
    cell_file=cell_webpage.read().decode() # read the URL contents 
    
    cell_list=cell_file.split("\n") # seperate sequence name from header
    cell_locations1=(str(cell_list[1])).replace("SUBCELLULAR LOCATION: ", "") # turn the locations to a string 
    
    pattern5=re.compile(r"Note=(.*?).+") # regular expression to remove notes found in the subcellular locations part of UniProt
    cellular_locations2=pattern5.sub("", cell_locations1) # remove notes found in the subcellular locations part of UniProt 
    
    cellular_locations3=str([loc for loc in sub_cell_locs_list if(loc in cellular_locations2)]) # if a part of the string contains
        # an item found in the subcellular locations list, then add to a list and turn it into a string 
    subcellular_location.append(cellular_locations3.replace("[", "").replace("]", "").replace("'", "")) # replace punctuation in string

df['Subcellular_Location']=subcellular_location # add the subcellular location list as a column to the dataframe


# __Step 3:__ Acquiring structure images from Protein Data Bank (PDB)

# __Step 3a:__ By going onto UniProt.org, and using the 'retrieve/ID mapping' function, it is possible to turn UniProt IDs into PDB IDs. To do this, use the 'uniprot_acc_nums.csv' file created earlier

# __Step 3b:__ Download the two lists: UniProt IDs that were mapped, and IDs that were not able to be mapped. Then load in the mapped list using pandas. As one UniProt ID maps to multiple PDB IDs, we can remove the duplicates and set up lists required to acquire the different information needed

# In[ ]:


pdb=pd.read_csv("uniprot_to_pdb.csv") # load in the list of uniprot IDs converted to PDB IDs

pdb_df = pdb.drop_duplicates(subset=['From'], keep='first') # remove the duplicates

final_pdb = pdb_df.rename(columns={'From': 'UniProt_ID', 'To': 'PDB_ID'}) # rename column headers

pdb_id=list(final_pdb.PDB_ID) # turn the PDB column into a list
uniprot=list(final_pdb.UniProt_ID) # turn the uniprot column into a list


# __Step 3c:__ Construct URLs to acquire the images, and extract the titles of the structures from the HTML source code

# In[ ]:


pdb_url=[]
pdb_title=[]

for z in range(len(pdb_id)):
    pdb_url.append("https://cdn.rcsb.org/pdb/images/"+pdb_id[z]+"_asym_r_500.jpg") # construct the url to acquire images and 
                                                                                   # append to url list
    
    structure_url="https://www.rcsb.org/structure/"+pdb_id[z] # retrieve webpage of kinase structures
    structure_webpage=urllib.request.urlopen(structure_url) # open the URL
    structure_file=structure_webpage.read().decode() # read the URL contents
    
    pattern=re.compile(r"content=\"(.*?): (.*?)\"") # regular expression to acquire the structure title
    match=pattern.search(structure_file) # match the title
    title=match.group(2) # assign the second match group to a new object
    
    title_string=(str(title)).replace("[", "").replace("]", "").replace("'", "").replace('"', "") # tidy up the title entry
    pdb_title.append(title_string) # append to the pdb title list


# __Step 3d:__ Open the unmapped UniProt IDs that were unable to be mapped using pandas, and append them to the uniprot list. Also add empty strings to the ID and URL lists, and in the title list add a message saying no PDB structure available

# In[ ]:


not_in_pdb=pd.read_csv("uniprot_not_mapped.csv") # load in the list of uniprot IDs that couldn't be mapped to pdb ids
not_in_pdb_list=list(not_in_pdb.not_mapped) # turn it into a python list

for j in range(len(not_in_pdb_list)):
    uniprot.append(not_in_pdb_list[j]) # append the uniprot IDs to list
    pdb_id.append("") # append nothing to the pdb id list
    pdb_url.append("") # append nothing to the pdb url list
    pdb_title.append("No Protein Data Bank structure available.") # append message to the pdb title list


# __Step 3e:__ Create a new table appending all the UniProt IDs, PDB IDs, URLs and titles

# In[ ]:


structure=pd.DataFrame({'UniProt_ID' : uniprot}) # make a new pandas table 
structure['PDB_ID']=pdb_id # add the following columns from the lists constructed above
structure['PDB_URL']=pdb_url
structure['PDB_title']=pdb_title


# __Step 4:__ Merge the two dataframes to acquire the final dataframe

# In[ ]:


final_df=structure.merge(df, on="UniProt_ID") # merge the two dataframes based on the uniprot IDs


# __Step 5:__ Export the pandas dataframe to a .csv file

# In[ ]:


final_df.to_csv("human_kinase_dataframe.csv", index=False) # export pandas table to .csv and exclude indexing values as a column

