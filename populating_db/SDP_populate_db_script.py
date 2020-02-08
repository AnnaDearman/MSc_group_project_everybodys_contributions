#!/usr/bin/env python
# coding: utf-8

# ### BIO727P - Bioinformatics Software Development Group Project (2019/20)

# #### AIM: To populate the database using SQLAlchemy

# In[ ]:


# Python version: Python 3.7.4

# import the required packages

get_ipython().system('pip install sqlalchemy')
get_ipython().system('pip install pandas')

from sqlalchemy import Column, Integer, String, ForeignKey
from sqlalchemy import create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship
from sqlalchemy.orm import sessionmaker
import pandas as pd


# __Step 1:__ Create an engine that allows us to comminicate with the database, and a session which is the middle ground to talk to our engine

# In[ ]:


# engine = create_engine('sqlite:///c:\\sqlite\\final_database.db', echo=True)

engine = create_engine('sqlite:///c__sqlite_final_database.db', echo=True)

Session = sessionmaker(bind=engine)
session = Session()


# __Step 2:__ The following returns a new base class from which all mapped classes should inherit

# In[ ]:


Base = declarative_base()


# __Step 3:__ Define the structures of all our tables within our database, including the primary and foreign keys

# In[ ]:


class HumanKinases(Base):
    __tablename__ = 'human_kinases'
    UniProt_ID = Column(String(15), primary_key=True)
    PDB_ID = Column(String(5))
    PDB_URL = Column(String(60))
    PDB_title = Column(String(250))
    Entry_name = Column(String(15))
    Primary_Protein_Name = Column(String(100))
    Alternative_Protein_Name = Column(String(350)) 
    Gene_Symbol = Column(String(15))
    Alternative_Gene_Name = Column(String(60)) 
    Families = Column(String(175))
    AA_Seq = Column(String(34400))
    Molecular_Mass = Column(String(10))
    Subcellular_Location = Column(String(350))


# In[ ]:


class Phosphosites(Base):
    __tablename__ = 'phosphosites'
    GENE = Column(String(20))
    PROTEIN = Column(String(20))
    ACC_ID = Column(String(19))
    HU_CHR_LOC = Column(String(26))
    MOD_RSD = Column(String(6))
    SITE_GRP_ID = Column(String(10))
    MW_kD = Column(Integer)
    DOMAIN = Column(String(30)) 
    SITE_7_AA = Column(String(15))
    LT_LIT = Column(Integer)
    MS_LIT = Column(Integer)
    MS_CST = Column(Integer)
    CST_CAT = Column(String(141))
    SOURCE = Column(String(66))
    PMID = Column(String(8))
    PHOS_ID5 = Column(String(24), primary_key=True)
    PHOS_ID = Column(String(31))
    PHOS_ID2 = Column(String(26))
    PHOS_ID3 = Column(String(32))
    PHOS_ID4 = Column(String(25))
    ISOFORM = Column(String(10))
    ID_PH = Column(String(9))


# In[ ]:


class Inhibitors(Base):
    __tablename__ = 'inhibitors'
    Inhibitor = Column(String(150), primary_key=True)
    Ki_nM = Column(Integer) 
    IC50_nM = Column(Integer) 
    Kd_nM = Column(Integer) 
    EC50_nM = Column(Integer) 
    POC = Column(Integer) 
    Source = Column(String(15))
    IMG_URL = Column(String(100))
    ID_IN = Column(String(10))


# In[ ]:


class KinasesPhosphosites(Base):
    __tablename__ = 'kinases_phosphosites'
    GENE = Column(String(13))
    KIN_ACC_ID = Column(String(9))
    SUB_ACC_ID = Column(String(17))
    IN_VIVO_RXN = Column(String(1))
    IN_VITRO_RXN = Column(String(1))
    CST_CAT = Column(String(141))
    SOURCE = Column(String(64))
    PMID = Column (String(8))
    PHOS_ID5 = Column(String(23), ForeignKey('phosphosites.PHOS_ID5'))
    PHOS_ID = Column(String(31))
    PHOS_ID2 = Column(String(23))
    PHOS_ID3 = Column(String(29))
    PHOS_ID4 = Column(String(25))
    KIN_ACC_ID_2 = Column(String(10), ForeignKey('human_kinases.UniProt_ID'))
    ID_KS = Column(String(9), primary_key=True)


# In[ ]:


class PhosphositesDiseases(Base):
    __tablename__ = 'phosphosites_diseases'
    DISEASE = Column(String(92))
    ALTERATION = Column(String(32))
    ACC_ID = Column(String(16))
    PMIDs = Column(String(20))
    LT_LIT = Column(String(20))
    MS_LIT = Column(String(20))
    MS_CST = Column(String(20))
    CST_CAT = Column(String(141))
    NOTES = Column(String(314))
    PHOS_ID = Column(String(22), ForeignKey('phosphosites.PHOS_ID5'))  # duplicates
    ID_DP = Column(String(9), primary_key=True)


# In[ ]:


class InhibKin(Base):
    __tablename__ = 'inhib_kin'
    Kinase = Column(String(8))
    Inhibitor = Column(String(134), ForeignKey('inhibitors.Inhibitor'))
    UniProt_ID = Column(String(6), ForeignKey('human_kinases.UniProt_ID'))
    ID_KI = Column(String(9), primary_key=True)


# __Step 4:__ Define the relationships between the tables within our database

# In[ ]:


kinases_phosphosites = relationship("KinasesPhosphosites", backref="human_kinases")
inhibitors_kinases = relationship("Inhibitors", backref="human_kinases")

phosphosites_kinases = relationship("KinasesPhosphosites", backref="phosphosites")
phosphosites_diseases = relationship("PhosphositesDiseases", backref="phosphosites")

inhib_kin = relationship("InhibKin", backref="inhibitors")


# __Step 5:__ The below creates the tables and relationships between the tables that form the basis of our database

# In[ ]:


Base.metadata.create_all(engine)


# __Step 6:__ Loading in the .csv files to populate the database with

# In[ ]:


kinases=pd.read_csv("human_kinase_dataframe.csv")

phosphosites=pd.read_csv("phosphosites.csv")

inhibs=pd.read_csv("inhibitors.csv")

phospho_kinases=pd.read_csv("kinases_phosphosites.csv")

phospho_diseases=pd.read_csv("phosphosites_diseases.csv")

inhib_kinases=pd.read_csv("inhib_kin.csv")


# __Step 7:__ Iterating over the columns of the .csv files, row by row to insert the information into the tables of the database

# In[ ]:


for k in range(len(kinases)):
    record = HumanKinases(**{
        "UniProt_ID" : kinases.iloc[k, 0],
        "PDB_ID" : kinases.iloc[k, 1],
        "PDB_URL" : kinases.iloc[k, 2],
        "PDB_title" : kinases.iloc[k, 3],
        "Entry_name" : kinases.iloc[k, 4],
        "Primary_Protein_Name" : kinases.iloc[k, 5],
        "Alternative_Protein_Name" : kinases.iloc[k, 6],
        "Gene_Symbol" : kinases.iloc[k, 7],
        "Alternative_Gene_Name" : kinases.iloc[k, 8],
        "Families" : kinases.iloc[k, 9],
        "AA_Seq" : kinases.iloc[k, 10],
        "Molecular_Mass" : kinases.iloc[k, 11],
        "Subcellular_Location" : kinases.iloc[k, 12]
    })
    session.add(record) # add all the records

session.commit() # commit all the records

session.close() # close the connection

print ("Finished")


# In[ ]:


for p in range(len(phosphosites)):
    record = Phosphosites(**{
        "GENE" : phosphosites.iloc[p, 0],
        "PROTEIN" : phosphosites.iloc[p, 1],
        "ACC_ID" : phosphosites.iloc[p, 2],
        "HU_CHR_LOC" : phosphosites.iloc[p, 3],
        "MOD_RSD" : phosphosites.iloc[p, 4],
        "SITE_GRP_ID" : phosphosites.iloc[p, 5],
        "MW_kD" : phosphosites.iloc[p, 6],
        "DOMAIN" : phosphosites.iloc[p, 7],
        "SITE_7_AA" : phosphosites.iloc[p, 8],
        "LT_LIT" : phosphosites.iloc[p, 9],
        "MS_LIT" : phosphosites.iloc[p, 10],
        "MS_CST" : phosphosites.iloc[p, 11],
        "CST_CAT" : phosphosites.iloc[p, 12],
        "SOURCE" : phosphosites.iloc[p, 13],
        "PMID" : str(phosphosites.iloc[p, 14]), # turning into a string to avoid the data being corrupted
        "PHOS_ID5" : phosphosites.iloc[p, 15],
        "PHOS_ID" : phosphosites.iloc[p, 16],
        "PHOS_ID2" : phosphosites.iloc[p, 17],
        "PHOS_ID3" : phosphosites.iloc[p, 18],
        "PHOS_ID4" : phosphosites.iloc[p, 19],
        "ISOFORM" : str(phosphosites.iloc[p, 20]),
        "ID_PH" : phosphosites.iloc[p, 21]
    })
    session.add(record) # add all the records

session.commit() # commit all the records

session.close() # close the connection

print ("Finished")


# In[ ]:


for i in range(len(inhibs)):
    record = Inhibitors(**{
        "Inhibitor" : inhibs.iloc[i, 0],
        "Ki_nM" : inhibs.iloc[i, 1],
        "IC50_nM" : inhibs.iloc[i, 2],
        "Kd_nM" : inhibs.iloc[i, 3],
        "EC50_nM" : inhibs.iloc[i, 4],
        "POC" : inhibs.iloc[i, 5],
        "Source" : inhibs.iloc[i, 6],
        "IMG_URL" : inhibs.iloc[i, 7],
        "ID_IN" : inhibs.iloc[i, 8]
    })
    session.add(record) # add all the records
    
session.commit() # commit all the records

session.close() # close the connection

print ("Finished")


# In[ ]:


for kp in range(len(phospho_kinases)):
    record = KinasesPhosphosites(**{
        "GENE" : phospho_kinases.iloc[kp, 0],
        "KIN_ACC_ID" : phospho_kinases.iloc[kp, 1],
        "SUB_ACC_ID" : phospho_kinases.iloc[kp, 2],
        "IN_VIVO_RXN" : phospho_kinases.iloc[kp, 3],
        "IN_VITRO_RXN" : phospho_kinases.iloc[kp, 4],
        "CST_CAT" : phospho_kinases.iloc[kp, 5],
        "SOURCE" : phospho_kinases.iloc[kp, 6],
        "PMID" : str(phospho_kinases.iloc[kp, 7]), # turning into a string to avoid the data being corrupted
        "PHOS_ID5" : phospho_kinases.iloc[kp, 8],
        "PHOS_ID" : phospho_kinases.iloc[kp, 9],
        "PHOS_ID2" : phospho_kinases.iloc[kp, 10],
        "PHOS_ID3" : phospho_kinases.iloc[kp, 11],
        "PHOS_ID4" : phospho_kinases.iloc[kp, 12],        
        "KIN_ACC_ID_2" : phospho_kinases.iloc[kp, 13],
        "ID_KS" : phospho_kinases.iloc[kp, 14]
    })
    session.add(record) # add all the records
    
session.commit() # commit all the records

session.close() # close the connection

print ("Finished")


# In[ ]:


for dp in range(len(phospho_diseases)):
    record = PhosphositesDiseases(**{
        "DISEASE" : phospho_diseases.iloc[dp, 0],
        "ALTERATION" : phospho_diseases.iloc[dp, 1],
        "ACC_ID" : phospho_diseases.iloc[dp, 2],
        "PMIDs" : str(phospho_diseases.iloc[dp, 3]), # turning into a string to avoid the data being corrupted
        "LT_LIT" : str(phospho_diseases.iloc[dp, 4]), # turning into a string to avoid the data being corrupted
        "MS_LIT" : str(phospho_diseases.iloc[dp, 5]), # turning into a string to avoid the data being corrupted
        "MS_CST" : str(phospho_diseases.iloc[dp, 6]), # turning into a string to avoid the data being corrupted
        "CST_CAT" : phospho_diseases.iloc[dp, 7],
        "NOTES" : phospho_diseases.iloc[dp, 8],
        "PHOS_ID" : phospho_diseases.iloc[dp, 9],
        "ID_DP" : phospho_diseases.iloc[dp, 10]
    })
    session.add(record) # add all the records
    
session.commit() # commit all the records

session.close() # close the connection

print ("Finished")


# In[ ]:


for ik in range(len(inhib_kinases)):
    record = InhibKin(**{
        "Kinase" : inhib_kinases.iloc[ik, 0],
        "Inhibitor" : inhib_kinases.iloc[ik, 1],
        "UniProt_ID" : inhib_kinases.iloc[ik, 2],
        "ID_KI" : inhib_kinases.iloc[ik, 3]
    })
    session.add(record) # add all the records
    
session.commit() # commit all the records

session.close() # close the connection

print ("Finished")

