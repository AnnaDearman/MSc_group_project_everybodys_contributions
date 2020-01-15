
# coding: utf-8

# In[16]:


from sqlalchemy import Column, Integer, String, ForeignKey
from sqlalchemy import create_engine
engine = create_engine('sqlite:///c:\\sqlite\\kinases_v1.db', echo=True)
from sqlalchemy.ext.declarative import declarative_base
Base = declarative_base()
from sqlalchemy.orm import relationship


class HumanKinases(Base):
    __tablename__ = 'human_kinases'
    Entry_name = Column(String (30))
    UniProt_ID = Column(String (30),primary_key=True) 
    Primary_Protein_Name = Column (String (100))
    Alternate_Protein_Name= Column("Alternate_Protein_Name(s)",String (200)) 
    Gene_Symbol = Column(String (20))
    Alternative_Gene_Name = Column("Alternative_Gene_Name(s)", String(50)) 
    Families= Column(String (100))
    AA_sequence = Column(String (1000))
    Molecular_Mass = Column("Molecular_Mass_(Da)",Integer)
    Subcellular_Location = Column(String (500))

    kinase_phosphosite = relationship("KinasePhosphosites", backref="human_kinases")
    inhibitors = relationship("Inhibitors", backref="human_kinases")
    # one-to many,
    # tablename_to_link = relationship("Class_of_tablename_to_link",backref="this_tablename" )


class Phosphosites(Base):
    __tablename__ = 'phosphosites'
    GENE = Column(String(20))
    PROTEIN = Column(String(20))
    ACC_ID = Column(String(10))
    HU_CHR_LOC = Column(String(20))
    MOD_RSD = Column(String(20))
    SITE_GRP_ID = Column(String(30)) # check the csv file, it seems integer
    MW_kD = Column(Integer)
    DOMAIN = Column(String(20)) # check the csv file, it would be a string
    SITE_7_AA = Column("SITE_+/-7_AA", String(30))
    LT_LIT = Column(Integer)
    MS_LIT = Column(Integer)
    MS_CST = Column(Integer)
    CST_CAT = Column("CST_CAT#", String(20))
    PHOS_ID = Column(String(50))# to link to phosphosites_diseases.csv and kinases_phosphosites.csv
    PHOS_ID2 = Column(String(30))
    PHOS_ID3 = Column(String(30))
    PHOS_ID4 = Column(String(30))
    ISOFORM = Column(Integer)
    ID_PH = Column(Integer, primary_key=True)
    
    kinase_phosphosites = relationship("KinasePhosphosites", backref="phosphosites")
    kinase_phosphosites = relationship("PhosphositesDiseases", backref="phosphosites")


class KinasePhosphosites(Base):
    __tablename__ = 'kinase_phosphosite'
    GENE = Column(String(20))
    IN_VIVO_RXN = Column(String(2))
    IN_VITRO_RXN = Column(String(2))
    CST_CAT = Column("CST_CAT#", String(150))
    PHOS_ID = Column(String(50), ForeignKey('phosphosites.PHOS_ID'))
    PHOS_ID2 = Column(String(30))
    PHOS_ID3 = Column(String(30))
    PHOS_ID4 = Column(String(30))
    HUMAN_KINASE = Column(String(30), ForeignKey('human_kinases.Entry_name'))
    ID_KS = Column(String(10), primary_key=True)


class PhosphositesDiseases(Base):
    __tablename__ = 'phosphosites_diseases'
    DISEASE = Column(String(100))
    ALTERATION = Column(String(40))
    PMIDs = Column(Integer)
    LT_LIT = Column(Integer)
    MS_LIT = Column(Integer)
    MS_CST = Column(Integer)
    CST_CAT = Column("CST_CAT#", String(150))
    NOTES = Column(String(350))
    PHOS_ID = Column(String(50), ForeignKey('phosphosites.PHOS_ID'))  # duplicates
    ID_PD = Column(String(10), primary_key=True)

    
class Inhibitors(Base):
    __tablename__ = 'inhibitors'
    inhibitor = Column (String (150))
    Ki_nM = Column (Integer) # does entry n/a affect?
    IC50_nM = Column (Integer)# does entry n/a affect?
    Kd_nM = Column (Integer)# does entry n/a affect?
    EC50_nM = Column (Integer)# does entry n/a affect?
    POC = Column (Integer)# does entry n/a affect?
    Source = Column (String(15))
    IMG_URL = Column (String (100))
    ID_IN = Column (String (10), primary_key=True)
    
    inhib_kin = relationship("InhibKin", backref="inhibitors")
    # tablename_to_link = relationship("Class_of_tablename_to_link",backref="this_tablename" )

class InhibKin(Base):
    __tablename__ = 'inhib_kin'
    kinase = Column(String(30), ForeignKey('human_kinases.Entry_name')) #this record MUST to have _HUMAN to match
    inhibitor = Column (String (150), ForeignKey('inhibitors.inhibitor')) 
    ID_KI = Column (String (10), primary_key=True)

Base.metadata.create_all(engine)

