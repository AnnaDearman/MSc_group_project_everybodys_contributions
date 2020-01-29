
# coding: utf-8

# In[29]:


from sqlalchemy import create_engine, Column, Integer, String, ForeignKey
engine = create_engine('sqlite:///c:\\sqlite\\mock_kinases_v2.1.2.db', echo=True)
from sqlalchemy.ext.declarative import declarative_base
Base = declarative_base()
from sqlalchemy.orm import relationship

class HumanKinases(Base):
    __tablename__ = 'human_kinases'
    Entry_name = Column(String (30),primary_key=True)
    UniProt_ID = Column(String (30)) 
    Primary_Protein_Name = Column (String (100))
    Alternative_Protein_Name= Column("Alternate_Protein_Name(s)",String (200)) 
    Gene_Symbol = Column(String (20))
    Alternative_Gene_Name = Column("Alternative_Gene_Name(s)", String(50)) 
    Families= Column(String (100))
    AA_seq = Column(String (1000))
    Molecular_Mass = Column("Molecular_Mass_(Da)",Integer)
    Subcellular_Location = Column(String (500))
    
    kinases_phosphosites = relationship("KinasesPhosphosites", backref="human_kinases")
    inhib_kin = relationship("InhibKin", backref="human_kinases")
    
class Phosphosites(Base):
    __tablename__ = 'phosphosites'
    GENE = Column(String(20))
    PROTEIN = Column(String(20))
    ACC_ID = Column(String(10))
    HU_CHR_LOC = Column(String(20))
    MOD_RSD = Column(String(20))
    SITE_GRP_ID = Column(String(30))
    MW_kD = Column(Integer)
    DOMAIN = Column(String(20)) 
    SITE_7_AA = Column("SITE_+/-7_AA", String(30))
    LT_LIT = Column(Integer)
    MS_LIT = Column(Integer)
    MS_CST = Column(Integer)
    CST_CAT = Column("CST_CAT#", String(20))
    PHOS_ID = Column(String(50), primary_key=True)
    PHOS_ID2 = Column(String(30))
    PHOS_ID3 = Column(String(30))
    PHOS_ID4 = Column(String(30))
    ISOFORM = Column(Integer)
    ID_PH = Column(String(9))  
    kinases_phosphosites = relationship("KinasesPhosphosites", backref="phosphosites")
    kinases_phosphosites = relationship("PhosphositesDiseases", backref="phosphosites")

class Inhibitors(Base):
    __tablename__ = 'inhibitors'
    Inhibitor = Column (String (150), primary_key=True)
    Ki_nM = Column (Integer) # does entry n/a affect?
    IC50_nM = Column (Integer)# does entry n/a affect?
    Kd_nM = Column (Integer)# does entry n/a affect?
    EC50_nM = Column (Integer)# does entry n/a affect?
    POC = Column (Integer)# does entry n/a affect?
    Source = Column (String(15))
    IMG_URL = Column (String (100))
    ID_IN = Column (String (10))
    inhib_kin = relationship("InhibKin", backref="inhibitors")
    #tablename_to_link = relationship("Class_of_tablename_to_link",backref="this_tablename" )


class KinasesPhosphosites(Base):
    __tablename__ = 'kinases_phosphosites'
    GENE = Column(String(20))
    IN_VIVO_RXN = Column(String(2))
    IN_VITRO_RXN = Column(String(2))
    CST_CAT = Column("CST_CAT#", String(150))
    PHOS_ID = Column(String(50), ForeignKey('phosphosites.PHOS_ID'))
    PHOS_ID2 = Column(String(30))
    PHOS_ID3 = Column(String(30))
    PHOS_ID4 = Column(String(30))
    HUMAN_KINASE = Column(String(30), ForeignKey('human_kinases.Entry_name'))
    ID_KS = Column(String(9), primary_key=True)


class InhibKin(Base):
    __tablename__ = 'inhib_kin'
    Kinase = Column(String(30), ForeignKey('human_kinases.Entry_name')) #this record MUST to have _HUMAN to match
    Inhibitor = Column (String (150), ForeignKey('inhibitors.Inhibitor'))
    ID_KI = Column (String (10), primary_key=True)

Base.metadata.create_all(engine)


#from sqlalchemy.orm import sessionmaker
#Session= sessionmaker(bind=engine)
#session= Session()



#This is a simple query to filter the KInases with names starting with AKT and give the info for the the name and the accession number

#first_search=session.query(HumanKinases).all()
#for x in first_search:
#    print("Kinase Name",x.Entry_name,"; Uniprot ID: ",x.UniProt_ID)

#first_search=session.query(HumanKinases).filter(HumanKinases.Entry_name.like('AKT%'))
#for x in first_search:
#    print("Kinase Name",x.Entry_name,"; Uniprot ID: ",x.UniProt_ID)

