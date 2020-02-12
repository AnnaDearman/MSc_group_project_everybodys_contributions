from sqlalchemy import create_engine, Column, Integer, String, ForeignKey
engine = create_engine('sqlite:///c:\\sqlite\\kinase_schema.db', echo=True)
from sqlalchemy.ext.declarative import declarative_base
Base = declarative_base()
from sqlalchemy.orm import relationship
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
    
class Inhibitors(Base):
    __tablename__ = 'inhibitors'
    BindingDB_ID = Column(String(30))
    chEMBL_ID = Column(String(30))
    Ki_nM = Column(String(20))
    IC50_nM = Column(String(20))
    Kd_nM = Column(String(20))
    EC50_nM = Column(String(20))
    Molecule_name = Column(String(100))
    Molecule_type = Column(String(50))
    Molecular_formula = Column(String(50))
    Molecular_weight = Column(String(20))
    Synonyms = Column(String(1000))
    IN_ID = Column(String(20), primary_key=True)
    chEMBL_URL = Column(String(100))

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
    PHOS_ID = Column(String(22), ForeignKey('phosphosites.PHOS_ID5'))
    ID_DP = Column(String(9), primary_key=True)

class InhibKin(Base):
    __tablename__ = 'inhib_kin'
    UniProt_ID = Column(String(20), ForeignKey('human_kinases.UniProt_ID'))
    BindingDB_ID = Column(String(30))
    chEMBL_ID = Column(String(30))
    Molecule_name = Column(String(100), ForeignKey('inhibitors.Molecule_name'))
    IN_KI = Column(String(20), primary_key=True)

kinases_phosphosites = relationship("KinasesPhosphosites", backref="human_kinases")
inhibitors_kinases = relationship("Inhibitors", backref="human_kinases")
phosphosites_kinases = relationship("KinasesPhosphosites", backref="phosphosites")
phosphosites_diseases = relationship("PhosphositesDiseases", backref="phosphosites")
inhib_kin = relationship("InhibKin", backref="inhibitors")

Base.metadata.create_all(engine)

