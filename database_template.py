from sqlalchemy import Column, Integer, String, ForeignKey
# specify the methods/object you want to add from the beginning

from sqlalchemy import create_engine

engine = create_engine('sqlite:///c:\\sqlite\\preliminary_kinases.db', echo=True)

from sqlalchemy.ext.declarative import declarative_base

Base = declarative_base()

from sqlalchemy.orm import relationship  # to define the related tables


class Human_kinases(Base):
    __tablename__ = 'human_kinases'

    Entry_name = Column(String(30))
    UniProt_ID = Column(String(30), primary_key=True)
    Primary_Protein_Name = Column(String(100))
    Alternate_Protein_Name = Column("Alternate_Protein_Name(s)", String(200))
    Gene_Symbol = Column(String(20))
    Alternative_Gene_Name = Column("Alternative_Gene_Name(s)", String(50))
    Families = Column(String(100))
    AA_sequence = Column(String(1000))
    Molecular_Mass = Column("Molecular_Mass_(Da)", Integer)
    Subcellular_Location = Column(String(500))
    # Fields not defined in the table yet
    Structure = Column(String(200))
    Targets = Column(String(50))
    Isoforms = Column(String(100))
    Inhibitors = Column(String(100))

    kinase_phosphosite = relationship("Kinase_phosphosites", backref="human_kinases")
    Inhibitors = relationship("inhibitors", backref="human_kinases")

    # To define relationships:
    # one-to many, foreign key will be on the table with multiple records
    # and this statement in the main table
    # ClassnameMultiple = relationship("tablenameMultiple",backref="maintablename" )


class Phosphosites(Base):
    __tablename__ = 'phosphosites'

    ID_P = Column(Integer, primary_key=True)
    GENE = Column(String(20), primary_key=True)
    HU_CHR_LOC = Column(String(20))
    MOD_RSD = Column(String(20))
    SITE_GRP_ID = Column(Integer)
    MW_kD = Column(Integer)
    DOMAIN = Column(String(20))
    SITE_7_AA = Column("SITE_+/-7_AA", String(30))
    LT_LIT = Column(Integer)
    MS_LIT = Column(Integer)
    MS_CST = Column(Integer)
    CST_CAT = Column("CST_CAT#", String(20))
    PHOS_ID = Column(String(30), primary_key=True)

    Kinase_phosphosites = relationship("kinase_phosphosites", backref="phosphosites")
    Kinase_phosphosites = relationship("phosphosites_diseases", backref="phosphosites")


class Kinase_phosphosites(Base):
    __tablename__ = 'kinase_phosphosite'

    ID_KP = Column(Integer, primary_key=True)
    GENE = Column(String(20))
    IN_VIVO_RXN = Column(String(10))
    IN_VITRO_RXN = Column(String(10))
    CST_CAT = Column("CST_CAT#", String(20))
    PHOS_ID = Column(String(30), ForeignKey('phosphosites.PHOS_ID'))
    HUMAN_KINASE = Column(String(30), ForeignKey('human_kinases.Primary_Protein_Name'))


class Phosphosites_diseases(Base):
    __tablename__ = 'phosphosites_diseases'

    ID_PD = Column(Integer, primary_key=True)
    DISEASE = Column(String(100))
    ALTERATION = Column(String(40))
    PMIDs = Column(Integer)
    LT_LIT = Column(Integer)
    MS_LIT = Column(Integer)
    MS_CST = Column(Integer)
    CST_CAT = Column("CST_CAT#", String(141))
    NOTES = Column(String(315))
    PHOS_ID = Column(String(32), ForeignKey('phosphosites.PHOS_ID'))  # duplicates


class Inhibitors(Base):
    __tablename__ = 'inhibitors'

# THIS IS PRELIMINARY BASED ON THE Database info file
    ID_I = Column(Integer, primary_key=True)
    Name = Column(String(20))
    Human_Kinase = Column(String(20), ForeignKey('human_kinases.Primary_Protein_Name'))
    Chemical_structure = Column(String(20))
    Phospho_sites_bind_to = Column(String(20), ForeignKey('phosphosites.PHOS_ID'))


Base.metadata.create_all(engine)

