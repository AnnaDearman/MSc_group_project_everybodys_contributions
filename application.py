
from flask_wtf import FlaskForm
from flask_wtf.file import FileField, FileRequired, FileAllowed


from sqlalchemy import Column, Integer, String, ForeignKey, create_engine, or_
engine = create_engine('sqlite:///c__sqlite_final_database.db', echo = True, connect_args={'check_same_thread': False})
from sqlalchemy.ext.declarative import declarative_base
Base = declarative_base()
from sqlalchemy.orm import sessionmaker, relationship
Session = sessionmaker(bind = engine)
session = Session()

#from database import engine, Session, session, Base, HumanKinases, Phosphosites, Inhibitors, KinasesPhosphosites, PhosphositesDiseases, InhibKin


from scipy.stats import norm
import math

from flask import Flask, flash, render_template, request, redirect, url_for, jsonify, send_from_directory, send_file
from wtforms import StringField, SelectField, SubmitField
from wtforms.validators import DataRequired
from werkzeug.utils import secure_filename
import pandas as pd
import numpy as np
import os
import re

#from datanalysis import Analysis_df1

#import KSEA as datanalysis

#--------------------------------------------------------------------------------------------------------------------------
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
    SEQUENCE = Column(String(3500))
    PMID = Column(String(8))
    PHOS_ID5 = Column(String(24), primary_key=True)
    PHOS_ID = Column(String(31))
    PHOS_ID2 = Column(String(26))
    PHOS_ID3 = Column(String(32))
    PHOS_ID4 = Column(String(25))
    ISOFORM = Column(Integer)
    ID_PH = Column(String(9))

class Inhibitors(Base):
    __tablename__ = 'inhibitors'
    Inhibitor = Column (String(150), primary_key=True)
    Ki_nM = Column (Integer) # does entry n/a affect?
    IC50_nM = Column (Integer) # does entry n/a affect?
    Kd_nM = Column (Integer) # does entry n/a affect?
    EC50_nM = Column (Integer) # does entry n/a affect?
    POC = Column (Integer) # does entry n/a affect?
    Source = Column (String(15))
    IMG_URL = Column (String(100))
    ID_IN = Column (String(10))

class KinasesPhosphosites(Base):
    __tablename__ = 'kinases_phosphosites'
    GENE = Column(String(13))
    KIN_ACC_ID = Column(String(9))
    SUB_ACC_ID = Column(String(17))
    IN_VIVO_RXN = Column(String(1))
    IN_VITRO_RXN = Column(String(1))
    CST_CAT = Column(String(141))
    SOURCE = Column(String(64))
    SEQUENCE = Column(String(8797))
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
    PMIDs = Column(String(8))
    LT_LIT = Column(Integer)
    MS_LIT = Column(Integer)
    MS_CST = Column(Integer)
    CST_CAT = Column(String(141))
    NOTES = Column(String(314))
    PHOS_ID = Column(String(22), ForeignKey('phosphosites.PHOS_ID5'))  # duplicates
    ID_DP = Column(String(9), primary_key=True)

class InhibKin(Base):
    __tablename__ = 'inhib_kin'
    Kinase = Column(String(30)) #ForeignKey('human_kinases.Entry_name')) 
    Inhibitor = Column (String(150), ForeignKey('inhibitors.Inhibitor'))
    ID_KI = Column (String(10), primary_key=True)
    UniProt_ID = Column(String(6), ForeignKey('human_kinases.UniProt_ID'))

kinases_phosphosites = relationship("KinasesPhosphosites", backref="human_kinases")
inhibitors_kinases = relationship("Inhibitors", backref="human_kinases")

phosphosites_kinases = relationship("KinasesPhosphosites", backref="phosphosites")
phosphosites_diseases = relationship("PhosphositesDiseases", backref="phosphosites")

inhib_kin = relationship("InhibKin", backref="inhibitors")

Base.metadata.create_all(engine)
#--------------------------------------------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------------------------------------------- 
#Flask application object
application = Flask(__name__)
application.config['SECRET_KEY'] = 'jacky'

UPLOAD_FOLDER =  os.path.dirname(os.path.abspath(__file__)) + '/uploads/'
application.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

ALLOWED_EXTENSIONS = {'tsv'}
def allowed_file(filename):
   return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

DOWNLOAD_FOLDER = os.path.dirname(os.path.abspath(__file__)) + '/downloads/'
application.config['DOWNLOAD_FOLDER'] = DOWNLOAD_FOLDER


#----------------------------------------------------------------------------------------------------------

def process_file(path, filename):
    #def Analysis_df1(path, filename):
   ### import sample data
    data_raw = pd.read_table(path, sep = "\t")
    

    data_raw = data_raw.loc[:, ~data_raw.columns.str.contains('^Unnamed')]

    data_raw.columns = ['Substrate', 'Control_mean', 'Treat_mean','Fold_change', 'p_value', 'ctrlCV', 'treatCV']
    ### drop any rows have a NAN
    global df1
    df1 = data_raw.copy()
    df1 = df1.dropna(how="any")

    ### split the first column
    df1[["Gene", "sites"]] = df1.Substrate.str.split("(", expand = True,)
    df1[["Phosphosite", "junk"]] = df1.sites.str.split(")", expand = True,)
    df1 = df1[['Substrate','Fold_change','p_value', 'Gene', 'Phosphosite']]

    ### drop rows contain None in Phosphosite column
    df1 = df1[df1.Phosphosite != "None"]

    ### transform Fold_change to Log2 and then drop original Fold_change column
    df1["Log_F"] = np.log2(df1['Fold_change'])
    df1 = df1.drop(['Fold_change'], axis = 1)

    ### transform p_value to -log10 p_value
    df1["-log10p"] = -np.log10(df1["p_value"])

    #########
    df_phos = pd.read_sql(session.query(KinasesPhosphosites).statement, session.query(KinasesPhosphosites).session.bind)
    df_merge = pd.merge(df1, df_phos, how = 'left', left_on = 'Substrate', right_on = 'PHOS_ID4')
    df_merge = df_merge[['Substrate','p_value','Gene','Phosphosite','Log_F','-log10p','GENE']]
    df_merge['GENE'].fillna('Not_Found', inplace=True)
    df_merge.columns = ['Substrate','p_value','Gene','Phosphosite','Log_F','-log10p','Kinase']
    df_merge = df_merge.drop_duplicates(subset = 'Substrate', keep = False)
    df1= df_merge
 
    df1.to_csv(os.path.join(application.config['DOWNLOAD_FOLDER'], r'table1_analysis.csv'))


    ### drop rows contain Not_Found in Kinase column
    global df_kinase
    df_kinase = df1[df1.Kinase != "Not_Found"]

    df_kinase.to_csv(os.path.join(application.config['DOWNLOAD_FOLDER'], r'table2_analysis.csv'))
    
##########################################################################
### Function to get the table with substrates without kinases
    ### select rows contain Not_Found in Kinase column
    global df_no_kinase
    df_no_kinase = df1[df1.Kinase == "Not_Found"]


    df_no_kinase.to_csv(os.path.join(application.config['DOWNLOAD_FOLDER'], r'table3_analysis.csv'))

    
    ### calculate the subgroup mean
    global df_submean
    df_submean = df_kinase.groupby(['Kinase']).mean()
    df_submean = df_submean.drop(['p_value'], axis = 1)
    df_submean = df_submean.drop(['-log10p'], axis = 1)

    ### calculate Z_score and convert Z_score to P_value
    def Z_score(mS):
        Z_score = (mS - df1['Log_F'].mean())*(len(df1['Log_F'])**(1/2))/np.std(df1['Log_F'])
        P_value = norm.sf(abs(Z_score))
    #Z_score = (mS - df1['Log_F'].mean())*(len(df1['Log_F'])**(1/2))/np.std(df1['Log_F'])
    #P_value = norm.sf(abs(Z_score))

    ### calculate Z_score and convert Z_score to P_value
    df_submean['p_value'] = df_submean['Log_F'].apply(Z_score)


    ### calculate the standard deviation of fold change in each subgroup
    df_std = df_kinase.groupby(['Kinase']).aggregate(np.std)
    df_submean['Std'] = df_std['Log_F']



    ### Calculate the relative kinase activity
    df_submean['Kinase_relative_activity_score'] = df_submean['Log_F']/df1['Log_F'].mean()

    ### add a kinase column to df_submean
    df_submean['Kinase'] = df_submean.index


    # Resets the index
    df_submean.reset_index(drop=True, inplace=True)


    ### set a count colum to show the number of substrates of each group
    df_submean['count'] = list(df_kinase['Kinase'].value_counts().sort_index())


    df_submean.to_csv(os.path.join(application.config['DOWNLOAD_FOLDER'], r'table4_analysis.csv'))


    Bon_cor = -math.log10(0.05/len(df1['p_value']))

    ### set a color colum to show significance of the p_value
    df1.loc[df1['-log10p'] < Bon_cor, 'color'] = 'Non-sig'
    df1.loc[df1['-log10p'] >= Bon_cor, 'color'] = 'Sig'

    ### volcano plot ### also need to install nbformat before running plotly
    import plotly.express as px
    global fig
    fig = px.scatter(df1, x = "Log_F", y="-log10p", hover_data=['Substrate'], color = 'color')
    fig.write_html("templates/Volcano.html") #Convert the figure to HTML, so it can be accessed on the web application

    ### bar plot of kinase relative activity scores
    global df_subgroub
    df_subgroub = df_submean
    df_subgroub.loc[df_subgroub['Kinase_relative_activity_score'] < 0, 'color'] = 'Treat_down'
    df_subgroub.loc[df_subgroub['Kinase_relative_activity_score'] >= 0, 'color'] = 'Treat_up'

    ### set a text colum to show significance of the p_value

    df_subgroub.loc[df_subgroub['p_value'] <= 0.05, 'Significance'] = '*'
    df_subgroub.loc[df_subgroub['p_value'] > 0.05, 'Significance'] = ''
    df_subgroub.loc[df_subgroub['p_value'] == 0, 'Significance'] = ' '


    ### show relative activity by a bar plot
    global fig2
    fig2 = px.bar(df_subgroub,
                 x = 'Kinase',
                 y = 'Kinase_relative_activity_score',
                 color = 'color',
                 text = 'Significance',
                 error_y = 'Std',
                 hover_data = ['p_value', 'count'])
    
    fig2.update_traces(textposition='outside')
    fig2.write_html("templates/Kinase_RKA_barplot.html") #Convert the figure to HTML, so it can be accessed on the web application


#--------------------------------------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------------------------------------





#-------------------------------------------------------------------------------------------------------------------------- 

#Search & Upload File forms
#--------------------------------------------------------------------------------------------------------------------------  
class SearchForm(FlaskForm):
    protein_name = StringField("", validators=[DataRequired()])
    submit = SubmitField('Submit')
    datanalysis = SubmitField('Download phosphoproteomics analysis')
#--------------------------------------------------------------------------------------------------------------------------  

#Main application script

#Homepage(s)
#--------------------------------------------------------------------------------------------------------------------------  
#init_db()

@application.route('/', methods=['GET', 'POST'])
#@app.route('/datanalysis', methods = ['GET', 'POST'])
@application.route('/home') 
def index():
    form = SearchForm()
    #form1 = DocumentUploadForm()
    if request.method == 'POST':
        if 'file' not in request.files:
           print('No file attached in request')
           return redirect(request.url)
        file = request.files['file']
        if file.filename == '':
           print('No file selected')
           return redirect(request.url)
        if file and allowed_file(file.filename):
           filename = secure_filename(file.filename)
           file.save(os.path.join(application.config['UPLOAD_FOLDER'], filename))
           process_file(os.path.join(application.config['UPLOAD_FOLDER'], filename), filename)
           return redirect(url_for('uploaded_file', filename=filename))
        else:
            return "The file uploaded is not in tsv format, please upload a file in this format"
    return render_template('homepage.html', form=form)

@application.route('/kin/', methods=['GET', 'POST']) 
def kin():
    form = SearchForm()
    protein_name = None
    if form.validate_on_submit():
        protein_name = form.protein_name.data
        return redirect(url_for('kinasesearch', kinase_search=protein_name))
    return render_template('homepagekin.html', form=form)

@application.route('/inh/', methods=['GET', 'POST']) 
def inh():
    form = SearchForm()
    protein_name = None
    if form.validate_on_submit():
        protein_name = form.protein_name.data
        return redirect(url_for('inhibitorsearch', inhibitor_search=protein_name))
    return render_template('homepageinh.html',form=form)

@application.route('/phos/', methods=['GET', 'POST']) #If Phosphosite selected on select search menu, redirects to homepage with an additional select box (genomic location & neighbouring sequence)
def phos():
    form = SearchForm()
    #protein_name = None
    return render_template('homepagephos.html', form=form)

@application.route('/nseq/', methods=['GET', 'POST']) #If Phosphosite selected on select search menu, redirects to homepage with an additional select box (genomic location & neighbouring sequence)
def phosnseq():
    form = SearchForm()
    protein_name = None
    if form.validate_on_submit():
        protein_name = form.protein_name.data
        return redirect(url_for('phosphositesearchnseq', phosphosite_search=protein_name))
    return render_template('homepagephos.html', form=form)

@application.route('/prot/', methods=['GET', 'POST']) #If Phosphosite selected on select search menu, redirects to homepage with an additional select box (genomic location & neighbouring sequence)
def phosprot():
    form = SearchForm()
    protein_name = None
    if form.validate_on_submit():
        protein_name = form.protein_name.data
        return redirect(url_for('phosphositesearchprot', phosphosite_search=protein_name))
    return render_template('homepagephos.html', form=form)

@application.route('/phoskin/', methods=['GET', 'POST']) #If Phosphosite selected on select search menu, redirects to homepage with an additional select box (genomic location & neighbouring sequence)
def phospkin():
    form = SearchForm()
    protein_name = None
    if form.validate_on_submit():
        protein_name = form.protein_name.data
        return redirect(url_for('phosphositesearchkin', phosphosite_search=protein_name))
    return render_template('homepagephos.html', form=form)



#The next routes are to redirect to the intended url
@application.route('/kin/inh/')
@application.route('/phos/inh/')
@application.route('/prot/inh/')
@application.route('/gen/inh/')
@application.route('/nseq/inh/')
@application.route('/phos/redirect/inh/')
@application.route('/kin/redirect/inh/')
@application.route('/uploads/inh/')
@application.route('/phoskin/inh/')
@application.route('/nseq/inh/')
@application.route('/prot/inh/')
def redihn():
    return redirect('/inh/')

@application.route('/inh/phos/')
@application.route('/kin/phos/')
@application.route('/inh/redirect/phos/')
@application.route('/kin/redirect/phos/')
@application.route('/uploads/phos/')
def redphos():
    return redirect('/phos/')

@application.route('/inh/kin/')
@application.route('/phos/kin/')
@application.route('/phoskin/kin/')
@application.route('/prot/kin/')
@application.route('/gen/kin/')
@application.route('/nseq/kin/')
@application.route('/inh/redirect/kin/')
@application.route('/phos/redirect/kin/')
@application.route('/uploads/kin/')
def redkin():
    return redirect('/kin/')

@application.route('/prot/nseq')
@application.route('/gen/nseq')
def rednseq():
    return redirect('/nseq/')

@application.route('/nseq/prot')
@application.route('/gen/prot')
def redprot():
    return redirect('/prot/')

@application.route('/nseq/gen')
@application.route('/prot/gen')
def redgen():
    return redirect('/gen/')


#------------------------------------------------------------------------------------------------------------------------------------------   

#Data analysis page
#-----------------------------------------------------------------------------------------------------------------------------------------------------
#Display the results of the analysis
@application.route('/uploads/<filename>')
def uploaded_file(filename):
    form = SearchForm()
    return render_template('datanalysis.html', filename=filename, tables1=[df1.to_html(classes='data')], titles1=df1.columns.values, tables2=[df_kinase.to_html(classes='data')], titles2=df_kinase.columns.values, tables3=[df_no_kinase.to_html(classes='data')], titles3=df_no_kinase.columns.values, tables4=[df_submean.to_html(classes='data')], titles4=df_submean.columns.values, fig=fig2, form=form)

#Download the phosphosites-kinases table with kinases match, and phosphosites with no kinase match
@application.route('/download/table1')
def download_table1():
    return send_from_directory(application.config['DOWNLOAD_FOLDER'], 'table1_analysis.csv' , as_attachment=True)

#Download the phosphosites-kinases table with only the phosphosites for which a kinase match was found
@application.route('/download/table2')
def download_table2():
    return send_from_directory(application.config['DOWNLOAD_FOLDER'], 'table2_analysis.csv' , as_attachment=True)

#Download the phosphosites-kinases table with the phosphosites for which no kinase match was found
@application.route('/download/table3')
def download_table3():
    return send_from_directory(application.config['DOWNLOAD_FOLDER'], 'table3_analysis.csv' , as_attachment=True)

#Download the relative kinase activity table
@application.route('/download/table4')
def download_table4():
    return send_from_directory(application.config['DOWNLOAD_FOLDER'], 'table4_analysis.csv' , as_attachment=True)

#Interactive version of the volcano plot
@application.route('/volcanoplot')
def volcanoplot():
    return render_template('Volcano.html')

#Interactive version of the relative kinase activity bar plot
@application.route('/RKAplot')
def RKAplot():
    return render_template('Kinase_RKA_barplot.html')
#---------------------------------------------------------------------------------------------------------------------------------------------------- 

#Hits pages
#------------------------------------------------------------------------------------------------------------------------------------------
#Inhibitor search results page
@application.route('/inh/<inhibitor_search>', methods=['GET', 'POST'])
def inhibitorsearch(inhibitor_search):
    form = SearchForm()
    inhibitor_search = inhibitor_search.upper()
    resultss = session.query(Inhibitors).filter(Inhibitors.Inhibitor.startswith(inhibitor_search)).all()
    protein_name = None
    if form.validate_on_submit():
        protein_name = form.protein_name.data
        return redirect(url_for('inhibitorsearch', inhibitor_search=protein_name))
    return render_template('inhibitorhits.html', resultss=resultss, form=form)

#Kinase search results page
@application.route('/kin/<kinase_search>', methods=['GET', 'POST'])
def kinasesearch(kinase_search):
    form = SearchForm()
    kinase_search = kinase_search.upper()
    resultss = session.query(HumanKinases).filter(or_(HumanKinases.Entry_name.startswith(kinase_search), HumanKinases.Gene_Symbol.startswith(kinase_search)) ).all()
    protein_name = None
    if form.validate_on_submit():
        protein_name = form.protein_name.data
        return redirect(url_for('kinasesearch', kinase_search=protein_name))
    return render_template('kinasehits.html', resultss=resultss, form=form)
 
#Phosphosite search by neighbouring sequence results page   
@application.route('/nseq/<phosphosite_search>', methods=['GET', 'POST'])
def phosphositesearchnseq(phosphosite_search):
    form = SearchForm()
    resultss = session.query(Phosphosites).filter(Phosphosites.SITE_7_AA.startswith(phosphosite_search)).all()
    protein_name = None
    if form.validate_on_submit():
        protein_name = form.protein_name.data
        return redirect(url_for('phosphositesearchnseq', phosphosite_search=protein_name))

    return render_template('phoshits.html', resultss=resultss, form=form)


#Phosphosite search by protein results page
@application.route('/prot/<phosphosite_search>', methods=['GET', 'POST'])
def phosphositesearchprot(phosphosite_search):
    form = SearchForm()
    resultss = session.query(Phosphosites).filter(Phosphosites.PROTEIN.startswith(phosphosite_search)).all()
    protein_name = None
    if form.validate_on_submit():
        protein_name = form.protein_name.data
        return redirect(url_for('phosphositesearchprot', phosphosite_search=protein_name))
    return render_template('phoshits.html', resultss=resultss, form=form)

#Phosphosite search by 'which kinase they're phosphorylated by' results page
@application.route('/phoskin/<phosphosite_search>', methods=['GET', 'POST'])
def phosphositesearchkin(phosphosite_search):
    form = SearchForm()
    resultss = session.query(Phosphosites).join(KinasesPhosphosites).join(HumanKinases).filter(or_(HumanKinases.Entry_name.startswith(phosphosite_search), HumanKinases.Gene_Symbol.startswith(phosphosite_search)) ).all() #This line has been changed
    protein_name = None
    if form.validate_on_submit():
        protein_name = form.protein_name.data
        return redirect(url_for('phosphositesearchkin', phosphosite_search=protein_name))
    return render_template('phoshits.html', resultss=resultss, form=form)

#Phosphosite browse by genomic location page
@application.route('/gen/', methods=['GET', 'POST'])
def phosphositesearchgenhome():
    form = SearchForm()
    resultss = session.query(Phosphosites).filter(Phosphosites.HU_CHR_LOC.startswith('t')).all()
    if form.validate_on_submit():
        protein_name = form.protein_name.data
        return redirect(url_for('phosphositesearchgen', gen_location=protein_name))
    return render_template('zphoshitsgen.html',form=form, resultss=resultss)

#Phosphosite browse by genomic location page, once you have selected a genomic location
@application.route('/gen/<gen_location>', methods=['GET', 'POST'])
def phosphositesearchgen(gen_location):
    form = SearchForm()
    resultss = session.query(Phosphosites).filter(Phosphosites.HU_CHR_LOC.startswith(gen_location)).all()
    protein_name = None
    if form.validate_on_submit():
        protein_name = form.protein_name.data
        gen_location=protein_name
        #return redirect(url_for('phosphositesearchgen', gen_location=protein_name))
        if gen_location== '1':
            p= session.query(Phosphosites).filter(Phosphosites.HU_CHR_LOC.startswith('1p')).all()
            q= session.query(Phosphosites).filter(Phosphosites.HU_CHR_LOC.startswith('1q')).all()
            return render_template('zphoshitsgenalt.html',form=form, p=p, q=q)
        if gen_location== '2':
            p2= session.query(Phosphosites).filter(Phosphosites.HU_CHR_LOC.startswith('2p')).all()
            q2= session.query(Phosphosites).filter(Phosphosites.HU_CHR_LOC.startswith('2q')).all()
            return render_template('zphoshitsgenalt.html',form=form, q=q2, p=p2)
        else:
            return render_template('zphoshitsgen.html',form=form, resultss=resultss)
    else:
        if gen_location== '1':
            p= session.query(Phosphosites).filter(Phosphosites.HU_CHR_LOC.startswith('1p')).all()
            q= session.query(Phosphosites).filter(Phosphosites.HU_CHR_LOC.startswith('1q')).all()
            return render_template('zphoshitsgenalt.html',form=form, p=p, q=q)
        if gen_location== '2':
            p2= session.query(Phosphosites).filter(Phosphosites.HU_CHR_LOC.startswith('2p')).all()
            q2= session.query(Phosphosites).filter(Phosphosites.HU_CHR_LOC.startswith('2q')).all()
            return render_template('zphoshitsgenalt.html',form=form, q=q2, p=p2)
        else:
            return render_template('zphoshitsgen.html',form=form, resultss=resultss)
    return render_template('zphoshitsgen.html',form=form, resultss=resultss)
#--------------------------------------------------------------------------------------------------------------------------  


#Kinase, Phosphosite and Inhibitor information pages
#------------------------------------------------------------------------------------------------------------------------------------------

#Kinase page
@application.route('/kin/redirect/<kinase_name>', methods=['GET', 'POST'])
def kinase(kinase_name):
    form = SearchForm()
    searchkin = session.query(HumanKinases).filter(HumanKinases.Entry_name.startswith(kinase_name)).first()
    searchkinphos = session.query(Phosphosites).join(KinasesPhosphosites).join(HumanKinases).filter(HumanKinases.Entry_name==kinase_name).all()
    searchphosdis = session.query(PhosphositesDiseases).join(Phosphosites).join(KinasesPhosphosites).join(HumanKinases).filter(HumanKinases.Entry_name==kinase_name).all()
    searchkininhibitors = session.query(Inhibitors).join(InhibKin).join(HumanKinases).filter(HumanKinases.Entry_name==kinase_name).all()
    protein_name = None
    if form.validate_on_submit():
        protein_name = form.protein_name.data
        return redirect(url_for('kinasesearch', kinase_search=protein_name))
    return render_template('kinasepage.html', kinase_name=kinase_name, UniProt_ID=searchkin.UniProt_ID, PDB_ID=searchkin.PDB_ID, PDB_URL=searchkin.PDB_URL, PDB_title=searchkin.PDB_title, Primary_Protein_Name=searchkin.Primary_Protein_Name, Alternate_Protein_Name=searchkin.Alternative_Protein_Name, Families=searchkin.Families, AA_sequence=searchkin.AA_Seq, Molecular_Mass=searchkin.Molecular_Mass, Subcellular_Location=searchkin.Subcellular_Location, Gene_Symbol=searchkin.Gene_Symbol, Alternative_Gene_Name=searchkin.Alternative_Gene_Name, searchkinphos=searchkinphos, searchkininhibitors=searchkininhibitors, searchphosdis=searchphosdis, form=form)

#Phosphosite page
@application.route('/phos/<phosphosite_search>/<phosphosite_name>', methods=['GET', 'POST'])
def phosphositepage(phosphosite_search, phosphosite_name):
    form = SearchForm()
    searchphos = session.query(Phosphosites).filter(Phosphosites.PHOS_ID==phosphosite_name).first()
    searchphoskin = session.query(HumanKinases).join(KinasesPhosphosites).filter(KinasesPhosphosites.PHOS_ID==phosphosite_name).all()
    searchphosdis = session.query(PhosphositesDiseases).join(Phosphosites).filter(Phosphosites.PHOS_ID==phosphosite_name).all() #This line has been changed
    protein_name = None
    if form.validate_on_submit():
        protein_name = form.protein_name.data
        return redirect(url_for('phosphositesearchprot', phosphosite_search=protein_name))
    return render_template('phosphositepage.html', form=form, ID_PH=searchphos.ID_PH, GENE=searchphos.GENE, PROTEIN=searchphos.PROTEIN, HU_CHR_LOC=searchphos.HU_CHR_LOC, MOD_RSD=searchphos.MOD_RSD, MW_kD=searchphos.MW_kD, SITE_7_AA=searchphos.SITE_7_AA, DOMAIN=searchphos.DOMAIN, SOURCE=searchphos.SOURCE, searchphoskin=searchphoskin, searchphosdis=searchphosdis )


#Inhibitor page
@application.route('/inh/redirect/<inhibitor_name>', methods=['GET', 'POST'])
def inhibitor(inhibitor_name):
    form = SearchForm()
    searchinh = session.query(Inhibitors).filter(Inhibitors.Inhibitor==inhibitor_name).first()
    searchinhkin = session.query(HumanKinases).join(InhibKin).join(Inhibitors).filter(Inhibitors.Inhibitor==inhibitor_name).all()
    protein_name = None
    if form.validate_on_submit():
        protein_name = form.protein_name.data
        return redirect(url_for('inhibitorsearch', inhibitor_search=protein_name))
    return render_template('inhibitorpage.html', Inhibitor=inhibitor_name, Ki_nM=searchinh.Ki_nM, IC50_nM=searchinh.IC50_nM, Kd_nM=searchinh.EC50_nM, Source=searchinh.Source, IMG_URL=searchinh.IMG_URL, searchinhkin=searchinhkin, form=form)

#--------------------------------------------------------------------------------------------------------------------------  


#--------------------------------------------------------------------------------------------------------------------------  
if __name__ == '__main__':
    application.run()


