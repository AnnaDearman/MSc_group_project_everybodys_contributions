#Import statements
#--------------------------------------------------------------------------------------------------------------------------
#Basic modules for a Flask app to work (set Flask app, return HTML template for specific @application.route, redirect to a different
#@application.route, generate URL to a given endpoint).
from flask import Flask, render_template, redirect, url_for

#Modules necessary to set the search engine 
from flask_wtf import FlaskForm
from wtforms import StringField, SelectField, SubmitField
from wtforms.validators import DataRequired

#Modules necessary to upload a phosphoproteomics data file to the website, and to be able to download the tables returned
from flask import request, send_from_directory, send_file
import os
from werkzeug.utils import secure_filename

#Modules necessary to perform the data analysis of the phosphoproteomics data file uploaded by the user
from scipy.stats import norm
import math
import pandas as pd
import numpy as np
import plotly.express as px
import re

#Modules necessary to set the structure of the database used, and make queries to said database
from sqlalchemy import Column, Integer, String, ForeignKey, create_engine, or_
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker, relationship

engine = create_engine('sqlite:///c__sqlite_final_database_v8.db', echo = True, connect_args={'check_same_thread': False}) #Connect database to Flask app
Base = declarative_base() #Construct a base class for declarative class definitions.
Session = sessionmaker(bind = engine) #Set Session, necessary to query the database (see any of the classes under #Hits pages for example)
session = Session()
#--------------------------------------------------------------------------------------------------------------------------


#Flask application object & paths to the uploads & downloads folders
#-------------------------------------------------------------------------------------------------------------------------- 
#Set Flask application object
application = Flask(__name__, static_url_path='/static') #Need to define static_url_path in order for images in the static folder to show on the HTML templates
application.config['SECRET_KEY'] = 'jacky' #SECRET_KEY,necessary for wtforms to work

#Set the path to the uploads folder, in which the phosphoproteomics data files uploaded by the user are stored
UPLOAD_FOLDER =  os.path.dirname(os.path.abspath(__file__)) + '/uploads/'
application.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

#Set the path to the downloads folder, in which files produced from the data analysis of the phosphoproteomics data (see def data_analysis) are stored.
DOWNLOAD_FOLDER = os.path.dirname(os.path.abspath(__file__)) + '/downloads/'
application.config['DOWNLOAD_FOLDER'] = DOWNLOAD_FOLDER

#Define allowed extensions for the phosphoproteomics data file uploaded by the user
ALLOWED_EXTENSIONS = {'tsv'}
def allowed_file(filename):
   return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS
#-------------------------------------------------------------------------------------------------------------------------- 

#Search form
#--------------------------------------------------------------------------------------------------------------------------  
class SearchForm(FlaskForm):
    protein_name = StringField("", validators=[DataRequired()]) #Box in which to write your search, requires data
    submit = SubmitField('Submit') #Submit button
#--------------------------------------------------------------------------------------------------------------------------



#Database structure
#--------------------------------------------------------------------------------------------------------------------------
#Database tables

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
    PHOS_ID = Column(String(22), ForeignKey('phosphosites.PHOS_ID5'))  # duplicates
    ID_DP = Column(String(9), primary_key=True)

class InhibKin(Base):
    __tablename__ = 'inhib_kin'
    UniProt_ID = Column(String(20), ForeignKey('human_kinases.UniProt_ID'))
    BindingDB_ID = Column(String(30))
    chEMBL_ID = Column(String(30))
    Molecule_name = Column(String(100), ForeignKey('inhibitors.Molecule_name'))
    IN_KI = Column(String(20), primary_key=True)
    
    
#Relationships between the database tables
kinases_phosphosites = relationship("KinasesPhosphosites", backref="human_kinases")
inhibitors_kinases = relationship("Inhibitors", backref="human_kinases")
phosphosites_kinases = relationship("KinasesPhosphosites", backref="phosphosites")
phosphosites_diseases = relationship("PhosphositesDiseases", backref="phosphosites")
inhib_kin = relationship("InhibKin", backref="inhibitors")

Base.metadata.create_all(engine) #Create the tables in the database 'sqlite:///c__sqlite_final_database.db'
#--------------------------------------------------------------------------------------------------------------------------


#Data Analysis of phosphoproteomics data (from file uploaded by the user)
#----------------------------------------------------------------------------------------------------------------------------
def data_analysis(path, filename):

    ### import sample data
    data_raw = pd.read_table(path, sep = "\t")

    ### drop unexpected Unnamed columns
    data_raw = data_raw.loc[:, ~data_raw.columns.str.contains('^Unnamed')]

    ### replace inf or -inf with NaN    
    data_raw = data_raw.replace([np.inf, -np.inf], np.nan)

    ## drop any rows that have a NAN
    global df1
    df1 = data_raw.copy()
    df1 = df1.dropna(how="any")
    
    ### drop first 2 columns to have df2
    global df2
    df2 = df1.iloc[:, 2:]

    ### set a STEP variable to show whether CV columns are included
    if df2.shape[1]%3 == 0:
        STEP = 3
    elif df2.shape[1]%5 == 0:
        STEP = 5

    ### extract the inhibitors from the df2
    global Inhibitor
    Inhibitor = list() 
    for x in list(df2.columns)[::STEP]:
        aa = x.split('_')[0]
        Inhibitor.append(aa)

    ### if the length of STEP equals the column length of the df2, just return df1
    if STEP == df2.shape[1]:
        # Reset the index
        df1.reset_index(drop=True, inplace=True)
        
    ### create a series used as the index of the transform table
        Inh_list = []
        for x in Inhibitor:
            Inh_list.extend([x]*df1.shape[0])
        Inh_series = pd.Series(Inh_list, name = "Inhibitor")

        ### concatenate df_multiple, df3 and Inh_series
        global df_transform
        df_transform = pd.concat([df1, Inh_series], axis = 1)
    
    elif STEP < df2.shape[1]:
      
        ### reshape the table by appending the repeated columns to the first columns
        global df3
        df3 = df2.iloc[:, 0:STEP]
        if STEP == 5: # If statement added by AD to include STEP == 3
            df3.columns = ['Treat_mean','Fold_change', 'p_value', 'ctrlCV', 'treatCV']
        else:
            df3.columns = ['Treat_mean','Fold_change', 'p_value']
        for i in range(STEP, df2.shape[1], STEP):
            df = df2.iloc[:, i:i+STEP]
            if STEP == 5: # If statement added by AD to include STEP == 3
                df.columns = ['Treat_mean','Fold_change', 'p_value', 'ctrlCV', 'treatCV']
            else:
                df.columns = ['Treat_mean','Fold_change', 'p_value']    
            df3 = df3.append(df, ignore_index = True)

        ### form a table with column Substrate and control_mean row binded to itself
        global df_scm
        df_scm = df1.iloc[:, 0:2]
        global df_multiple
        df_multiple = df_scm
        i = 0
        while i < len(Inhibitor)-1:
            i += 1
            df_add = df_scm
            df_multiple = df_multiple.append(df_add, ignore_index = True)

        ### create a series used as the index of the transform table
        Inh_list = []
        for x in Inhibitor:
            Inh_list.extend([x]*df1.shape[0])
        Inh_series = pd.Series(Inh_list, name = "Inhibitor")
    
        ### concatenate df_multiple, df3 and Inh_series
        df_transform = pd.concat([df_multiple, df3, Inh_series], axis = 1)

        ################################
    global inh_len
    inh_len = len(Inhibitor)
    global df_submeans
    df_submeans = []
    ### import sample data
    for i in range(inh_len):
        global inh
        inh = Inhibitor[i]
        df1 = df_transform[df_transform.Inhibitor == inh] 
        df1 = df1.drop(['Inhibitor'], axis = 1)
        if len(df1.columns) == 7:
            df1.columns = ['Substrate', 'Control_mean', 'Treat_mean','Fold_change', 'p_value', 'ctrlCV', 'treatCV']
        elif len(df1.columns) == 5:
            df1.columns = ['Substrate', 'Control_mean', 'Treat_mean','Fold_change', 'p_value']
        
        ### split the first column
        df1[["Gene", "sites"]] = df1.Substrate.str.split("(", expand = True,)
        df1[["Residue", "junk"]] = df1.sites.str.split(")", expand = True,)
        df1 = df1[['Substrate','Fold_change','p_value', 'Gene', 'Residue']]

        ### drop rows contain None in Residue column
        df1 = df1[df1.Residue != "None"]

        ### drop rows contain 0 in Fold_change column
        df1 = df1[df1.Fold_change != 0]

        ### transform Fold_change to Log2 and then drop original Fold_change column
        df1["Log_Fold_Change"] = np.log2(df1['Fold_change'])
        df1 = df1.drop(['Fold_change'], axis = 1)

        ### transform p_value to -log10 p_value for the volcano plot
        df1["-log10p"] = -np.log10(df1["p_value"])

        ######### Merge user submitted table with kinase-substrate table
        df_phos = pd.read_sql(session.query(KinasesPhosphosites).statement, session.query(KinasesPhosphosites).session.bind)
        df_merge = pd.merge(df1, df_phos, how = 'left', left_on = 'Substrate', right_on = 'PHOS_ID4')
        df_merge = df_merge[['Substrate','p_value','-log10p','Gene','Residue','Log_Fold_Change','GENE']]
        df_merge['GENE'].fillna('Not_Found', inplace=True)
        df_merge.columns = ['Substrate','p_value','-log10p','Gene','Residue','Log_Fold_Change','Kinase']
        df_merge = df_merge.drop_duplicates(subset = 'Substrate', keep = False)
        df1 = df_merge
        df8 = df1.copy()
        df8 = df8[['Substrate','Log_Fold_Change','p_value','Kinase']]
 
        df8.to_csv(os.path.join(application.config['DOWNLOAD_FOLDER'], f'table1{inh}_analysis.csv')) #Convert table to csv and save on the downloads folder

        ### drop rows contain Not_Found in Kinase column
        global df_kinase
        df_kinase = df1[df1.Kinase != "Not_Found"]
        global df_kinase2
        df_kinase2 = df_kinase.copy()
        df_kinase2 = df_kinase2[['Substrate','Kinase','p_value','Log_Fold_Change']]

        df_kinase2.to_csv(os.path.join(application.config['DOWNLOAD_FOLDER'], f'table2{inh}_analysis.csv')) #Convert table to csv and save on the downloads folder

    
##########################################################################
### Function to get the table with substrates without kinases
        ### select rows contain Not_Found in Kinase column
        global df_no_kinase
        df_no_kinase = df1[df1.Kinase == "Not_Found"]
        global df_no_kinase2
        df_no_kinase2 = df_no_kinase.copy()
        df_no_kinase2 = df_no_kinase2[['Substrate','Kinase','p_value','Log_Fold_Change']]

        df_no_kinase2.to_csv(os.path.join(application.config['DOWNLOAD_FOLDER'], f'table3{inh}_analysis.csv')) #Convert table to csv and save on the downloads folder

    ### calculate the subgroup mean
        global df_submean
        df_submean = df_kinase.groupby(['Kinase']).mean()
        df_submean = df_submean.drop(['p_value'], axis = 1)
        df_submean = df_submean.drop(['-log10p'], axis = 1)
    ### calculate Z_score and convert Z_score to P_value
        def Z_score(mS):
            Z_score = (mS - df1['Log_Fold_Change'].mean())*(len(df1['Log_Fold_Change'])**(1/2))/np.std(df1['Log_Fold_Change'])
            P_value = norm.sf(abs(Z_score))
            return(P_value)
    #Z_score = (mS - df1['Log_Fold_Change'].mean())*(len(df1['Log_Fold_Change'])**(1/2))/np.std(df1['Log_Fold_Change'])
    #P_value = norm.sf(abs(Z_score))

    ### calculate Z_score and convert Z_score to P_value
        df_submean['p_value'] = df_submean['Log_Fold_Change'].apply(Z_score)
    ### calculate the standard deviation of fold change in each subgroup
        df_std = df_kinase.groupby(['Kinase']).aggregate(np.std)
        df_submean['Std_Dev'] = df_std['Log_Fold_Change']

    ### Calculate the relative kinase activity
        df_submean['Kinase_RAS'] = df_submean['Log_Fold_Change']/df1['Log_Fold_Change'].mean()

    ### add a kinase column to df_submean
        df_submean['Kinase'] = df_submean.index

    # Resets the index
        df_submean.reset_index(drop=True, inplace=True)

    ### set a count colum to show the number of substrates of each group
        df_submean['Number_of_substrates'] = list(df_kinase['Kinase'].value_counts().sort_index())
        global df_submean2
        df_submean2 = df_submean.copy()
        df_submean2 = df_submean2[['Kinase','Kinase_RAS','Log_Fold_Change','p_value','Number_of_substrates','Std_Dev']]
        df_submean2.to_csv(os.path.join(application.config['DOWNLOAD_FOLDER'], f'table4{inh}_analysis.csv'))

        global df_submean2_top10
        global df_submean2_bot10
        global df_submean2_20
        global RKA_cols
        df_submean2_top10 = df_submean2.sort_values(by=['Kinase_RAS'],ascending=False)
        df_submean2_top10 = df_submean2_top10.reset_index(drop=True)
        df_submean2_top10 = df_submean2_top10.loc[0:10]
        df_submean2_bot10 = df_submean2.sort_values(by=['Kinase_RAS'],ascending=True)
        df_submean2_bot10 = df_submean2_bot10.reset_index(drop=True)
        df_submean2_bot10 = df_submean2_bot10.loc[0:10]
        df_submean2_bot10 = df_submean2_bot10.sort_values(by=['Kinase_RAS'],ascending=False)
        df_submean2_bot10 = df_submean2_bot10.reset_index(drop=True)
        df_submean2_20 = pd.concat([df_submean2_top10,df_submean2_bot10],axis=0)
        df_submean2_20 = df_submean2_20.reset_index(drop=True)
        RKA_cols = df_submean2_20.columns.values
        df_submean2_20 = df_submean2_20.to_html(classes="data")
        df_submeans.append(df_submean2_20)

        Bon_cor = -math.log10(0.05/len(df1['p_value']))

    ### set a Significant? colum to show significance of the p_value
        df1.loc[df1['-log10p'] < Bon_cor, 'Significant?'] = 'Non-sig'
        df1.loc[df1['-log10p'] >= Bon_cor, 'Significant?'] = 'Sig'

    ### volcano plot ### also need to install nbformat before running plotly
        global fig
        fig = px.scatter(df1, x = "Log_Fold_Change", y="-log10p", hover_data=['Substrate'], color = 'Significant?')
        fig.write_html(f"templates/Volcano{inh}.html") #Convert the figure to HTML, save it on the templates folder

    ### bar plot of kinase relative activity scores
        global df_subgroub
        df_subgroub = df_submean
        df_subgroub.loc[df_subgroub['Kinase_RAS'] < 0, 'Effect'] = 'Treat_down'
        df_subgroub.loc[df_subgroub['Kinase_RAS'] >= 0, 'Effect'] = 'Treat_up'

    ### set a text colum to show significance of the p_value
        df_subgroub.loc[df_subgroub['p_value'] <= 0.05, 'Significance'] = '*'
        df_subgroub.loc[df_subgroub['p_value'] > 0.05, 'Significance'] = ''
        df_subgroub.loc[df_subgroub['p_value'] == 0, 'Significance'] = ' '

    ### show relative activity by a bar plot
        global fig2
        fig2 = px.bar(df_subgroub.sort_values(by=['Kinase_RAS']),
            x = 'Kinase',
            y = 'Kinase_RAS',
            color = 'Effect',
            text = 'Significance',
            error_y = 'Std_Dev',
            hover_data = ['p_value', 'Number_of_substrates'])
    
        fig2.update_traces(textposition='outside')
        fig.update_layout(title_text=inh)
        fig2.write_html(f"templates/Kinase_RKA_barplot{inh}.html") #Convert the figure to HTML, so it can be accessed on the web application
#--------------------------------------------------------------------------------------------------------------------



#Main application script

#Homepage(s)
#--------------------------------------------------------------------------------------------------------------------------  
#Main homepage
@application.route('/', methods=['GET', 'POST'])
@application.route('/home', methods=['GET', 'POST']) 
def index():
    form = SearchForm() 
    #If the user has submitted form (in this clase, clicking on the 'Upload' button of the homepage)
    if request.method == 'POST': 
        if 'file' not in request.files: # Check that the post request has the file part
           print('No file attached in request')
           return redirect(request.url)
        file = request.files['file'] 
        #Check that the user has selected a file, if they haven't, return a message asking them to return to the homepage and select one
        if file.filename == '': 
           return "You have not selected a file, please return to the homepage, click on 'Choose File', select a phosphoproteomics data file (tsv extension) to upload, and the click on 'Upload'."
        #If a file has been selected, and it is in the allowed extension (tsv)
        if file and allowed_file(file.filename): 
           filename = secure_filename(file.filename)
           file.save(os.path.join(application.config['UPLOAD_FOLDER'], filename)) #Save uploaded file to the uploads folder
           data_analysis(os.path.join(application.config['UPLOAD_FOLDER'], filename), filename) #Perform data analysis of the uploaded folder (def data_analysis function above)
           return redirect(url_for('uploaded_file', filename=filename)) #Return url showing the results of the data_analysis function
        #If the file selected is not in the allowed extension, return message asking to upload file with tsv extension
        else: 
            return "The file uploaded is not in tsv format, please upload a file in this format"
    return render_template('homepage.html', form=form)

#From homepage.html redirects to this URL if you select 'Kinase' on the 'Search' dropdown select menu
@application.route('/kin/', methods=['GET', 'POST']) 
def kin():
    form = SearchForm()
    protein_name = None
    if form.validate_on_submit(): #If the user writes in the searchbox (StringField) and clicks Submit/press Enter
        protein_name = form.protein_name.data #What the user has written in the searchbox
        return redirect(url_for('kinasesearch', kinase_search=protein_name)) #Redirect to def kinasesearch(kinase_search)
    return render_template('homepagekin.html', form=form) #Same appearance as the homepage but the dropdown select menu displays the word 'Kinase' instead of 'Search'


#From homepage.html redirects to this URL if you select 'Inhibitor' on the 'Search' dropdown select menu
@application.route('/inh/', methods=['GET', 'POST']) 
def inh():
    form = SearchForm()
    protein_name = None
    if form.validate_on_submit(): #If the user writes in the searchbox (StringField) and clicks Submit/press Enter
        protein_name = form.protein_name.data #What the user has written in the searchbox
        return redirect(url_for('inhibitorsearch', inhibitor_search=protein_name))#Redirect to def inhibitorsearch(inhibitor_search)
    return render_template('homepageinh.html',form=form) #Same appearance as the homepage.html but the dropdown select menu displays the word 'Inhibitor' instead of 'Search'

#From homepage.html redirects to this URL if you select 'Phosphosite' on the 'Search' dropdown select menu
@application.route('/phos/', methods=['GET', 'POST']) 
def phos():
    form = SearchForm()
    return render_template('homepagephos.html', form=form) #Same appearance as homepage.html but the dropdown select menu displays the word 'Phosphosite', and there is an additional dropdown select menu (genomic location & neighbouring sequence) compared to the homepage.html

#From homepagephos.html redirects to this URL if you select 'Phosphosite' on the 'Search' dropdown select menu and then 'Neighbouring Sequence' on the second dropdown select menu
@application.route('/nseq/', methods=['GET', 'POST']) 
def phosnseq():
    form = SearchForm()
    protein_name = None
    if form.validate_on_submit(): #If the user writes in the searchbox (StringField) and clicks Submit/press Enter
        protein_name = form.protein_name.data #What the user has written in the searchbox
        return redirect(url_for('phosphositesearchnseq', phosphosite_search=protein_name)) #Redirect to def phosphositesearchnseq(phosphosite_search)
    return render_template('homepagephos.html', form=form) #Same appearance as homepage.html but the dropdown select menu displays the word 'Phosphosite', and there is an additional dropdown select menu (genomic location & neighbouring sequence) compared to the homepage.html 

#From homepagephos.html redirects to this URL if you select 'Phosphosite' on the 'Search' dropdown select menu and then 'Protein' on the second dropdown select menu
@application.route('/prot/', methods=['GET', 'POST']) 
def phosprot():
    form = SearchForm()
    protein_name = None
    if form.validate_on_submit(): #If the user writes in the searchbox (StringField) and clicks Submit/press Enter
        protein_name = form.protein_name.data #What the user has written in the searchbox
        return redirect(url_for('phosphositesearchprot', phosphosite_search=protein_name)) #Redirect to def phosphositesearchprot(phosphosite_search)
    return render_template('homepagephos.html', form=form) #Same appearance as homepage.html but the dropdown select menu displays the word 'Phosphosite', and there is an additional dropdown select menu (genomic location & neighbouring sequence) compared to the homepage.html

#From homepagephos.html redirects to this URL if you select 'Phosphosite' on the 'Search' dropdown select menu and then 'Protein' on the second dropdown select menu
@application.route('/phoskin/', methods=['GET', 'POST']) 
def phospkin():
    form = SearchForm()
    protein_name = None
    if form.validate_on_submit(): #If the user writes in the searchbox (StringField) and clicks Submit/press Enter
        protein_name = form.protein_name.data #What the user has written in the searchbox
        return redirect(url_for('phosphositesearchkin', phosphosite_search=protein_name)) #Redirect to def phosphositesearchkin(phosphosite_search)
    return render_template('homepagephos.html', form=form) #Same appearance as homepage.html but the dropdown select menu displays the word 'Phosphosite', and there is an additional dropdown select menu (genomic location & neighbouring sequence) compared to the homepage.html



#The next routes are to redirect to the intended url

#Return to the @application.route('/inh/') if you select 'Inhibitor' on the first select dropdown menu from any other page on the website
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

#Return to the @application.route('/phos/') if you select 'Phosphosite' on the first select dropdown menu from any other page on the website
@application.route('/inh/phos/')
@application.route('/kin/phos/')
@application.route('/inh/redirect/phos/')
@application.route('/kin/redirect/phos/')
@application.route('/uploads/phos/')
def redphos():
    return redirect('/phos/')

#Return to the @application.route('/kin/') if you select 'Kinase' on the select dropdown menu from any other page on the website
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

#Return to the @application.route('/nseq/') if you select 'Neighbouring Sequence' on the second select dropdown menu from a different URL
@application.route('/prot/nseq')
@application.route('/gen/nseq')
def rednseq():
    return redirect('/nseq/')

#Return to the @application.route('/prot/') if you select 'Protein' on the second select dropdown menu from a different URL
@application.route('/nseq/prot')
@application.route('/gen/prot')
def redprot():
    return redirect('/prot/')

#Return to the @application.route('/gen/') if you select 'Genomic Location' on the second select dropdown menu from a different URL
@application.route('/nseq/gen')
@application.route('/prot/gen')
def redgen():
    return redirect('/gen/')
#------------------------------------------------------------------------------------------------------------------------------------------   



#Data analysis page
#-----------------------------------------------------------------------------------------------------------------------------------------------------
#Display the results of the data analysis
@application.route('/uploads/<filename>')
def uploaded_file(filename):
    form = SearchForm()
    return render_template('datanalysis.html', len= len(Inhibitor), Inhibitor= Inhibitor, tables4=[df_submeans], titles4=RKA_cols, form=form)

#Download the phosphosites-kinases table with kinases match, and phosphosites with no kinase match (from the downloads folder)
@application.route('/download/table1/<inhibitor>')
def download_table1(inhibitor):
    return send_from_directory(application.config['DOWNLOAD_FOLDER'], filename='table1'+inhibitor+'_analysis.csv' , as_attachment=True)

#Download the phosphosites-kinases table with only the phosphosites for which a kinase match was found (from the downloads folder)
@application.route('/download/table2/<inhibitor>')
def download_table2(inhibitor):
    return send_from_directory(application.config['DOWNLOAD_FOLDER'], filename='table2'+inhibitor+'_analysis.csv' , as_attachment=True)

#Download the phosphosites-kinases table with the phosphosites for which no kinase match was found (from the downloads folder)
@application.route('/download/table3/<inhibitor>')
def download_table3(inhibitor):
    return send_from_directory(application.config['DOWNLOAD_FOLDER'], filename='table3'+inhibitor+'_analysis.csv'  , as_attachment=True)

#Download the relative kinase activity table (from the downloads folder)
@application.route('/download/table4/<inhibitor>')
def download_table4(inhibitor):
    return send_from_directory(application.config['DOWNLOAD_FOLDER'], filename='table4'+inhibitor+'_analysis.csv' , as_attachment=True)

#Show Interactive version of the volcano plot
@application.route('/volcanoplot/<inhibitor>')
def volcanoplot(inhibitor):
    template= 'Volcano'+inhibitor+'.html'
    return render_template(template)

#Show Interactive version of the relative kinase activity bar plot
@application.route('/RKAplot/<inhibitor>')
def RKAplot(inhibitor):
    template= 'Kinase_RKA_barplot' + inhibitor+'.html'
    return render_template(template)
#---------------------------------------------------------------------------------------------------------------------------------------------------- 

#Hits pages
#------------------------------------------------------------------------------------------------------------------------------------------
#Inhibitor search results page
@application.route('/inh/<inhibitor_search>', methods=['GET', 'POST'])
def inhibitorsearch(inhibitor_search):
    form = SearchForm()
    inhibitor_search = inhibitor_search.upper() #Make sure that the inhibitor entered by the user is in uppercase
    # resultss = session.query(Inhibitors).filter(Inhibitors.Inhibitor.startswith(inhibitor_search)).all() #Query the Inhibitors table on the database for entries on the Inhibitor column that start with the inhibitor entered by the user
    resultss = session.query(Inhibitors).filter(or_(Inhibitors.BindingDB_ID.startswith(inhibitor_search), Inhibitors.chEMBL_ID.startswith(inhibitor_search), Inhibitors.Molecule_name.startswith(inhibitor_search), Inhibitors.Synonyms.contains(inhibitor_search)) ).all()
    protein_name = None
    if form.validate_on_submit(): #If the user searches for an inhibitor while on this page
        protein_name = form.protein_name.data
        return redirect(url_for('inhibitorsearch', inhibitor_search=protein_name))
    return render_template('inhibitorhits.html', resultss=resultss, form=form)

#Kinase search results page
@application.route('/kin/<kinase_search>', methods=['GET', 'POST'])
def kinasesearch(kinase_search):
    form = SearchForm()
    kinase_search = kinase_search.upper() #Make sure that the kinase entered by the user is in uppercase
    resultss = session.query(HumanKinases).filter(or_(HumanKinases.Entry_name.startswith(kinase_search), HumanKinases.Gene_Symbol.startswith(kinase_search)) ).all() #Query the HumanKinases table on the database for entries on the Entry_name or Gene_Symbol column that start with the kinase entered by the user
    protein_name = None
    if form.validate_on_submit(): #If the user searches for a kinase while being on this page
        protein_name = form.protein_name.data
        return redirect(url_for('kinasesearch', kinase_search=protein_name))
    return render_template('kinasehits.html', resultss=resultss, form=form)
 
#Phosphosite search by neighbouring sequence results page   
@application.route('/nseq/<phosphosite_search>', methods=['GET', 'POST'])
def phosphositesearchnseq(phosphosite_search):
    form = SearchForm()
    resultss = session.query(Phosphosites).filter(Phosphosites.SITE_7_AA.startswith(phosphosite_search)).all() #Query the Phosphosites table on the database for entries on the SITE_7_AA column that start with the sequence entered by the user
    protein_name = None
    if form.validate_on_submit(): #If the user searches a phosphosite by neighbouring sequence while being on this page
        protein_name = form.protein_name.data
        return redirect(url_for('phosphositesearchnseq', phosphosite_search=protein_name))

    return render_template('phoshits.html', resultss=resultss, form=form)


#Phosphosite search by protein results page
@application.route('/prot/<phosphosite_search>', methods=['GET', 'POST'])
def phosphositesearchprot(phosphosite_search):
    form = SearchForm()
    resultss = session.query(Phosphosites).filter(or_(Phosphosites.PROTEIN.startswith(phosphosite_search),Phosphosites.GENE.startswith(phosphosite_search))).all() #Query the Phosphosites table on the database for entries on the PROTEIN column that start with the protein entered by the user
    protein_name = None
    if form.validate_on_submit(): #If the user searches a phosphosite by protein while being on this page
        protein_name = form.protein_name.data
        return redirect(url_for('phosphositesearchprot', phosphosite_search=protein_name))
    return render_template('phoshits.html', resultss=resultss, form=form)

#Phosphosite search by 'which kinase they're phosphorylated by' results page
@application.route('/phoskin/<phosphosite_search>', methods=['GET', 'POST'])
def phosphositesearchkin(phosphosite_search):
    form = SearchForm()
    resultss = session.query(Phosphosites).join(KinasesPhosphosites).join(HumanKinases).filter(or_(HumanKinases.Entry_name.startswith(phosphosite_search), HumanKinases.Gene_Symbol.startswith(phosphosite_search)) ).all() #Join three tables from the database (Phosphosites, KinasesPhosphosites & HumanKinases), query the joint table for entries on the Entry_name or Gene_Symbol column that start with the kinase entered by the user
    protein_name = None
    if form.validate_on_submit(): #If the user searches a phosphosite by kinase that phosphorylates it while being on this page
        protein_name = form.protein_name.data
        return redirect(url_for('phosphositesearchkin', phosphosite_search=protein_name))
    return render_template('phoshits.html', resultss=resultss, form=form)

#Phosphosite browse by genomic location page
@application.route('/gen/', methods=['GET', 'POST'])
def phosphositesearchgenhome():
    form = SearchForm()
    resultss = session.query(Phosphosites).filter(Phosphosites.HU_CHR_LOC.startswith('t')).all() #Empty query, we just want to return the headers of resultss on phoshitsgen.html
    if form.validate_on_submit(): #If the user types a genomic location on the genomic browser
        protein_name = form.protein_name.data
        return redirect(url_for('phosphositesearchgen', gen_location=protein_name))
    return render_template('zphoshitsgen.html',form=form, resultss=resultss)

#Phosphosite browse by genomic location page, once you have selected or typed in a genomic location
@application.route('/gen/<gen_location>', methods=['GET', 'POST'])
def phosphositesearchgen(gen_location):
    form = SearchForm()
    resultss = session.query(Phosphosites).filter(Phosphosites.HU_CHR_LOC.startswith(gen_location)).all() #Query the Phosphosites table on the database for entries on the HU_CHR_LOC column that start with the chromosomal location entered/selected by the user 
    protein_name = None
    if form.validate_on_submit(): #If the user types in the chromosomal location 
        protein_name = form.protein_name.data
        gen_location=protein_name
        # Since we are using SQLAlchemy queries with 'startwith', if the user selects 1 (for Chromosome 1) it will return not
        #only phosphosites on Chromosome 1 but also Chromosome 10,11,etc, and if the user selects 2 (for Chromosome 2) it will 
        #also return phosphosites on Chromosomes 21 and 22, therefore:
        if gen_location== '1': #If the user selects/types in 1, return only the phosphosites on Chromosome 1
            p= session.query(Phosphosites).filter(Phosphosites.HU_CHR_LOC.startswith('1p')).all()
            q= session.query(Phosphosites).filter(Phosphosites.HU_CHR_LOC.startswith('1q')).all()
            return render_template('zphoshitsgenalt.html',form=form, p=p, q=q)
        if gen_location== '2': #If the user selects/types in 2, return only the phosphosites on Chromosome 2
            p2= session.query(Phosphosites).filter(Phosphosites.HU_CHR_LOC.startswith('2p')).all()
            q2= session.query(Phosphosites).filter(Phosphosites.HU_CHR_LOC.startswith('2q')).all()
            return render_template('zphoshitsgenalt.html',form=form, q=q2, p=p2)
        else: #If the user selects/types in a chromosomal location other than 1 or 2, simply return the resultss query on zphoshitsgen.html
            return render_template('zphoshitsgen.html',form=form, resultss=resultss)
    else: #If the user selects the chromosomal location on the HTML image map 
        if gen_location== '1': #If the user selects/types in 1, return only the phosphosites on Chromosome 1
            p= session.query(Phosphosites).filter(Phosphosites.HU_CHR_LOC.startswith('1p')).all()
            q= session.query(Phosphosites).filter(Phosphosites.HU_CHR_LOC.startswith('1q')).all()
            return render_template('zphoshitsgenalt.html',form=form, p=p, q=q)
        if gen_location== '2': #If the user selects/types in 2, return only the phosphosites on Chromosome 2
            p2= session.query(Phosphosites).filter(Phosphosites.HU_CHR_LOC.startswith('2p')).all()
            q2= session.query(Phosphosites).filter(Phosphosites.HU_CHR_LOC.startswith('2q')).all()
            return render_template('zphoshitsgenalt.html',form=form, q=q2, p=p2)
        else: #If the user selects/types in a chromosomal location other than 1 or 2, simply return the resultss query on zphoshitsgen.html
            return render_template('zphoshitsgen.html',form=form, resultss=resultss)
    return render_template('zphoshitsgen.html',form=form, resultss=resultss)
#--------------------------------------------------------------------------------------------------------------------------  


#Kinase, Phosphosite and Inhibitor information pages
#------------------------------------------------------------------------------------------------------------------------------------------

#Kinase page
@application.route('/kin/redirect/<kinase_name>', methods=['GET', 'POST'])
def kinase(kinase_name):
    form = SearchForm()
    searchkin = session.query(HumanKinases).filter(HumanKinases.Entry_name.startswith(kinase_name)).first() #Query the HumanKinases table on the database for entries on the Entry_name column that start with the kinase selected on the kinase hits page
    searchkinphos = session.query(Phosphosites).join(KinasesPhosphosites).join(HumanKinases).filter(HumanKinases.Entry_name==kinase_name).all() ##Join three tables from the database (Phosphosites, KinasesPhosphosites & HumanKinases), query the joint table for entries on the Entry_name column that start with the kinase selected on the kinase hits page
    searchphosdis = session.query(PhosphositesDiseases).join(Phosphosites).join(KinasesPhosphosites).join(HumanKinases).filter(HumanKinases.Entry_name==kinase_name).all() #Join four tables from the database (PhosphositesDiseases, Phosphosites, KinasesPhosphosites & HumanKinases), query the joint table for entries on the Entry_name column that start with the kinase selected on the kinase hits page, query the joint table for entries on the Entry_name column that start with the kinase selected on the kinase hits page
    searchkininhibitors = session.query(Inhibitors).join(InhibKin).join(HumanKinases).filter(HumanKinases.Entry_name==kinase_name).all() #Join three tables from the database (Inhibitors, InhibKin & HumanKinases), query the joint table for entries on the Entry_name column that start with the kinase selected on the kinase hits page 
    protein_name = None
    if form.validate_on_submit(): #If the user searches a Kinase from this page
        protein_name = form.protein_name.data
        return redirect(url_for('kinasesearch', kinase_search=protein_name))
    return render_template('kinasepage.html', kinase_name=kinase_name, UniProt_ID=searchkin.UniProt_ID, PDB_ID=searchkin.PDB_ID, PDB_URL=searchkin.PDB_URL, PDB_title=searchkin.PDB_title, Primary_Protein_Name=searchkin.Primary_Protein_Name, Alternate_Protein_Name=searchkin.Alternative_Protein_Name, Families=searchkin.Families, AA_sequence=searchkin.AA_Seq, Molecular_Mass=searchkin.Molecular_Mass, Subcellular_Location=searchkin.Subcellular_Location, Gene_Symbol=searchkin.Gene_Symbol, Alternative_Gene_Name=searchkin.Alternative_Gene_Name, searchkinphos=searchkinphos, searchkininhibitors=searchkininhibitors, searchphosdis=searchphosdis, form=form)

#Phosphosite page
@application.route('/phos/<phosphosite_search>/<phosphosite_name>', methods=['GET', 'POST'])
def phosphositepage(phosphosite_search, phosphosite_name):
    form = SearchForm()
    searchphos = session.query(Phosphosites).filter(Phosphosites.PHOS_ID5==phosphosite_name).first() #Query the Phosphosites table on the database for entries on the PHOS_ID column that start with the phosphosite selected on the phosphosite hits page
    searchphoskin = session.query(HumanKinases).join(KinasesPhosphosites).filter(KinasesPhosphosites.PHOS_ID5==phosphosite_name).all() #Join two tables from the database (HumanKinases, KinasesPhosphosites) to get the entries on the HumanKinases table for kinases which in the KinasesPhosphosites table are in the same row as the phosphosite (PHOS_ID) selected on the phosphosite hits page
    searchphosdis = session.query(PhosphositesDiseases).join(Phosphosites).filter(Phosphosites.PHOS_ID5==phosphosite_name).all() #Join two tables from the database (PhosphositesDiseases, Phosphosites) to isolate entries on the PhosphositesDiseases table for the phosphosite selected on the phosphosite hits page
    return render_template('phosphositepage.html', form=form, ACC=searchphos.ACC_ID, ID_PH=searchphos.ID_PH, GENE=searchphos.GENE, PROTEIN=searchphos.PROTEIN, HU_CHR_LOC=searchphos.HU_CHR_LOC, MOD_RSD=searchphos.MOD_RSD, MW_kD=searchphos.MW_kD, SITE_7_AA=searchphos.SITE_7_AA, DOMAIN=searchphos.DOMAIN, SOURCE=searchphos.SOURCE, searchphoskin=searchphoskin, searchphosdis=searchphosdis )

#Inhibitor page
@application.route('/inh/redirect/<inhibitor_name>', methods=['GET', 'POST'])
def inhibitor(inhibitor_name):
    form = SearchForm()
    searchinh = session.query(Inhibitors).filter(Inhibitors.Molecule_name==inhibitor_name).first() #Query the Inhibitor table on the database for entries on the Inhibitor columns that start with the inhibitor selected on the inhibitor hits page
    searchinhkin = session.query(HumanKinases).join(InhibKin).join(Inhibitors).filter(Inhibitors.Molecule_name==inhibitor_name).all() #Join three tables from the database (HumanKinases, InhibKin, Inhibitors) to isolate the entries on the HumanKinase table that are for kinases that on the InhibKin table are on the same row as the inhibitor selected on the inhibitor hits page
    protein_name = None
    if form.validate_on_submit(): #If the user searches an Inhibitor from this page
        protein_name = form.protein_name.data
        return redirect(url_for('inhibitorsearch', inhibitor_search=protein_name))
    return render_template('inhibitorpage.html', Molecular_weight=searchinh.Molecular_weight, Molecular_formula=searchinh.Molecular_formula, Molecule_type=searchinh.Molecule_type, Synonyms=searchinh.Synonyms, chEMBL_URL=searchinh.chEMBL_URL, BindingDB_ID=searchinh.BindingDB_ID, Molecule_name=inhibitor_name, Ki_nM=searchinh.Ki_nM, EC50_nM=searchinh.EC50_nM, IC50_nM=searchinh.IC50_nM, Kd_nM=searchinh.EC50_nM, chEMBL_ID=searchinh.chEMBL_ID, searchinhkin=searchinhkin, form=form)

#--------------------------------------------------------------------------------------------------------------------------  


#--------------------------------------------------------------------------------------------------------------------------  
if __name__ == '__main__':
    application.run()