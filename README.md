## JACKY: Software Development Group Project 2020
**JACKY** is a new web-based software that unifies kinase, phosphosite, disease and inhibitor information collected form different sources (UniProt, PDB, phosphosite.org, phospho.ELM and KIDFamMap) in a single tool. 

The user, individuals interesting in kinases in an academic or research setting, not olny will find any required information on the topic available in the internet in this platform but also they will be allowed to upload their experimental quantitative phosphoproteomics results to our data analysis tool to compute estimated relative kinase activity scores.

Visit JACKY at http://jacky-03.ehym3crjpy.eu-west-2.elasticbeanstalk.com/

## Getting Started
The following instructions will guide you on how to install JACKY on your local machine for development and testing purposes.

### Prerequisites
JACKY was developed using Ancaonda distribution and Python 3.6. To run JACKY in you personal computer please firstly create a virtual environment using Anaconda prompt as follows:
```
conda create -n jacky_env python=3.6
```
Then install the packages listed below using

```
conda create -n jacky_env name-of-the-package
```



- certifi, 2019.11.28
- click, 7.0
- flask, 1.1.1
- flask-wtf, 0.14.2
- itsdangerous, 1.1.0
- jinja2, 2.11.0
- jquery, 3.4.1
- markupsafe, 1.1.1
- numpy, 1.18.1
- pandas, 0.25.3
- pip, 20.0.2
- plotly, 4.5.0
- python, 3.6.10
- python-dateutil, 2.8.1
- pytz, 2019.3
- retrying, 1.3.3
- scipy, 1.4.1
- setuptools, 45.1.0
- six, 1.14.0
- sqlalchemy, 1.3.12
- sqlite, 3.30.1
- vc, 14.1
- vs2015_runtime, 14.16.27012
- werkzeug, 0.16.1
- wheel, 0.34.1
- wincertstore, 0.2
- wtforms, 2.2.1

To activate the environment and start running the software please use:

```
conda activate jacky_env
```


### Installing
The files to run this software are in the repository https://github.com/celiaccb/Software-Development-Group-Project-2020
Clone the repository using git:
```
git clone https://github.com/celiaccb/Software-Development-Group-Project-2020

```
The following folder should be donwloaded in your computer:
 * "Generating phosphosite and inhibitor files" and "Kinase dataframe", these two folder contains explanations and the codes with the steps to extract the data from external websites and store it in csv files.  
 * "database schema", which contains a the script to create the schema of the database and a diagram withg the tables.
 * "population_db", which contains the csv files with the data, the final database with all the information transferred and the code to transfer the data from csv files to the database.
 * "application", this folder contains all the files to deploy the software in a live system. The structure of the folder should follow the organisation below:
```
.
├── .ebextensions/
│ └── environmentvariables.config
├── downloads/
├── static/
│ └── aesthetic.css
│ └── Chromosomes1.jpg
│ └── Chromosomes12.jpg
│ └── file_structure.jpg
│ └── JACKY.jpg
│ └── RKA_table.jpg
├── templates/
│ └── datanalysis.html
│ └── homepage.html
│ └── homepageinh.html
│ └── homepagekin.html
│ └── homepagephos.html
│ └── inhibitorhits.html
│ └── inhibitorpage.html
│ └── kinasehits.html
│ └── kinasepage.html
│ └── phoshits.html
│ └── phoshitsgen.html
│ └── phosphositepage.html
│ └── zphoshitsgen.html
│ └── zphoshitsgenalt.html
├── uploads/
├── application.py
├── c__sqlite_final_database.db
├── requirements.txt
```
### Executing the software

Once the virtual environment is activated, then move to the folder where the software has been installed and type:
```
python application.py runserver -d
```
The following messages will be displayed in the screen

```
2020-02-09 09:48:39,506 INFO sqlalchemy.engine.base.Engine SELECT CAST('test plain returns' AS VARCHAR(60)) AS anon_1
2020-02-09 09:48:39,508 INFO sqlalchemy.engine.base.Engine ()
2020-02-09 09:48:39,512 INFO sqlalchemy.engine.base.Engine SELECT CAST('test unicode returns' AS VARCHAR(60)) AS anon_1
2020-02-09 09:48:39,514 INFO sqlalchemy.engine.base.Engine ()
2020-02-09 09:48:39,516 INFO sqlalchemy.engine.base.Engine PRAGMA main.table_info("human_kinases")
2020-02-09 09:48:39,516 INFO sqlalchemy.engine.base.Engine ()
2020-02-09 09:48:39,521 INFO sqlalchemy.engine.base.Engine PRAGMA main.table_info("phosphosites")
2020-02-09 09:48:39,522 INFO sqlalchemy.engine.base.Engine ()
2020-02-09 09:48:39,524 INFO sqlalchemy.engine.base.Engine PRAGMA main.table_info("inhibitors")
2020-02-09 09:48:39,525 INFO sqlalchemy.engine.base.Engine ()
2020-02-09 09:48:39,527 INFO sqlalchemy.engine.base.Engine PRAGMA main.table_info("kinases_phosphosites")
2020-02-09 09:48:39,527 INFO sqlalchemy.engine.base.Engine ()
2020-02-09 09:48:39,530 INFO sqlalchemy.engine.base.Engine PRAGMA main.table_info("phosphosites_diseases")
2020-02-09 09:48:39,531 INFO sqlalchemy.engine.base.Engine ()
2020-02-09 09:48:39,532 INFO sqlalchemy.engine.base.Engine PRAGMA main.table_info("inhib_kin")
2020-02-09 09:48:39,533 INFO sqlalchemy.engine.base.Engine ()
 * Serving Flask app "application" (lazy loading)
 * Environment: production
   WARNING: This is a development server. Do not use it in a production deployment.
   Use a production WSGI server instead.
 * Debug mode: off
 * Running on http://127.0.0.1:5000/ (Press CTRL+C to quit)
```

The http address will give you access to the software in a web framework from your local terminal. 

If you want to test the scripts to extract data from the external sources on kinases and transfer those data to a database, follow the instructions and explantions from the corresponding folders in section **Installing**.

### Built With
1. **Flask**  - The web framework used
2. **SQLite** - The Relational Database Management System
3. **KSEA (Kinase-Substrate Enrichment Analysis)** - method to calculate the relative activity of the kinases based on the phosphoproteomics data
4. **Amazon Web Services** - to deplot the software on a live system

### Version History
* 0.1, Initial release (Feburary 14th 2020)

### Authors and contacts
The developers of this software were the students of the Msc in Bioinformatics at Queen Mary Univeristy of London:

* **Yutang Chen**, y.chen@se19.qmul.ac.uk
* **Celia De Los Angeles Colomina Basanta**, c.colominabasanta@se16.qmul.ac.uk
* **Anna Dearman**, a.dearman@se19.qmul.ac.uk
* **Juan Ledesma Moreno**, j.ledesmamoreno@se19.qmul.ac.uk
* **Katie Skinner**, k.skinner@se19.qmul.ac.uk

### License
This project is totally free. Please cite our repository if you are using JACKY as a tool for your research in case of publication. 

### Acknowledgments
JACKY team would like to thank all the people who have helped and supported us with their comments and help during the development of the software.

Juan Ledesma Moreno wants to thank to the National Infection Service from Public Health England (PHE) for the financial support, to Dr Jean Lutamyo Mbisa, Dr Eileen Gallagher, Dr David Williams and Dr Ana Penedos from PHE for their unconditional support, and, in special to my partner Diana for her love.  
