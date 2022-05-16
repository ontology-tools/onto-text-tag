# Copyright 2018 Google LLC
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
from logging import error
import pyhornedowl
from flask import Flask, request, redirect, url_for, session, Response, stream_with_context, current_app
from flask.templating import render_template
# from ontotagtext import ExtractorComponent
from ontotagtext import MultiExtractorComponent
from ontotagtext import PREFIXES
from ontotagtext import RDFSLABEL
import spacy
import en_core_web_sm
from Bio import Entrez
import requests
from urllib.request import urlopen
import json
import re
import io
import traceback

import pprint

#circle graph generator
import hv_generate
from hv_generate import hv_generator

from build_ontotermentions import buildit

from bokeh.embed import json_item 
from jinja2 import Template


from bokeh.resources import CDN

from bokeh.plotting import figure
from bokeh.sampledata.iris import flowers

import os
import csv #for writing test_terms.csv (once only)
import pickle
import shelve

#OGER:
import oger
from oger.ctrl.router import Router, PipelineServer
conf = Router(termlist_path='static/test_terms.tsv')
pl = PipelineServer(conf)

import inflect

development = False

pp = pprint.PrettyPrinter(depth=4)

app = Flask(__name__)


app.config.from_object('config')
#idName = "ID"
# python -m spacy download en_core_web_md
# or: en_core_web_md or en_core_web_lg
# nlp = spacy.load('en_core_web_md')
nlp = en_core_web_sm.load()

location = f"https://raw.githubusercontent.com/addicto-org/addiction-ontology/master/addicto-merged.owx"
location2 = f"https://raw.githubusercontent.com/HumanBehaviourChangeProject/ontologies/master/Upper%20Level%20BCIO/bcio-merged.owx"

print("Fetching release file from", location)
ontol1 = pyhornedowl.open_ontology(urlopen(location).read().decode('utf-8'))
print("Fetching release file from", location2)
ontol2 = pyhornedowl.open_ontology(urlopen(location2).read().decode('utf-8'))

# pickle_in = open("allAbstracts.pkl","rb")
# abstract_associations = pickle.load(pickle_in)
development = (os.environ.get("FLASK_ENV")=='development')
#if development:
abstract_ass_db = shelve.open('static/allAbstracts.db', flag='r', writeback=False)
print("loaded abstract associations db")
all_titles_db = shelve.open('static/allTitles.db', flag='r', writeback=False)
print("loaded abstract titles db")
all_dates_db = shelve.open('static/allDates.db', flag='r', writeback=False)
print("loaded abstract dates db")
all_authors_db = shelve.open('static/allAuthorAffils.db', flag='r', writeback=False)
print("loaded abstract authors db")

for prefix in PREFIXES:
    ontol1.add_prefix_mapping(prefix[0], prefix[1])
    ontol2.add_prefix_mapping(prefix[0], prefix[1])

# combined test
# populated with {"label1": name1, ontofile1}, {"label2": ...}
ontoDict = {
    "ontologies": [       
        {
            "label": "BCIO",
            "name": "BCIO",
            "ontology": ontol2
        },
        
        {
            "label": "AddictO",
             "name": "AddictO",
            "ontology": ontol1
        },
    ]
}
    
def get_all_descendents(id_list):   
    descendent_ids = []

    #todo: refactor below:

    for entry in id_list:
        entryIri = ontol1.get_iri_for_id(entry.replace("_", ":"))
        if entryIri:
            descs = pyhornedowl.get_descendants(ontol1, entryIri)
            if len(descs) > 0:
                for d in descs:
                    try:
                        add_id = ontol1.get_id_for_iri(d).replace(":", "_")
                        descendent_ids.append(add_id.replace("_", ":"))
                    except:
                        print("error when getting descendents of ",d)
                
    for entry in id_list:
        entryIri = ontol2.get_iri_for_id(entry.replace("_", ":"))
        if entryIri:
            descs = pyhornedowl.get_descendants(ontol2, entryIri)
            if len(descs) > 0:
                for d in descs:
                    try:
                        add_id = ontol1.get_id_for_iri(d).replace(":", "_")
                        descendent_ids.append(add_id.replace("_", ":"))
                    except:
                        print("error getting descendents of ",d)
                
    if len(descendent_ids) == 0:
        return id_list
    else:
        #remove duplicates from descendent_ids:
        return_list = id_list + list(set(descendent_ids))
        return return_list
        # return descendent_ids #test only descendents

onto_extractor3 = MultiExtractorComponent(
    nlp,
    ontoDict=ontoDict
    )
nlp.add_pipe(onto_extractor3, after="ner")


# Interaction with PubMed: get detailed results for a list of IDs
def fetch_details(id_list):
   results = ""
   ids = ','.join(id_list)
   Entrez.email = 'janna.hastings@ucl.ac.uk'
   handle = Entrez.efetch(db='pubmed',
                          retmode='xml',
                          id=ids)
   results = Entrez.read(handle)
   return results

# parse the title, authors and date published
# try return separate values for year, day, month, AuthourList and ArticleTitle


def get_article_details(result):
    # articleDetails = "DATE;TITLE;AUTHORS"
    articleDetails = ""
    dayCompleted = ""
    monthCompleted = ""
    yearCompleted = ""
    for detail in result:
        if 'MedlineCitation' in detail:
            if 'DateCompleted' in detail['MedlineCitation']:  
                dayCompleted = str(
                    detail['MedlineCitation']['DateCompleted']['Day'])
                monthCompleted = str(
                    detail['MedlineCitation']['DateCompleted']['Month'])
                yearCompleted = str(
                    detail['MedlineCitation']['DateCompleted']['Year'])
            else:
                yearCompleted = ""
                monthCompleted = ""
                # todo: we still end up with // here if no data returned..
                dayCompleted = ""

            if 'ArticleTitle' in detail['MedlineCitation']['Article']:
                articleDetails = dayCompleted + "/" + monthCompleted + "/" + \
                    yearCompleted + ";" + \
                    str(detail['MedlineCitation']['Article']['ArticleTitle'])
            else:
                # print(f"ArticleTitle not found")
                articleDetails = "/ / ;"

            # this works, need to refine though
            if 'AuthorList' in detail['MedlineCitation']['Article']:
                articleDetails += ";Authors: "
                for s in range(len(detail['MedlineCitation']['Article']['AuthorList'])):
                    articleDetails += detail['MedlineCitation']['Article']['AuthorList'][s]['LastName']
                    if(s == (len(detail['MedlineCitation']['Article']['AuthorList'])-1)):
                        articleDetails += "."
                    else:
                        articleDetails += ", "
            else:
                # print(f"AuthorList not found")
                articleDetails += ";Authors: "
    return articleDetails

# Parse the PubMed result to get the abstract text if it is there


def get_abstract_text(result):
    abstractText = None
    fixed_abstractText = ""
    fixed = []
    for detail in result:
        if 'MedlineCitation' in detail:
            if 'Article' in detail['MedlineCitation']:
                if 'Abstract' in detail['MedlineCitation']['Article']:
                    # print("Article is: ", detail['MedlineCitation']['Article'])
                    if 'AbstractText' in detail['MedlineCitation']['Article']['Abstract']:
                        abstractText = str(detail['MedlineCitation']['Article']['Abstract']['AbstractText'])
                        # fixed_abstractText = abstractText
                        if("StringElement" in abstractText):
                            fixed = re.findall(r'StringElement\((.+?)attributes',abstractText)
                            fixed_abstractText = "".join(fixed)
                        else:
                            fixed_abstractText = abstractText.strip('[]') # remove "[]"
                        
    return fixed_abstractText

def get_ids(ontol_list):
    # print("get_ids running here")
    checklist = []
    for ontol in [ontol1,ontol2]:
        for classIri in ontol.get_classes():
            # print("for classIri running")        
            classId = ontol.get_id_for_iri(classIri)
            label = ontol.get_annotation(classIri, RDFSLABEL)
            # label = ontol.get_annotation(classId, RDFSLABEL)
            if classId:
                # print("got classId and labels") 
                checklist.append(classId + "|"+ label)                   
                # print(classId)
                # print(label)
                    
    return checklist

#strip html: 
from io import StringIO
from html.parser import HTMLParser
class MLStripper(HTMLParser):
    def __init__(self):
        super().__init__()
        self.reset()
        self.strict = False
        self.convert_charrefs= True
        self.text = StringIO()
    def handle_data(self, d):
        self.text.write(d)
    def get_data(self):
        return self.text.getvalue()

def strip_tags(html):
    s = MLStripper()
    s.feed(html)
    returnString = s.get_data()
    returnString = re.sub(r'\\x..', '', returnString) #from StackOverflow, works
    returnString = re.sub(r'\\u....', '', returnString) #
    return returnString


@app.route('/build')
@app.route('/build')
def build():
    try:
        # todo: use RQ scheduler instead of below, need better system to indicate success
        # todo: relative path to addiction-ontology
        os.chdir(os.path.dirname(os.path.realpath(__file__)))
        os.system("python -W ignore build_ontotermentions.py --path /home/tom/addiction-ontology")
        os.chdir(os.path.dirname(os.path.realpath(__file__)))
        os.system("touch success.txt")
        return "Build Successful"
    except: 
        os.chdir(os.path.dirname(os.path.realpath(__file__)))
        os.system("touch failed.txt")
        return "Build call FAILED"
    


# Pages for the app
@app.route('/')
@app.route('/home')
def home():
    development = (os.environ.get("FLASK_ENV")=='development')
    return render_template('index.html', development = development)



@app.route('/associations', methods=['GET', 'POST'])
def associations():
    # todo: filter by get_ids
    ontologies = ontoDict["ontologies"]
    # print(ontologies)
    labels = get_ids(ontologies) 
    label_list={'labels': labels}
    json.dumps(label_list)
    return render_template('associations.html', label_list = label_list)



@app.route('/visualise_associations', methods=['POST'])
def visualise_associations():  
    ontology_id_list = json.loads(request.form.get('ontology_id_list')) 
    # print("got ontology id list: ", ontology_id_list)
    include_descendents = request.form.get("include_descendent_classes")
    session['saved_ontology_id_list'] = ontology_id_list
    #todo: saved_ontology_id_list not being saved to session..??
    print("saved_ontology_id_list should be: ", session['saved_ontology_id_list'])
    session['get_descendents'] = include_descendents
    iframe = url_for('chordout')
    return render_template("chord.html", iframe=iframe)


#generate chord plot as html, to insert into iframe:
@app.route('/chordout')
def chordout():  
    if 'saved_ontology_id_list' in session:
        saved_ontology_id_list = session['saved_ontology_id_list']
        get_descendents = session['get_descendents']
        print("saved_ontology_id_list in /chordout ", saved_ontology_id_list)
        if get_descendents == "true":
            ontology_id_list = get_all_descendents(saved_ontology_id_list)
            # print("looking for: ", ontology_id_list)
            html = json.loads(hv_generator(ontology_id_list))
            return render_template('chordout2.html', html=html)
        else:
            html = json.loads(hv_generator(saved_ontology_id_list))
            return render_template('chordout2.html', html=html)
    # return("No ontology ID List error") #todo: error message in index.html here

@app.route('/visualise_similarities', methods=['POST'])
def visualise_similarities():  
    ontology_id_list = json.loads(request.form.get('ontology_id_list')) 
    include_descendent_classes = request.form.get('include_descendent_classes')
    print("time for a similarity visual!")
    # hv_generator(ontology_id_list) #todo: something different here?
    # return ( json.dumps({"message":"Success"}), 200 )
    return render_template('similarity.html')

@app.route('/pubmed', methods=['POST'])
def pubmed():
    # print("pubmed")
    id = request.form.get('pubmed_id')
    if id:
        fixed_id = id.strip()
    else:
        return render_template('index.html',
            error_msg=f"Error: cannot retrieve PubMed for empty ID.",
            development=development)
    dateA = ""
    titleA = ""
    authorsA = "" 

#This try/except fails on any issue with missing pubmed ID. Deal with each case individually. 
    try:
        if fixed_id is not None and fixed_id in abstract_ass_db.keys():
            fixed_abstractText = abstract_ass_db[fixed_id]
            articleDetails = id
            if all_dates_db[fixed_id] is not None:
                dateA = all_dates_db[fixed_id]
            if all_titles_db[fixed_id] is not None: 
                titleA = all_titles_db[fixed_id]
            if all_authors_db[fixed_id] is not None: 
                authorsA = all_authors_db[fixed_id]
        else:
            #not in pre-downloaded abstracts, get pubmed info direct from pubmed here:
            #todo: refactor this duplicate block
            try:
                detailResults = fetch_details([id])
                # print("got fetched_details: ", detailResults)
                articleDetails = id
                for result in detailResults:
                    resultDetail = detailResults[result]
                    for detail in resultDetail:
                        if 'MedlineCitation' in detail:
                            PMID = str(detail['MedlineCitation']['PMID'])
                            if 'Article' in detail['MedlineCitation']:
                                #title: 
                                if 'ArticleTitle' in detail['MedlineCitation']['Article']:
                                    articleTitle = str(detail['MedlineCitation']['Article']['ArticleTitle'])
                                    titleA = articleTitle
                                #date: 
                                if 'ArticleDate' in detail['MedlineCitation']['Article']:
                                    try: 
                                        year = detail['MedlineCitation']['Article']['ArticleDate'][0]['Year']
                                        month = detail['MedlineCitation']['Article']['ArticleDate'][0]['Month']
                                        day = detail['MedlineCitation']['Article']['ArticleDate'][0]['Day']
                                    except: 
                                        year = ""
                                        month = ""
                                        day = ""
                                    articleDate = day + "/" + month + "/" + year
                                    dateA = articleDate
                                if 'Abstract' in detail['MedlineCitation']['Article']:
                                    if 'AbstractText' in detail['MedlineCitation']['Article']['Abstract']: 
                                        res = str(detail['MedlineCitation']['Article']['Abstract']['AbstractText']).strip('][').split(', ')
                                        abstractText = ""
                                        for r in res: 
                                            if "StringElement" in r: 
                                                r = r.replace("StringElement('", "")
                                                r = r.replace("StringElement(\"", "")
                                            if not "})" in r:         
                                                if not "attributes" in r: 
                                                    abstractText = abstractText.strip() + " " + r
                                        abstractText = strip_tags(abstractText) #strip html tags
                                        fixed_abstractText = abstractText
                                        
                                authorDetails = ""
                                if 'AuthorList' in detail['MedlineCitation']['Article']:
                                    try: 
                                        for s in range(len(detail['MedlineCitation']['Article']['AuthorList'])):
                                            authorDetails += detail['MedlineCitation']['Article']['AuthorList'][s]['LastName']
                                            if(s == (len(detail['MedlineCitation']['Article']['AuthorList'])-1)):
                                                authorDetails += "."
                                            else:
                                                authorDetails += ", "
                                        authorsA = authorDetails
                                    except:
                                        authorDetails += ""
                                        AuthorsA = authorDetails
                                        
                                else:
                                    authorDetails += ""
                                    authorsA = authorDetails
                #returning details from fetch_details here:
                # print("should get tag here") 
                # # try:
                # return render_template('tag.html', data={
                #                         "inputDetails": articleDetails,
                # "inputText": fixed_abstractText, "dateDetails": dateA,
                # "titleDetails": titleA, "authorsDetails": authorsA})

                # r = requests.post(url_for("tag", _external=True), data={
                #                         "inputDetails": articleDetails,
                # "inputText": fixed_abstractText, "dateDetails": dateA,
                # "titleDetails": titleA, "authorsDetails": authorsA})
                # print("got result: ", r.text, r.status_code, r.header)
                # return r.text, r.status_code, r.headers.items() #this line is giving issue..
                # except: 
                #     return render_template('index.html',
                #                        error_msg=f"Error tagging {id}",
                #                        development=development)
            except Exception as exe:
                #no pubmed ID found, return error message
                print(exe)
                traceback.print_exc()
                return render_template('index.html',
                                       error_msg=f"The PubMed ID {id} was not indexed  - try pasting in the abstract text instead",)
        # get from all_abstracts db:
        # todo: below if statements causing errors - is this still needed? 
        # if all_dates_db[fixed_id] is not None:
        #     dateA = all_dates_db[fixed_id]
        # if all_titles_db[fixed_id] is not None: 
        #     titleA = all_titles_db[fixed_id]
        # if all_authors_db[fixed_id] is not None: 
        #     authorsA = all_authors_db[fixed_id]
        r = requests.post(url_for("tag", _external=True), data={
                                    "inputDetails": articleDetails,
            "inputText": fixed_abstractText, "dateDetails": dateA,
            "titleDetails": titleA, "authorsDetails": authorsA})
        print("pubmed return")
        return r.text, r.status_code, r.headers.items()
    except Exception as e:
        print(e)
        traceback.print_exc()
        try:
            #get pubmed info direct from pubmed here
            fetched_details = fetch_details([id])
            #todo: refactor this duplicate block
            # print("got fetched_details: ", fetched_details)
            articleDetails = id
            for result in detailResults:
                    resultDetail = detailResults[result]
                    for detail in resultDetail:
                        if 'MedlineCitation' in detail:
                            PMID = str(detail['MedlineCitation']['PMID'])
                            if 'Article' in detail['MedlineCitation']:
                                #title: 
                                if 'ArticleTitle' in detail['MedlineCitation']['Article']:
                                    articleTitle = str(detail['MedlineCitation']['Article']['ArticleTitle'])
                                    titleA = articleTitle
                                #date: 
                                if 'ArticleDate' in detail['MedlineCitation']['Article']:
                                    try: 
                                        year = detail['MedlineCitation']['Article']['ArticleDate'][0]['Year']
                                        month = detail['MedlineCitation']['Article']['ArticleDate'][0]['Month']
                                        day = detail['MedlineCitation']['Article']['ArticleDate'][0]['Day']
                                    except: 
                                        year = ""
                                        month = ""
                                        day = ""
                                    articleDate = day + "/" + month + "/" + year
                                    dateA = articleDate
                                if 'Abstract' in detail['MedlineCitation']['Article']:
                                    if 'AbstractText' in detail['MedlineCitation']['Article']['Abstract']: 
                                        res = str(detail['MedlineCitation']['Article']['Abstract']['AbstractText']).strip('][').split(', ')
                                        abstractText = ""
                                        for r in res: 
                                            if "StringElement" in r: 
                                                r = r.replace("StringElement('", "")
                                                r = r.replace("StringElement(\"", "")
                                            if not "})" in r:         
                                                if not "attributes" in r: 
                                                    abstractText = abstractText.strip() + " " + r
                                        abstractText = strip_tags(abstractText) #strip html tags
                                        fixed_abstractText = abstractText
                                        
                                authorDetails = ""
                                if 'AuthorList' in detail['MedlineCitation']['Article']:
                                    try: 
                                        for s in range(len(detail['MedlineCitation']['Article']['AuthorList'])):
                                            authorDetails += detail['MedlineCitation']['Article']['AuthorList'][s]['LastName']
                                            if(s == (len(detail['MedlineCitation']['Article']['AuthorList'])-1)):
                                                authorDetails += "."
                                            else:
                                                authorDetails += ", "
                                        authorsA = authorDetails
                                    except:
                                        authorDetails += ""
                                        AuthorsA = authorDetails
                                        
                                else:
                                    authorDetails += ""
                                    authorsA = authorDetails

            #returning details from fetch_details here: 
            print("should get tag here2")
            # try:
            r = requests.post(url_for("tag", _external=True), data={
                                    "inputDetails": articleDetails,
            "inputText": fixed_abstractText, "dateDetails": dateA,
            "titleDetails": titleA, "authorsDetails": authorsA})
            print("got result: ", r.text, r.status_code, r.header)
            # except: 
            #     return render_template('index.html',
            #                         error_msg=f"Error tagging {id}",
            #                         development=development)
            print("pubmed return2")
            return r.text, r.status_code, r.headers.items()
        except Exception as exe:
            #no pubmed found, return error message
            print(exe)
            traceback.print_exc()
            return render_template('index.html',
                                   error_msg=f"This PubMed ID {id} was not indexed  - try pasting in the abstract text instead")

    
@ app.route('/tag', methods=['POST'])
def tag():
    print("/tag")
    # development = (os.environ.get("FLASK_ENV")=='development')
    development="development"
    text = request.form['inputText']
    # print("got text", text)
    details = request.form.get('inputDetails') #pmid
    # print("should have id: ", details)
    date = request.form.get('dateDetails')
    title = request.form.get('titleDetails')
    authors = request.form.get('authorsDetails')
    if details is None:
        details = ""
    if date is None:
        date = ""
    if title is None:
        title = ""
    if authors is None:
        authors = ""

    # process the text
    tag_results = []
    
    # NOTE: build_terms is for re-building test_terms.tsv
    # to do this, set build_terms to True and delete the test_terms.tsv, and test_terms.tsv.pickle (generated by OGER)
    # only necessary if the ontotermentions.csv file has been updated. 
    # todo: re-factor to build terms when ontotermentions.csv is updated

    #new idea: if test_terms.tsv is not in the /static directory, then build it.
    build_terms = False

    if not os.path.isfile(os.path.join(app.root_path, 'static/test_terms.tsv.pickle')):
        build_terms = True

    engine = inflect.engine()

           
    if build_terms: 
        print("Building terms, please wait...")
        # stop words, don't try to match these
        stopwords = nlp.Defaults.stop_words
        stopwords.add("ands")
        stopwords.add("ends")
        stopwords.add("ci")
        #test build test_terms.tsv from onto_extractor3:
        mydict = []
        for f in onto_extractor3.terms:
            l = onto_extractor3.get_label(f)
            if l is not None:
                term=onto_extractor3.get_term(l['id'])
                if term:                        
                    ont = term['id'][0:term['id'].index(":")]
                    #checking different ontologies:
                    # if ont == "BCIO":
                    #     print("BCIO")
                    # else: 
                    #     print(ont)
                    if term['id'] == "BCIO:010055":
                        continue
                    else:
                        # order: 'a', 'ont', 'id', 'alt_name', 'name', 'definition'
                        sing = {'a': '', 'ont': ont, 'id': term['id'], 'alt_name': term['name'], 'name': term['name'], 'definition': term['definition']}
                        mydict.append(sing)
                        
                        #plurals:                    
                        try:
                            plural = engine.plural(term['name'].strip())
                            plur = {'a': '', 'ont': ont, 'id': term['id'], 'alt_name': plural, 'name': term['name'], 'definition': term['definition']}
                            mydict.append(plur)
                        except: 
                            # print("Problem getting plural of ", term['name'].strip())
                            continue        
                else:
                    print("No term for label: ", l)
        # use mydict id to check for synonyms and plurals of synonyms:    
        
        for ontol in onto_extractor3.ontols:
            for termid in ontol.get_classes():                   
                SYN = "http://purl.obolibrary.org/obo/IAO_0000118"
                DEFINITION = "http://purl.obolibrary.org/obo/IAO_0000115"
                synonyms = ontol.get_annotations(termid, SYN)
                label = ontol.get_annotation(termid, RDFSLABEL)
                definition = ontol.get_annotation(termid, DEFINITION)
                termshortid = ontol.get_id_for_iri(termid)
                if termshortid: 
                    # termshortid from last /:
                    ontol_namespace = termshortid[termshortid.rfind("/")+1:].strip()
                    #ontol_namespace string until :
                    ontol_namespace = ontol_namespace[0:ontol_namespace.index(":")]
                else:
                    ontol_namespace = "" 
                
                for s in synonyms:
                    # print("adding SYNONYM: ", s)
                    # print("got ontol_namespace: ", ontol_namespace)
                    if s.strip().lower() not in stopwords:                            
                        syn1 = {'a': '', 'ont': ontol_namespace, 'id': termshortid, 'alt_name': s, 'name': label, 'definition': definition}
                        mydict.append(syn1)
                        try:
                            plural2 = engine.plural(s.strip())
                            plur2 = {'a': '', 'ont': ontol_namespace, 'id': termshortid, 'alt_name': plural2, 'name': label, 'definition': definition}
                            mydict.append(plur2)
                        except:
                            print("Problem getting plural of ",s)
                            pass

                          
        filename = 'static/test_terms.tsv'
        fields = ['a', 'ont', 'id', 'alt_name', 'name', 'definition'] 
        with open(filename, 'w') as tsvfile: 
            # creating a csv dict writer object 
            writer = csv.DictWriter(tsvfile, delimiter='\t', fieldnames=fields) 
            writer.writerows(mydict) 
        print("done creating test_terms.tsv")
    
    #file_path = os.path.join(current_app.root_path,'static')
    # file_path = url_for('static', filename = 'text.txt')

    #with open(file_path + '/text.txt', "w") as textfile:
    #    textfile.write(text)
    #    textfile.close()

    textfile = io.StringIO(text)


    # coll_pmid = []    
    # coll_pmid.append(idName) #idName is the pubmed id
    # pass
    # print("coll_pmid = ", coll_pmid)
    #load_file = file_path + '/text.txt'
    coll = pl.load_one(textfile, fmt='txt')#, iter_mode='document')
    pl.process(coll)
    
    
    # note: the entity.names below are just to fit in with OGER's un-related column naming. 
    for entity in coll[0].iter_entities():
        span_text = entity.text.strip()
        ontol_id = entity.cid.strip() 
        link_prefix_http = ""
        if 'BCIO' in ontol_id.strip():
            link_prefix_http = "https://bciovocab.org/"
        else: 
            link_prefix_http = "http://addictovocab.org/"
        ontol_label = entity.pref.strip()
        ontol_def = entity.type.strip()
        ontol_namespace = entity.db.strip()
        tag_results.append({"ontol_id": ontol_id,
                                "span_text": span_text,
                                "ontol_label": ontol_label,
                                "ontol_def": ontol_def,
                                "ontol_namespace": ontol_namespace,
                                "ontol_link": link_prefix_http+ontol_id,
                                "match_index": ontol_id})
    print("should render index.html here (/tag)")
    return render_template('index.html',
                           text=text,
                           details=details,
                           date=date,
                           title=title,
                           authors=authors,
                           id='',
                           tag_results=tag_results,
                           development=development)



if __name__ == "__main__":
    app.run(host='0.0.0.0')
