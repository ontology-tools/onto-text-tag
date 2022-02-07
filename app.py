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

import pprint

#circle graph generator
import hv_generate
from hv_generate import hv_generator

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
idName = "ID"
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
abstract_ass_db = shelve.open('allAbstracts.db')
print("loaded abstract associations db")
all_titles_db = shelve.open('allTitles.db')
print("loaded abstract titles db")
all_dates_db = shelve.open('allDates.db')
print("loaded abstract dates db")
all_authors_db = shelve.open('allAuthors.db')
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
    # ids = ','.join(id_list)
    # Entrez.email = 'janna.hastings@ucl.ac.uk'
    # handle = Entrez.efetch(db='pubmed',
    #                        retmode='xml',
    #                        id=ids)
    # results = Entrez.read(handle)
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



# Pages for the app
@app.route('/')
@app.route('/home')
def home():
    development = (os.environ.get("FLASK_ENV")=='development')
    return render_template('index.html', development = development)



@app.route('/associations', methods=['GET', 'POST'])
def associations():
    ontologies = ontoDict["ontologies"]
    labels = get_ids(ontologies) 
    label_list={'labels': labels}
    json.dumps(label_list)
    return render_template('associations.html', label_list = label_list)



@app.route('/visualise_associations', methods=['POST'])
def visualise_associations():  
    ontology_id_list = json.loads(request.form.get('ontology_id_list')) 
    include_descendents = request.form.get("include_descendent_classes")
    session['saved_ontology_id_list'] = ontology_id_list
    session['get_descendents'] = include_descendents
    iframe = url_for('chordout')
    return render_template("chord.html", iframe=iframe)


#generate chord plot as html, to insert into iframe:
@app.route('/chordout')
def chordout():    
    if 'saved_ontology_id_list' in session:
        saved_ontology_id_list = session['saved_ontology_id_list']
        get_descendents = session['get_descendents']
        if get_descendents == "true":
            html = json.loads(hv_generator(saved_ontology_id_list, True))
        else:
            html = json.loads(hv_generator(saved_ontology_id_list, False))    
    return render_template('chordout2.html', html=html)

@app.route('/visualise_similarities', methods=['POST'])
def visualise_similarities():  
    ontology_id_list = json.loads(request.form.get('ontology_id_list')) 
    include_descendent_classes = request.form.get('include_descendent_classes')
    print("time for a similarity visual!")
    # hv_generator(ontology_id_list) #todo: something different here?
    # return ( json.dumps({"message":"Success"}), 200 )
    return render_template('similarity.html')

@app.route('/pubmed', methods=['POST', 'GET'])
def pubmed():
    development = (os.environ.get("FLASK_ENV")=='development')
    id = request.form.get('pubmed_id')
    fixed_id = ""
    if id:
        fixed_id = id.strip()
    global idName
    articleDetails = ""
    idName = ""
    dateA = ""
    titleA = ""
    authorsA = "" 

    try: 
        idName = fixed_id
        fixed_abstractText = abstract_ass_db[fixed_id]
        articleDetails = id

        if all_dates_db[fixed_id] is not None: 
            dateA = all_dates_db[fixed_id]
        if all_titles_db[fixed_id] is not None: 
            titleA = all_titles_db[fixed_id]
        #todo: authors incorrect - showing country instead?!
        if all_authors_db[fixed_id] is not None: 
            authorsA = all_authors_db[fixed_id]
    # if fixed_id in abstract_associations:
    #     one_abstract = abstract_associations[fixed_id]         
    #     if("StringElement" in one_abstract):
    #         fixed = re.findall(r'StringElement\((.+?)attributes',one_abstract)
    #         fixed_abstractText = "".join(fixed)
    #     else:
    #         fixed_abstractText = one_abstract.strip('[]') # remove "[]"
    #     print(fixed_abstractText)
        r = requests.post(url_for("tag", _external=True), data={
                                    "inputDetails": articleDetails, "inputText": fixed_abstractText, "dateDetails": dateA, "titleDetails": titleA, "authorsDetails": authorsA})
        return r.text, r.status_code, r.headers.items()
    except: 
        return render_template('index.html', error_msg=f"This PubMed ID {id} was not indexed  - try pasting in the abstract text instead", development = development)

    # else:     
    #     return render_template('index.html', error_msg=f"No abstract found for PubMed ID {id}", development = development)

    #old method using entrez:
    # if id:
    #     idName = f"{id}"
    #     try:
    #         results = fetch_details([id])
    #         for result in results:
    #             resultDetail = results[result]
    #             abstractText = get_abstract_text(resultDetail)
    #             articleDetails = get_article_details(resultDetail)
    #             try:
    #                 dateA, titleA, authorsA = articleDetails.split(';')
    #             except:
    #                 pass
    #             if abstractText:
    #                 r = requests.post(url_for("tag", _external=True), data={
    #                                   "inputDetails": articleDetails, "inputText": abstractText, "dateDetails": dateA, "titleDetails": titleA, "authorsDetails": authorsA})
    #                 return r.text, r.status_code, r.headers.items()
    #     except Exception as err:  # 400 bad request handling, also if no internet connection
    #         print(err)
    # return render_template('index.html', error_msg=f"No abstract found for PubMed ID {id}", development = development)


@ app.route('/tag', methods=['POST'])
def tag():    
    development = (os.environ.get("FLASK_ENV")=='development')

    text = request.form['inputText']
    # print("got text", text)
    details = request.form.get('inputDetails') #pmid
    print("should have id: ", details)
    date = request.form.get('dateDetails')
    title = request.form.get('titleDetails')
    authors = request.form.get('authorsDetails')
    if details == None:
        details = ""
    if date == None:
        date = ""
    if title == None:
        title = ""
    if authors == None:
        authors = ""

    # process the text
    tag_results = []
    
   

    # NOTE: build_terms is for re-building test_terms.tsv
    # to do this, set build_terms to True and delete the test_terms.tsv, and test_terms.tsv.pickle (generated by OGER)
    # only necessary if the ontotermentions.csv file has been updated. 
    # todo: re-factor to build terms when ontotermentions.csv is updated
    build_terms = False

    engine = inflect.engine()

           
    if build_terms:
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
                            print("Problem getting plural of ", term['name'].strip())
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
                    print("adding SYNONYM: ", s)
                    print("got ontol_namespace: ", ontol_namespace)
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
    

    # replaced nlp with OGER:
    # fields list for entity is here: https://github.com/OntoGene/OGER/blob/f23cf9bec70ba51f85605f26f3de2df72f7c4d5a/oger/doc/document.py
    # if (len(idName.strip()) > 0): #todo: fix this or it won't work when text is in the text field
    #     coll_pmid = []    
    #     coll_pmid.append(idName) #idName is the pubmed id
        
    #     # print("coll_pmid = ", coll_pmid)
    #     coll = pl.load_one(coll_pmid, fmt='pubmed')
    #     pl.process(coll)
        
    #     # note: the entity.names below are just to fit in with OGER's un-related column naming. 
    #     for entity in coll[0].iter_entities():
    #         span_text = entity.text.strip()
    #         ontol_id = entity.cid.strip() 
    #         ontol_label = entity.pref.strip()
    #         ontol_def = entity.type.strip()
    #         ontol_namespace = entity.db.strip()
    #         tag_results.append({"ontol_id": ontol_id,
    #                                 "span_text": span_text,
    #                                 "ontol_label": ontol_label,
    #                                 "ontol_def": ontol_def,
    #                                 "ontol_namespace": ontol_namespace,
    #                                 "ontol_link": "http://addictovocab.org/"+ontol_id,
    #                                 "match_index": ontol_id})
    # else:
        #todo: can we load text directly without saving? 
    file_path = os.path.join(current_app.root_path,'static') 
    # file_path = url_for('static', filename = 'text.txt')
    with open(file_path + 'text.txt', "w") as textfile:
        textfile.write(text)
        textfile.close()


    # coll_pmid = []    
    # coll_pmid.append(idName) #idName is the pubmed id
    # pass
    # print("coll_pmid = ", coll_pmid)
    load_file = file_path + 'text.txt'
    coll = pl.load_one(load_file, fmt='txt')#, iter_mode='document')
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
    # else: #nlp/spacy
    #     doc3 = nlp(text)
    #     # get ontology IDs identified
    #     for token in doc3:
    #         if token._.is_ontol_term:
    #             term=onto_extractor3.get_term(token._.ontol_id)
    #             if term:
    #                 ontol_label = term['name']
    #                 ontol_def = str(term['definition'])                    
    #                 ontol_namespace = term['id'][0:term['id'].index(":")]
    
    #             tag_results.append({"ontol_id": token._.ontol_id,
    #                                 "span_text": token.text,
    #                                 "ontol_label": ontol_label,
    #                                 "ontol_def": ontol_def,
    #                                 "ontol_namespace": ontol_namespace,
    #                                 "ontol_link": "http://addictovocab.org/"+token._.ontol_id,
    #                                 "match_index": token.idx})

    return render_template('index.html',
                           text=text,
                           details=details,
                           date=date,
                           title=title,
                           authors=authors,
                           id=idName,
                           tag_results=tag_results,
                           development=development)



if __name__ == "__main__":
    app.run()
