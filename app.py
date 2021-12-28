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
from flask import Flask, request, redirect, url_for, session, Response, stream_with_context 
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


#test OGER:
import oger
from oger.ctrl.router import Router, PipelineServer
conf = Router(termlist_path='static/test_terms.tsv')
pl = PipelineServer(conf)

development = False

pp = pprint.PrettyPrinter(depth=4)

app = Flask(__name__)

page = Template("""
<!DOCTYPE html>
<html lang="en">
<head>
  {{ resources }}
</head>
<body>
  <div id="myplot"></div>
  <script>
  fetch('/plot')
    .then(function(response) { return response.json(); })
    .then(function(item) { return Bokeh.embed.embed_item(item); })
  </script>
</body>
""")

colormap = {'setosa': 'red', 'versicolor': 'green', 'virginica': 'blue'}
colors = [colormap[x] for x in flowers['species']]

def make_plot(x, y):
    p = figure(title = "Iris Morphology", sizing_mode="fixed", width=400, height=400)
    p.xaxis.axis_label = x
    p.yaxis.axis_label = y
    p.circle(flowers[x], flowers[y], color=colors, fill_alpha=0.2, size=10)
    return p

app.config.from_object('config')
idName = "ID"
# python -m spacy download en_core_web_md
# or: en_core_web_md or en_core_web_lg
# nlp = spacy.load('en_core_web_md')
nlp = en_core_web_sm.load()

location = f"https://raw.githubusercontent.com/addicto-org/addiction-ontology/master/addicto-merged.owx"
location2 = f"https://raw.githubusercontent.com/HumanBehaviourChangeProject/ontologies/master/Upper%20Level%20BCIO/bcio-merged.owx"
# location = f"https://raw.githubusercontent.com/HumanBehaviourChangeProject/ontologies/master/Upper%20Level%20BCIO/bcio-merged.owx"

print("Fetching release file from", location)
ontol1 = pyhornedowl.open_ontology(urlopen(location).read().decode('utf-8'))
print("Fetching release file from", location2)
ontol2 = pyhornedowl.open_ontology(urlopen(location2).read().decode('utf-8'))

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

    #print("should be getting descendants here")
    #todo: refactor below:
    for entry in id_list:
        #print("looking at entry: ", entry)
        entryIri = ontol1.get_iri_for_id(entry.replace("_", ":"))
        if entryIri:
            #print("looking at entryIri: ", entryIri)
            descs = pyhornedowl.get_descendants(ontol1, entryIri)
            if len(descs) > 0:
                for d in descs:
                    try:
                        add_id = ontol1.get_id_for_iri(d).replace(":", "_")
                        descendent_ids.append(add_id.replace("_", ":"))
                    except:
                        print("error when getting descendents of ",d)
                # if add_id:
                #     if add_id not in descendent_ids:
                #         print("adding id: ", add_id)
                #         descendent_ids.append(add_id)
                # id_list.append(repo1.get_id_for_iri(d).replace(":", "_")) #todo: does adding this to same array cause issues? 
    for entry in id_list:
        #print("looking at entry: ", entry)
        entryIri = ontol2.get_iri_for_id(entry.replace("_", ":"))
        if entryIri:
            #print("looking at entryIri: ", entryIri)
            descs = pyhornedowl.get_descendants(ontol2, entryIri)
            if len(descs) > 0:
                for d in descs:
                    try:
                        add_id = ontol1.get_id_for_iri(d).replace(":", "_")
                        #print("add_id is: ", add_id)
                        descendent_ids.append(add_id.replace("_", ":"))
                    except:
                        print("error getting descendents of ",d)
                #todo: remove duplicates
                # if add_id:
                #     if add_id not in descendent_ids:
                #         print("adding id: ", add_id)
                        # descendent_ids.append(add_id)
                # id_list.append(repo2.get_id_for_iri(d).replace(":", "_"))
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



# Pages for the app
@app.route('/')
@app.route('/home')
def home():
    if os.environ.get("FLASK_ENV")=='development':
        development = True
    # return redirect(url_for('associations')) #todo: temporary testing redirect - remove this
    return render_template('index.html', development = development)



@app.route('/associations', methods=['GET', 'POST'])
def associations():
    ontologies = ontoDict["ontologies"]
    labels = get_ids(ontologies) 
    label_list={'labels': labels}
    #test values:
    # label_list = {'labels': ["ADDICTO:123457|smoking", "ADDICTO:123456|vaping", "ADDICTO:123458|human being", "BCIO:123456|addiction", "BCIO:223456|long term ", "ADDICTO:133456|ontology" ]}
    json.dumps(label_list)
    return render_template('associations.html', label_list = label_list)



@app.route('/visualise_associations', methods=['POST'])
def visualise_associations():  
    ontology_id_list = json.loads(request.form.get('ontology_id_list')) 
    print("got ontology_id_list: ", ontology_id_list)
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
    # print("ontology_id_list is: ", ontology_id_list)
    include_descendent_classes = request.form.get('include_descendent_classes')
    # print("checkbox says: ", include_descendent_classes)
    print("time for a similarity visual!")
    # hv_generator(ontology_id_list) #todo: something different here?
    # return ( json.dumps({"message":"Success"}), 200 )
    return render_template('similarity.html')

#todo: below no longer needed?
# @app.route('/similarity')
# def similarity():    
#     # iframe = url_for('chordout')
#     return render_template("similarity.html") #, iframe=iframe) 

@app.route('/pubmed', methods=['POST', 'GET'])
def pubmed():
    #test OGER:
    coll = pl.load_one(['29148565'], fmt='pubmed')
    print(coll[0][0].text) # title
    pl.process(coll) 
    # entity = next(coll[0].iter_entities())
    # print(entity.info)
    for entity in coll[0].iter_entities():
        print(entity.start, entity.text, entity.end)

    if os.environ.get("FLASK_ENV")=='development':
        development = True
    id = request.form.get('pubmed_id')
    global idName
    articleDetails = ""
    idName = ""
    #print(f"Pubmed id {id}")
    if id:
        #print(f"Got it {id}")
        idName = f"{id}"
        try:
            results = fetch_details([id])
            for result in results:
                resultDetail = results[result]
                abstractText = get_abstract_text(resultDetail)
                # print(f"Got abstract text {abstractText}")
                articleDetails = get_article_details(resultDetail)
                # print("Got articleDetails: ", articleDetails) #when we get the right details...
                try:
                    dateA, titleA, authorsA = articleDetails.split(';')
                except:
                    pass
                if abstractText:
                    r = requests.post(url_for("tag", _external=True), data={
                                      "inputDetails": articleDetails, "inputText": abstractText, "dateDetails": dateA, "titleDetails": titleA, "authorsDetails": authorsA})
                    return r.text, r.status_code, r.headers.items()
        except Exception as err:  # 400 bad request handling, also if no internet connection
            print(err)
    return render_template('index.html', error_msg=f"No abstract found for PubMed ID {id}", development = development)


# Text tagging app

@ app.route('/tag', methods=['POST'])
def tag():
    if os.environ.get("FLASK_ENV")=='development':
        development=True
    text = request.form['inputText']
    details = request.form.get('inputDetails')
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

    doc3 = nlp(text)
    # get ontology IDs identified
    for token in doc3:
        if token._.is_ontol_term:
            # print("token details: ", token._.ontol_id, token.text, token.idx)
            term=onto_extractor3.get_term(token._.ontol_id)
            # print("ontol_id is: ", token._.ontol_id)
            # print("term is: ", term)
            if term:
                ontol_label = term['name']
                # print("ontol_label: ", ontol_label)
                ontol_def = str(term['definition'])
                # print("ontol_def: ", ontol_def)
                ontol_namespace = term['id'][0:term['id'].index(":")]
                # print("ontol_namespace: ", ontol_namespace)
            else:
                ontol_label = token.idx
                ontol_def = token.text
                ontol_namespace = ""
                # print("ontol_namespace not found")
            tag_results.append({"ontol_id": token._.ontol_id,
                                "span_text": token.text,
                                "ontol_label": ontol_label,
                                "ontol_def": ontol_def,
                                "ontol_namespace": ontol_namespace,
                                "ontol_link": "http://addictovocab.org/"+token._.ontol_id,
                                "match_index": token.idx})

    # print(f"Got tag results {tag_results}")

    return render_template('index.html',
                           text=text,
                           details=details,
                           date=date,
                           title=title,
                           authors=authors,
                           id=idName,
                           tag_results=tag_results,
                           development=development)

# get id from label here:
# @app.route('/get_ids', methods=['GET', 'POST'])
# def get_ids():
#     ontology_id = request.form.get('ontology_id') #get ontology_label?
#     print("got ontology_id: ", ontology_id)
#     return ( json.dumps({"message":"Success",
#                              "response": ontology_id}), 200 )

if __name__ == "__main__":
    app.run()
