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

from flask import Flask, request, redirect, url_for
from flask.templating import render_template
# from ontotagtext import ExtractorComponent
from ontotagtext import MultiExtractorComponent
import spacy
from Bio import Entrez
import requests
from urllib.request import urlopen

import re

import pprint
pp = pprint.PrettyPrinter(depth=4)

app = Flask(__name__)
app.config.from_object('config')
idName = "ID"
# python -m spacy download en_core_web_md
# or: en_core_web_sm or en_core_web_lg
nlp = spacy.load('en_core_web_md')
nlp2 = spacy.load('en_core_web_md')
nlp3 = spacy.load('en_core_web_md')


location = f"https://raw.githubusercontent.com/addicto-org/addiction-ontology/master/addicto-merged.owx"
location2 = f"https://raw.githubusercontent.com/HumanBehaviourChangeProject/ontologies/master/Upper%20Level%20BCIO/bcio-merged.owx"
# location = f"https://raw.githubusercontent.com/HumanBehaviourChangeProject/ontologies/master/Upper%20Level%20BCIO/bcio-merged.owx"

print("Fetching release file from", location)
data = urlopen(location).read()  # bytes
print("Fetching release file from", location2)
data2 = urlopen(location2).read()  # bytes

ontofile1 = data.decode('utf-8')
ontofile2 = data2.decode('utf-8')


# combined test
# todo: create a dictionary which can be populated with {"label1": name1, ontofile1}, {"label2": ...}
ontoDict = {
    "ontologies": [
        {
            "label": "AddictO",
             "name": "AddictO",
            "ontologyfile": ontofile1
        },
        {
            "label": "BCIO",
            "name": "BCIO",
            "ontologyfile": ontofile1
        }
    ]
}
    

onto_extractor3 = MultiExtractorComponent(
    nlp3,
    ontoDict=ontoDict
    )
nlp3.add_pipe(onto_extractor3, after="ner")


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


# Pages for the app
@app.route('/')
@app.route('/home')
def home():
    return redirect(url_for('associations')) #todo: temporary testing redirect - remove this
    # return render_template('index.html')

@app.route('/associations')
def associations():
    return render_template('associations.html')

@app.route('/pubmed', methods=['POST', 'GET'])
def pubmed():
    id = request.form.get('pubmed_id')
    global idName
    articleDetails = ""
    idName = ""
    print(f"Pubmed id {id}")
    if id:
        print(f"Got it {id}")
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
    return render_template('index.html', error_msg=f"No abstract found for PubMed ID {id}")


# Text tagging app

@ app.route('/tag', methods=['POST'])
def tag():
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

    doc3 = nlp3(text)
    # get ontology IDs identified
    for token in doc3:
        if token._.is_ontol_term:
            # print("token details: ", token._.ontol_id, token.text, token.idx)
            term=onto_extractor3.get_term(token._.ontol_id)

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
                           tag_results=tag_results)


if __name__ == "__main__":
    app.run()
