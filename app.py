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
from ontotagtext import ExtractorComponent
import spacy
from Bio import Entrez
import requests

import pprint
pp = pprint.PrettyPrinter(depth=4)

app = Flask(__name__)
app.config.from_object('config')

# python -m spacy download en_core_web_md
# or: en_core_web_sm or en_core_web_lg
nlp = spacy.load('en_core_web_md')


onto_extractor = ExtractorComponent(
    nlp,
    name="ADDICTO",
    label="ADDICTO",
    ontologyfile="static/addicto.obo")
nlp.add_pipe(onto_extractor, after="ner")


# Interaction with PubMed: get detailed results for a list of IDs
def fetch_details(id_list):
    ids = ','.join(id_list)
    Entrez.email = 'janna.hastings@ucl.ac.uk'
    handle = Entrez.efetch(db='pubmed',
                           retmode='xml',
                           id=ids)
    results = Entrez.read(handle)
    # print(f"")
    # print(f"results are: " + results)
    # print(f"")
    return results
    
#parse the title, authors and date published 
def get_article_details(result):
    articleDetails = None
    for detail in result:
        if 'MedlineCitation' in detail:
            # print(f"")
            print(f"MedlineCitation is: ")
            pp.pprint(detail['MedlineCitation']) #pretty print
            # print(str(detail['MedlineCitation']))
            if 'DateCompleted' in detail['MedlineCitation']: #this works
                print(f"")
                print(f"DateCompleted is: ")
                print(str(detail['MedlineCitation']['DateCompleted']))
            if 'AuthorList' in detail['MedlineCitation']['Article']: #this works, need to refine though
                print(f"")
                print(f"AuthorList is: ")
                pp.pprint(detail['MedlineCitation']['Article']['AuthorList']) 
            else:
                print(f"AuthorList not found")
            if 'ArticleTitle' in detail['MedlineCitation']['Article']: 
                print(f"")
                print(f"ArticleTitle is: ") #got title!
                pp.pprint(detail['MedlineCitation']['Article']['ArticleTitle'])
                # articleDetails = str(detail['MedlineCitation']['Article']['Abstract']['ArticleTitle'])
            else:
                print(f"ArticleTitle not found")

    return articleDetails

# Parse the PubMed result to get the abstract text if it is there
def get_abstract_text(result):
    abstractText = None
    for detail in result:
        if 'MedlineCitation' in detail:
            if 'Article' in detail['MedlineCitation']:
                if 'Abstract' in detail['MedlineCitation']['Article']:
                    if 'AbstractText' in detail['MedlineCitation']['Article']['Abstract']:
                        abstractText = str(detail['MedlineCitation']['Article']['Abstract']['AbstractText'])

    return abstractText


# Pages for the app
@app.route('/')
@app.route('/home')
def home():
    return render_template('index.html')


@app.route('/pubmed', methods=['POST', 'GET'])
def pubmed():
    id = request.form.get('pubmed_id')
    # id = request.get_json()
    print(f"Pubmed id {id}")
    if id:
        print(f"Got it {id}")
        try:
            results = fetch_details([id]) 
            # print(f"results are: {results}") 
            for result in results:
                resultDetail = results[result]
                abstractText = get_abstract_text(resultDetail)
                # print(f"Got abstract text {abstractText}")
                articleDetails = get_article_details(resultDetail)
                # print(f"Got abstract text {articleDetails}") #when we get the right details...
                if abstractText:
                    r = requests.post(url_for("tag", _external=True), data={"inputText":abstractText})
                    return r.text, r.status_code, r.headers.items()
        except Exception as err: #400 bad request handling, also if no internet connection
            print(err)
    return render_template('index.html', error_msg = f"No abstract found for PubMed ID {id}")            


# Text tagging app

@ app.route('/tag', methods=['POST'])
def tag():
    text=request.form['inputText']
    # print(f"Got input text {text}")
    # process the text
    doc=nlp(text)
    tag_results=[]
    # get ontology IDs identified
    for token in doc:
        if token._.is_ontol_term:
            # print(token._.ontol_id, token.text, token.idx)
            term=onto_extractor.get_term(token._.ontol_id)
            if term:
                ontol_label=term.name
                ontol_def=str(term.definition)
                ontol_namespace=term.namespace
                if ontol_namespace is None:
                    ontol_namespace=term.id[0:term.id.index(":")]
            else:
                ontol_label=""
                ontol_def=""
                ontol_namespace=""
            tag_results.append({"ontol_id": token._.ontol_id,
                                "span_text": token.text,
                                "ontol_label": ontol_label,
                                "ontol_def": ontol_def,
                                "ontol_namespace": ontol_namespace,
                                "ontol_link": "http://addictovocab.org/"+token._.ontol_id,
                                "match_index": token.idx})
    # print(f"Got tag results {tag_results}")

    return render_template('index.html',
                           text = text,
                           tag_results = tag_results)


if __name__ == "__main__":
    app.run()
