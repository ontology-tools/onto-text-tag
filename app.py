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

from flask import Flask, request
from flask.templating import render_template
from ontotagtext import ExtractorComponent
import spacy

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


# Pages for the app
@app.route('/')
@app.route('/home')
def home():
    return render_template('index.html')


@app.route('/pubmed', methods=['POST'])
def pubmed():
    id = request.form.get('pubmed_id')
    # id = request.get_json()
    print(f"Pubmed id {id}")
    if id:
        print(f"Got it {id}")
        return render_template('index.html')
    else:
        print(f"Got nothing")
        return render_template('index.html')
    #     text = id)
    #    tag_results=tag_results)


# Text tagging app

@ app.route('/tag', methods=['POST'])
def tag():
    text=request.form['inputText']
    print(f"Got input text {text}")
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
    print(f"Got tag results {tag_results}")

    return render_template('index.html',
                           text = text,
                           tag_results = tag_results)


if __name__ == "__main__":
    app.run()
