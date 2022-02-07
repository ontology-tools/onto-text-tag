# Building - notes:

addiction-ontology branch: updateIndex

onto-text-tag branch: shelve_update

environment: python 3.8

## 1.
In onto-text-tag root folder:
build ontotermentions.csv, as well as allTitles.pkl, allDates.pkl, allAuthors.pkl, allAbstracts.pkl:

> python build_ontotermentions.py --path "path/to/addiction-ontology"

*note this takes 20+ hours on my machine, and requires a stable internet connection for the first few hours at least (fetching abstracts part)*

## 2.
Copy ontotermentions.csv across:
Manually copy ontotermentions.csv from addiction-ontology to onto-text-tag /static folder

## 3. 
Create shelve db's from allTitles.pkl, allDates.pkl, allAuthors.pkl, allAbstracts.pkl:

> python generate_abstract_associations_shelve.py –path “/path/to/addiction-ontology”

*each shelve db consists of 3 files - .bak .dat and .dir - but addressed in code as .db*

## 4.
Create ontotermentions.pkl from ontotermentions.csv, then turn it into a shelve db:

> python load_ontotermentions.py

## 5.
Build test_terms.tsv for use by OGER:
First delete old static/test_terms.tsv and static/test_terms.tsv.pickle

> rm static/test_terms.tsv static/test_terms.tsv.pickle

In `app.py` change build_terms to True
 
> build_terms = True

Run the app: 

> flask run

On the front end (127.0.0.1:5000) go through tagging a PMID and "visualise associations" for any PMID
This will initialise the static/test_terms.tsv and static/test_terms.tsv.pickle
In `app.py` change build_terms to False again for normal use (test_terms.tsv also not needed anymore - it is ignored by OGER if the static/test_terms.tsv.pickle is present)

> build_terms = False
