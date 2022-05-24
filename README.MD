# Building - notes:

addiction-ontology branch: updateIndex

onto-text-tag branch: development

environment: python 3.8

## 1.
In onto-text-tag root folder:
build ontotermentions.csv, as well as allTitles.pkl, allDates.pkl, allAuthors.pkl, allAbstracts.pkl:

> python build_ontotermentions.py --path "path/to/addiction-ontology"

*note this takes 20+ hours on my machine, and requires a stable internet connection for the first few hours at least (fetching abstracts part)*

## 2. 
Run the app: 

> flask run

On the front end (127.0.0.1:5000) go through tagging a PMID and "visualise associations" for any PMID (find one from static/allAbstracts.db.bak). This will initialise the static/test_terms.tsv and static/test_terms.tsv.pickle

