
import pickle
import os
import csv
import shelve
import pyhornedowl
from urllib.request import urlopen
from ontotagtext import MultiExtractorComponent
from ontotagtext import PREFIXES
from ontotagtext import RDFSLABEL
import en_core_web_sm
import spacy

# todo: pickle file, shelf db is getting put into incorrect format here
# correct format for reference: ontoterminology_old.pkl
# old ontoterminology file: ontotermentionsJAN2022.csv
# how to add dictionary to shelf? 
nlp = en_core_web_sm.load()

location = f"https://raw.githubusercontent.com/addicto-org/addiction-ontology/master/addicto-merged.owx"
location2 = f"https://raw.githubusercontent.com/HumanBehaviourChangeProject/ontologies/master/Upper%20Level%20BCIO/bcio-merged.owx"

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

onto_extractor3 = MultiExtractorComponent(
    nlp,
    ontoDict=ontoDict
    )
nlp.add_pipe(onto_extractor3, after="ner")

def load_ontotermentions(): 
    # all_titles_db = shelve.open('allTitles.db') #for labels
    # print("loaded abstract associations db")
    """
    Loads the ontoterminology.csv file and returns a dictionary of the form:
    ID: [pmID, pmID, pmID]
    1259047,ENVO:01000838,smoke,22241744
    """
    ontoterminology = {}
    #build ontotermentions.csv path from current directory:
    ontopath = os.path.dirname(os.path.realpath(__file__))
    ontoterminology_path = ontopath + "/static/ontotermmentions.csv"
    with open(ontoterminology_path, mode = 'r') as f:
        csvFile = csv.reader(f, delimiter = ',')
        count = 0
        for line in csvFile:
            count = count + 1
            print(count)
            number, ID, name, pmID = line[0], line[1], line[2], line[3]
            split_ID = ID.rsplit("/", 1)
            try: 
                ID = split_ID[1] #new version of ontotermentions is full path now, split out ID
                ID = ID.replace("_", ":")
            except: 
                ID = ID.replace("_", ":")
            #this is a workaround to get missing names:
            
                
            try: 
                name = onto_extractor3.get_term.get_term(ID)['name']


                # for f in onto_extractor3.terms:
                #     l = onto_extractor3.get_label(f)
                #     if l is not None:
                #         # print("looking at label: ", l)
                #         term=onto_extractor3.get_term(l['id'])
                #         if term:                        
                #             if term['id'] == ID:
                #                 name = onto_extractor3.get_term(l['name'])
                #                 print("found a name: ", name)
                # name = ontol1.get_iri_for_id(ID).label #todo: did get_iri_for_id work??
                # if name is None or "": 
                #     name = ontol2.get_iri_for_id(ID).label
            except: 
                name = name #probably ""
            if ID in ontoterminology:                
                if pmID not in ontoterminology[ID]:
                    ontoterminology[ID]['PMID'].append(pmID)
                    #got name already!
                    # print("added pmID: ", pmID)
            else: 
                ontoterminology[ID] = {'PMID': [pmID], 'NAME': name} 
                # print("added pmID: ", pmID)
    # print("got full ontoterminology: ", ontoterminology)  
    # ontoterminology to pickle:
    with open('ontoterminology.pkl', 'wb') as f:
        pickle.dump(ontoterminology, f)
    

load_ontotermentions()


def ontotermentions_to_shelve():
    
    terms_in = open("ontoterminology.pkl","rb")
    terms = pickle.load(terms_in)
    
    with shelve.open('ontoterminology.db', "c",) as db: #will overwrite
        for id in terms: 
            # print("adding: ", terms[id], "with id: ", id)
            one_term = terms[id]   
            #add to shelf:
            db[id]=one_term
            # print("added id: ", id, ", one_term: ", one_term)
            # print("added: ", id, " -> ", db[id])
        db.close()
    
    # terms_db = shelve.open('ontoterminology')
    # #test inside function:
    # all_keys = list(terms_db.keys())
    # print(list)
    # terms_db.close()

ontotermentions_to_shelve()

def check_shelve_db():
    terms_db = shelve.open('ontoterminology.db')
    print("loaded terms db")
    selectedID = 'ADDICTO:0000990'
    selectedData = terms_db[selectedID]
    print("selectedData: ", selectedData)
    #todo: get selectedID from shelve db terms_db: 
    # all_keys = list(terms_db.keys())
    # all_items = list(terms_db.items())
    # all_values = list(terms_db.values())
    # print(all_keys)
    # print(all_items)
    # print(all_values)


    # print(terms_db[selectedID]) #todo: fix KeyError
    # if terms_db[selectedID] is not None: 
    #     print(terms_db[selectedID])
    # else: 
    #     print("no match found")

    # for key in terms_db.keys():
    #     print("checking ", key)
    #     if selectedID == key: #found a match.
    #         print("NAME: ", terms_db[key]['NAME'], " PMID: ", terms_db[key]['PMID'])
    #     else: 
    #         print("no match found")

    # with shelve.open('ontoterminology.db') as terms: 
    #     print(terms['ADDICTO:0000692'])
        # ontology_id_list =  ['http://addictovocab.org/ADDICTO_0000468', 'http://humanbehaviourchange.org/ontology/BCIO_038000', 'http://purl.obolibrary.org/obo/IAO_0000007', 'http://addictovocab.org/ADDICTO_0000372']
        # for term in terms:
        #     print(term.items())
            

check_shelve_db()

def load_check_pickle():
    print("CHECKING PICKLE")
    with open('ontoterminology.pkl', 'rb') as f:
        check_dict = pickle.load(f)
    #test new ontotermentions limit to selected ID's: 
    ontology_id_list =  ['ADDICTO:0000692']
    for ont in ontology_id_list:
        for key in check_dict:
            # print(key)
            if ont == key: #found a match. 
                print("key: ", key)
                print(check_dict[key])
    # print(check_dict)

# load_check_pickle()

