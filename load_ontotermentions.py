
import pickle
import os
import csv
import shelve

# todo: pickle file, shelf db is getting put into incorrect format here
# correct format for reference: ontoterminology_old.pkl
# old ontoterminology file: ontotermentionsJAN2022.csv
# how to add dictionary to shelf? 

def load_ontotermentions(): 
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
    

# load_ontotermentions()


def ontotermentions_to_shelve():
    
    terms_in = open("ontoterminology.pkl","rb")
    terms = pickle.load(terms_in)
    
    with shelve.open('ontoterminlology.db', "c") as db: #will overwrite
        for id in terms: 
            print("adding: ", terms[id], "with id: ", id)
            one_term = terms[id]   
            #add to shelf:
            db[id]=one_term
            # print("db[id]: ", id, ", one_term: ", one_term)
        db.close()

ontotermentions_to_shelve()

def check_shelve_db():
    with shelve.open('ontoterminology.db') as terms: 
        print(terms['ADDICTO:0000717'])
        # ontology_id_list =  ['http://addictovocab.org/ADDICTO_0000468', 'http://humanbehaviourchange.org/ontology/BCIO_038000', 'http://purl.obolibrary.org/obo/IAO_0000007', 'http://addictovocab.org/ADDICTO_0000372']
        # for term in terms:
        #     print(term.items())
            

check_shelve_db()

def load_check_pickle():
    print("CHECKING PICKLE")
    with open('ontoterminology.pkl', 'rb') as f:
        check_dict = pickle.load(f)
    #test new ontotermentions limit to selected ID's: 
    ontology_id_list =  ['ADDICTO:0000692', 'BCIO:038000', 'IAO:0000007', 'ADDICTO:0000372']
    for ont in ontology_id_list:
        for key in list(check_dict.keys()):
            # print(key)
            if ont == key: #found a match. 
                print("key: ", key)
                print(check_dict[key])
    # print(check_dict)

# load_check_pickle()

