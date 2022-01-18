
import pickle
import os
import csv

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
            if ID in ontoterminology:
                if str(pmID) not in ontoterminology[ID]:
                    ontoterminology[ID].append(str(pmID))
                    # print("added pmID: ", pmID)
            else: 
                ontoterminology[ID] = [str(pmID)]  
                # print("added pmID: ", pmID)
    print("got full ontoterminology: ", ontoterminology)  
    # ontoterminology to pickle:
    with open('ontoterminology.pkl', 'wb') as f:
        pickle.dump(ontoterminology, f)
    

# load_ontotermentions()
# print("done")
   
def load_check_pickle():
    with open('ontoterminology.pkl', 'rb') as f:
        check_dict = pickle.load(f)
    print(check_dict)

# load_check_pickle()

