
import pickle
import re
import shelve

#the allAbstracts.pkl file is generated when we create ontotermentions
def associate():
    #convert pickle to shelve db: 
    pickle_in = open("/home/tom/Documents/PROGRAMMING/Python/addiction-ontology/scripts/allAbstracts.pkl","rb")
    abstract_associations = pickle.load(pickle_in)
    with shelve.open('allAbstracts.db', "c") as db: 
        for id in abstract_associations: 
            # print("adding: ", id)
            one_abstract = abstract_associations[id]    
            #clean up the string:     
            if("StringElement" in one_abstract):
                fixed = re.findall(r'StringElement\((.+?)attributes',one_abstract)
                fixed_abstractText = "".join(fixed)
            else:
                fixed_abstractText = one_abstract.strip('[]') # remove "[]"
            # print(fixed_abstractText)

            #add to shelf:
            db[id]=fixed_abstractText
        db.close()

    #test that it worked:
    with shelve.open('allAbstracts.db') as dbread: 
        print(dbread['35021268'])
    


associate()


