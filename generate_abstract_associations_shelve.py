
import pickle
import re
import shelve

#the allAbstracts.pkl file is generated when we create ontotermentions
def associate():
    #convert pickle to shelve db: 
    pickle_in = open("/home/tom/Documents/PROGRAMMING/Python/addiction-ontology/scripts/allAbstracts.pkl","rb")
    abstract_associations = pickle.load(pickle_in)
    with shelve.open('allAbstracts.db', "c") as db: #will overwrite
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

def other_values():
    titles_in = open("/home/tom/Documents/PROGRAMMING/Python/addiction-ontology/scripts/allTitles.pkl","rb")
    titles = pickle.load(titles_in)
    
    with shelve.open('allTitles.db', "c") as db: #will overwrite
        for id in titles: 
            # print("adding: ", id)
            one_title = titles[id]    
            #add to shelf:
            db[id]=one_title
        db.close()

    dates_in = open("/home/tom/Documents/PROGRAMMING/Python/addiction-ontology/scripts/allDates.pkl","rb")
    dates = pickle.load(dates_in)
    with shelve.open('allDates.db', "c") as db: #will overwrite
        for id in dates: 
            # print("adding: ", id)
            one_date = dates[id]    
            #add to shelf:
            db[id]=one_date
        db.close()  

    authors_in = open("/home/tom/Documents/PROGRAMMING/Python/addiction-ontology/scripts/allAuthorAffiliations.pkl","rb")  
    authors = pickle.load(authors_in)
    with shelve.open('allAuthors.db', "c") as db: #will overwrite
        for id in authors: 
            # print("adding: ", id)
            one_author = authors[id]    
            #add to shelf:
            db[id]=one_author
        db.close()


    
    
    

#associations:
associate()
#date, title and authors:
other_values()

#test that it worked:
with shelve.open('allTitles.db') as db_titles: 
    print("title: ", db_titles['35021268'])
with shelve.open('allDates.db') as db_dates: 
    print("date: ", db_dates['35021268'])
with shelve.open('allAuthors.db') as db_authors: 
    print("authors: ", db_authors['35021268'])
with shelve.open('allAbstracts.db') as dbread: 
    print("abstract: ", dbread['35021268'])


