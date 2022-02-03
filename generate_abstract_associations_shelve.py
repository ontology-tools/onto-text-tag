#usage: 
#python generate_abstract_associations_shelve.py --path "/home/tom/Documents/PROGRAMMING/Python/addiction-ontology"


import pickle
import re
import shelve
import os
import argparse
import sys

#the allAbstracts.pkl file is generated when we create ontotermentions
def associate(addictopath):
    #convert pickle to shelve db: 
    pickle_in = open(addictopath + "allAbstracts.pkl","rb")
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

def other_values(addictopath):
    titles_in = open(addictopath + "allTitles.pkl","rb")
    titles = pickle.load(titles_in)
    
    with shelve.open('allTitles.db', "c") as db: #will overwrite
        for id in titles: 
            # print("adding: ", id)
            one_title = titles[id]    
            #add to shelf:
            db[id]=one_title
        db.close()

    dates_in = open(addictopath + "allDates.pkl","rb")
    dates = pickle.load(dates_in)
    with shelve.open('allDates.db', "c") as db: #will overwrite
        for id in dates: 
            # print("adding: ", id)
            one_date = dates[id]    
            #add to shelf:
            db[id]=one_date
        db.close()  

    authors_in = open(addictopath + "allAuthorAffiliations.pkl","rb")  
    authors = pickle.load(authors_in)
    with shelve.open('allAuthors.db', "c") as db: #will overwrite
        for id in authors: 
            # print("adding: ", id)
            one_author = authors[id]    
            #add to shelf:
            db[id]=one_author
        db.close()


    
    
    


if __name__ == '__main__':
    print("File one executed when ran directly")
    parser=argparse.ArgumentParser()
    parser.add_argument('--path', '-i',help='Name of the path')
    
    args=parser.parse_args()
    path = args.path
    
    
    if args is None :
        parser.print_help()
        sys.exit('Not enough arguments. Expected at least -i "Path to ontology folder" ')
    else:
        print(path)
    
    if path is None: 
        parser.print_help()
        sys.exit('Not enough arguments. Expected at least -i "Path to ontology folder" ')
    else:
        print("got it: ", path)
        #todo: get ontopath from current directory
        ontopath = os.path.dirname(os.path.realpath(__file__))
        print("ontopath is: ", ontopath)

        addictopath = path + "/scripts/"

        #associations:
        associate(addictopath)
        #date, title and authors:
        other_values(addictopath)

        #optional - test that it worked:
        with shelve.open('allTitles.db') as db_titles: 
            print("title: ", db_titles['35100637'])
        with shelve.open('allDates.db') as db_dates: 
            print("date: ", db_dates['35100637'])
        with shelve.open('allAuthors.db') as db_authors: 
            print("authors: ", db_authors['35100637'])
        with shelve.open('allAbstracts.db') as dbread: 
            print("abstract: ", dbread['35100637'])








