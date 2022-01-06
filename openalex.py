

import requests
from urllib.request import urlopen
import json


def get_publications(page):
    counter = 0
    #exit recursion:
    print("page: ", page)
    if page > 20: 
        return page
    
    #todo: per-page doesn't seem to do anything? per_page also same... default 25 per page
    location = f"https://api.openalex.org/works?page=" + str(page) + "?per-page=50?filter=publication_year:2021"
    # print("Fetching api json data from query: ", location)
    
    response = requests.get(location)
    alexjson = response.json() if response and response.status_code == 200 else None
    
    if alexjson:
        alex_meta = alexjson['meta']
        print("meta: ", alex_meta)
        # print("Received", alex_meta['count'], "publications")        

        #test:
        # print("alexjson is: ", alexjson)
        # print("alexjson type is: ", type(alexjson))
        # print("alexjson keys are: ", alexjson.keys())
        # print("datajson values are: ", alexjson.values())
        # print("datajson items are: ", alexjson.items())
        
        for result in alexjson['results']:
            counter = counter + 1
            print("checking result: ", counter)
            # print("result: ", result['id'])
            # print("name: ", result['display_name'])
            if 'abstract_inverted_index' in result:                 
                if result['abstract_inverted_index'] == None:
                    continue
                else:
                    print("abstract_inverted_index is: ", result['abstract_inverted_index'])
                    return "got one finally!"
        #recurse: 
        page = page+1
        get_publications(page)

        
    else:
        print("no response?")


pubs = get_publications(1)
print(pubs)
