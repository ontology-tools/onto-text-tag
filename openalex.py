

import requests
from urllib.request import urlopen
import json

def get_publications(page):

    #exit recursion:
    print("page: ", page)
    if page > 1: 
        return page

    location = f"https://api.openalex.org/works?page=" + str(page) + "?per-page=50?filter=publication_year:2021"
    print("Fetching release file from", location)
    # data = urlopen(location).read()  # bytes
    # datajson = json.load(data)
    response = requests.get(location)
    alexjson = response.json() if response and response.status_code == 200 else None
    
    if alexjson:
        alex_meta = alexjson['meta']
        print("Received", alex_meta['count'], "publications")
        # print("alexjson is: ", alexjson)
        # print("alexjson type is: ", type(alexjson))
        # print("alexjson keys are: ", alexjson.keys())

        #test:
        # print("datajson values are: ", alexjson.values())
        # print("datajson items are: ", alexjson.items())
        
        for result in alexjson['results']:
            print("result: ", result['id'])
            if 'abstract_inverted_index' in result:                 
                if result['abstract_inverted_index'] == None:
                    continue
                else:
                    print("abstract_inverted_index is: ", result['abstract_inverted_index'])
        #recurse: 
        page = page+1
        get_publications(page)

        
    else:
        print("no response?")


pubs = get_publications(1)
print(pubs)
