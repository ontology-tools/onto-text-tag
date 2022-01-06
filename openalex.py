

import requests
from urllib.request import urlopen
import json



def get_one_with_id(id):
    location = f"https://api.openalex.org/works/" + str(id)+ ","    
    headers = {
        'User-Agent': 'Python-requests/2.25.1',
        'mailto': 'tom@devsoft.co.za'  
    }
    response = requests.get(location, headers=headers)
    alexjson = response.json() if response and response.status_code == 200 else None
    print(alexjson)
    

def get_publications(page, end, alex_id_list):
    counter = 0
    #exit recursion:
    print("page: ", page)
    if page > end:
        # print("ales_id_list is: ", alex_id_list)
        return (alex_id_list)
    
    #todo: per-page doesn't seem to do anything? per_page also same... default 25 per page
    location = f"https://api.openalex.org/works?page=" + str(page) + "?per-page=50?filter=publication_year:2021"
    # print("Fetching api json data from query: ", location)
    
    response = requests.get(location)
    alexjson = response.json() if response and response.status_code == 200 else None
    
    if alexjson:
        # alex_meta = alexjson['meta']
        # print("meta: ", alex_meta)
        # print("Received", alex_meta['count'], "publications")        

        #test:
        # print("alexjson is: ", alexjson)
        # print("alexjson type is: ", type(alexjson))
        # print("alexjson keys are: ", alexjson.keys())
        # print("datajson values are: ", alexjson.values())
        # print("datajson items are: ", alexjson.items())
        
        for result in alexjson['results']:
            counter = counter + 1
            if 'ids' in result:
                # print("ids: ", result['ids'])
                if 'openalex' in result['ids']:
                    print("openalex: ", result['ids']['openalex'])
                    alex_id_list.append(result['ids']['openalex'])

            # print("checking result: ", counter)
            # print("result: ", result['id'])
            # print("name: ", result['display_name'])
            if 'abstract_inverted_index' in result:                 
                if result['abstract_inverted_index'] == None:
                    continue
                else:
                    print("abstract_inverted_index is: ", result['abstract_inverted_index'])
                    return "got one finally!" #none from this type of request, todo: remove this check
        #recurse: 
        page = page+1
        return_alex_id_list = get_publications(page, end, alex_id_list)
        return return_alex_id_list

        
    else:
        print("no response?")


# test getting data from alex api: 
alex_id_list = []
alex_id_url_list = get_publications(1, 1, [])
# print("alex_id_list is: ", alex_id_list)
#get alex after last "/" in url:
for alex in alex_id_url_list:
    a = alex.rsplit("/")
    print("a: ", a[-1])
    alex_id_list.append(a[-1])
print("done")


one = get_one_with_id(alex_id_list[0])
print("got one: ", one)

# one = get_one_with_id("W2741809807")
# print("got one: ", one)
