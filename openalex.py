

import requests
import json
import collections

def unpack_inverted_index(input):
    all_dict = {}
    for ind in input:
        for i in range(len(input[ind])):
            key = input[ind][i]
            all_dict[key] = ind
    sorted_dict = collections.OrderedDict(sorted(all_dict.items()))

    full_string = ""
    for item in sorted_dict.values():
        full_string += item + " "

    # print(full_string)
    return full_string

def get_one_with_id(id):
    location = f"https://api.openalex.org/works/" + str(id)+ ","    
    #header is for the "polite pool": https://docs.openalex.org/api#the-polite-pool
    headers = {
        'User-Agent': 'Python-requests/2.25.1',
        'mailto': 'tom@devsoft.co.za'  
    }
    response = requests.get(location, headers=headers)
    alexjson = response.json() if response and response.status_code == 200 else None
    # print(alexjson)
    if alexjson:
        #title:
        if 'display_name' in alexjson: 
            print("title: ", alexjson['display_name'])

        #authors:
        all_authors = []
        if 'authorships' in alexjson:
            for auth in alexjson['authorships']:
                if 'author' in auth:
                    if 'display_name' in auth['author']:
                        all_authors.append(auth['author']['display_name'])
        print("all_authors: ", all_authors)
        
        #abstracts:
        abstracts = {}
        if 'abstract_inverted_index' in alexjson:                 
            if alexjson['abstract_inverted_index'] == None:
                pass
            else:
                # print("abstract_inverted_index is: ", alexjson['abstract_inverted_index'])
                abstracts = alexjson['abstract_inverted_index']   
                full_string = unpack_inverted_index(abstracts)
                print("full abstract: ", full_string)





def get_publications(page, end, alex_id_list):
    print("page: ", page)
    #exit recursion:
    if page > end:
        return (alex_id_list)
    
    #todo: per-page doesn't seem to do anything? per_page also same... default 25 per page
    location = f"https://api.openalex.org/works?page=" + str(page) + "?per-page=50?filter=publication_year:2021"
    headers = {
        'User-Agent': 'Python-requests/2.25.1',
        'mailto': 'tom@devsoft.co.za'  
    }
    # print("Fetching api json data from query: ", location)
    
    response = requests.get(location, headers=headers)
    alexjson = response.json() if response and response.status_code == 200 else None
    
    if alexjson:
        for result in alexjson['results']:
            if 'ids' in result:
                if 'openalex' in result['ids']:
                    alex_id_list.append(result['ids']['openalex'])
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
    # print("a: ", a[-1])
    alex_id_list.append(a[-1])
print("done")

one = get_one_with_id(alex_id_list[0]) #test only one
print("got one: ", one)

# one = get_one_with_id("W2741809807")
# print("got one: ", one)
