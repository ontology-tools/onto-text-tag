

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
    full_string = ""
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
    #just the string for test: 
    return full_string





def get_publications(page, end, alex_id_list):
    print("page: ", page)
    #exit recursion:
    if page > end:
        return (alex_id_list)
    
    #todo: per-page doesn't seem to do anything? per_page also same... default 25 per page
    # location = f"https://api.openalex.org/works?page=" + str(page) + "?per-page=50?filter=publication_year:2021"
    location = f"https://api.openalex.org/works?filter=publication_year:2021&page=" + str(page) + "&per-page=50" #todo: test if this is the correct format

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

def get_alex():
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
    return one

#test get_alex(): 
# get_alex()

def alexquery(query, years):
    # location = f"https://api.openalex.org/works?filter=title.search:cubist&page=2&per-page=50" #working

    #test get paginated list of records > 10 000 - so page > 200. Doesn't return anything. 
    location = f"https://api.openalex.org/works?filter=publication_year:2021&page=199&per-page=50" #working query, nothing over page 200

    # location = f"https://api.openalex.org/works?filter=title.search:addiction,smoking,vaping,behaviour&filter=publication_year:2011,2012,2013,2014,2015,2016,2017,2018,2019,2020,2021" #not working
    
    # location = f"https://api.openalex.org/works?page=" + str(page) + "?per-page=50?filter=publication_year:2021" #per-page doesn't seem to do anything? per_page also same... default 25 per page
    headers = {
        'User-Agent': 'Python-requests/2.25.1',
        'mailto': 'tom@devsoft.co.za'  
    }
    # print("Fetching api json data from query: ", location)
    
    response = requests.get(location, headers=headers)
    alexjson = response.json() if response and response.status_code == 200 else None
    # print("alexjson: ", alexjson)
    if alexjson:
        print("meta: ", alexjson['meta'])
        num_results = alexjson['meta']['count']
        num_pages = alexjson['meta']['count']/25
        print("page: ", alexjson['meta']['page'], " of ", num_pages, ' in ', num_results, ' results.')
        test_id_list = []
        for result in alexjson['results']:
            if 'ids' in result:
                if 'openalex' in result['ids']:
                    a = result['ids']['openalex'].rsplit("/")
                    test_id_list.append(a[-1])
        print("test_id_list: ", test_id_list)

#test alexquery():
alexquery("addiction", "2011")