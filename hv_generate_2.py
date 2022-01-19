import pandas as pd
import holoviews as hv
from holoviews import opts, dim
hv.extension('bokeh')
hv.renderer('bokeh')
hv.output(size=200)
from bokeh.embed import json_item
from bokeh.embed import file_html
import pprint as pp
import pyhornedowl
import requests
from urllib.request import urlopen
import json
from memory_profiler import profile
from timeit import default_timer as timer
import pickle 
import io

#todo: memory usage testing: changing to category

# @profile
def hv_generator(ontology_id_input, should_get_descendents):
    start_time = timer()
    from app import get_all_descendents #todo: refactor to avoid this circular import
    #load ontotermilology.pkl:
    with open('ontoterminology.pkl', 'rb') as f:
        ontoterminology = pickle.load(f)
    #todo: idea here, what about creating a new ontotermentions from the ontology_id_input, with only the ones we need?

    # load csv:
    df2 = pd.read_csv("static/ontotermmentions.csv", na_filter = False, delimiter = ",", index_col=0, engine='c')    
    # load csv with dtype specified (not working): 
    # df2 = pd.read_csv("static/ontotermmentions.csv", delimiter = ",", engine="c", na_filter = False, index_col=0, dtype={'ADDICTOID': 'object' ,'LABEL': 'object', 'PMID': 'int64'}) 
    # df2 = pd.read_csv("static/ontotermmentions.csv", delimiter = ",", engine="c", na_filter = False, index_col=0, dtype={'ADDICTOID': 'category' ,'LABEL': 'category', 'PMID': 'int64'}) 

    print("read csv into dataframe")
    # change dtype to category on the fly (not working):
    # df2["LABEL"] = df2["LABEL"].astype("category")
    # df2["ADDICTOID"] = df2["ADDICTOID"].astype("category")   

    # print("dtypes: ", df2.dtypes)

    # get descendants?
    if should_get_descendents == True:
        ontology_id_list = get_all_descendents(ontology_id_input)
    else:
        ontology_id_list = ontology_id_input 
    
    df2 = df2.drop(df2[~df2.ADDICTOID.isin(ontology_id_list)].index)

    # This creates a table of pairs of terms in the same abstract
    dcp = pd.merge(df2,df2,on="PMID",how="inner")

    dcp = dcp.drop(dcp[dcp.LABEL_x == dcp.LABEL_y].index)    

    # We filter the table just to the ones in the ID list we provided as input 
    
    # We filter the table so that pairs are only represented in one direction, i.e. if we have both (smoking, children) and (children, smoking) for the same PMID we drop the second one

    # solution 2 from Stack Overflow - replaces iterrows():
    dcp['ADDICTOID'] = dcp[['ADDICTOID_x', 'ADDICTOID_y']].apply(sorted, axis=1).apply(tuple)
    dcp = dcp.drop_duplicates(subset=['ADDICTOID', 'PMID'], keep='first')
    dcp = dcp.drop(['ADDICTOID'], axis=1)

    # Now we count the distinct numbers of abstracts this combination appeared in    
    data_chord_plot = dcp.groupby(['LABEL_x', 'LABEL_y'], as_index=False)[['PMID']].count()
    data_chord_plot.columns = ['source','target','value']           

    # print("Final chord plot: ", data_chord_plot)

    #new method using ontoterminology: 
    mentions = {}
    for selectedID in ontology_id_list:
        for key in list(ontoterminology.keys()):
            if selectedID == key: #found a match. 
                # print("key: ", key)
                mentions[ontoterminology[key]['NAME']] = ontoterminology[key]['PMID']
    # print("mentions: ", list(mentions.keys()))

    #todo: below is taking too long... 
    chn_list = []
    for source in list(mentions.keys()):
        for target in list(mentions.keys()): 
            if source.strip() == "" or target.strip() == "":
                pass
            else:
                intersection = list(set(mentions[source]).intersection(set(mentions[target])))
                if len(intersection) > 0: 
                    chn = {"source": source, "target": target, "PMID": len(intersection)}
                    # print("adding chn: ", chn)
                    chn_list.append(chn)
    #todo: join the inverse mirror entries: 
    de_dup_chn_list = []
    for i in chn_list:
        for j in chn_list: 
            if (i['source'] + i['target']) == (j['target'] + j['source']): #todo: finds all, should be half
                # print("found inverse: ", i['source'] +", " + i['target'])
                de_dup_chn_list.append(i) #todo: need to still add the PMID numbers together here!
            

    # print(de_dup_chn_list)
    print("length: ", len(chn_list))
    print("de-dup len: ", len(de_dup_chn_list))

    # Build the data table expected by the visualisation library
    links = data_chord_plot
    node_names = links.source.append(links.target)
    node_names = node_names.unique()
    node_info = {"index":node_names,"name":node_names,"group":[1]*len(node_names)}

    nodes = hv.Dataset(pd.DataFrame(node_info), 'index')
    nodes.data.head()

    chord = hv.Chord((links, nodes)).select(value=(5, None)) #todo: was value=5 - changed to 0 now it works for more? 

    chord.opts(
        opts.Chord(cmap='Category20', edge_cmap='Category20', edge_color=dim('source').str(),
                labels='name', node_color=dim('index').str()))
    print("The time difference is :", timer() - start_time)
    start_time_2 = timer()
    #error message html if no chord plot to show:
    if dcp.empty:
        print('empty dataframe, should create an error message chordout.html here')
        html_error_message = "<!doctype html><div><h4>ERROR CREATING TABLE - no associations found, or possibly some of the ID's were incorrect?</h4></div></html>"
        return(json.dumps(html_error_message))
    else:
        renderer = hv.renderer('bokeh')
        hvplot = renderer.get_plot(chord)
        html = renderer.static_html(hvplot)
        print("Time taken for preparing render :", timer() - start_time_2)
        return json.dumps(html)
        