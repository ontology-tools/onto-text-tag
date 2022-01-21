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

# @profile
def hv_generator(ontology_id_input, should_get_descendents):
    try:
        start_time = timer() #test code
        from app import get_all_descendents #todo: refactor to avoid this circular import
        #load ontotermilology.pkl:
        with open('ontoterminology.pkl', 'rb') as f:
            ontoterminology = pickle.load(f)
        #todo: idea here, what about creating a new ontotermentions from the ontology_id_input, with only the ones we need?

        # get descendants
        if should_get_descendents == True:
            print("getting descendents now")
            ontology_id_list = get_all_descendents(ontology_id_input)
            print("got them")
        else:
            ontology_id_list = ontology_id_input 

        #new method using ontoterminology:         
        mentions = {}
        for selectedID in ontology_id_list:
            for key in list(ontoterminology.keys()):
                if selectedID == key: #found a match. 
                    mentions[ontoterminology[key]['NAME']] = ontoterminology[key]['PMID']
        print("loaded mentions")
        chn_list = []
        for source in list(mentions.keys()):
            for target in list(mentions.keys()): 
                if source.strip() == "" or target.strip() == "":
                    pass
                elif source.strip() == target.strip():
                    pass
                else:
                    intersection = list(set(mentions[source]).intersection(set(mentions[target])))
                    if len(intersection) > 0: 
                        chn = {"source": source, "target": target, "PMID": len(intersection)}
                        #attempt to do inverse checking here instead of separately: 
                        add_item = True
                        for k in chn_list: #todo: another whole loop here, any faster way to do this?
                            if source + target == k['target'] + k['source']:
                                add_item = False
                        if add_item:
                            chn_list.append(chn)
        print("finished checking for inverse duplicates..")        
        print("length: ", len(chn_list))

        # Build the data table expected by the visualisation library
        links = pd.DataFrame.from_dict(chn_list)        
        node_names = links.source.append(links.target)
        node_names = node_names.unique()
        node_info = {"index":node_names,"name":node_names,"group":[1]*len(node_names)}

        nodes = hv.Dataset(pd.DataFrame(node_info), 'index')
        nodes.data.head()

        chord = hv.Chord((links, nodes)).select(value=(5, None)) #todo: was value=5 - changed to 0 now it works for more? 

        chord.opts(
            opts.Chord(cmap='Category20', edge_cmap='Category20', edge_color=dim('source').str(),
                    labels='name', node_color=dim('index').str()))
        print("Time taken to find mentions :", timer() - start_time) #test code
        start_time_2 = timer() #test code
        renderer = hv.renderer('bokeh')
        hvplot = renderer.get_plot(chord)
        html = renderer.static_html(hvplot)
        print("Time taken for preparing render :", timer() - start_time_2) #test code
        return json.dumps(html)
    except: 
        html_error_message = "<!doctype html><div><h4>ERROR CREATING TABLE - no associations found, or possibly some of the ID's were incorrect?</h4></div></html>"
        return(json.dumps(html_error_message))
        