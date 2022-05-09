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
import pickle 
import io
# import shelve
import traceback
# from memory_profiler import profile
import shelve 
# import dbm.gnu

# def gdbm_shelve(filename, flag="c", writeback=False): 
#     return shelve.Shelf(dbm.open(filename, flag, writeback)) 

# ontoterminology = gdbm_shelve.open('static/ontoterminology.db', flag='r', writeback=False) 
# ontoterminology = gdbm_shelve('static/ontoterminology.db', flag='r', writeback=False) 

ontoterminology = shelve.open('static/ontoterminology.db', flag='r', writeback=False) #works

print("loaded terms db")

for ont in ontoterminology:
    print("ontoterminology", ont)

def hv_generator(ontology_id_list):
    try:
        #from app import get_all_descendents #todo: refactor to avoid this circular import
        #load ontotermilology.pkl:
        #todo: check that replacing .pkl with shelve .db works, then delete below load .pkl
        # with open('ontoterminology.pkl', 'rb') as f:
        #     ontoterminology = pickle.load(f)

        # get descendants
        #if should_get_descendents == True:
        #    # print("getting descendents now")
        #    ontology_id_list = get_all_descendents(ontology_id_input)
        #    # print("got them")
        #else:
        #    ontology_id_list = ontology_id_input

        #new method using ontoterminology:         
        mentions = {}
        for selectedID in ontology_id_list:
            if selectedID in ontoterminology.keys():
                mentions[ontoterminology[selectedID]['NAME']] = set(ontoterminology[selectedID]['PMID'])
                print("got one: ", mentions[ontoterminology[selectedID]['NAME']])
            else:
                print("No mentions found for ",selectedID)
        # print("loaded mentions")
        chn_list = []
        for source in mentions:
            for target in mentions: 
                if source.strip() == "" or target.strip() == "":
                    print("blank source or target")
                    pass
                elif source.strip() == target.strip():
                    print("intersection: ", source.strip())
                    pass
                else:
                    intersection = mentions[source].intersection(mentions[target])
                    if len(intersection) > 0: 
                        chn = {"source": source, "target": target, "PMID": len(intersection)}
                        #inverse duplicate checking here: 
                        add_item = True
                        for k in chn_list: #todo: another whole loop here, any faster way to do this? I can't see it.
                            if source + target == k['target'] + k['source']:
                                add_item = False
                        if add_item:
                            chn_list.append(chn)
        print("finished checking for inverse duplicates..")        
        print("length of intersection list: ", len(chn_list))
        print(chn_list)

        # Build the data table expected by the visualisation library
        links = pd.DataFrame.from_dict(chn_list)      
        node_names = links.source.append(links.target)
        node_names = node_names.unique()
        print(node_names)
        node_info = {"index":node_names,"name":node_names,"group":[1]*len(node_names)}
        # node_info = {"index":node_names,"name":node_names,"group":node_names}
        print(node_info)
        nodes = hv.Dataset(pd.DataFrame(node_info), 'index')
        nodes.data.head()

        chord = hv.Chord((links, nodes)).select(value=(0, None)) # value=5 - changing to 0 works for more?

        chord.opts(
            opts.Chord(cmap='Category20', edge_cmap='Category20', edge_color='source',
                    labels='name', node_color='index'))
        renderer = hv.renderer('bokeh')
        hvplot = renderer.get_plot(chord)
        html = renderer.static_html(hvplot)
        return json.dumps(html)
    except Exception as e:
        print(e)
        traceback.print_exc()
        html_error_message = "<!doctype html><div><h4>ERROR CREATING TABLE - no associations found, or possibly some of the ID's were incorrect?</h4></div></html>"
        return(json.dumps(html_error_message))
        