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
import traceback
import shelve 


ontoterminology = shelve.open('static/ontoterminology.db', flag='r', writeback=False) #works

print("loaded terms db")


def hv_generator(ontology_id_list):
    try:      
        mentions = {}
        for selectedID in ontology_id_list:
            if selectedID in ontoterminology.keys():
                # print("ontoterminology selectedID 'NAME'", ontoterminology[selectedID]['NAME'])
                # print("set:  ", set(ontoterminology[selectedID]['PMID']))
                mentions[ontoterminology[selectedID]['NAME']] = set(ontoterminology[selectedID]['PMID'])
                # print("got one: ", mentions[ontoterminology[selectedID]['NAME']])
            else:
                # print("No mentions found for ",selectedID)
                pass
        # print("loaded mentions", mentions)
        chn_list = []
        for source in mentions:
            # print("plain source: ", source)
            for target in mentions: 
                if source.strip() == "" or target.strip() == "":
                    # print("blank source or target")
                    pass
                elif source.strip() == target.strip():
                    # print("intersection: ", source.strip())
                    pass
                else:
                    intersection = mentions[source].intersection(mentions[target])
                    if len(intersection) > 0: 
                        chn = {"source": source, "target": target, "PMID": len(intersection)}
                        #inverse duplicate checking here: 
                        add_item = True
                        for k in chn_list: 
                            if source + target == k['target'] + k['source']:
                                add_item = False
                        if add_item:
                            chn_list.append(chn)
        print("finished checking for inverse duplicates..")        
        # print("length of intersection list: ", len(chn_list))
        # print(chn_list)

        # Build the data table expected by the visualisation library
        links = pd.DataFrame.from_dict(chn_list)      
        node_names = links.source.append(links.target)
        node_names = node_names.unique()
        # print(node_names)
        node_info = {"index":node_names,"name":node_names,"group":[1]*len(node_names)}
        # print(node_info)
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
        