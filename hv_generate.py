import pandas as pd
import holoviews as hv
from holoviews import opts, dim
hv.extension('bokeh')
hv.renderer('bokeh')
hv.output(size=200)
# from hv.extension.bokeh.embed import json_item
from bokeh.embed import json_item
from bokeh.embed import file_html
# from hv.extension.bokeh.
import pprint as pp
import pyhornedowl
import requests
from urllib.request import urlopen
import json
from memory_profiler import profile
from timeit import default_timer as timer

import io

#todo: memory usage testing: changing to category

@profile
def hv_generator(ontology_id_input, should_get_descendents):
    start_time = timer()
    from app import get_all_descendents #todo: refactor to avoid this circular import
    # load csv:
    df2 = pd.read_csv("static/ontotermmentions.csv", delimiter = ",", index_col=0, engine='c')
    
    # load csv with dtype specified (not working): 
    # df2 = pd.read_csv("static/ontotermmentions.csv", index_col=0, dtype={'ADDICTOID': 'category' ,'LABEL': 'category', 'PMID': 'int64'}) 
    print("read csv into dataframe")
    # change dtype to category on the fly (not working):
    # df2["LABEL"] = df2["LABEL"].astype("category")
    # df2["ADDICTOID"] = df2["ADDICTOID"].astype("category")   

    print("dtypes: ", df2.dtypes)

    # get descendants?
    if should_get_descendents == True:
        ontology_id_list = get_all_descendents(ontology_id_input)
    else:
        ontology_id_list = ontology_id_input 
    #test: 
    # ontology_id_list =  ['ADDICTO:0000386', 'ADDICTO:0000803', 'ADDICTO:0000828', 'ADDICTO:0000678']
    # ontology_id_list = ['ADDICTO:0000311', 'ADDICTO:0000303', 'ADDICTO:0000139', 'ADDICTO:0000138', 'ADDICTO:0000649', 'BCIO:037000', 'ADDICTO:0000941', 'ADDICTO:0000410', 'MF:0000016', 'ADDICTO:0000311', 'ADDICTO:0000303', 'ADDICTO:0000523', 'BCIO:037000', 'ADDICTO:0000352', 'ADDICTO:0000523', 'ADDICTO:0000941', 'ADDICTO:0000410', 'MF:0000016', 'ADDICTO:0000311', 'ADDICTO:0000303', 'ADDICTO:0000311', 'ADDICTO:0000303', 'ADDICTO:0000349', 'IAO:0000088', 'ADDICTO:0000941', 'ADDICTO:0000410', 'MF:0000016', 'ADDICTO:0000311', 'ADDICTO:0000303', 'ADDICTO:0000405', 'ADDICTO:0000836', 'ADDICTO:0000370', 'ADDICTO:0000773', 'ADDICTO:0000405', 'ADDICTO:0000788', 'BCIO:037000', 'ADDICTO:0000139', 'ADDICTO:0000649', 'ADDICTO:0000697', 'BCIO:037000', 'ADDICTO:0000149', 'IAO:0000027', 'BCIO:037000', 'BCIO:037000', 'ADDICTO:0000149', 'ADDICTO:0000941', 'ADDICTO:0000410', 'MF:0000016', 'BCIO:037000', 'BCIO:037000', 'BCIO:037000', 'BCIO:037000', 'CHEBI:18723', 'ADDICTO:0000263', 'CHEBI:3219', 'BCIO:037000', 'ADDICTO:0000147', 'BCIO:037000', 'BCIO:037000', 'ADDICTO:0000147', 'CHEBI:18723', 'ADDICTO:0000262', 'CHEBI:3219', 'http://purl.obolibrary.org/obo/OAE_0000001', 'ADDICTO:0000500', 'http://purl.obolibrary.org/obo/OAE_0000001', 'BCIO:037000', 'BCIO:037000', 'ADDICTO:0000139', 'ADDICTO:0000649', 'ADDICTO:0000697', 'ADDICTO:0000752', 'ADDICTO:0000941', 'ADDICTO:0000410', 'MF:0000016', 'BCIO:037000', 'BCIO:037000', 'BCIO:037000', 'ADDICTO:0000370', 'ADDICTO:0000773', 'ADDICTO:0000405', 'ADDICTO:0000788']
    # print("ontology_id_list is now: ", ontology_id_list)
        # print("get_descendents is: ", should_get_descendents)
    #todo: trying filtering by ontology_id_list before merge - works
    df2 = df2.drop(df2[~df2.ADDICTOID.isin(ontology_id_list)].index)
    # This creates a table of pairs of terms in the same abstract
    # print("df2 after drop: ", df2)

    dcp = pd.merge(df2,df2,on="PMID",how="inner")
    # print("dcp after merge: ", dcp)
    # change nan values to "":  - doesn't seem to help any
    # dcp['LABEL_y'].fillna('', inplace=True) 
    # dcp['LABEL_x'].fillna('', inplace=True) 

    dcp = dcp.drop(dcp[dcp.LABEL_x == dcp.LABEL_y].index) 
    # print("dcp after drop Label_x == Label_y ", dcp)

    # We filter the table just to the ones in the ID list we provided as input 
    # print("about to filter dcp to correct values from ", dcp)
    
    # We filter the table so that pairs are only represented in one direction, i.e. if we have both (smoking, children) and (children, smoking) for the same PMID we drop the second one
    # print("about to drop duplicates")

    # solution 2 from Stack Overflow - replaces iterrows():
    dcp['ADDICTOID'] = dcp[['ADDICTOID_x', 'ADDICTOID_y']].apply(sorted, axis=1).apply(tuple)
    dcp = dcp.drop_duplicates(subset=['ADDICTOID', 'PMID'], keep='first')
    # drop "ADDICTOID" column - this column is not needed anymore:  
    dcp = dcp.drop(['ADDICTOID'], axis=1)

    # print("final dcp is: ", dcp)

    # Now we count the distinct numbers of abstracts this combination appeared in    
    data_chord_plot = dcp.groupby(['LABEL_x', 'LABEL_y'], as_index=False)[['PMID']].count()
    data_chord_plot.columns = ['source','target','value']

    # print("Final chord plot: ", data_chord_plot)
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
    #error message html if no chord plot to show:
    if dcp.empty:
        print('empty dataframe, should create an error message chordout.html here')
        html_error_message = "<!doctype html><div><h4>ERROR CREATING TABLE - no associations found, or possibly some of the ID's were incorrect?</h4></div></html>"
        return(json.dumps(html_error_message))
    else:
        renderer = hv.renderer('bokeh')
        hvplot = renderer.get_plot(chord)
        html = renderer.static_html(hvplot)
        return json.dumps(html)
        