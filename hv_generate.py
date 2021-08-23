import pandas as pd
import holoviews as hv
from holoviews import opts, dim
hv.extension('bokeh')
hv.output(size=200)
# import argparse
import pprint as pp


def hv_generator(ontology_id_input):
    
    df2 = pd.read_csv("static/ontotermmentions.csv",index_col=0)
   
    # print("ontology_id_input type is: ", type(ontology_id_input))
    print("ontology_id_input values: ", ontology_id_input)
    ontology_id_list = ontology_id_input #todo: uncomment, commented for testing - using values below
    #test values which work:
    # ontology_id_list = ["BFO:0000023", "ADDICTO:0000349", "MF:0000016", "ADDICTO:0000632", "ADDICTO:0000904", "ADDICTO:0000491","ADDICTO:0000872" ]
    # ontology_id_list = ["BFO:0000023", "ADDICTO:0000349", "ADDICTO:0000175", "ADDICTO:0000717", "ADDICTO:0000687"]
    #test Ontology IDs to use: 
    # BFO:0000023|role,ADDICTO:0000349|addiction,MF:0000016|human being,ADDICTO:0000632|cannabis use,ADDICTO:0000904|sales,ADDICTO:0000491|opioid agonist treatment,ADDICTO:0000872|Food and Drug Administration
    #test values which don't work:
    # ontology_id_list = ["ADDICTO:0000687"]
    # ontology_id_list = ["ADDICTO:0000175", "ADDICTO:0000717", "ADDICTO:0000687"]
    
    # This creates a table of pairs of terms in the same abstract

    dcp = pd.merge(df2,df2,on="PMID")
    print("about to filter dcp to correct values from ", dcp)
    # We filter the table just to the ones in the ID list we provided as input 
    dcp = dcp.drop(dcp[dcp.LABEL_x == dcp.LABEL_y].index)
    print("dcp after drop Labelx == Labely ", dcp)    
    dcp = dcp.drop(dcp[~dcp.ADDICTOID_y.isin(ontology_id_list)].index) #todo: this one is causing empty dataframe most
    pp.pprint(dcp)
    # dcp2 = dcp
    dcp = dcp.drop(dcp[~dcp.ADDICTOID_x.isin(ontology_id_list)].index)
    # dcp = dcp.drop(dcp2[~dcp2.ADDICTOID_x.isin(ontology_id_list)].index)
    print("dcp after dropping all: ", dcp)
    # We filter the table so that pairs are only represented in one direction, i.e. if we have both (smoking, children) and (children, smoking) for the same PMID we drop the second one
    print("about to drop duplicates")
    # solution 2 from Stack Overflow - replaces iterrows():
    dcp['ADDICTOID'] = dcp[['ADDICTOID_x', 'ADDICTOID_y']].apply(sorted, axis=1).apply(tuple)
    dcp = dcp.drop_duplicates(subset=['ADDICTOID', 'PMID'], keep='first')

    # for index, row in dcp.iterrows():  # THIS IS SLOW
    #     if index % 100 == 0:
    #         print(".",index)
    #     if ((dcp['ADDICTOID_x'] == row['ADDICTOID_y'])
    #         & (dcp['ADDICTOID_y'] == row['ADDICTOID_x'])
    #         & (dcp['PMID'] == row['PMID'])).any():  # Does the inverse of this row exist in the table?
    #         dcp.drop(index, inplace=True)

    # Now we count the distinct numbers of abstracts this combination appeared in
    print("dcp is: ", dcp)
    data_chord_plot_2 = dcp.groupby(['LABEL_x', 'LABEL_y'], as_index=False)[['PMID']].count()
    data_chord_plot_2.columns = ['source','target','value']

    # Build the data table expected by the visualisation library

    links = data_chord_plot_2
    node_names = links.source.append(links.target)
    node_names = node_names.unique()
    node_info = {"index":node_names,"name":node_names,"group":[1]*len(node_names)}

    nodes = hv.Dataset(pd.DataFrame(node_info), 'index')
    nodes.data.head()

    chord = hv.Chord((links, nodes)).select(value=(5, None))

    chord.opts(
        opts.Chord(cmap='Category20', edge_cmap='Category20', edge_color=dim('source').str(),
                labels='name', node_color=dim('index').str()))
    #todo: save error message html if no data returned
    if dcp.empty:
        print('empty dataframe, should create an error message chordout.html here')
        html_error_message = "<!doctype html><div><h4>ERROR CREATING TABLE - possibly some of the ID's were incorrect?</h4></div></html>"
        html_chord_error = open("templates/chordout.html", 'w')
        html_chord_error.write(html_error_message)
    else:
        hv.save(chord, 'templates/chordout.html') 