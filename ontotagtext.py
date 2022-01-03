
import pyhornedowl
import spacy
from spacy import displacy
from spacy.tokens import Doc, Span, Token
from spacy.lang.en import English
from spacy.matcher import PhraseMatcher
from spacy.util import filter_spans
import inflect

RDFSLABEL = "http://www.w3.org/2000/01/rdf-schema#label"
SYN = "http://purl.obolibrary.org/obo/IAO_0000118"
DEFINITION = "http://purl.obolibrary.org/obo/IAO_0000115"
PREFIXES = [ ["ADDICTO","http://addictovocab.org/ADDICTO_"],
             ["BFO","http://purl.obolibrary.org/obo/BFO_"],
             ["CHEBI","http://purl.obolibrary.org/obo/CHEBI_"],
             ["UBERON","http://purl.obolibrary.org/obo/UBERON_"],
             ["PATO","http://purl.obolibrary.org/obo/PATO_"],
             ["BCIO","http://humanbehaviourchange.org/ontology/BCIO_"],
             ["SEPIO","http://purl.obolibrary.org/obo/SEPIO_"],
             ["OMRSE","http://purl.obolibrary.org/obo/OMRSE_"],
             ["OBCS","http://purl.obolibrary.org/obo/OBCS_"],
             ["OGMS","http://purl.obolibrary.org/obo/OGMS_"],
             ["ENVO","http://purl.obolibrary.org/obo/ENVO_"],
             ["OBI", "http://purl.obolibrary.org/obo/OBI_"],
             ["MFOEM","http://purl.obolibrary.org/obo/MFOEM_"],
             ["MF","http://purl.obolibrary.org/obo/MF_"],
             ["CHMO","http://purl.obolibrary.org/obo/CHMO_"],
             ["DOID","http://purl.obolibrary.org/obo/DOID_"],
             ["IAO","http://purl.obolibrary.org/obo/IAO_"],
             ["ERO","http://purl.obolibrary.org/obo/ERO_"],
             ["PO","http://purl.obolibrary.org/obo/PO_"],
             ["RO","http://purl.obolibrary.org/obo/RO_"],
             ["APOLLO_SV","http://purl.obolibrary.org/obo/APOLLO_SV_"],
             ["PDRO","http://purl.obolibrary.org/obo/PDRO_"],
             ["GAZ","http://purl.obolibrary.org/obo/GAZ_"],
             ["GSSO","http://purl.obolibrary.org/obo/GSSO_"]
           ]



class MultiExtractorComponent(object):
    def __init__(self, nlp, ontoDict):
        # add ontology and label from ontoDict
        self.ontoDict = ontoDict
        self.all_labels = ""

        # stop words, don't try to match these
        stopwords = nlp.Defaults.stop_words
        stopwords.add("ands")
        stopwords.add("ends")
        stopwords.add("ci")

        self.ontols = []

        ontologies = ontoDict["ontologies"]
        for ontology in ontologies:
            for key, value in ontology.items():
                if(key == "label"):
                    self.all_labels = self.all_labels + value
                if (key == "ontology"):
                    self.ontols.append(value)

        # print("all_labels = ", self.all_labels)
        
        # for making plural forms of labels for text matching
        engine = inflect.engine()

        # init terms and patterns
        self.terms = {}
        patterns = []

        #todo: get this to work with "ontols":
        # nr_terms = len(self.ontol2.get_classes())+len(self.ontol.get_classes())
        # print(f"Importing {nr_terms} terms")

        #build unified table of all ID, IRI, Label and Synonyms:
        for k, ontol in [self.ontols]: #should be all ontols in 
            # print("checking ontol: ", ontol)
            for termid in ontol.get_classes():
                # print("k is: ", k)
                termshortid = ontol.get_id_for_iri(termid)            
                label = ontol.get_annotation(termid, RDFSLABEL)
                definition = ontol.get_annotation(termid, DEFINITION)
                if label: 
                    if label.strip().lower() == "bupropion":
                            print("got bupropion")
                    if label.strip().lower() == "intervention":
                            print("got intervention")                   
                    term_entry = {'id': termid if termshortid is None else termshortid,
                                'name': label.strip(),
                                'definition': definition}
                if label is not None and label.strip().lower() not in stopwords:
                    self.terms[label.strip().lower()] = term_entry
                    patterns.append(nlp.make_doc(label.strip().lower()))
                    plural = engine.plural(label.strip())
                    self.terms[plural.lower()] = term_entry
                    patterns.append(nlp.make_doc(plural.lower()))
                synonyms = ontol.get_annotations(termid, SYN)
                for s in synonyms:
                    if s.strip().lower() not in stopwords:
                        if s.strip().lower() == "tobacco":
                            print("got tobacco")
                        if s.strip().lower() == "intervention":
                            print("got intervention")
                        self.terms[s.strip().lower()] = term_entry
                        patterns.append(nlp.make_doc(s.strip().lower()))
                        try:
                            plural = engine.plural(s.strip().lower())
                            self.terms[plural.lower()] = term_entry
                            patterns.append(nlp.make_doc(plural.lower()))
                        except:
                            print("Problem getting plural of ",s)
                            continue
        
        # initialize matcher and add patterns
        self.matcher = PhraseMatcher(nlp.vocab, attr='LOWER')        
        self.matcher.add(self.all_labels, None, *patterns)

        # set extensions to tokens, spans and docs
        Token.set_extension("is_ontol_term", default=False, force=True)
        Token.set_extension("ontol_id", default=False, force=True)
        Token.set_extension("merged_concept", default=False, force=True)
        Doc.set_extension("has_ontols", getter=self.has_ontols, force=True)
        Doc.set_extension("ontols", default=[], force=True)
        Span.set_extension("has_ontols", getter=self.has_ontols, force=True)

    def __call__(self, doc):
        matches = self.matcher(doc)
        spans = [Span(doc, match[1], match[2], label=self.all_labels) for match in matches]
        for i, span in enumerate(spans):
          span._.set("has_ontols", True)
          for token in span:
                if span.text.lower() in self.terms:
                    token._.set("is_ontol_term", True)
                    token._.set("ontol_id", self.terms[span.text.lower()]["id"])
                else:
                    print("Term not found: ",span.text.lower())

        with doc.retokenize() as retokenizer:
            for span in filter_spans(spans):
                retokenizer.merge(span, attrs={"_": {"merged_concept": True}})
                doc._.ontols = list(doc._.ontols) + [span]

        return doc

    # getter function for doc level
    def has_ontols(self, tokens):
        return any([t._.get("is_ontol_term") for t in tokens])

    def get_term(self, term_id): 
        if term_id.strip() in [ v['id'] for v in self.terms.values()]: 
            keys = [k for k, v in self.terms.items() if v['id'] == term_id]
            return self.terms[keys[0]]
        else:
            return None             

    def get_label(self, label): 
        if label in [ v['name'] for v in self.terms.values()]: 
            keys = [k for k, v in self.terms.items() if v['name'].strip().lower() == label.strip().lower()]
            return self.terms[keys[0]]
        else:
            return None

# Testing
#if __name__ == "__main__":
#    import requests
#    from urllib.request import urlopen
#    location = f"https://raw.githubusercontent.com/addicto-org/addiction-ontology/master/addicto-merged.owx"
#    location2 = f"https://raw.githubusercontent.com/HumanBehaviourChangeProject/ontologies/master/Upper%20Level%20BCIO/bcio-merged.owx"
#    # location = f"https://raw.githubusercontent.com/HumanBehaviourChangeProject/ontologies/master/Upper%20Level%20BCIO/bcio-merged.owx"

#    print("Fetching release file from", location)
#    data = urlopen(location).read()  # bytes
#    print("Fetching release file from", location2)
#    data2 = urlopen(location2).read()  # bytes

#    ontofile1 = data.decode('utf-8')
#    ontofile2 = data2.decode('utf-8')


#    ontoDict = {
#        "ontologies": [
#            {
#                "label": "BCIO",
#                "name": "BCIO",
#                "ontologyfile": ontofile2, #todo: why ontofile2 not working here if BCIO added after AddictO?
#            },
#            {
#                "label": "AddictO",
#                "name": "AddictO",
#                "ontologyfile": ontofile1,
#            },
#        ]
#    }
    # python -m spacy download en_core_web_md
    # or: en_core_web_sm or en_core_web_lg
#    nlp = spacy.load('en_core_web_md')

#    onto_extractor = MultiExtractorComponent(
#        nlp,
#        ontoDict)
#        # name="ADDICTO",
#        # label="ADDICTO",
#        # ontologyfile="/home/tom/Documents/PROGRAMMING/Python/addiction-ontology/addicto-merged.owx")
#    nlp.add_pipe(onto_extractor, after="parser")

    # test = '''
    # The promotion of the London Smoking Cessation Transformation Programme during September 2017 was associated with a significant increase in quit attempts compared with the rest of England. The results were inconclusive regarding an effect on quit success among those who tried.
    # Smokers were more successful than non-smokers and this was good. Those who smoked were associated with ...
    # '''

    # doc = nlp(test)

    # # print ontology IDs identified
    # for token in doc:
    #     if token._.is_ontol_term:
    #         print(token._.ontol_id, token.text, token.idx)


    # from spacy import displacy

    # svg = displacy.serve(doc, style='ent')

    # from spacy.matcher import PhraseMatcher

    # matcher = PhraseMatcher(nlp.vocab, attr='LEMMA')

    # terms = ['smoking cessation', 'smoker', 'non-smoker', 'smoke', 'effect']
    # patterns = [nlp(text) for text in terms]
    # matcher.add("TerminologyList", patterns)

    # matches = matcher(doc)
    # print(matches)

    # for (match_id, start, end) in matches:
    #     print(nlp.vocab.strings[match_id], doc[start:end])


    ## Another approach using spacylookup

    #nlp = English()
    #from spacy_lookup import Entity


    #addicto = Entity(keywords_list=['smoking cessation','smoker','effect'])
    #nlp.add_pipe(addicto, last=False)

    #doc = nlp(test)

    #print([(token.text, token._.canonical) for token in doc if token._.is_entity])


    #svg = spacy.displacy.render(doc, style="dep")
    #output_path = Path(os.path.join("./", "sentence.svg"))
    #output_path.open('w', encoding="utf-8").write(svg)


    # Test using EntityRuler rather than PatternMatcher
#    from spacy.lang.en import English
#    from spacy.pipeline import EntityRuler

#    nlp = English()
#    ruler = EntityRuler(nlp)
#    patterns = [{"label": "ORG", "pattern": "Apple", "id": "apple"},
#                {"label": "GPE", "pattern": [{"LOWER": "san"}, {"LOWER": "francisco"}],
#                 "id": "san-francisco"},
#                {"label": "GPE", "pattern": [{"LOWER": "san"}, {"LOWER": "fran"}],
#                 "id": "san-francisco"}]
#    ruler.add_patterns(patterns)
#    nlp.add_pipe(ruler)

#    doc1 = nlp("Apple is opening its first big office in San Francisco.")
#    print([(ent.text, ent.label_, ent.ent_id_) for ent in doc1.ents])

#    doc2 = nlp("Apple is opening its first big office in San Fran.")
#    print([(ent.text, ent.label_, ent.ent_id_) for ent in doc2.ents])

#    '''
#    {"label": "ORG", "pattern": "Apple"}
#    {"label": "GPE", "pattern": [{"LOWER": "san"}, {"LOWER": "francisco"}]}
#    '''

#    ruler.to_disk("./patterns.jsonl")
#    new_ruler = EntityRuler(nlp).from_disk("./patterns.jsonl")


