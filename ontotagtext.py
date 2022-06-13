
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
        # print("self.ontols: ", self.ontols)
        # for x in self.ontols: 
        #     print("got x: ", x)
        # print("all_labels = ", self.all_labels)
        
        # for making plural forms of labels for text matching
        engine = inflect.engine()

        # init terms and patterns
        self.terms = {}
        patterns = []

        #build unified table of all ID, IRI, Label and Synonyms:
        for ontol in self.ontols: #should be all ontols in 
            # print("checking ontol: ", ontol)
            for termid in ontol.get_classes():
                # print("k is: ", k)
                termshortid = ontol.get_id_for_iri(termid)    
                        
                label = ontol.get_annotation(termid, RDFSLABEL)
                definition = ontol.get_annotation(termid, DEFINITION)
                if label:                    
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
                    # print("adding SYNONYM in ontotagtext: ", s)
                    if s.strip().lower() not in stopwords:
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
        # print("getting term")
        if term_id in [ v['id'] for v in self.terms.values()]: 
            keys = [k for k, v in self.terms.items() if v['id'].strip() == term_id.strip()]
            return self.terms[keys[0]]
        else:
            return None             

    def get_label(self, label): 
        # print("getting label")
        if label.strip().lower() in [ v['name'].strip().lower() for v in self.terms.values()]: 
            keys = [k for k, v in self.terms.items() if v['name'].strip().lower() == label.strip().lower()]
            return self.terms[keys[0]]
        else:
            return None

