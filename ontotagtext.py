
from pronto import Ontology
import spacy
from spacy import displacy
from spacy.tokens import Doc, Span, Token
from spacy.lang.en import English
from spacy.matcher import PhraseMatcher
from spacy.util import filter_spans

class ExtractorComponent(object):
    def __init__(self, nlp, name, label, ontologyfile):
        # label that is applied to the matches
        self.label = label
        self.name = name

        # load ontology
        print("Loading ontology")
        ontol = Ontology(ontologyfile)

        # init terms and patterns
        self.terms = {}
        patterns = []

        i = 0
        nr_terms = len(ontol.terms())
        # init progress bar as loading terms takes long
        print("Importing terms")

        # iterate over terms in ontology
        for term in ontol.terms():
          # if term has a name
          if term.name is not None:
            self.terms[term.name.lower()] = {'id': term.id}
            patterns.append(nlp(term.name))
          for s in term.synonyms:
            self.terms[s.description.lower()] = {'id': term.id}
            patterns.append(nlp(s.description))
          extra_synonyms = [ s.literal for s in term.annotations if s.property=="IAO:0000118"]
          for s in extra_synonyms:
            self.terms[s.lower()] = {'id': term.id}
            patterns.append(nlp(s))
          i += 1

        # initialize matcher and add patterns
        self.matcher = PhraseMatcher(nlp.vocab, attr='LOWER')
        self.matcher.add(label, None, *patterns)

        # set extensions to tokens, spans and docs
        Token.set_extension("is_ontol_term", default=False, force=True)
        Token.set_extension("ontol_id", default=False, force=True)
        Token.set_extension("merged_concept", default=False, force=True)
        Doc.set_extension("has_ontols", getter=self.has_ontols, force=True)
        Doc.set_extension("ontols", default=[], force=True)
        Span.set_extension("has_ontols", getter=self.has_ontols, force=True)

    def __call__(self, doc):
        matches = self.matcher(doc)
        spans = [Span(doc, match[1], match[2], label=self.label) for match in matches]
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
    # setter function for doc level
    def has_ontols(self, tokens):
        return any([t._.get("is_ontol_term") for t in tokens])



# Testing
if __name__ == "__main__":

    # python -m spacy download en_core_web_md
    # or: en_core_web_sm or en_core_web_lg
    nlp = spacy.load('en_core_web_md')

    onto_extractor = ExtractorComponent(
        nlp,
        name="ADDICTO",
        label="ADDICTO",
        ontologyfile="/Users/hastingj/Work/Onto/addiction-ontology/addicto.obo")
    nlp.add_pipe(onto_extractor, after="ner")

    test = '''
    The promotion of the London Smoking Cessation Transformation Programme during September 2017 was associated with a significant increase in quit attempts compared with the rest of England. The results were inconclusive regarding an effect on quit success among those who tried.
    '''

    doc = nlp(test)

    # print ontology IDs identified
    for token in doc:
        if token._.is_ontol_term:
            print(token._.ontol_id, token.text, token.idx)
