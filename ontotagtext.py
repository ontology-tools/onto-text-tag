
from pronto import Ontology
import spacy
from spacy import displacy
from spacy.tokens import Doc, Span, Token
from spacy.lang.en import English
from spacy.matcher import PhraseMatcher
from spacy.util import filter_spans
import inflect

class ExtractorComponent(object):
    def __init__(self, nlp, name, label, ontologyfile):

        # label that is applied to the matches
        self.label = label
        self.name = name
        # stop words, don't try to match these
        stopwords = nlp.Defaults.stop_words
        stopwords.add("ands")
        stopwords.add("ends")

        # load ontology
        print("Loading ontology")
        self.ontol = Ontology(ontologyfile)
        self.ontol_ids = [t.id for t in self.ontol.terms()]

        # for making plural forms of labels for text matching
        engine = inflect.engine()

        # init terms and patterns
        self.terms = {}
        patterns = []

        i = 0
        nr_terms = len(self.ontol.terms())
        print(f"Importing {nr_terms} terms")

        # iterate over terms in ontology
        for term in self.ontol.terms():
          # if term has a name
          if term.name is not None and term.name.strip().lower() not in stopwords:
            self.terms[term.name.strip().lower()] = {'id': term.id}
            patterns.append(nlp.make_doc(term.name.strip()))
            plural = engine.plural(term.name.strip())
            self.terms[plural.lower()] = {'id': term.id}
            patterns.append(nlp.make_doc(plural))
          for s in term.synonyms:
              if s.description.strip().lower() not in stopwords:
                self.terms[s.description.strip().lower()] = {'id': term.id}
                patterns.append(nlp.make_doc(s.description.strip()))
                plural = engine.plural(s.description.strip())
                self.terms[plural.lower()] = {'id': term.id}
                patterns.append(nlp.make_doc(plural))
          extra_synonyms = [ s.literal for s in term.annotations if s.property=="IAO:0000118"]
          for s in extra_synonyms:
              if s.strip().lower() not in stopwords:
                self.terms[s.strip().lower()] = {'id': term.id}
                patterns.append(nlp.make_doc(s.strip()))
                plural = engine.plural(s.strip())
                self.terms[plural.lower()] = {'id': term.id}
                patterns.append(nlp.make_doc(plural))
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

    # getter function for doc level
    def has_ontols(self, tokens):
        return any([t._.get("is_ontol_term") for t in tokens])

    def get_term(self, term_id):
        if term_id in self.ontol_ids:
            return self.ontol.get_term(term_id)
        else:
            return None


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
    nlp.add_pipe(onto_extractor, after="parser")

    test = '''
    The promotion of the London Smoking Cessation Transformation Programme during September 2017 was associated with a significant increase in quit attempts compared with the rest of England. The results were inconclusive regarding an effect on quit success among those who tried.
    Smokers were more successful than non-smokers and this was good. Those who smoked were associated with ...
    '''

    doc = nlp(test)

    # print ontology IDs identified
    for token in doc:
        if token._.is_ontol_term:
            print(token._.ontol_id, token.text, token.idx)


    from spacy import displacy

    svg = displacy.serve(doc, style='ent')

    from spacy.matcher import PhraseMatcher

    matcher = PhraseMatcher(nlp.vocab, attr='LEMMA')

    terms = ['smoking cessation', 'smoker', 'non-smoker', 'smoke', 'effect']
    patterns = [nlp(text) for text in terms]
    matcher.add("TerminologyList", patterns)

    matches = matcher(doc)
    print(matches)

    for (match_id, start, end) in matches:
        print(nlp.vocab.strings[match_id], doc[start:end])


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


