'''
Created on Jul 12, 2012, amended April 2015

Module for managing Gene Ontology data, in particular gene:terms
annotations and term definitions

It can be used on files you can download from geneontology.org.
The class GO is constructed from:
- annotation file which is (usually) specific to the species of interest
- OBO file which defines the GO terms and their relationships
  e.g.
  > go = GO('gene_association.goa_ref_human', 'go-basic.obo')
Internal data structures are created so that you can query
- what are the terms of my gene (or genes)? Use getTerms
- what are the genes of my term? Use getGenes
- what terms occur amongst my genes, ranked by their absolute count? Use getGOReport without background
- what terms are statistically enriched in my genes, relative a background set of genes? Use getGOReport with background

The class BinGO works with a compact (memory saving) binary format that aggregates information from an annotation
file and an OBO file. Therefore, you first need to construct this binary file, using writeBitFile.
Subsequently you can construct instances of BinGO and query terms and genes, roughly in the manner identified above for GO.

@author: mikael
'''

from struct import pack, unpack, calcsize, error
import operator
import time
import os
import stats

# Character codes used by binary format to identify ontology
onto_codes = {
              'P': 'Biological process',
              'F': 'Molecular function',
              'C': 'Cellular component'}

# Labels for edges in the ontology graph, index is used in binary format
onto_rel = ['is_a', 'isect', 'part_of', 'has_part', 'regulates']

# Evidence codes assigned to annotations, an index is assigned when creating binary file and is stored in its header
evid_codes = { # Experimental Evidence Codes
    'EXP': 'Inferred from Experiment',
    'IDA': 'Inferred from Direct Assay',
    'IPI': 'Inferred from Physical Interaction',
    'IMP': 'Inferred from Mutant Phenotype',
    'IGI': 'Inferred from Genetic Interaction',
    'IEP': 'Inferred from Expression Pattern',
               #Computational Analysis Evidence Codes
    'ISS': 'Inferred from Sequence or Structural Similarity',
    'ISO': 'Inferred from Sequence Orthology',
    'ISA': 'Inferred from Sequence Alignment',
    'ISM': 'Inferred from Sequence Model',
    'IGC': 'Inferred from Genomic Context',
    'IBA': 'Inferred from Biological aspect of Ancestor',
    'IBD': 'Inferred from Biological aspect of Descendant',
    'IKR': 'Inferred from Key Residues',
    'IRD': 'Inferred from Rapid Divergence',
    'RCA': 'inferred from Reviewed Computational Analysis',
    'TAS': 'Traceable Author Statement',
    'NAS': 'Non-traceable Author Statement',
               #Curator Statement Evidence Codes
    'IC': 'Inferred by Curator',
    'ND': 'No biological Data available',
               #Automatically-assigned Evidence Codes
    'IEA': 'Inferred from Electronic Annotation',
               #Obsolete Evidence Codes
    'NR': 'Not Recorded'}

class GO():
    """ Classical interface for working with GO terms usually within the same species and when memory is not a major issue.
        Implementations are relatively efficient (for Python at least).
        Major functions:
        __init__: construct instance of GO session from an annotation file and an OBO file (geneontology.org)
        getTerms: get GO terms from gene or genes (transitively or not)
        getGenes: get genes that are annotated with given term or terms
        getGOReport: perform basic gene set enrichment
        """

    # Structures to hold all data relevant to session
    annots = {}     # annotations: annots[gene] = (taxa, terms[term] = (evid, T/F))
    termdefs = {}   # definitions: termdefs[term] = (onto, set((term, rel)), name)
    children = {}   # redundant, parent-to-child structure: children[term] = set((term, rel))

    def __init__(self, annotFile, obofile, annotfile_columns = (1,2,3,4,6,8)):
        """ Start GO session with specified data loaded:
        annotfile: name of annotation file, e.g.'gene_association.tair'
        OBO file: name of gene ontology definition file, e.g. 'gene_ontology_ext.obo'
        Optionally, specify what columns in the annotation file that contains in order:
        gene, symb, qual, term, evid, onto. Note that index starts at 0 NOT 1.
        (The default seems to work for most annotation files, but sometime if you wish to cross reference
        say gene names, you need to point to an alternate column, e.g. 9 for TAIR's A. thaliana annotations:
        go = GO('gene_association.tair', 'gene_ontology_ext.obo', (9,2,3,4,6,8))
        """
        print(("Started at", time.asctime()))
        # Get GO definitions
        terms = readOBOFile(obofile)
        for term in terms:
            (term_name, term_onto, term_is) = terms[term]
            self.termdefs[term] = (term_onto, term_is, term_name)
            self.children[term] = set()

        for term in self.termdefs:
            (term_onto, term_is, term_name) = self.termdefs[term]
            for (parent, prel) in term_is:
                try:
                    cset = self.children[parent]
                    cset.add((term, prel))
                except KeyError:
                    pass
        print(("Read %d GO definitions" % len(terms)))
        # open annotation file to analyse and index data
        src = open(annotFile, 'r')
        gene_cnt = 0
        cnt = 0
        for line in src:
            cnt += 1
            if line.startswith('!'):
                continue
            (gene, symb, qual, term, evid, onto, taxa) = _extractAnnotFields(line, annotfile_columns)
            #print(gene, symb, qual, term, evid, onto, taxa)
            try:
                (taxa_q, terms_map) = self.annots[gene]
                terms_map[term] = (evid, qual != 'NOT')
            except KeyError: # not a previously encountered gene
                gene_cnt += 1
                terms_map = {term: (evid, qual != 'NOT')}
                self.annots[gene] = (taxa, terms_map)
        src.close()
        print(("Read annotations for %d genes" % gene_cnt))

    def _makeIntoList(self, id_or_ids):
        if type(id_or_ids) != list and type(id_or_ids) != set and type(id_or_ids) != tuple:
            return [id_or_ids]
        return id_or_ids

    def getTerms(self, genes_or_gene, evid = None, onto = None, include_more_general = True):
        """ Retrieve all terms for a gene or a set/list/tuple of genes.
            If evid(ence) is specified the method returns only entries with that specific evidence code (see header of file for codes).
            If onto(logy) is specified the method includes only entries from specified ontology ('P', 'F' or 'C').
            If include_more_general is true, terms that are transitively related are included.
            With multiple genes provided in query, the result is a map, keyed by gene (each identifying a set of terms).
            When only one gene is provided, the result is simply a set of terms.
        """
        if type(genes_or_gene) != list and type(genes_or_gene) != set and type(genes_or_gene) != tuple:
            return self.getTerms4Gene(genes_or_gene, evid, onto, include_more_general)
        else:
            return self.getTerms4Genes(genes_or_gene, evid, onto, include_more_general)

    def getTerms4Genes(self, genes, evid = None, onto = None, include_more_general = True):
        """ Retrieve all GO terms for a given set/list/tuple of genes.
            If evid(ence) is specified the method returns only entries with that specific evidence code (see header of file for codes).
            If onto(logy) is specified the method includes only entries from specified ontology ('P', 'F' or 'C').
            If include_more_general is True (default) then transitively related terms are included.
            With multiple genes provided in query, the result is a map, keyed by gene (each identifying a set of terms).
        """
        gomap = {} # gene to GO terms map
        genes = self._makeIntoList(genes)
        for gene in genes:
            gomap[gene] = self.getTerms4Gene(gene, evid, onto, include_more_general)
        return gomap

    def getTerms4Gene(self, gene, evid = None, onto = None, include_more_general = True):
        """ Retrieve all GO terms for a given (single) gene.
            If evid(ence) is specified the method returns only entries with that specific evidence code (see header of file for codes).
            If onto(logy) is specified the method includes only entries from specified ontology ('P', 'F' or 'C').
            If include_more_general is True (default) then transitively related terms are included
            When only one gene is provided, the result is simply a set of terms.
        """
        direct = set()
        # STEP 1: Find all terms directly associated with specified genes
        try:
            (taxa, terms_map) = self.annots[gene]
            for term in terms_map:
                (term_evid, term_qual) = terms_map[term]
                if (evid == None or evid == term_evid) and term_qual:
                    direct.add(term)
        except KeyError:
            return set() # gene was not found, hence no annotations for it
        # STEP 2: Find terms associated with (indirect) parents of terms from STEP 1
        indirect = set()
        if include_more_general:
            for term in direct:
                parents = self.getParents(term, include_more_general)
                for parent in parents:
                    indirect.add(parent)
        return direct.union(indirect)

    def getGenes(self, terms_or_term, evid = None, taxa = None, rel = None, include_more_specific = False):
        """ Retrieve all genes that are annotated with specified term or terms,
            qualified by evidence, taxa and relation type, e.g. "is_a".
            If multiple terms are provided, a map is returned keyed by term (each identifying set of genes).
            With a single term provided, a set of genes is returned.
        """
        if type(terms_or_term) != list and type(terms_or_term) != set and type(terms_or_term) != tuple:
            return self.getGenes4Term(terms_or_term, evid, taxa, rel, include_more_specific)
        else:
            return self.getGenes4Terms(terms_or_term, evid, taxa, rel, include_more_specific)

    def getGenes4Terms(self, terms, evid = None, taxa = None, rel = None, include_more_specific = False):
        """ Retrieve all genes that are annotated with specified terms,
            qualified by evidence, taxa and relation type, e.g. "is_a".
            Since multiple terms are provided, a map is returned keyed by term (each identifying set of genes).
        """
        gomap = {} # term to genes map
        terms = self._makeIntoList(terms)
        for term in terms:
            gomap[term] = self.getGenes4Term(term, evid, taxa, rel, include_more_specific)
        return gomap

    def getGenes4Term(self, term, evid = None, taxa = None, rel = None, include_more_specific = False):
        """ Retrieve all genes that are annotated with specified term or terms,
            qualified by evidence, taxa and relation type, e.g. "is_a".
            With a single term provided, a set of genes is returned.
        """
        genes = self._getGenes4Term(term, evid, taxa, rel)
        if include_more_specific:
            terms = self.getChildren(term, rel, True) # not recursive yet
            for t in terms:
                tgenes = self._getGenes4Term(t, evid, taxa, rel)
                for g in tgenes:
                    genes.add(g)
        return genes

    def _getGenes4Term(self, term, evid = None, taxa = None, rel = None):
        """ Retrieve all genes that are annotated with specified term, and qualified by evidence, taxa etc. """
        genes = set()
        # Scour through all genes
        for gene in self.annots: # annotations: annots[gene] = (taxa, terms[term] = (evid, T/F))
            (qtaxa, qterms) = self.annots[gene]
            if taxa == None or taxa == qtaxa:
                for qterm in qterms:
                    if qterm != term:
                        continue
                    (qevid, qqual) = qterms[term]
                    if (evid == None or evid == qevid) and qqual:
                        genes.add(gene)
                        break
        return genes

    def getChildren(self, parent_term_id_or_ids, rel = None, include_more_specific = False):
        """ Retrieve all direct children of the given (parent) term.
        """
        parent_terms = self._makeIntoList(parent_term_id_or_ids)
        cset = set()
        for parent in parent_terms:
            # definitions: children[term] = set((term, relation), ...)
            current = self.children[parent]
            for (child_term, child_rel) in current:
                if rel == None or rel == child_rel:
                    cset.add(child_term)
        if len(cset) > 0 and include_more_specific:
            grandkids = self.getChildren(cset, rel, True)
            for grandkid in grandkids:
                cset.add(grandkid)
        return cset

    def getParents(self, child_term_id, include_more_general = True):
        """ Retrieve all parents of the given term, transitively or not.
        """
        direct = set() # all GO terms which are parents to given term
        try:
            (onto_ch, terms_ch, name_ch) = self.termdefs[child_term_id]
            for (parent_id, parent_rel) in terms_ch:
                (onto_pa, terms_pa, name_pa) = self.termdefs[parent_id]
                direct.add(parent_id)
                if (include_more_general):
                    parents = self.getParents(parent_id, True)
                    for parent in parents:
                        direct.add(parent)
        except KeyError:
            pass # term was not found, possibly throw error?
        return direct

    def getTermdef(self, term_id):
        """ Retrieve information about a given term:
            ontology, parent terms, and name as a tuple.
        """
        try:
            (onto_ch, terms_set, term_name) = self.termdefs[term_id]
            return (onto_ch, terms_set, term_name)
        except KeyError:
            return ('Unknown', 'Unknown', 'Unknown')

    def getAllAnnots(self):
        """ Retrieve all annotated gene products """
        return list(self.annots.keys())

    def getAllBackground(self, positives = [], taxa = None, evid = None, include_more_general = False):
        """ Retrieve all genes and terms that are annotated but not in a list of positives (gene products).
        """
        # (taxa, terms[term] = (evid, T/F))
        bg_genes = set()
        bg_list = []
        for gene in self.annots:
            if not gene in positives:
                bg_genes.add(gene)
                (qtaxa, qterms) = self.annots[gene]
                if taxa == None or qtaxa == taxa:
                    for t in qterms:
                        (qevid, qqual) = qterms[t]
                        if (evid == None or qevid == evid) and qqual:
                            bg_list.append(t)
                            if include_more_general:
                                for parent in self.getParents(t, True):
                                    bg_list.append(parent)
        return (bg_genes, bg_list)

    def getCountReport(self, positives, threshold = None, include_more_general = True):
        """ For a set of named gene products (positives) this method determines the counts of GO terms.
            Returns a list of tuples (GO_Term_ID[str], Foreground_no[int], Term_description[str]) sorted by count.
            positives: names of gene products
            threshold: the count that must be reached for term to be reported (default is 0)
            If evid(ence) is specified the method returns only entries with that specific evidence code (see header of file for codes).
            include_more_general: if True, include also more general GO terms annotated to gene products (default is True)
            """
        fg_list = [] # all terms, with multiple copies for counting
        fg_map = self.getTerms4Genes(positives, include_more_general = include_more_general) #
        for id in fg_map:
            for t in fg_map[id]:
                fg_list.append(t)
        term_set = set(fg_list)
        term_cnt = {}

        nPos = len(positives)
        if threshold == None:
            threshold = 0 # include all terms
        for t in term_set:
            cnt = fg_list.count(t)
            if cnt >= threshold:
                term_cnt[t] = cnt
        sorted_cnt = sorted(list(term_cnt.items()), key=lambda v: v[1], reverse=True)
        ret = []
        for t in sorted_cnt:
            defin = self.getTermdef(t[0])
            if defin == None:
                print(('Could not find definition of %s' % t[0]))
            else:
                ret.append((t[0], t[1], defin[2], defin[0]))
        return ret

    def getEnrichmentReport(self, positives, background = None, evid = None, threshold = None, include_more_general = True):
        """ For a set of named gene products (positives) this method determines the enrichment of GO terms.
            Each GO term is also assigned an enrichment p-value (on basis of provided background, or on basis of all annotated genes, if not provided).
            Note that to use the full set as background can be computationally expensive, so to speed up subsequent runs, the results are cached.
            Returns a list of tuples (GO_Term_ID[str], E-value[float], Foreground_no[int], Background_no[int], Term_description[str]).
            E-value is a Bonferroni-corrected p-value.
            positives: names of gene products
            background: names of gene products (or None if all annotated gene products should be used; default)
            threshold: E-value that must be reached for term to be reported (default is 0.05)
            If evid(ence) is specified the method returns only entries with that specific evidence code (see header of file for codes).
            include_more_general: if True, include also more general GO terms annotated to gene products (default is True)
            """
        # Process foreground: find terms of genes
        fg_list = [] # all terms, with multiple copies for counting
        fg_map = self.getTerms4Genes(positives, evid = evid, include_more_general = include_more_general) #
        for fg_gene in fg_map:
            for t in fg_map[fg_gene]:
                fg_list.append(t)
        nPos = len(positives)
        # Process background: find terms of genes
        bg_list = []
        if background == None: # need to use the full set
            background = list(self.annots.keys())
        negatives = set(background).difference(set(positives)) # remove the positives from the background to create genuine negatives
        nNeg = len(negatives)
        bg_map = self.getTerms4Genes(negatives, evid = evid, include_more_general = include_more_general)
        for bg_gene in bg_map:
            for t in bg_map[bg_gene]:
                bg_list.append(t)

        term_set = set(fg_list)
        term_cnt = {}

        if threshold == None:
            threshold = 0.05

        for t in term_set:
            fg_hit = fg_list.count(t) # number of foreground genes WITH GO term (number of terms in the list for the collective set of foreground genes)
            bg_hit = bg_list.count(t) # number of background genes WITH GO term (number of terms in the list for the collective set of background genes)
            fg_nohit = nPos - fg_hit  # total number of genes in foreground minus that number of hits
            bg_nohit = nNeg - bg_hit  # total number of genes in background minus that number of hits
            pvalue = stats.getFETpval(fg_hit, bg_hit, fg_nohit, bg_nohit, False) # one-tailed FET
            evalue = pvalue * len(term_set) # Bonferroni correction
            if evalue <= threshold: # check if significance req is fulfilled
                term_cnt[t] = (fg_hit, fg_hit + bg_hit, evalue)
        sorted_cnt = sorted(list(term_cnt.items()), key=lambda v: v[1][2], reverse=False)
        ret = []
        for t in sorted_cnt:
            defin = self.getTermdef(t[0])
            if defin == None:
                print(('Could not find definition of %s' % t[0]))
            else:
                ret.append((t[0], t[1][2], t[1][0], t[1][1], defin[2], defin[0]))
        return ret

class BinGO():

    # Structures to hold all data relevant to session, all keys are "encoded"
    annots = {}     # annotations: annots[gene] = (taxa, terms[term] = (evid, T/F))
    termdefs = {}   # definitions: termdefs[term] = (onto, terms[term] = relation, name)
    # Codes for encoding and decoding
    gene_code = None
    term_code = None
    evid_code = None
    # indices
    annot_index = {}
    # Files
    f = None

    def __init__(self, filename, taxa = None):
        """ The binary file contains all the data and will initialise
            gene annotations (annots) and term definitions (termdefs)
            and the encoding/decoding keys. """
        self.f = self._readBitFile(filename, taxa = taxa)

    def _decodeGeneIDs(self, gene_codes):
        if type(gene_codes) != list and type(gene_codes) != set and type(gene_codes) != tuple:
            gene_codes = [gene_codes]
        ids = []
        for i in gene_codes:
            s = decode(i, self.gene_code)
            ids.append(s)
        return ids

    def _encodeGeneIDs(self, gene_names):
        if type(gene_names) != list and type(gene_names) != set and type(gene_names) != tuple:
            gene_names = [gene_names]
        ids = []
        for i in gene_names:
            y = encode(i, self.gene_code)
            ids.append(y)
        return ids

    def _getGeneEntry(self, gene):
        peek = self.annot_index[gene]
        self.f.seek(peek, 0)
        buf = self.f.read(calcsize('IIH'))
        (gene_int, taxa_int, nterms) = unpack('IIH', buf)
        buf = self.f.read(nterms * calcsize('?BI'))
        terms_dict = {}
        for pos in range(0, len(buf) - 1, calcsize('?BI')):
            (qual_bool, evid_int, term_int) = unpack('?BI', buf[pos:pos+calcsize('?BI')])
            terms_dict[term_int] = (evid_int, qual_bool)
        return (taxa_int, terms_dict)

    def _getSuperTerms(self, term, rel = None):
        """ Recursively compute the transitive closure. """
        found = set()
        try:
            (_, closure, _) = self.termdefs[term]
            for (t, r) in list(closure.items()):
                if (not rel) or r == rel:
                    found.add(t)
                    found.update(self._getSuperTerms(t, rel))
        except KeyError:
            print(('Could not find GO:%s' % (''.join(decode(term, self.term_code)))))
        return found

    def _getChildTerms(self, term, rel = None):
        found = set()
        for (child, termdef) in list(self.termdefs.items()):
            (_, parents_dict, _) = termdef
            try:
                myrel = parents_dict[term]
                if rel == myrel or not rel: found.add(child)
            except KeyError:
                pass
        return found

    def _getSpecificTerms(self, term, rel = None):
        direct = self._getChildTerms(term, rel)
        found = set()
        for t in direct:
            found.add(t)
            found.update(self._getSpecificTerms(t, rel))
        return found

    def getTerms(self, genes, evid = None, onto = None, include_more_general = True):
        """
        Retrieve all GO terms for a given set of genes (or single gene).
        The result is given as a map (key=gene name, value=list of unique terms) OR
        in the case of a single gene as a list of unique terms.
        If include_more_general is True (default) then transitively related terms are included
        """
        mymap = dict()
        # STEP 1: Find all terms directly associated with specified genes
        direct = set() # all GO terms (encoded)
        ids = self._encodeGeneIDs(genes)
        for i in ids:
            gene_name = ''.join(decode(i, self.gene_code))
            mymap[gene_name] = set()
            try:
                (taxa, terms) = self._getGeneEntry(i)
                for (term, evid_and_qual) in list(terms.items()):
                    if evid_and_qual[1] and not evid: # if True and no evidence is specified
                        direct.add(term)
                        mymap[gene_name].add(term)
                    elif self.evid_code[evid_and_qual[0]] == evid:
                        direct.add(term)
                        mymap[gene_name].add(term)
            except KeyError:
                pass
                #print 'Failed to find annotations for gene %s' % gene_name
        if include_more_general:
            # STEP 2: Find the transitive closure of each term identified, store as a dictionary
            indirect = {}
            for t in direct:
                if not t in indirect:
                    indirect[t] = set(self._getSuperTerms(t))
        # STEP 3: compile and return results
        for gene in mymap:
            term_ids = mymap[gene]
            all_ids = set(term_ids)
            if include_more_general:
                for term_id in term_ids:
                    all_ids.update(indirect[term_id])
            mymap[gene] = set()
            for term_enc in all_ids:
                mymap[gene].add('GO:'+''.join(decode(term_enc, self.term_code)))
        return mymap

    def getAllGenes(self):
        names = []
        for g in self._decodeGeneIDs(list(self.annot_index.keys())):
            names.append(''.join(g))
        return names

    def getGenes(self, terms, evid = None, taxa = None, rel = None, include_more_specific = True):
        """ Retrieve all genes that are annotated with specified terms, and qualified by evidence, taxa etc. """
        """ TODO: Debug--suspect this implementation is incorrect. """
        term_ids = set()
        for t in terms:
            term_ids.add(encode(t[3:], self.term_code))
        # STEP 1 (optional): determine more specific terms to be included in query
        if include_more_specific:
            myterms = set()
            for t in term_ids:
                myterms.add(t)
                children = self._getSpecificTerms(t, rel)
                myterms.update(children)
            term_ids = myterms
        # STEP 2: identify genes with those terms
        found = {}
        for g in self.annot_index:
            gene_name = decode(g, self.gene_code)
            (mytaxa, tdict) = self._getGeneEntry(g)
            if not taxa or taxa == mytaxa:
                for annot_term in list(tdict.keys()):
                    if tdict[annot_term] == evid:
                        if annot_term in terms:
                            try:
                                added = found[gene_name]
                                added.add(annot_term)
                            except KeyError:
                                found[gene_name] = set([annot_term])
        # STEP 3: compile and return results
        for gene in found:
            term_ids = found[gene]
            all_ids = set(term_ids)
            found[gene] = set()
            for term_enc in all_ids:
                found[gene].add('GO:'+''.join(decode(term_enc, self.term_code)))
        return found

    def getTermdef(self, term):
        term_id = encode(term[3:], self.term_code)
        try:
            (onto_ch, terms_dict, name_peek) = self.termdefs[term_id]
            self.f.seek(name_peek, 0)
            term_name = self.f.readline()
            return (onto_codes[onto_ch], terms_dict, term_name)
        except KeyError:
            return ('Unknown', 'Unknown', 'Unknown')

    def _readBitFile(self, filename, taxa, termnames = False):
        f = open(filename, 'r')
        # STEP 1: header info
        ngene_code = None
        nterm_code = None
        nevid_code = None
        ngene_cnt = 0
        nterm_cnt = 0
        nevid_cnt = 0
        header = True
        total_gene_cnt = None
        current_gene_cnt = 0
        current_terms_cnt = 0
        annot_offset = 0
        obo_offset = 0
        while f:
            if not ngene_code:
                line = f.readline()
                fields = line.split()
                total_gene_cnt = int(fields[0])
                total_terms_cnt = int(fields[1])
                ngene_code = int(fields[2])
                nterm_code = int(fields[3])
                nevid_code = int(fields[4])
                self.gene_code = ['' for _ in range(ngene_code)]
                self.term_code = ['' for _ in range(nterm_code)]
                self.evid_code = ['' for _ in range(nevid_code)]
            elif ngene_cnt < ngene_code:
                line = f.readline()
                self.gene_code[ngene_cnt] = line.strip()
                ngene_cnt += 1
            elif nterm_cnt < nterm_code:
                line = f.readline()
                self.term_code[nterm_cnt] = line.strip()
                nterm_cnt += 1
            elif nevid_cnt < nevid_code:
                line = f.readline()
                self.evid_code[nevid_cnt] = line.strip()
                nevid_cnt += 1
            else: # we're not in the header
                if header: offset = f.tell()
                header = False
                try:
                    if current_gene_cnt < total_gene_cnt: # we are reading gene:terms annotations
                        peek = f.tell()
                        buf = f.read(calcsize('IIH'))
                        (gene_int, taxa_int, nterms) = unpack('IIH', buf)
                        current_gene_cnt += 1
                        if (not taxa) or (taxa_int == taxa or taxa_int in taxa):
                            self.annot_index[gene_int] = peek
                        bufsize = calcsize('?BI')
                        f.read(nterms * bufsize)
                    elif current_terms_cnt < total_terms_cnt: # we are reading term definitions (term is_a term, term, term, ...)
                        buf = f.read(calcsize('IcH'))
                        (term_int, onto_ch, nterms) = unpack('IcH', buf)
                        current_terms_cnt += 1
                        bufsize = calcsize('BI')
                        buf = f.read(nterms * bufsize)
                        terms_dict = {}
                        for pos in range(0, len(buf) - 1, bufsize):
                            (rel_ndx, sup_int) = unpack('BI', buf[pos:pos+bufsize])
                            terms_dict[sup_int] = rel_ndx
                        name_peek = f.tell()
                        f.readline() # skip putting name in memory, instead refer to the position in the file
                        self.termdefs[term_int] = (onto_ch, terms_dict, name_peek)
                    else:
                        buf = f.read(calcsize('II'))
                        (annot_offset, obo_offset) = unpack('II', buf)
                        break
                except error as inst:
                    print(("Problem reading binary file: ", inst, "at gene ", current_gene_cnt, "at definition ", current_terms_cnt, "at", f.tell()))
                    exit(3)
        print(("Read %d genes and %d term definitions" % (current_gene_cnt, current_terms_cnt)))
        print(("Annotations start at", annot_offset, "\nDefinitions start at", obo_offset))
        return f

    #FIXME: write code to perform test of taxa enrichment

    def getGOReport_byScore(self, gene_score_map, negatives_score_map = {}, include_more_general = True, descending_order = True):
        """ Generate a complete GO term report for a set of genes with associated scores.
            Uses the Wilcoxon Ranksum test for each GO term to assign a p-value,
            indicating the enrichment of term to "top" genes in descending order by score (by default).
        """
        fg_map = self.getTerms(list(gene_score_map.keys()), include_more_general = include_more_general)
        fg_list = []
        for id in fg_map:
            for t in fg_map[id]:
                fg_list.append(t)
        term_set = set(fg_list)
        term_pval = {}
        if len(negatives_score_map) > 0:
            bg_map = self.getTerms(list(negatives_score_map.keys()), include_more_general = include_more_general)
        for t in term_set:
            pos = []
            neg = []
            for gene in gene_score_map:
                annot = fg_map[gene]
                if not annot == None:
                    if t in annot:
                        pos.append(gene_score_map[gene])
                    else:
                        neg.append(gene_score_map[gene])
            if len(pos) > 0 and len(neg) > 0:
                if descending_order:
                    p = stats.getRSpval(neg, pos)
                else:
                    p = stats.getRSpval(pos, neg)
            if len(negatives_score_map) > 0 and p <= 0.05:
                mpos = pos # scores of foreground genes with matching GO term
                mneg = [] # scores of background genes with matching GO terms
                for gene in negatives_score_map:
                    annot = bg_map[gene]
                    if not annot == None:
                        if t in annot:
                            mneg.append(negatives_score_map[gene])
                if len(mneg) > 0:
                    if descending_order:
                        p2 = stats.getRSpval(mneg, mpos)
                    else:
                        p2 = stats.getRSpval(mpos, mneg)
                else:
                    p2 = 0.0
                term_pval[t] = (p, p2)
            else:
                term_pval[t] = (p, 1.0)

        sorted_pval = sorted(list(term_pval.items()), key=lambda v: v[1][0], reverse=False)

        ret = []
        for t in sorted_pval:
            defin = self.getTermdef(t[0])
            if defin == None:
                print(('Could not find definition of %s' % t[0]))
            else:
                ret.append((t[0], t[1][0], t[1][1], defin[2].strip(), defin[0]))
        return ret

    def getGOReport(self, positives, background = None, taxa = None, include_more_general = True):
        """ Generate a complete GO term report for a set of genes (positives).
            Each GO term is also assigned an enrichment p-value (on basis of background, if provided).
            Returns a list of tuples (GO_Term_ID[str], Foreground_no[int], Term_description[str]) with no background, OR
            (GO_Term_ID[str], E-value[float], Foreground_no[int], Background_no[int], Term_description[str]).
            E-value is a Bonferroni-corrected p-value.
            """
        pos = set(positives)
        fg_map = self.getTerms(pos, include_more_general = include_more_general)
        fg_list = []
        for id in fg_map:
            for t in fg_map[id]:
                fg_list.append(t)
        bg_map = {}
        bg_list = []
        neg = set()
        if background != None:
            neg = set(background).difference(pos)
            bg_map = self.getTerms(neg, include_more_general = include_more_general)
            for id in bg_map:
                for t in bg_map[id]:
                    bg_list.append(t)
        term_set = set(fg_list)
        term_cnt = {}

        nPos = len(pos)
        nNeg = len(neg)
        if background == None:
            for t in term_set:
                term_cnt[t] = fg_list.count(t)
            sorted_cnt = sorted(list(term_cnt.items()), key=lambda v: v[1], reverse=True)
        else: # a background is provided
            for t in term_set:
                fg_hit = fg_list.count(t)
                bg_hit = bg_list.count(t)
                fg_nohit = nPos - fg_hit
                bg_nohit = nNeg - bg_hit
                term_cnt[t] = (fg_hit, fg_hit + bg_hit, stats.getFETpval(fg_hit, bg_hit, fg_nohit, bg_nohit, False))
            sorted_cnt = sorted(list(term_cnt.items()), key=lambda v: v[1][2], reverse=False)

        ret = []
        for t in sorted_cnt:
            defin = self.getTermdef(t[0])
            if defin == None:
                print(('Could not find definition of %s' % t[0]))
            else:
                if background != None:
                    ret.append((t[0], t[1][2] * len(term_set), t[1][0], t[1][1], defin[2], defin[0]))
                else:
                    ret.append((t[0], t[1], defin[2], defin[0]))
        return ret

def encode(code_me, encode_strings):
    code = 0
    accum = 1
    try:
        for pos in range(len(code_me)):
            codelen = len(encode_strings[pos])
            for i in range(codelen):
                if encode_strings[pos][i] == code_me[pos]:
                    code += accum * i
                    accum *= codelen
                    break
    except IndexError as e:
        print((e, code_me))
    return code

def decode(code, encode_strings):
    npos    = len(encode_strings)
    accum   = [1 for _ in range(npos)]
    try:
        for pos in range(1, npos): accum[pos] = accum[pos - 1] * len(encode_strings[pos - 1])
        indices = [-1 for _ in range(npos)]
        for pos in range(npos - 1, -1, -1): # go backwards, start at last (most significant) position
            indices[pos] = code / accum[pos]
            code -= accum[pos] * indices[pos]
        string = [encode_strings[pos][indices[pos]] for pos in range(len(encode_strings))]
    except IndexError as e:
        print((e, code))
    return string

def _extractAnnotFields(line, columns = (1,2,3,4,6,8)):
    """ Extract appropriate details from annotation file. This typically follows:
        1.  DB
            Database from which entry has been taken.
            Example: PDB
        2.  DB_Object_ID
            A unique identifier in the DB for the item being annotated.
            Here: PDB ID and chain ID of the PDB entry.
            Example: 2EKB_A
        3.  DB_Object_Symbol
            Here: PDB ID and chain ID of the PDB entry.
            Example:EKB_A
        4.  Qualifiers
            This column is used for flags that modify the interpretation of an annotation.
            This field may be equal to: NOT, colocalizes_with, contributes_to,
            NOT|contributes_to, NOT|colocalizes_with
            Example: NOT
        5.  GO Identifier
            The GO identifier for the term attributed to the DB_Object_ID.
            Example: GO:0005625
        6.  DB:Reference
            A single reference cited to support an annotation.
            Where an annotation cannot reference a paper, this field will contain
            a GO_REF identifier. See section 8 and
            http://www.geneontology.org/doc/GO.references
            for an explanation of the reference types used.
            Example: PMID:9058808
        7.  Evidence
            One of either EXP, IMP, IC, IGI, IPI, ISS, IDA, IEP, IEA, TAS, NAS,
            NR, ND or RCA.
            Example: TAS
        9.  Aspect
            One of the three ontologies: P (biological process), F (molecular function)
            or C (cellular component).
            Example: P

        In columns specify index (starts with 0 NOT 1) for gene, symb, qual, term, evid, onto
        """
    fields = line.strip().split('\t')
    gene = fields[columns[0]]
    symb = fields[columns[1]]
    qual = fields[columns[2]]
    term = fields[columns[3]]
    if not term.startswith('GO:'):
        term = None
        raise "No GO term on line: " + line
    evid = fields[columns[4]]
    if not evid in evid_codes:
        evid = None
    onto = fields[columns[5]]
    if not onto in onto_codes:
        onto = None
    taxa_idx = line.find('taxon:')
    if taxa_idx == -1:
        taxa = None
    else:
        taxa = line[taxa_idx:]
        taxa = taxa.split('\t')
        taxa_spec = taxa[0].split(':')
        taxa = int(taxa_spec[len(taxa_spec) - 1]) # pick last taxon ID
    return (gene, symb, qual, term, evid, onto, taxa)

def readOBOFile(obofile):
    """
    http://www.geneontology.org/GO.format.obo-1_2.shtml
    """
    src = open(obofile, 'r')
    terms = {}
    in_term_def = False
    in_type_def = False
    for line in src:
        if in_term_def:
            if line.startswith('id: '):
                term_id = line[4:14]
                term_is = set()
            elif line.startswith('name: '):
                term_name = line[6:].strip()
            elif line.startswith('def: '):
                # Note this is a multi-line field, delimited by "'s
                pass
            elif line.startswith('namespace: '):
                if   line[11] == 'b': term_onto = 'P'
                elif line[11] == 'm': term_onto = 'F'
                elif line[11] == 'c': term_onto = 'C'
            elif line.startswith('is_a: '):
                term_is.add((line[6:16], 'is_a'))
            elif line.startswith('relationship: '):
                fields = line.split()
                term_is.add((fields[2], fields[1]))
            elif line.startswith('intersection_of: '):
                fields = line.split()
                if fields[1].startswith('GO:'):
                    term_is.add((fields[1], 'isect'))
                else:
                    term_is.add((fields[2], fields[1]))
            elif line.startswith('is_obsolete: '):
                in_term_def = False # ignore this entry
        if line.startswith('[Term]'):
            if in_term_def: # already defining one, stash it before moving on to the next...
                terms[term_id] = (term_name, term_onto, term_is)
            elif in_type_def:
                in_type_def = False
            in_term_def = True
        if line.startswith('[Typedef]'):
            if in_term_def: # already defining one, stash it before moving on to the next...
                in_term_def= False
            in_type_def = True
    if in_term_def: #  defining one, stash it
        terms[term_id] = (term_name, term_onto, term_is)
    return terms

def writeBitFile(annotFile, obofile, destFile, taxas = None):
    print(("Started at", time.asctime()))
    # open annotation file to analyse and index data
    src = open(annotFile, 'r')
    gene_index = [{} for _ in range(6)]  # count different characters in different positions
    term_index = [{} for _ in range(7)]  # count different characters in different positions
    evid_index = {}
    gene_cnt = 0
    cnt = 0
    prev_gene = None
    for line in src:
        cnt += 1
        #if cnt > 100000:
        #    break
        if line.startswith('!'):
            continue
        (gene, symb, qual, term, evid, onto, taxa) = _extractAnnotFields(line)
        if not (taxas and ((taxa == taxas) or (taxa in taxas))):  # The gene does NOT belong to a nominated taxon
            continue
        if not (gene == prev_gene): # not the same gene
            gene_cnt += 1
        try:
            evid_index[evid]
        except:
            evid_index[evid] = len(evid_index)
        pos = 0
        for ch in gene[0:6]:
            try:
                gene_index[pos][ch]
            except KeyError: # no match
                gene_index[pos][ch] = len(gene_index[pos])
            pos += 1
        pos = 0
        for ch in term[3:10]:
            try:
                term_index[pos][ch]
            except KeyError: # no match
                term_index[pos][ch] = len(term_index[pos])
            pos += 1
        prev_gene = gene
    src.close()
    print(("Read annotations for %d genes" % gene_cnt))

    gene_code = ['' for _ in range(6)]
    term_code = ['' for _ in range(7)]
    for d in range(len(gene_index)):
        arr = ['?' for _ in gene_index[d]]
        for (ch, index) in list(gene_index[d].items()):
            arr[index] = ch
        gene_code[d] = ''.join(arr)
    for d in range(len(term_index)):
        arr = ['?' for _ in term_index[d]]
        for (ch, index) in term_index[d].items():
            arr[index] = ch
        term_code[d] = ''.join(arr)
    evid_code = ['' for _ in range(len(evid_index))]
    for (e, ndx) in list(evid_index.items()):
        evid_code[ndx] = e

    # Get GO definitions
    terms = readOBOFile(obofile)
    print(("Read %d GO definitions" % len(terms)))

    # re-open, now with the aim of copying info
    src = open(annotFile, 'r')
    dst = open(destFile, 'w')
    # STEP 1: header info
    dst.write("%d\t%d\t%d\t%d\t%d\n" % (gene_cnt, len(terms), len(gene_code), len(term_code), len(evid_index)))
    for code_str in gene_code:
        dst.write(code_str+"\n")
    for code_str in term_code:
        dst.write(code_str+"\n")
    for e_str in evid_code:
        dst.write(e_str+'\n')
    print(("Wrote header %d\t%d\t%d\t%d\t%d, now at @%d" % (gene_cnt, len(terms), len(gene_code), len(term_code), len(evid_index), dst.tell())))

    # STEP 2: write annotations
    annot_offset = dst.tell()
    prev_gene = None
    concat_terms = {}
    cnt = 0
    for line in src:
        cnt += 1
        #if cnt > 100000:
        #    break
        if line.startswith('!'):
            continue
        (gene, symb, qual, term, evid, onto, taxa) = _extractAnnotFields(line)
        if not (taxas and ((taxa == taxas) or (taxa in taxas))): # The gene does NOT belong to a nominated taxon
            continue
        if gene != prev_gene: # new gene is found
            if prev_gene != None:
                # write prev data
                s = pack('IIH', encode(prev_gene, gene_code), taxa, len(concat_terms))
                dst.write(s)
                for t in concat_terms:
                    (o, q, e) = concat_terms[t]
                    s = pack('?BI', q, evid_index[e], encode(t, term_code))
                    dst.write(s)
            # re-init
            prev_gene = gene
            concat_terms = {}
        concat_terms[term[3:]] = (onto, qual, evid)
    if len(concat_terms) > 0:
        # write data in buffer
        s = pack('IIH', encode(prev_gene, gene_code), taxa, len(concat_terms))
        dst.write(s)
        for t in concat_terms:
            (o, q, e) = concat_terms[t]
            s = pack('?BI', q, evid_index[e], encode(t, term_code))
            dst.write(s)
    print(("Wrote GO annotations, now at @%d" % dst.tell()))

    # Next, the ontology definition...
    obo_offset = dst.tell() # remember the position where the OBO starts
    sorted_terms = sorted(iter(terms.items()), key=operator.itemgetter(0))
    for [t, _] in sorted_terms:
        (term_name, term_onto, term_is) = terms[t]
        s = pack('IcH', encode(t[3:], term_code), term_onto, len(term_is))
        dst.write(s)
        for (sup_term, sup_rel) in term_is:
            try:
                index = onto_rel.index(sup_rel)
            except ValueError:
                index = 9
            s = pack('BI', index, encode(sup_term[3:], term_code))
            dst.write(s)
        dst.write(term_name + '\n')
    print(("Wrote %d GO definitions, now at @%d" % (len(sorted_terms), dst.tell())))

    # Finally, write the offsets to quickly access annotations and definitions, resp
    dst.write(pack('II', annot_offset, obo_offset))
    # done, close
    dst.close()
    print(("Completed at", time.asctime()))
