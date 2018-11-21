import re
import locale
from collections import defaultdict
from flib.core.gmt import GMT
from flib.core.url import URLResource
from six import iteritems,itervalues
import logging
logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


class OBO:

    def __init__(self, obo_file=None):
        """Initialize with optional obo file"""
        self.heads = []
        self.go_terms = {}
        self.go_obsolete = {}
        self.alt_id2std_id = {}
        self.name2synonyms = {}
        self.populated = False
        self._meta = {}

        if obo_file:
            self.load_obo(obo_file)

    def load_obo(self, obo_file, remote_location=False, timeout=5):

        if remote_location:
            lines = URLResource(obo_file).get_lines()
        else:
            obo = open(obo_file, 'r')
            lines = obo.readlines()

        inside = False
        gterm = None
        for line in lines:
            fields = line.rstrip().split()
            if len(fields) < 1:
                # Blank line
                continue
            elif len(fields) < 2 and not fields[0].startswith('['):
                # All other lines should be key:value pairs
                logger.debug('Skipping unrecognized line: %s', line)
                continue

            # Load the meta data commonly included at the
            # start of the obo file.
            # Example keys: format-version, data-version, etc
            elif not inside and not len(self.go_terms.keys()) and len(fields) > 1:
                key = fields[0]
                if key.endswith(':'):
                    key = key[:-1]
                self._meta[key] = fields[1]

            # Term definition block
            elif fields[0] == '[Term]':
                if gterm:
                    if gterm.head:
                        self.heads.append(gterm)
                inside = True

            # Relationship definition block
            elif fields[0] == '[Typedef]':
                if gterm:
                    if gterm.head:
                        self.heads.append(gterm)
                inside = False

            # Term identifier field (e.g. GO:00008150)
            elif inside and fields[0] == 'id:':
                if fields[1] in self.go_terms:
                    gterm = self.go_terms[fields[1]]
                else:
                    gterm = GOTerm(fields[1])
                    self.go_terms[gterm.get_id()] = gterm

            # Term name field (e.g. biological process)
            elif inside and fields[0] == 'name:':
                # Term name is underscore delimited; term fullname is
                # space delimited
                fields.pop(0)
                gterm.fullname = ' '.join(fields)
                name = '_'.join(fields)
                name = name.replace('\'', '')
                name = re.sub('[^\w\s_-]', '_', name).strip().lower()
                name = re.sub('[-\s_]+', '_', name)
                gterm.name = name

            # Ontology namespace field (e.g. molecular function)
            elif inside and fields[0] == 'namespace:':
                gterm.namespace = fields[1]

            # Term definition field (e.g. "The maintenance of the...")
            elif inside and fields[0] == 'def:':
                gterm.desc = ' '.join(fields[1:]).split('\"')[1]

            # Alternative identifiers, likely historical, now obsolete terms
            elif inside and fields[0] == 'alt_id:':
                gterm.alt_id.append(fields[1])
                self.alt_id2std_id[fields[1]] = gterm.get_id()

            # is_a relationship field
            elif inside and fields[0] == 'is_a:':
                # If term has a parent, it can't be a head (root) term
                gterm.head = False

                fields.pop(0)
                pgo_id = fields.pop(0)
                if pgo_id not in self.go_terms:
                    self.go_terms[pgo_id] = GOTerm(pgo_id)

                gterm.is_a.append(self.go_terms[pgo_id])
                self.go_terms[pgo_id].parent_of.add(gterm)
                gterm.child_of.add(self.go_terms[pgo_id])

            # Relationship field defined in Typedef
            # e.g. "relationship: regulates GO:0000XXX"
            elif inside and fields[0] == 'relationship:':
                if fields[1].find('has_part') != -1:
                    # has part is not a parental relationship -- it is actually
                    # for children.
                    continue

                # If term has a parent, it can't be a head (root) term
                gterm.head = False

                pgo_id = fields[2]
                if pgo_id not in self.go_terms:
                    self.go_terms[pgo_id] = GOTerm(pgo_id)

                # Check which relationship you are with this parent go term
                if fields[1] == 'regulates' or \
                        fields[1] == 'positively_regulates' or \
                        fields[1] == 'negatively_regulates':
                    gterm.relationship_regulates.append(self.go_terms[pgo_id])
                elif fields[1] == 'part_of':
                    gterm.relationship_part_of.append(self.go_terms[pgo_id])
                else:
                    logger.info(
                        "Unknown relationship %s",
                        self.go_terms[pgo_id].name)
                    continue
                self.go_terms[pgo_id].parent_of.add(gterm)
                gterm.child_of.add(self.go_terms[pgo_id])

            # Term obsolete flag
            elif inside and fields[0] == 'is_obsolete:':
                # If term is obsolete, it can't be a root term
                gterm.head = False

                # Only keep current terms
                del self.go_terms[gterm.get_id()]

                gterm.obsolete = True
                self.go_obsolete[gterm.get_id()] = gterm

            # Term name synonyms
            elif inside and fields[0] == 'synonym:':
                syn = ' '.join(fields[1:]).split('\"')[1]
                syn = syn.replace('lineage name: ', '')
                gterm.synonyms.append(syn)
                if gterm.name in self.name2synonyms:
                    self.name2synonyms[gterm.name].append(syn)
                else:
                    self.name2synonyms[gterm.name] = [syn]

            # Term id mapping to another resources (e.g. OMIM, Wikipedia)
            elif inside and fields[0] == 'xref:':
                tok = fields[1].split(':')
                if len(tok) > 1:
                    (xrefdb, xrefid) = tok[:2]
                    gterm.xrefs.setdefault(xrefdb, set()).add(xrefid)
        if not remote_location:
            obo.close()
        return True

    def propagate(self):
        """Propagate all gene annotations"""
        logger.info("Propagate gene annotations")
        for head_gterm in self.heads:
            logger.info("Propagating %s", head_gterm.name)
            self._propagate_recurse(head_gterm)

    def _propagate_recurse(self, gterm):
        if not len(gterm.parent_of):
            logger.debug("Base case with term %s", gterm.name)
            return

        for child_term in gterm.parent_of:
            self._propagate_recurse(child_term)
            new_annotations = set()

            regulates_relation = (gterm in child_term.relationship_regulates)
            part_of_relation = (gterm in child_term.relationship_part_of)

            for annotation in child_term.annotations:
                copied_annotation = None
                # If this relation with child is a regulates (and its sub class)
                # filter annotations
                if regulates_relation:
                    # Only add annotations that didn't come from a part_of or
                    # regulates relationship. "ready_regulates_cutoff" indicates
                    # the annotation was propagated from a regulates or part_of
                    # relationship already.
                    # See: http://www.geneontology.org/page/ontology-relations##reg_reas
                    if annotation.ready_regulates_cutoff:
                        continue
                    else:
                        copied_annotation = annotation.prop_copy(
                            ready_regulates_cutoff=True)
                elif part_of_relation:
                    copied_annotation = annotation.prop_copy(
                        ready_regulates_cutoff=True)
                else:
                    copied_annotation = annotation.prop_copy()

                new_annotations.add(copied_annotation)
            gterm.annotations = gterm.annotations | new_annotations

    def get_term(self, tid):
        """Return GOTerm object corresponding with id=tid"""
        logger.debug('get_term: %s', tid)
        term = None
        try:
            term = self.go_terms[tid]
        except KeyError:
            try:
                term = self.go_terms[self.alt_id2std_id[tid]]
            except KeyError:
                logger.warning('Term name does not exist: %s', tid)
        return term

    def get_meta_data(self, key):
        """Return metadata in obo corresponding to key"""
        if key in self._meta:
            return self._meta[key]
        else:
            return None

    def get_termobject_list(self, terms=None, p_namespace=None):
        """Return list of all GOTerms"""
        logger.info('get_termobject_list')
        if terms is None:
            terms = self.go_terms.keys()
        reterms = []
        for tid in terms:
            obo_term = self.get_term(tid)
            if obo_term is None:
                continue
            if p_namespace is not None and obo_term.namespace != p_namespace:
                continue
            reterms.append(obo_term)
        return reterms

    def get_obsolete_terms(self):
        """Return list of all obsolete GOTerms"""
        logger.info('get_obsolete_list')
        return self.go_obsolete.values()

    def get_xref_mapping(self, prefix):
        """Return dict of terms mappings to external database ids"""
        xrefs = defaultdict(set)
        for term in self.get_termobject_list():
            ids = term.get_xrefs(prefix)
            if ids:
                for xref in ids:
                    xrefs[xref].add(term.go_id)
        return xrefs

    def as_gmt(self):
        """Return gene annotations as GMT object"""
        gmt = GMT()
        tlist = sorted(self.get_termobject_list(),key=cmp_to_key(go_term_id_comparison))
        for term in tlist:
            if len(term.annotations):
                gmt.add_geneset(gsid=term.go_id, name=term.name)
            for annotation in term.annotations:
                gmt.add_gene(term.go_id, annotation.gid)
        return gmt

    def map_genes(self, id_name, xdb_prefixed=False):
        """Map gene names using the idmap object id_name"""
        for go_term in itervalues(self.go_terms):
            go_term.map_genes(id_name, xdb_prefixed=xdb_prefixed)

    def filter_annotations(self, evidence_codes):
        """Filter out gene annotations by their annotation evidence"""
        for go_term in itervalues(go_terms):
            go_term.filter_annotations(evidence_codes)

    def populate_annotations(self, annotation_file, remote_location=False, xdb_col=0,
                             gene_col=1, term_col=4, ref_col=5, ev_col=6, date_col=13):
        """Populate the ontology with gene annotations from an association file"""
        logger.info('Populate gene annotations: %s', annotation_file)

        if remote_location:
            lines = URLResource(annotation_file).get_lines()
            #lines = ass_file.readlines()
        else:
            ass_file = open(annotation_file, 'r')
            lines = ass_file.readlines()
            ass_file.close()

        details_col = 3
        for line in lines:
            if line.startswith('!'):
                continue

            fields = line.rstrip('\n').split('\t')
            if len(fields) == 1:
                continue

            xdb = fields[xdb_col]
            gene = fields[gene_col]
            go_id = fields[term_col]

            try:
                ref = fields[ref_col]
            except IndexError:
                ref = None
            try:
                ev = fields[ev_col]
            except IndexError:
                ev = None
            try:
                date = fields[date_col]
            except IndexError:
                date = None

            if date_col < len(fields):
                date = fields[date_col]
            else:
                date = None

            try:
                details = fields[details_col]
                if details == 'NOT':
                    continue
            except IndexError:
                pass
            go_term = self.get_term(go_id)
            if go_term is None:
                continue
            logger.debug('Gene %s and term %s', gene, go_term.go_id)
            annotation = Annotation(
                xdb=xdb,
                gid=gene,
                ref=ref,
                evidence=ev,
                date=date,
                direct=True)
            go_term.annotations.add(annotation)

        self.populated = True

    def populate_annotations_from_gmt(self, gmt):
        """Populate the ontology with gene annotations from a GMT file"""
        for (gsid, genes) in iteritems(gmt.genesets):
            term = self.get_term(gsid)
            if term:
                for gid in genes:
                    term.add_annotation(gid)

    def add_annotation(self, go_id, gid, ref, direct):
        """Add a gene annotation to a term
        Args:
            go_id:  term identifier
            gid:    gene identifier
            ref:    publication reference (e.g. pubmed id)
            direct: boolean indicating direct or propagated annotation

        Returns:
            True for succes, False otherwise
        """
        go_term = self.get_term(go_id)
        if not go_term:
            return False
        annot = Annotation(xdb=None, gid=gid, direct=direct, ref=ref)
        go_term.annotations.add(annot)
        return True

    def get_descendents(self, gterm):
        """Return propagated descendents of term"""
        if gterm not in self.go_terms:
            return set()
        term = self.go_terms[gterm]

        if len(term.parent_of) == 0:
            return set()

        child_terms = set()
        for child_term in term.parent_of:
            if child_term.namespace != term.namespace:
                logger.info("Parent and child terms are different namespaces: %s and %s",
                        child_term, term)
                continue
            child_terms.add(child_term.go_id)
            child_terms = child_terms | self.get_descendents(child_term.go_id)

        return child_terms

    def get_ancestors(self, gterm):
        """Return propagated ancestors of term"""
        if (gterm in self.go_terms) is False:
            return set()
        term = self.go_terms[gterm]

        if len(term.child_of) == 0:
            return set()

        parent_terms = set()
        for parent_term in term.child_of:
            if parent_term.namespace != term.namespace:
                logger.info("Parent and child terms are different namespaces: %s and %s",
                        parent_term, term)
                continue
            parent_terms.add(parent_term.go_id)
            parent_terms = parent_terms | self.get_ancestors(parent_term.go_id)

        return parent_terms

    def get_leaves(self, namespace='biological_process', min_annot=10):
        """Return a set of leaf terms from the ontology"""
        leaves, bottom = set(), set()
        for term in self.go_terms.values():
            if len(term.parent_of) == 0 and term.namespace == namespace and len(
                    term.annotations) >= min_annot:
                leaves.add(term)
        return leaves

    def print_to_dir(self, out_dir, terms=None, p_namespace=None):
        """Writes to out_dir each term and its gene annotations in individual files"""
        logger.info('print_terms')
        tlist = self.get_termobject_list(terms=terms, p_namespace=p_namespace)
        # print terms
        for term in tlist:
            id_set = set(term.get_annotated_genes())
            if len(id_set) == 0:
                continue
            output_fh = open(out_dir + '/' + term.name, 'w')
            # keep previous behavior w/ newline at end
            output_fh.write('\n'.join(id_set) + '\n')
            output_fh.close()

    def print_to_single_file(self, out_file, terms=None,
                             p_namespace=None, gene_asso_format=False):
        logger.info('print_to_single_file')
        tlist = sorted(
            self.get_termobject_list(
                terms=terms,
                p_namespace=p_namespace),key=cmp_to_key(go_term_id_comparison))
        f = open(out_file, 'w')
        for term in tlist:
            for annotation in term.annotations:
                if gene_asso_format:
                    to_print = [annotation.xdb if annotation.xdb else '',
                                annotation.gid if annotation.gid else '',
                                '', '',  # Gene Symbol, NOT/''
                                term.go_id if term.go_id else '',
                                annotation.ref if annotation.ref else '',
                                annotation.evidence if annotation.evidence else '',
                                annotation.date if annotation.date else '',
                                str(annotation.direct),
                                # Direct is added in to indicate prop status
                                # cross annotated is added in to indicate cross
                                # status
                                str(annotation.cross_annotated),
                                # if cross annotated, where the annotation is
                                # from
                                annotation.origin if annotation.cross_annotated else '',
                                # if cross annotated, then the evidence of the
                                # cross_annotation (e.g. bootstrap value,
                                # p-value)
                                str(annotation.ortho_evidence) if annotation.ortho_evidence else '', '', '']
                    #print >> f, '\t'.join([str(x) for x in to_print])
                    line = '\t'.join([str(x) for x in to_print])+'\n'
                    f.write(line)
                else:
                    #print >> f, term.go_id + '\t' + term.name + '\t' + annotation.gid
                    f.write(term.go_id + '\t' + term.name + '\t' + annotation.gid+'\n')

        f.close()

    def print_to_gmt_file(self, out_file, terms=None, p_namespace=None):
        logger.info('print_to_gmt_file')
        tlist = sorted(self.get_termobject_list(terms=terms,p_namespace=p_namespace),key=cmp_to_key(go_term_id_comparison))
        f = open(out_file, 'w')
        for term in tlist:
            genes = set()
            for annotation in term.annotations:
                genes.add(annotation.gid)
            if len(genes) > 0:
                line = "%s\t%s\t" % (term.go_id,term.name)
                genes_str = "\t".join(genes)
                line+=genes_str
                line+="\n"
                f.write(line)
        f.close()

    def print_to_mat_file(self, out_file, terms=None, p_namespace=None):
        logger.info('print_to_mat_file')
        tlist = sorted(
            self.get_termobject_list(
                terms=terms,
                p_namespace=p_namespace),key=cmp_to_key(go_term_id_comparison))
        f = open(out_file, 'w')

        allgenes = set()
        genedict = defaultdict(set)
        termlist = []
        for term in tlist:
            if len(term.annotations) == 0:
                continue

            termlist.append(term.go_id)

            for annotation in term.annotations:
                allgenes.add(annotation.gid)
                genedict[annotation.gid].add(term.go_id)

        #print >> f, '\t' + '\t'.join(termlist)
        line = '\t' + '\t'.join(termlist) + '\n'
        f.write(line)

        for g in list(allgenes):
            row = []
            row.append(g)
            for termid in termlist:
                row.append('1' if termid in genedict[g] else '0')
                line = '\t'.join(row) + '\n'
                f.write(line)
        f.close()


class Annotation(object):

    def __init__(self, xdb=None, gid=None, ref=None, evidence=None, date=None, direct=False,
                 cross_annotated=False, origin=None, ortho_evidence=None, ready_regulates_cutoff=False):
        # Annotation source
        super(Annotation, self).__setattr__('xdb', xdb)
        # Gene identifier
        super(Annotation, self).__setattr__('gid', gid)

        # Publication reference
        super(Annotation, self).__setattr__('ref', ref)

        # Evidence code
        super(Annotation, self).__setattr__('evidence', evidence)

        # Date of annotation
        super(Annotation, self).__setattr__('date', date)

        # Direct annotation or possibly propagated
        super(Annotation, self).__setattr__('direct', direct)

        # Annotated from another organism
        super(Annotation, self).__setattr__('cross_annotated', cross_annotated)
        super(Annotation, self).__setattr__('origin', origin)
        super(Annotation, self).__setattr__('ortho_evidence', ortho_evidence)

        # Boolean indicating whether the annotation can be propagated along a
        # regulates relationship.
        # See: http://www.geneontology.org/page/ontology-relations##reg_reas
        super(
            Annotation,
            self).__setattr__(
            'ready_regulates_cutoff',
            ready_regulates_cutoff)

    def prop_copy(self, ready_regulates_cutoff=None):
        if ready_regulates_cutoff is None:
            ready_regulates_cutoff = self.ready_regulates_cutoff

        return Annotation(xdb=self.xdb, gid=self.gid, ref=self.ref,
                          evidence=self.evidence, date=self.date, direct=False, cross_annotated=False,
                          ortho_evidence=self.ortho_evidence, ready_regulates_cutoff=ready_regulates_cutoff)

    def __hash__(self):
        return hash((self.xdb, self.gid, self.ref, self.evidence, self.date,
                     self.direct, self.cross_annotated, self.ortho_evidence,
                     self.ready_regulates_cutoff, self.origin))

    def __eq__(self, other):
        return (self.xdb, self.gid, self.ref, self.evidence, self.date,
                self.direct, self.cross_annotated, self.ortho_evidence,
                self.ready_regulates_cutoff, self.origin).__eq__((other.xdb,
                                                                  other.gid, other.ref, other.evidence, other.date,
                                                                  other.direct, other.cross_annotated, other.ortho_evidence,
                                                                  other.ready_regulates_cutoff, other.origin))

    def __setattr__(self, *args):
        raise TypeError("Attempt to modify immutable object.")
    __delattr__ = __setattr__


class GOTerm:

    def __init__(self, go_id):
        # Indicator of whether the term is a root node
        self.head = True

        # Term identifier
        self.go_id = go_id

        # Set of gene annotations
        self.annotations = set([])

        # List of is_a parents
        self.is_a = []

        # List of regulates parents
        # Note: if A regulates B, B is A's parent in gene ontology
        self.relationship_regulates = []

        # List of part_of parents
        self.relationship_part_of = []

        # All parent terms
        self.parent_of = set()

        # All child terms
        self.child_of = set()

        # Alternative IDs, likely to be obsolete
        self.alt_id = []

        # Namespace of the term
        self.namespace = None

        # Term description
        self.desc = None

        # Term name, delimited by underscores
        self.name = None

        # Official term name, unadulterated
        self.fullname = None

        # Term name synonyms
        self.synonyms = []

        # Term ID mappings to other resources
        self.xrefs = {}

        # Boolean indicated whether the term is now obsolete
        self.obsolete = False

        # As far as I can tell, no longer used (7/13/2017)
        # self.cross_annotated_genes = set([])
        # self.included_in_all = True
        # self.valid_go_term = True
        # self.base_counts = None
        # self.counts = None
        # self.votes = set([])

    # def __cmp__(self, other):
    #     return cmp(self.go_id, other.go_id)


    def __eq__(self, other):
        """Override the default Equals behavior"""
        return self.go_id == other.go_id

    def __hash__(self):
        return(self.go_id.__hash__())

    def __repr__(self):
        return(self.go_id)
        #+ ': ' + self.name)

    def __str__(self):
        return(":"+self.go_id)

    def get_id(self):
        return self.go_id

    def map_genes(self, id_name, xdb_prefixed=False):
        """Map gene ids"""
        mapped_annotations_set = set([])
        for annotation in self.annotations:
            if xdb_prefixed:
                mapped_genes = id_name.get(annotation.xdb + ':' + annotation.gid)
            else:
                mapped_genes = id_name.get(annotation.gid)

            if mapped_genes is None and 'CELE_' in annotation.gid:
                mapped_genes = id_name.get(
                    annotation.gid[5:len(annotation.gid)])

            if mapped_genes is None:
                logger.warning('No matching gene id: %s', annotation.gid)
                continue

            for mgene in mapped_genes:
                mapped_annotations_set.add(Annotation(xdb=None, gid=mgene,
                                                      direct=annotation.direct,
                                                      ref=annotation.ref,
                                                      evidence=annotation.evidence,
                                                      date=annotation.date,
                                                      cross_annotated=annotation.cross_annotated))
        self.annotations = mapped_annotations_set

    def filter_annotations(self, evidence_codes):
        """Filter gene annotations by their evidence code"""
        for annotation in list(self.annotations):
            if annotation.evidence not in evidence_codes:
                self.remove_annotation(annotation)

    def get_annotated_genes(self, include_cross_annotated=True):
        genes = []
        for annotation in self.annotations:
            if (not include_cross_annotated) and annotation.cross_annotated:
                continue
            genes.append(annotation.gid)
        return genes

    def remove_annotation(self, annot):
        try:
            self.annotations.remove(annot)
        except KeyError:
            return

    def add_annotation(self, gid, ref=None, cross_annotated=False,
                       allow_duplicate_gid=True, origin=None, ortho_evidence=None):
        if not allow_duplicate_gid:
            for annotated in self.annotations:
                if annotated.gid == gid:
                    return
        self.annotations.add(
            Annotation(
                gid=gid,
                ref=ref,
                cross_annotated=cross_annotated,
                origin=origin,
                ortho_evidence=ortho_evidence))

    def get_annotation_size(self):
        return len(self.annotations)

    def get_namespace(self):
        return self.namespace

    def get_xrefs(self, dbid):
        if dbid in self.xrefs:
            return self.xrefs[dbid]
        else:
            return None



def go_term_id_comparison(go_term_x, go_term_y):
    '''Compares two strings according to the current LC_COLLATE setting.
    As any other compare function, returns a negative, or a positive value, or 0,
    depending on whether string1 collates before or after string2 or is equal to it.
    '''
    return locale.strcoll(go_term_x.go_id, go_term_y.go_id)


def cmp_to_key(mycmp):
    'Convert a cmp= function into a key= function'
    class K:
        def __init__(self, obj, *args):
            self.obj = obj
        def __lt__(self, other):
            return mycmp(self.obj, other.obj) < 0
        def __gt__(self, other):
            return mycmp(self.obj, other.obj) > 0
        def __eq__(self, other):
            return mycmp(self.obj, other.obj) == 0
        def __le__(self, other):
            return mycmp(self.obj, other.obj) <= 0
        def __ge__(self, other):
            return mycmp(self.obj, other.obj) >= 0
        def __ne__(self, other):
            return mycmp(self.obj, other.obj) != 0
    return K



if __name__ == '__main__':
    from argparse import ArgumentParser
    import sys

    usage = "usage: %(prog)s [options]"
    parser = ArgumentParser(prog=usage)
    parser.add_argument('-o',
        '--obo_file',
        dest='obo_file',
        help='Disease Ontology obo file')
    parser.add_argument("-v", "--verbose", dest="verbose", action='store_true',
                        help="output debug loglevel")
    parser.add_argument('-V', '--version', action='version',
                        version="%(prog)s dev-unreleased")
    args = parser.parse_args()

    if args.verbose:  # Setup logging at desired level
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.WARNING)

    logger.debug("Args: %s", args)

    if args.obo_file is None:
        sys.stderr.write("--obo_file file is required.\n")
        sys.exit()
    obo_file=args.obo_file
    logger.info('Loading disease obo from %s', obo_file)
    do = OBO()
    do.load_obo(obo_file=obo_file)
    print("Loaded")
    do.populated = True  # mark annotated
    do.propagate()  # prop annotations
    print("Populated")
