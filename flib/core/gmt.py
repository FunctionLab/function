import sys
import os
from collections import defaultdict


class GMT:

    """
    Pass the filename of gmt file.
    """

    def __init__(self, filename=None):
        self._genesets = defaultdict(set)
        self._setnames = {}
        self._genes = set()

        if filename:
            gmtfile = open(filename)

            for line in gmtfile:
                tok = line.strip().split('\t')
                (gsid, name, genes) = tok[0], tok[1], tok[2:]
                self._genesets[gsid] = set(genes)
                self._setnames[gsid] = name
                self._genes |= self._genesets[gsid]

    def ids(self):
        return self.genesets.keys()

    @property
    def genesets(self):
        return self._genesets

    @property
    def genes(self):
        return self._genes

    @property
    def setnames(self):
        return self._setnames

    def get_genes(self, gsid):
        return self._genesets[gsid]

    def add_geneset(self, gsid=None, name=None):
        self._setnames[gsid] = name
        self._genesets[gsid] = set()

    def add_gene(self, gsid, gene):
        self._genesets[gsid].add(gene)

    def add(self, gmt):
        for gsid, genes in gmt.genesets.iteritems():
            self._genesets[gsid] |= genes
        for gsid, name in gmt.setnames.iteritems():
            if gsid not in self._setnames:
                self._setnames[gsid] = gmt.setnames[gsid]

    def write(self, outfile):
        outf = open(outfile, 'w')
        for gsid, genes in self._genesets.iteritems():
            outf.write(
                gsid +
                '\t' +
                self._setnames[gsid] +
                '\t' +
                '\t'.join(
                    list(genes)) +
                '\n')
        outf.close()

    def __repr__(self):
        return self._genesets.__repr__()


if __name__ == '__main__':
    from optparse import OptionParser

    usage = "usage: %prog [options]"
    parser = OptionParser(usage, version="%prog dev-unreleased")
    parser.add_option(
        "-i",
        "--gmt-file",
        dest="gmt",
        help="GMT file",
        metavar="FILE")
    parser.add_option(
        "-d",
        "--directory",
        dest="dir",
        help="directory of gene files",
        metavar="FILE")
    parser.add_option(
        "-c",
        "--gmt-comp-file",
        dest="gmt2",
        help="GMT file to combine",
        metavar="FILE")
    parser.add_option(
        "-o",
        "--output-gmt-file",
        dest="output",
        help="GMT file to output",
        metavar="FILE")
    parser.add_option(
        "-g",
        "--gene-file",
        dest="genef",
        help="gene file for comparision",
        metavar="FILE")

    (options, args) = parser.parse_args()

    if options.dir:
        gs = GMT()
        for f in os.listdir(options.dir):
            lines = open(options.dir + '/' + f).readlines()
            for l in lines:
                gs.genesets[f].add(l.strip())
            gs.setnames[f] = f
        gs.write(options.output)

    else:
        # lines = open(options.genef).readlines()
        # genes = set()
        # for l in lines:
        #    genes.add(l.strip())

        gs = GMT(options.gmt)
        for gname, gset in gs.overlap().iteritems():
            for gname2, ovlp in gset.iteritems():
                print(("%s %s %s") % (gname, gname2, ovlp))
