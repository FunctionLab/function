#!/usr/bin/python

from __future__ import division
from __future__ import print_function
#from builtins import str
#from builtins import range
#from builtins import object

import sys
import array
import logging
import math

logger = logging.getLogger(__name__)


class Dab(object):

    def __init__(self, filename):
        self.gene_list = []
        self.gene_table = {}
        if filename.endswith('.qdab'):
            self.open_file(filename, qdab = True)
        else:
            self.open_file(filename)
        self.gene_index = {}
        for i in range(len(self.gene_list)):
            self.gene_index[self.gene_list[i]] = i
        logger.debug("Got %s genes.", len(self.gene_list))

    def open_file(self, filename, qdab=False):
        logger.debug("Opening %s", filename)
        dab_file = open(filename, 'rb')

        # get number of genes
        a = array.array('I')
        a.fromfile(dab_file, 1)
        size = a[0]
        logger.debug("Expecting %s genes.", size)

        # get gene names
        start = 4
        end = 4
        count = 0
        while count < size:
            dab_file.seek(end)
            chrs = dab_file.read(2)
            if chrs == b'\x00\x00':
                dab_file.seek(start)

                gene = dab_file.read(end - start).decode().strip()
                gene = gene.replace('\x00','')

                self.gene_list.append(gene)
                self.gene_table[gene] = count

                start = end + 2
                count += 1
                end += 1
            if not count % 5000 and start == end:  # Read 5k and on to next set
                logger.debug("Read %s gene names.", count)
            end += 1

        if qdab:
            # get number of bins
            dab_file.seek(start)
            a = array.array('B')
            a.fromfile(dab_file, 1)
            nbins = a[0]
            logger.debug("Number of bins: %s .",nbins)
            start=start+1

            # get the bin boundaries
            a = array.array('f')
            a.fromfile(dab_file, nbins)
            boundaries = a
            logger.debug("Bin boundaries: %s .",boundaries)
            start=start+4*len(boundaries)

            # get number of bits (+1 for NaN)
            nbits = int(math.ceil(math.log(nbins+1,2)))
            nan_val = math.pow(2,nbits)-1
            logger.debug("Number of bits for each value: %s .",nbits)


            # get half matrix values
            total = size * (size - 1) // 2

            a = array.array('B')
            a.fromfile(dab_file,1)
            bufferA = a[0]
            a = array.array('B')
            a.fromfile(dab_file,1)
            bufferB = a[0]

            iTotal = 0
            self.dat = array.array('f')
            for i in range(size-1):
                for j in range(0,size-i-1):
                    try:
                        iPos = (iTotal *nbits) % 8
                        if iPos + nbits > 8:
                            btmpb = (bufferA << iPos)
                            btmpf = ((bufferB >> (16 - nbits - iPos)) << (8-nbits))
                            self.dat.append((((btmpb | btmpf)& 0x000000FF) >> (8 - nbits)))
                            bufferA = bufferB
                            a = array.array('B')
                            a.fromfile(dab_file,1)
                            bufferB = a[0]
                        else:
                            self.dat.append((((bufferA << iPos) & 0x000000FF) >> (8 - nbits)))
                            if iPos + nbits == 8:
                                bufferA = bufferB
                                a = array.array('B')
                                a.fromfile(dab_file,1)
                                bufferB = a[0]
                    except:
                        #check we are reaching the boundary of the file
                        assert iTotal - len(self.dat) <= 1 + 8 // nbits

                    iTotal = iTotal + 1
                    if self.dat[-1] ==  nan_val:
                        self.dat[-1] = float('nan')
        else:
            # get half matrix values
            total = size * (size - 1) // 2
            dab_file.seek(start)
            self.dat = array.array('f')
            self.dat.fromfile(dab_file, total)

        assert len(self.dat) == total

    def get_size(self):
        return len(self.gene_list)

    def get_gene(self, id):
        return self.gene_list[id]

    def get_value_genestr(self, gene1, gene2):
        g1 = self.get_index(gene1)
        g2 = self.get_index(gene2)
        if g1 is None or g2 is None:
            return None
        else:
            return self.get_value(g1, g2)

    def get_value(self, gene1, gene2):
        g1 = min(gene1, gene2)
        g2 = max(gene1, gene2)

        # index of first id
        start = self.arith_sum((len(self.gene_list)) - g1,
                               (len(self.gene_list) - 1))
        start += (g2 - g1) - 1  # index of second id
        try:
            v = self.dat[int(start)]
        except IndexError:
            print('Error: ', start, gene1, gene2)
            exit()

        return v

    def get_scaled_value(self, gene1, gene2, prior_new, prior_old):
        r = prior_new / prior_old
        r_diff = (1 - prior_new) / (1 - prior_old)
        weight = self.get_value(gene1, gene2)
        return weight * r / (weight * r + (1 - weight) * r_diff)

    def get_index(self, gene):
        try:
            return self.gene_index[gene]
        except KeyError:
            return None

    def arith_sum(self, x, y):
        return .5 * (y - x + 1) * (x + y)

    def print_table(self, out_file=sys.stdout):
        cols = ['GENE']
        cols.extend(self.gene_list)
        print("\t".join(cols), file=out_file)

        for i in range(0, self.get_size()):
            line = []
            line.append(self.gene_list[i])
            for j in range(0, i):
                v = self.get_value(i, j)
                line.append(str(v))
            line.append("1")
            for j in range(i + 1, self.get_size()):
                v = self.get_value(i, j)
                line.append(str(v))

            print("\t".join(line), file=out_file)

    def print_flat(self, out_file=sys.stdout):
        for i in range(0, self.get_size()):
            for j in range(i + 1, self.get_size()):
                print(self.gene_list[i] + '\t' +
                      self.gene_list[j] + '\t' + str(self.get_value(i, j)),
                      file=out_file)
    '''
    def get_neighbors(self, gene_str, cutoff):
        neighbors = set()
        gene_id = self.get_index(gene_str)
        if gene_id is None:
            return neighbors
        for i in range(0, len(self.gene_list)):
            if self.get_value(gene_id, i) > cutoff:
                neighbors.add(self.gene_list[i])
        return neighbors

    def get_all_neighbor_vals(self, gene_id):
        vals = list()
        # gene_id = self.get_index(gene_str)
        if gene_id is None:
            return vals

        for i in range(gene_id):
            vals.append(self.get_value(gene_id, i))
        for i in range(gene_id + 1, self.get_size()):
            vals.append(self.get_value(gene_id, i))

        return vals

    def get_all_scaled_neighbor_vals(self, gene_id, prior_new, prior_old):
        vals = list()
        if gene_id is None:
            return vals

        for i in range(gene_id):
            vals.append(self.get_scaled_value(gene_id, i, prior_new,
                                              prior_old))
        for i in range(gene_id + 1, self.get_size()):
            vals.append(self.get_scaled_value(gene_id, i, prior_new,
                                              prior_old))

        return vals

    def get_all_neighbor_val_dict(self, gene_id):
        n_vals = dict()

        if gene_id is None:
            return n_vals

        for i in range(gene_id):
            n_vals[i] = self.get_value(gene_id, i)
        for i in range(gene_id + 1, self.get_size()):
            n_vals[i] = self.get_value(gene_id, i)

        return n_vals
    '''
    def get(self, gene_str):
        vals = []
        idx = self.get_index(gene_str)
        if idx is None:
            return vals

        for i in range(0, idx):
            # Get values from 0 to idx (not including idx)
            v = self.get_value(i, idx)
            vals.append(v)

        # Append self interaction value
        vals.append(1)

        start = self.arith_sum((len(self.gene_list)) - idx,
                               (len(self.gene_list) - 1))
        vals += self.dat[int(start):int(start) + len(self.gene_list) - (idx+1)]

        return vals



if __name__ == '__main__':
    from argparse import ArgumentParser

    usage = "usage: %(prog)s [options]"
    parser = ArgumentParser(prog=usage)
    parser.add_argument("-i", "--dab-file", dest="dab", help="DAB file",
                        metavar="FILE")
    parser.add_argument("-o", "--output-file", dest="out",
                        help="Output file (DAT or PCL)", metavar="FILE")
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

    if args.dab is None:
        sys.stderr.write("--dab file is required.\n")
        sys.exit()
    pcl_out = args.out.endswith('.pcl')
    dat_out = args.out.endswith('.dat')
    if args.out is not None and not pcl_out and not dat_out:
        sys.stderr.write("Unknown file format for: " + args.out + "\n")
        sys.exit()

    dab = dat(args.dab)

    if args.out is None:
        dab.print_table()
    else:
        ofile = open(args.out, 'w')
        if pcl_out:
            dab.print_table(ofile)
        elif dat_out:
            dab.print_flat(ofile)
        ofile.close()
