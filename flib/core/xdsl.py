'''
Created on Feb 4, 2013

@author: atadych
'''

from xml.etree import ElementTree
import sys
import numpy

def isodd(num):
    return num & 1 and True or False

class CptNodesHolder:
    '''
    Parses xdsl file and holds all the nodes as dictionary
    '''
    def __init__(self, filename=None):
        #dictionary of CptNode objects
        self.nodes = {}
        if filename:
            self.parse(filename)        
        
    def __str__(self):
        return self.nodes
    
    def __repr__(self):
        return self.__str__()
    
    def parse(self,filename):
        '''Parses the file and fills out the dictionary nodes'''
        try:
            tree = ElementTree.parse(filename)
            #smile entry
            root = tree.getroot()
            for entry in root:
                if entry.tag == 'nodes':
                    for elem in entry:
                        if elem.tag == 'cpt' and 'id' in elem.attrib:
                            node_id = elem.attrib['id']
                            node = CptNode(node_id=node_id)
                            for child in elem:
                                #add states
                                if child.tag == 'state' and 'id' in child.attrib:
                                    node.states.append(child.attrib['id'])
                                #add parents
                                elif child.tag == 'parents':
                                    node.parents.append(child.text)
                                elif child.tag == 'probabilities':                                    
                                    vals = self._parse_probabilities(child.text)
                                    node.neg_probabilities = vals[0]
                                    node.pos_probabilities = vals[1]
                            self.nodes[node_id] = node
                            del node
        except IOError:
            print 'Can\'t open the file: %s.\nPlease check if file exists and if it has right permission.'  % filename
                  
    def _parse_probabilities(self,prob_str):
        '''
        Splits the prob_str values into 2 lists of float values
        '''
        def_vals = [None,None]
        if prob_str == None or len(prob_str.strip())==0:
            return def_vals
        vals = prob_str.split()
        strlen = len(vals)
        if isodd(strlen):
            print 'Incorrect # of probabilities for this str: %s' % prob_str 
            return def_vals
        #Convert to float
        float_vals = [float(s) for s in vals]
        return [float_vals[0:strlen/2],float_vals[strlen/2:]]
            
    def get_node(self,node_id):
        ''' Returns CptNode obj '''
        if node_id in self.nodes:
            return self.nodes[node_id]

    def get_probabilities(self,node_id):
        ''' Returns list of probabilities for given node_id  '''        
        if node_id in self.nodes:
            return self.nodes[node_id].get_probabilities() 
    
    def get_nodes_ids(self):
        ''' Returns all nodes ids in the current holder '''        
        return self.nodes.keys()
    
    
    
class CptNode:
    '''
    Single node with states and probabilities
    '''
    def __init__(self,node_id=None):
        self.node_id = node_id
        self.parents = []
        #list of state ids
        self.states = []
        self.pos_probabilities = []
        self.neg_probabilities = []

    def __str__(self):
        string = '%s, states:%s, pos:%s, neg:%s' % (self.node_id, self.states, self.pos_probabilities, self.neg_probabilities)
        return string          

        
    def __repr__(self):
        return self.__str__()          

    def get_probabilities(self):
        return [self.neg_probabilities,self.pos_probabilities]

    def get_logratios(self):
        bineffects = []
        for (pneg, ppos) in zip(self.neg_probabilities, self.pos_probabilities):
            if ppos == 0 or pneg == 0:
                bineffects.append(0)
            else:
                bineffects.append(numpy.log(pneg/ppos))

        return bineffects


if __name__ == '__main__':
    from optparse import OptionParser

    usage = "usage: %prog [options]"
    parser = OptionParser(usage, version="%prog dev-unreleased")
    parser.add_option("-i", "--xdsl-file", dest="xdsl", help="XDSL file", metavar="FILE")
    parser.add_option("-o", "--out-file", dest="out", help="Output file", metavar="FILE")

    (options, args) = parser.parse_args()

    if options.xdsl is None:
        sys.stderr.write("--xdsl file is required.\n")
        sys.exit()   
        
    filename = options.xdsl     
    nodes = CptNodesHolder(filename=filename)
    
    if options.out:
        out_file = open(options.out, 'w')
        # such bad code...
        dsets = nodes.get_nodes_ids()
        for i in xrange(-1,7):
            for dataset in dsets:
                if i==-1:
                    out_file.write(dataset + '_neg' + '\t')
                    out_file.write(dataset + '_pos' + '\t')
                else:
                    node = nodes.get_node(dataset)
                    probs = node.get_probabilities()
                    if i < len(probs[1]):
                        out_file.write(str(probs[0][i]) + '\t')
                        out_file.write(str(probs[1][i]) + '\t')
                    else:
                        out_file.write('NA\tNA\t')

            out_file.write('\n')
        out_file.close()

    for s in nodes.get_nodes_ids():
        node = nodes.get_node(s)
        befs = node.get_logratios()

        llsum = 0
        for i in range(len(befs)):
            llsum += (befs[i])

        if max(befs) > 0:
            print s, llsum, [x for x in befs if x > 0]

