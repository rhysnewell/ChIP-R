'''
Module with methods and classes for phylogeny.
@author: mikael
'''
import sequence

class PhyloTree:
    """ Rooted, binary (bifurcating) tree for representing phylogenetic relationships.
        Functionality includes labelling and traversing nodes; reading and writing to Newick format;
        association with sequence alignment; maximum parsimony inference of ancestral sequence;
        generation of single, bifurcating rooted tree by UPGMA.
        Known issues: Binary only; Parsimony does not handle gaps in alignment.
        Programmers should note that almost all functionality is implemented through recursion. """

    def __init__(self, root):
        """ Create a tree from a node that is "root" in the tree."""
        self.root = root

    def putAlignment(self, aln):
        """ Associate the tree with a set of sequences/alignment.
            Involves assigning the sequence to the leaf nodes. """
        self.aln = aln
        self.root._assignAlignment(aln)

    def __str__(self):
        """ Produce a printable representation of the tree, specifically the root of the tree. """
        return str(self.root)

    def strSequences(self, start = None, end = None):
        """ Produce a sequence representation of the tree, specifically the root of the tree.
            Specify the start and end positions in the alignment for the sequence to be printed
            (if None the min and max positions will be used). """
        if self.aln != None:
            my_start = start or 0
            my_end = end or self.aln.alignlen
            return self.root._printSequences(my_start, my_end)

    def findLabel(self, label):
        """ Retrieve/return the node with the specified label.
            Returns None if not found."""
        return self.root._findLabel(label)

    def getDescendantsOf(self, node, transitive = False):
        """ Retrieve and return the (list of) descendants (children) of a specified node.
            Node can be the label or the instance.
            transitive indicates if only the direct descendants (False) or if all descendants
            should be returned.
            If node does not exist, None is returned.
            If node has no descendants, an empty list will be returned."""
        if not isinstance(node, PhyloNode):
            node = self.findLabel(node)
        if node:
            return node.getDescendants(transitive)
        return None

    def getAncestorsOf(self, node, transitive = False):
        """ Retrieve and return the ancestor (transitive=False) or
            ancestors (transitive=True) of a specified node.
            Node can be the label or the instance.
            If node does not exist, None is returned.
            If node is the root of the tree, None is returned."""
        if not isinstance(node, PhyloNode):
            node = self.findLabel(node)
        if node:
            myroot = self.root
            found = False
            branching = []
            while not found and myroot != None:
                branching.append(myroot)
                # check if "myroot" is a leaf node, i.e. does not have children
                if myroot.left == node or myroot.right == node:
                    found = True
                    break
                if myroot.left != None: # myroot has a "left" child
                    # check if the "left" child of "myroot" is the ancestor of "node"
                    if myroot.left.isAncestorOf(node, transitive = True): # if yes,
                        myroot = myroot.left    # move to the "left" child
                    else:                       # if not,
                        myroot = myroot.right   # move to the "right" child
                else: # myroot does NOT have a "left" child, so let's move "right"
                    myroot = myroot.right
            if found and transitive:
                return branching
            elif found and len(branching) > 0:
                return branching[len(branching)-1]
            return None

    def parsimony(self):
        """ Solve the "small parsimony problem",
            i.e. find the sequences on each of the internal nodes.
            See Jones and Pevzner, p. 368 and onwards, for details. """
        self.root._forwardParsimony(self.aln)  # setup and compute scores for all nodes
        self.root._backwardParsimony(self.aln) # use scores to determine sequences
        return self.root.getSequence() # return the sequence found at the root

    def canonise(self):
        self.root._canonise()

class PhyloNode:
    """ A class for a node in a rooted, binary (bifurcating) tree.
        Contains pointers to descendants/daughters (left and right),
        optional fields include data, label, sequence and dist.
        If parsimony is used scores and traceback pointers are available.
        A number of methods are named with a _ prefix. These can be, but
        are not intended to be used from outside the class. """

    def __init__(self, label = ''):
        """ Initialise an initially unlinked node.
            Populate fields left and right to link it with other nodes.
            Set label to name it.
            Use field data for any type of information associated with node.
            Use dist to indicate the distance to its parent (if any).
            Other fields are used internally, including sequence for associated alignment,
            seqscores, backleft and backright for maximum parsimony. """
        self.left = None
        self.right = None
        self.data = None
        self.label = label
        self.dist = None
        self.sequence = None # The sequence after an alignment have been mapped (leaf) or the most parsimonous sequence (ancestral)
        self.seqscores = None # The scores propagated from leaves via children
        self.backleft = None # Pointers back to left child: what symbol rendered current/parent symbols
        self.backright = None # Pointers back to right child: what symbol rendered current/parent symbols

    def __str__(self):
        """ Returns string with node (incl descendants) in a Newick style. """
        left = right = label = dist = ''
        if self.left:
            left = str(self.left)
        if self.right:
            right = str(self.right)
        if self.dist or self.dist == 0.0:
            dist = ':' + str(self.dist)
        if self.label != None:
            label = str(self.label)
            if not self.left and not self.right:
                return label + dist
            else:
                return '(' + left + ',' + right + ')' + label + dist
        else: # there is no label
            if not self.left and self.right:
                return ','+right
            elif self.left and not self.right:
                return left+','
            elif self.left and self.right:
                return '(' + left + ',' + right + ')' + dist
                
    def __le__(self, other):
        """ Returns indication of less than other node. """
        return other and self.__hash__() <= other.__hash__()
    	
    def __eq__(self, other):
        """ Returns indication of equivalence to other node. """
        return other and self.__hash__() == other.__hash__()

    def __hash__(self):
        """ Returns hash of object. """
        return hash((self.label, self.dist, self.sequence))
        
    def _printSequences(self, start, end):
        """ Returns string with node (incl descendants) in a Newick style. """
        left = right = label = dist = ''
        if self.left:
            left = self.left._printSequences(start, end)
        if self.right:
            right = self.right._printSequences(start, end)
        if self.dist:
            dist = ':' + str(self.dist)
        if self.sequence != None:
            label = "".join(self.sequence[start:end]) + ""
            if not self.left and not self.right:
                return label + dist
            else:
                return '(' + left + ',' + right + ')' + label + dist
        else: # there is no label
            if not self.left and self.right:
                return ','+right
            elif self.left and not self.right:
                return left+','
            elif self.left and self.right:
                return '(' + left + ',' + right + ')' + dist

    def _findLabel(self, label):
        """ Find a node by label at this node or in any descendants (recursively). """
        if self.label == label:
            return self
        else:
            if self.left:
                foundLeft = self.left._findLabel(label)
                if foundLeft:
                    return foundLeft
            if self.right:
                return self.right._findLabel(label)
            return None

    def _propagateDistance(self, parent_dist):
        """ Convert absolute distances to relative.
            The only parameter is the absolute distance to the parent of this node. """
        travelled = self.dist               # absolute distance to this node
        self.dist = parent_dist - self.dist # relative distance to this node
        if self.left != None:               # if there is a child node...
            self.left._propagateDistance(travelled) # pass absolute distance to this node
        if self.right != None:
            self.right._propagateDistance(travelled)

    def _assignAlignment(self, aln):
        """ Assign an alignment to the node, which implies assigning a sequence to it if one is
            available in the alignment. """
        self.sequence = None
        if self.left != None:
            self.left._assignAlignment(aln)
        if self.right != None:
            self.right._assignAlignment(aln)
        for seq in aln.seqs:
            if seq.name == self.label:
                self.sequence = seq
                break

    def _canonise(self):
        if self.left == None and self.right == None: # at leaf
            return self.label
        myleft = self.left._canonise()
        myright = self.right._canonise();
        if myleft > myright:
            tmpnode = self.left
            self.left = self.right
            self.right = tmpnode
            return myright
        return myleft

    def _forwardParsimony(self, aln):
        """ Internal function that operates recursively to first initialise each node (forward),
            stopping only once a sequence has been assigned to the node,
            then to propagate scores from sequence assigned nodes to root (backward). """
        if self.sequence == None: # no sequence has been assigned
            if self.left == None and self.right == None:    # no children, so terminal, cannot propagate scores
                raise RuntimeError("No sequence assigned to leaf node:", self.label)
            scoresleft = scoresright = None
            if self.left != None:
                scoresleft = self.left._forwardParsimony(aln)
            if self.right != None:
                scoresright = self.right._forwardParsimony(aln)
            # for each position in the alignment,
            # introduce (initially zero) score for each symbol in alphabet
	#Project "Substitution weights" should focus on this line of code
            self.seqscores = [[0 for _ in aln.alphabet] for col in range(aln.alignlen)]
            # for each position in the alignment,
            # allocate a position to put the left child symbol from which each current node symbol score was determined
            self.backleft = [[None for _ in aln.alphabet] for _ in range(aln.alignlen)]
            # allocate a position to put the right child symbol from which each current node symbol score was determined
            self.backright = [[None for _ in aln.alphabet] for _ in range(aln.alignlen)]
            for col in range(aln.alignlen):
                for a_parent in range(len(aln.alphabet)):
                    best_score_left = +9999999
                    best_score_right = +9999999
                    best_symb_left = 0
                    best_symb_right = 0
                    for a_left in range(len(aln.alphabet)):
                        score = (scoresleft[col][a_left] + (1 if a_left != a_parent else 0)) # if we want to weight scores, this would need to change
                        if score < best_score_left:
                            best_symb_left = a_left
                            best_score_left = score
                    for a_right in range(len(aln.alphabet)):
                        score = (scoresright[col][a_right] + (1 if a_right != a_parent else 0)) # if we want to weight scores, this would need to change
                        if score < best_score_right:
                            best_symb_right = a_right
                            best_score_right = score
                    self.seqscores[col][a_parent] = best_score_left + best_score_right
                    self.backleft[col][a_parent] = best_symb_left
                    self.backright[col][a_parent] = best_symb_right
        else:
            self.seqscores = [[0 if a==sym else 999999 for a in aln.alphabet] for sym in self.sequence] # if we want to weight scores, this would need to change
        return self.seqscores

    def _backwardParsimony(self, aln, seq = None):
        """ Internal function that operates recursively to inspect scores to determine
            most parsimonious sequence, from root to leaves. """
        if self.sequence == None: # no sequence has been assigned
            leftbuf = []
            rightbuf = []
            if self.left == None and self.right == None:    # no children, so terminal, cannot propagate scores
                raise RuntimeError("No sequence assigned to leaf node:", self.label)
            if seq == None: # Only root can do this, no parents to consider, so we pick the lowest scoring symbol
                currbuf = []
                for col in range(aln.alignlen):
                    min_score = 999999
                    min_symb = None
                    left_symb = None
                    right_symb = None
                    for a_parent in range(len(aln.alphabet)):
                        if self.seqscores[col][a_parent] < min_score:
                            min_score = self.seqscores[col][a_parent]
                            min_symb = a_parent
                            left_symb = self.backleft[col][a_parent]
                            right_symb = self.backright[col][a_parent]
                    currbuf.append(aln.alphabet[min_symb])
                    leftbuf.append(aln.alphabet[left_symb])
                    rightbuf.append(aln.alphabet[right_symb])
                self.sequence = sequence.Sequence(currbuf, aln.alphabet, self.label, gappy = True)
            else: # Non-root, but not leaf
                self.sequence = seq
                col = 0
                for sym_parent in self.sequence:
                    a_parent = aln.alphabet.index(sym_parent)
                    left_symb = self.backleft[col][a_parent]
                    right_symb = self.backright[col][a_parent]
                    leftbuf.append(aln.alphabet[left_symb])
                    rightbuf.append(aln.alphabet[right_symb])
                    col += 1
            self.left._backwardParsimony(aln, sequence.Sequence(leftbuf, aln.alphabet, self.label, gappy = True))
            self.right._backwardParsimony(aln, sequence.Sequence(rightbuf, aln.alphabet, self.label, gappy = True))
        return self.sequence

    def getSequence(self):
        """ Get the sequence for the node. Return None if no sequence is assigned.
            Requires that an alignment is associated with the tree, and that sequence names match node labels.
            If the explored node is not a leaf, the sequence can be determined by parsimony. """
        if self.sequence != None: # a sequence has been assigned
            return self.sequence
        elif self.seqscores != None: # inferred by parsimony but not yet assigned
            return None # determine most parsimonous sequence, not yet implemented

    def isAncestorOf(self, node, transitive = True):
        """ Decide if this node is the ancestor of specified node.
            If transitive is True (default), all descendants are included.
            If transitive is False, only direct descendants are included. """
        if node == self.left or node == self.right:
            return True
        elif transitive:
            if self.left:
                statusLeft = self.left.isAncestorOf(node, transitive)
                if statusLeft: return True
            if self.right:
                return self.right.isAncestorOf(node, transitive)
        else:
            return False

    def getDescendants(self, transitive = False):
        """ Retrieve and return (list of) nodes descendant of this.
            If transitive is False (default), only direct descendants are included.
            If transitive is True, all descendants are (recursively) included. """
        children = []
        if self.left:
            children.append(self.left)
        if self.right:
            children.append(self.right)
        if not transitive:
            return children
        else:
            grandchildren = []
            for c in children:
                d = c.getDescendants(transitive)
                if d:
                    grandchildren.extend(d)
            children.extend(grandchildren)
            return children

""" ----------------------------------------------------------------------------------------
    Methods for generating a single tree by clustering, here UPGMA Zvelebil and Baum p. 278
    ----------------------------------------------------------------------------------------"""

def runUPGMA(aln, measure, absoluteDistances=False):
    """ Generate an ultra-metric, bifurcating, rooted tree from an alignment based on pairwise distances.
        Use specified distance metric (see sequence.calcDistances).
        If absoluteDistances is True, the tree will be assigned the total distance from provided species.
        Otherwise, the relative addition at each path will be assigned."""
    D = {}
    N = {}  # The number of sequences in each node
    M = aln.calcDistances(measure)  # determine all pairwise distances
    nodes = [PhyloNode(seq.name) for seq in aln.seqs]  # construct all leaf nodes
    """ For each node-pair, assign the distance between them. """
    for i in range(len(nodes)):
        nodes[i].sequence = aln.seqs[i]
        nodes[i].dist = 0.0
        N[nodes[i]] = 1  # each cluster contains a single sequence
        for j in range(0, i):
            D[frozenset([nodes[i], nodes[j]])] = M[i, j]
    """ Now: treat each node as a cluster,
        until there is only one cluster left,
        find the *closest* pair of clusters, and
        merge that pair into a new cluster (to replace the two that merged).
        In each case, the new cluster is represented by the (phylo)node that is formed. """
    while len(N) > 1:  # N will contain all "live" clusters, to be reduced to a signle below
        closest_pair = (None, None)  # The two nodes that are closest to one another according to supplied metric
        closest_dist = None  # The distance between them
        for pair in D:  # check all pairs which should be merged
            dist = D[pair]
            if closest_dist == None or dist < closest_dist:
                closest_dist = dist
                closest_pair = list(pair)
        # So we know the closest, now we need to merge...
        x = closest_pair[0]  # See Zvelebil and Baum p. 278 for notation
        y = closest_pair[1]
        z = PhyloNode()  # create a new node for the cluster z
        z.dist = D.pop(frozenset([x, y])) / 2.0  # assign the absolute distance, travelled so far, note: this will change to relative distance later
        Nx = N.pop(x)  # find number of sequences in x, remove the cluster from list N
        Ny = N.pop(y)  # find number of sequences in y, remove the cluster from list N
        dz = {}  # new distances to cluster z
        for w in N:  # for each node w ...
            # we will merge x and y into a new cluster z, so need to consider w (which is not x or y)
            dxw = D.pop(frozenset([x, w]))  # retrieve and remove distance from D: x to w
            dyw = D.pop(frozenset([y, w]))  # retrieve and remove distance from D: y to w
            dz[w] = (Nx * dxw + Ny * dyw) / (Nx + Ny)  # distance: z to w
        N[z] = Nx + Ny  # total number of sequences in new cluster, insert new cluster in list N
        for w in dz:  # we have to run through the nodes again, now not including the removed x and y
            D[frozenset([z, w])] = dz[w]  # for each "other" cluster, update distance per EQ8.16 (Z&B p. 278)
        z.left = x  # link the phylogenetic tree
        z.right = y
        nodes.append(z)
    if not absoluteDistances:
        x._propagateDistance(z.dist)  # convert absolute distances to relative by recursing down left path
        y._propagateDistance(z.dist)  # convert absolute distances to relative by recursing down right path
        z.dist = 0.0  # root z is at distance 0 from merged x and y
    return PhyloTree(z)  # make it to tree, return

""" ----------------------------------------------------------------------------------------
    Methods for processing files of trees on the Newick format
    ----------------------------------------------------------------------------------------"""

def _findComma(string, level = 0):
    """ Find first comma at specified level of embedding """
    mylevel = 0
    for i in range(len(string)):
        if string[i] == '(':
            mylevel += 1
        elif string[i] == ')':
            mylevel -= 1
        elif string[i] == ',' and mylevel == level:
            return i
    return -1

def parseNewickNode(string):
    """ Utility function that recursively parses embedded string using Newick format. """
    first = string.find('(')
    last = string[::-1].find(')') # look from the back
    if first == -1 and last == -1: # we are at leaf
        y = string.split(':')
        node = PhyloNode(y[0])
        if len(y) >= 2:
            node.dist = float(y[1])
        return node
    elif first >= 0 and last >= 0:
        # remove parentheses
        last = len(string) - last - 1 # correct index to refer from start instead of end of string
        embed = string[first + 1:last]
        tail = string[last + 1:]
        # find where corresp comma is
        comma = _findComma(embed)
        if comma == -1:
            raise RuntimeError('Invalid format: invalid placement of "," in sub-string "' + embed + '"')
        left = embed[0:comma].strip()
        right = embed[comma + 1:].strip()
        y = tail.split(':')
        node = PhyloNode(y[0]) #node is an instance of the PhyloNode() class
        if len(y) >= 2:
            node.dist = float(y[1])
        node.left = parseNewickNode(left)
        node.right = parseNewickNode(right)
        return node
    else:
        raise RuntimeError('Invalid format: unbalanced parentheses in sub-string "' + string + '"')

def parseNewick(string):
    """ Main method for parsing a Newick string into a (phylogenetic) tree.
        Handles labels (on both leaves and internal nodes), and includes distances (if provided).
        Returns an instance of a PhyloTree. """
    if string.find(';') != -1:
        string = string[:string.find(';')]
    return PhyloTree(parseNewickNode(string))

def readNewick(filename):
    """ Read file on Newick format.
        Returns an instance of a PhyloTree."""
    f = open(filename)
    string = ''.join(f)
    return parseNewick(string)

def writeNewickFile(filename, my_tree):
    with open(filename, 'w') as fh:
        print(my_tree, end="", file=fh)
