import math
import numpy as np

class GeneExpression:

    dataset = ''            # name of data set (if any)
    genes = {}              # a dictionary of gene names to profile matrix index
    matrix = None           # a numpy two-dim array holding all expression values
    headers = []            # the names of the samples/experiments, e.g. GSM123
    default_value_if_null = None # Default value to use if entry is not set (e.g. addSamples may not add values for all genes)
    
    def __init__(self, datasetname='', headerlist=[], genedict={}):
        """ Create a gene expression data set.
            The class will store gene names and associated profiles (in which values correspond to "samples"). 
            It also stores headers (which names correspond to samples, i.e. experiments). 
            Data should be provided as 
            (0) a name of the set
            (1) a list of sample names (headerlist; must agree with the number of values in each gene profile)
            (2) a gene name dictionary where values contain the expression profile (genedict; profile is an iterable with the same number of elements)
            For example
            >>> g = GeneExpression("MySet", ['Sample1', 'Sample2'], {'G1': [0.13, 1.23], 'G2': [4.1, -0.9], 'G3': [2.1, -2.1]})
        """
        self.dataset = datasetname
        self.genes = {}
        ndx = 0
        for gene in genedict:
            self.genes[gene] = ndx
            ndx += 1 
        self.matrix = self._createMatrix(genedict)
        if len(self.matrix) == 0:
            nsamples = 0
        else:
            nsamples = len(self.matrix[0])  
        if isinstance(headerlist, str):
            headerlist = [headerlist]
        if len(headerlist) != nsamples:
            raise RuntimeError("The number of headers (%d) is not equal to the number of samples (%d)" % (len(headerlist), nsamples))
        self.headers = headerlist or ['S%d' % cnt for cnt in range(nsamples)]
        
    def _createMatrix(self, genedict):
        """ Internal method for constructing a numpy matrix from a gene-profile dictionary. """
        ngenes = len(self.genes)
        allow_new_genes = False
        if ngenes == 0:  # if instance is empty, include all genes in dict
            ngenes = len(genedict)
            allow_new_genes = True 
        nsamples = 0
        for gene in genedict:
            profile = genedict[gene]
            try:
                actual = len(profile)
            except TypeError:
                actual = 1
                genedict[gene] = [profile]
            if nsamples == 0:
                nsamples = actual
            elif nsamples != actual:
                raise RuntimeError("Each gene must have the same number of samples (see %s)" % gene)
        matrix = np.empty((ngenes, nsamples))
        matrix[:] = self.default_value_if_null
        ndx = 0
        for gene in genedict:
            try:
                ndx = self.genes[gene]
                matrix[ndx] = genedict[gene]
            except: # no match in current gene list
                if allow_new_genes:
                    matrix[ndx] = genedict[gene]
                    self.genes[gene] = ndx
                    ndx += 1
        return matrix
     
    def getHeaders(self, indices = None):
        """ Retrieve headers (names of experiments/samples).
            If indices is None (default), all headers are returned, e.g.
            >>> g.getHeaders()
            ['Sample1', 'Sample2']
            If indices is a single integer, the header for the corresponding entry is returned, e.g.
            >>> g.getHeaders(1)
            'Sample2'
            If indices is an iterable of integers (multiple indices), the list of corresponding headers is returned, e.g.
            >>> g.getHeaders([1,0])
            ['Sample2', 'Sample1']
        """
        if indices == None: 
            return self.headers
        elif isinstance(indices, int) or indices is slice:
            return self.headers[indices]
        else:
            ret = []
            for index in indices:
                ret.append(self.headers[index]) 
            return ret

    def getGenes(self, names = None):
        """ Retrieve applicable gene-profile entries.
            If names is None (default), all gene names are returned, e.g.
            >>> g.getGenes()
            ['G1', 'G2', 'G3']
            If names is a single string, the profile for the corresponding entry is returned, e.g. 
            >>> g.getGenes('G2')
            array([ 4.1, -0.9])
            If names is an iterable of strings (multiple gene names), a dictionary with gene name as key and profile as value is returned.
            >>> g.getGenes(['G3','G2'])
            {'G2': array([ 4.1, -0.9]), 'G3': array([ 2.1, -2.1])}
            """
        if names == None:
            return list(self.genes.keys())
        elif isinstance(names, str):
            return self.matrix[self.genes[names],:]
        else:
            ret = {}
            for name in names:
                ret[name] = self.matrix[self.genes[name],:] 
            return ret

    def __getitem__(self, ndx):
        """ Retrieve a specified sample (or a "slice" of samples) for all genes, e.g.
        >>> g[0:2]
            array([[ 2.1 , -2.1 ],
                   [ 4.1 , -0.9 ],
                   [ 0.13,  1.23]])
            Note that the order of rows/genes is NOT necessarily the same as that used for inserting the data. 
        """
        return self.matrix[:,ndx]
    
    def getHeaderIndex(self, headers):
        """ Find the index of the named experiment.
            Raises a ValueError if not in list. """
        if isinstance(headers, str):
            return self.headers.index(headers)
        else:
            return [self.headers.index(header) for header in headers]
    
    def getSamples(self, samples):
        """Construct a gene dictionary including only samples in specified indices, e.g.
        >>> g.getSamples(0)
            {'G1': 0.13, 'G2': 4.0999999999999996, 'G3': 2.1000000000000001}
        >>> g.getSamples('Sample2')
            {'G1': 1.23, 'G2': -0.90000000000000002, 'G3': -2.1000000000000001}
        >>> g.getSamples(['Sample2','Sample1'])
            {'G1': array([ 1.23,  0.13]),
             'G2': array([-0.9,  4.1]),
             'G3': array([-2.1,  2.1])}
        """
        try:
            index = self.getHeaderIndex(samples)
        except:
            index = samples
        mygenes = {}
        for (name, ndx) in list(self.genes.items()):
            mygenes[name] = self.matrix[ndx, index]
        return mygenes

    def sort(self, sample, descending=True):
        """Get a list of gene names, sorted by order of value in specified sample, e.g.
        >>> g.sort(0)
            ['G2', 'G3', 'G1']
        Then retrieve actual genes using e.g. 
        >>> g.getGenes('G2')
            array([-0.9,  4.1])
        """
        try:
            index = self.getHeaderIndex(sample)
            sort_ndx = np.nan_to_num(self.matrix[:,index]).argsort()
        except:
            sort_ndx = np.nan_to_num(self.matrix[:,sample]).argsort()
        name_tuples = sorted(list(self.genes.items()), key=lambda v: v[1]) # put all gene names in order of the matrix of profiles
        names = []
        if descending:
            for (name, index) in [name_tuples[index] for index in sort_ndx[::-1]]: # reverse the order 
                names.append(name)
        else:
            for (name, index) in [name_tuples[index] for index in sort_ndx]: # maintain order 
                names.append(name)
        return names

    def addSamples(self, headerlist, genedict):
        """Add a sample or multiple samples to the current data set.
           genedict is a dictionary with the same keys as the current gene set.
           Only values for genes in the current set will be added (others are ignored).
           >>> g.addSamples('Sample3', {'G1': 3.4, 'G2': -3.0})
        """
        newmat = self._createMatrix(genedict)
        nsamples = len(newmat[0])  
        if headerlist != None:
            if isinstance(headerlist, str):
                headerlist = [headerlist]
            if len(headerlist) != nsamples:
                raise RuntimeError("The number of headers (%d) is not equal to the number of samples (%d)" % (len(headerlist), nsamples))
        if len(self.matrix) == 0:
            self.matrix = newmat
        else:
            self.matrix = np.hstack((self.matrix, newmat))
        self.headers.extend(headerlist or ['S%d' % cnt+len(self.headers) for cnt in range(nsamples)])

    def getRatio(self, index1, index2):
        """ Get the ratio of two samples in the data set (index1 and index2).
            Creates and returns a gene dictionary with the corresponding ratios. """
        mygenes = {}
        mdiv = self.matrix[:, index1] / self.matrix[:, index2]
        for (name, ndx) in list(self.genes.items()):
            mygenes[name] = mdiv[ndx]
        return mygenes

    def getLogRatio(self, index1, index2):
        """ Get the log2-transformed ratio of two samples (index1 and index2)
            Creates and returns a gene dictionary with the corresponding log-ratios. """
        mygenes = {}
        mlr = np.log2(self.matrix[:, index1] / self.matrix[:, index2])
        for (name, ndx) in list(self.genes.items()):
            mygenes[name] = mlr[ndx]
        return mygenes
    
    def getPearson(self, probeID):
        """ Given a probe identifier, returns a gene/probe dictionary:
        identifiers to correlation coefficients with the specified probe. """
        index = self.genes[probeID]
        profile = self.matrix[index, :]
        mygenes = {}
        for (name, ndx) in list(self.genes.items()):
            other = self.matrix[ndx, :]
            mygenes[name] = pearson(profile, other)
        return mygenes
    
    def writeGEOFile(self, filename):
        """ Save data as a truncated GEO SOFT file named filename. """
        line = '^DATASET = ' + self.dataset + '\n'
        line += '!dataset_table_begin\nID_REF\tIDENTIFIER\t'
        for header in self.headers:
            line += header + '\t'
        line += '\n'
        for gene in self.genes:
            line += gene + '\t' + gene + '\t'
            index = self.genes[gene]
            for value in self.matrix[index, :]:
                line += format(value, '5.3f') + '\t'
            line += '\n'
        line += '!dataset_table_end\n'
        fh = open(filename, 'w')
        fh.write(line)
        fh.close()
    
    def getZScore(self, index):
        """ Get the Z-score of each expression value.
            index can be a list of indices (for which the z-score is computed independently).
            Important: assumes that values are normally distributed.
            For example use log-transformed ratios. """
        # Calculate mean and standard deviation of the list of values
        mu = np.mean(self.matrix[:, index], axis=0)
        sd = np.std(self.matrix[:, index], axis=0)
        # Calculate Z-score for the given column for each gene
        zscore = (self.matrix[:, index] - mu) / sd
        mygenes = {}
        for (name, ndx) in list(self.genes.items()):
            try:
                mygenes[name] = zscore[ndx, :]
            except IndexError:
                mygenes[name] = zscore[ndx]
        # Return the dictionary of Z-scores
        return mygenes

# Utility functions

def readGEOFile(filename, id_column=0):
    """Read a Gene Expression Omnibus file; return a GeneExpression instance.
    id_column indicates what field of each row that should be taken as
    gene identifier.
    """
    fh = open(filename, "rU")
    manylines = fh.read()
    fh.close()
    # If True, ignore genes with null samples; if False, use default value
    ignore_gene_if_null = False
    default_value_if_null = None # Default value to use if entry is null and not ignored
    # Indicates whether we're reading the data section or metadata
    data_rows = False
    cnt_data = 0
    cnt_null = 0
    dataset = ''    # name of dataset
    headers = []    # list of headers
    genes = {}      # dict with gene-name as key, expression profile as a list of floats
    for line in manylines.splitlines():
        if line.startswith('^DATASET'):
            dataset = line.split('= ')[1]
            continue
        if line.startswith('!dataset_table_begin'):
            data_rows = True
            continue
        if line.startswith('!dataset_table_end'):
            data_rows = False
            continue
        if line.startswith('!') or line.startswith('#') \
            or line.startswith('^'):
                continue
        if len(line.strip()) == 0:
            continue
        if data_rows:
            cnt_data += 1
            ignore = False
            name = line.split('\t')[id_column]
            # Ignore control probes
            if name.startswith("AFFX"):
                continue
            if (cnt_data == 1):  # First line contains the headers
                headers = line.split('\t')
            else:
                values = []
                cnt_word = 0
                for word in line.split('\t'):
                    cnt_word += 1
                    if cnt_word <= (id_column + 1):
                        continue
                    if word == 'null':
                        cnt_null += 1
                        if ignore_gene_if_null:
                            ignore = True
                            break
                        else:
                            word = default_value_if_null
                    try:
                        if word == None:
                            values.append(None)
                        else:
                            values.append(float(word))
                    except:  # ignore values that are not "float"
                        continue
                if ignore:
                    pass
                elif not name in genes:
                    genes[name] = values
    if len(genes) == 0:
        raise RuntimeError('No data in file')
    print('Data set %s contains %d genes' % (dataset, len(genes)))
    if cnt_null > 0:
        print('Data set has %d null-values' % (cnt_null))
    return GeneExpression(dataset, headers[2:], genes)


# ------------------ Helpful Extra Functions ------------------

def pearson(X, Y):
    """ Pearson correlation coefficient (r).
    Note that we are using the standard deviation of the sample, NOT the sample
    standard deviation (see http://en.wikipedia.org/wiki/Standard_deviation). """
    Xmu = np.mean(X)
    Xvar = np.var(X)
    Ymu = np.mean(Y)
    Yvar = np.var(Y)
    if len(X) != len(Y):
        raise RuntimeError('vectors are of uneven length')
    n = len(X)
    sum = 0.0
    for i in range(n):
        sum += (X[i] * Y[i])
    if n == 0 or Xvar == 0 or Yvar == 0:
        return 0
    return (sum - n * (Xmu * Ymu)) / (n * math.sqrt(Xvar) * math.sqrt(Yvar))

# ------------------- Example (basically exercise 7 in prac 9)---------------------

if __name__=='__main__':

    g = readGEOFile('GDS3198.soft', id_column = 1)
    meanfold = {}
    for gene in g.genes:
        profile = g.getGenes(gene)
        meanfold[gene] = (np.log2(profile[0] / profile[3]) + np.log2(profile[1] / profile[4]) + np.log2(profile[2] / profile[5])) / 3
    
    import matplotlib.pyplot as plt
    scores = [y for y in list(meanfold.values()) if not np.isnan(y)]
    hist, bins = np.histogram(scores, bins=50)
    width = 0.7 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2
    plt.bar(center, hist, align='center', width=width)
    plt.show()

    result = sorted(list(meanfold.items()), key=lambda v: v[1])
    print('========== Wildtype may down-regulate ==========')
    for r in result[0:100]:
        print(r[0], r[1])
    print('========== Wildtype may up-regulate ==========')
    for r in result[-1:-100:-1]:
        print(r[0], r[1])
    
