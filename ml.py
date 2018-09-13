import numpy
import numpy.random
import math
import random

class NN():
    """
    A basic implementation of a standard, multi-layer, feed-forward neural network
    and back-propagation learning.
    """
    def __init__(self, nInput, nHidden, nOutput):
        """ Constructs a neural network and initializes its weights to small random values.
            nInput  Number of input nodes
            nHidden Number of hidden nodes
            nOutput Number of output nodes
        """
        self.ninput = nInput
        self.hidden = numpy.empty(nHidden)                  # hidden nodes
        self.output = numpy.empty(nOutput)                  # output nodes
        self.w_hid  = numpy.random.randn(nHidden, nInput)   # weights in -> hid
        self.b_hid  = numpy.random.randn(nHidden)           # biases hidden layer
        self.w_out  = numpy.random.randn(nOutput, nHidden)  # weights hid -> out
        self.b_out  = numpy.random.randn(nOutput)           # biases output layer
        print("Constructed NN with %d inputs, %d hidden and %d output nodes." % (self.ninput, len(self.hidden), len(self.output)))

    def writeFile(self, filename):
        """ Save NN to a file. """
        f = open(filename, 'w')
        f.write(str(self.ninput)+'\n')
        f.write(str(len(self.hidden))+'\n')
        f.write(str(len(self.output))+'\n')
        for row in self.w_hid:
            for w in row:
                f.write(str(w)+'\n')
        for b in self.b_hid:
            f.write(str(b)+'\n')
        for row in self.w_out:
            for w in row:
                f.write(str(w)+'\n')
        for b in self.b_out:
            f.write(str(b)+'\n')
        f.close()

    def _fLogistic(self, net):
        """ The logistic output function.
            Computes the output value of a node given the summed incoming activation,
            values bounded between 0 and 1.
            net: The summed incoming activation. """
        return 1.0 / (1.0 + numpy.exp(-net))

    def _fSoftmax(self, net):
        """ The softmax output function.
            Computes the output value of a node given the summed incoming activation,
            values bounded between 0 and 1, where all add to 1.0.
            net: The summed incoming activation for each output (must be the full layer). """
        tmp = numpy.exp(net)
        sum = numpy.sum(tmp)
        out = tmp / sum
        return out

    def _fprimeLogistic(self, y):
        """ The derivative of the logistic output function.
            y: The value by which the gradient is determined.
            returns the gradient at output y. """
        return y * (1.0 - y)

    def feedforward(self, input):
        """ Computes the output values of the output nodes in the network given input values.
            input: the one-dim array of input values
            returns the one-dim array of computed output values. """
        # compute the activation of each hidden node (depends on supplied input values)
        self.hidden = self._fLogistic(self.w_hid.dot(input) + self.b_hid)
        # compute the activation of each output node (depends on hidden node activations computed above)
        if len(self.output) == 1:
            self.output = self._fLogistic(self.w_out.dot(self.hidden) + self.b_out)
        else:
            self.output = self._fSoftmax(self.w_out.dot(self.hidden) + self.b_out)
        return self.output

    def test(self, inputs, targets):
        """ Create a confusion matrix for all predictions with known target classes. """
        cm = numpy.zeros((len(self.output), len(self.output))) # confusion matrix
        for p in range(len(inputs)):
            input   = inputs[p]
            target  = targets[p]
            # present the input and calculate the outputs
            output = self.feedforward(input)
            # which class?
            c_targ = maxIndex(target)
            c_pred = maxIndex(output)
            cm[c_targ, c_pred] += 1
        return cm

    def train(self, input, target, eta = 0.1, niter = 1, shuffle = True):
        """ Adapts weights in the network given the values that should appear at the output (target)
            when the input has been presented. The procedure is known as error back-propagation.
            This implementation is "online" rather than "batched", that is, the change is not based
            on the gradient of the golbal error, merely the local, pattern-specific error.
            target:  The desired output values
            eta:     The learning rate, always between 0 and 1, typically a small value (default 0.1)
            shuffle: If true, input rows are shuffled before training (reduces bias imposed by order
            in online training)
            returns an error value (the root-mean-squared-error). """
        try:
            len(input[0])
            multi_input = input
            multi_targ  = target
        except TypeError:
            multi_input = [ input ]
            multi_targ  = [ target ]
        for i in range(niter):
            mse = 0.0
            entries = list(range(len(multi_input)))
            if shuffle:
                random.shuffle(entries)
            for p in entries:
                input = multi_input[p]
                target  = multi_targ[p]
                # present the input and calculate the outputs
                self.feedforward(input)
                # compute the error of output nodes (explicit target is available -- so quite simple)
                # also, calculate the root-mean-squared-error to indicate progress
                dif_out = (target - self.output)
                if len(self.output) == 1:
                    err_out = dif_out * self._fprimeLogistic(self.output)
                else:
                    err_out = dif_out #* self._fprimeSoftmax(self.output)
                # compute the error of hidden nodes (indirect contribution to error at output layer)
                err_hid = self.w_out.T.dot(err_out) * self._fprimeLogistic(self.hidden)
                # change weights according to errors
                self.w_out += numpy.outer(err_out, self.hidden) * eta
                self.b_out += err_out * eta
                self.w_hid += numpy.outer(err_hid, input) * eta
                self.b_hid += err_hid * eta
                if i == niter - 1: # last round
                    mse += float(numpy.mean(numpy.square(dif_out)))
        return math.sqrt(mse / len(entries)) # Root of mean squared error (RMSE)

def readNNFile(filename):
    """ Load a NN from a file. """
    f = open(filename, 'r')
    nInput = int(f.readline())
    nHidden = int(f.readline())
    nOutput = int(f.readline())
    nn = NN(nInput, nHidden, nOutput)
    for i in range(nHidden):
        for j in range(nInput):
            nn.w_hid[i, j] = float(f.readline())
    for i in range(nHidden):
        nn.b_hid[i] = float(f.readline())
    for i in range(nOutput):
        for j in range(nHidden):
            nn.w_out[i, j] = float(f.readline())
    for i in range(nOutput):
        nn.b_out[i] = float(f.readline())
    f.close()
    return nn

def maxIndex(output):
    """ Figure out the index of the largest value in the specified array/list. """
    if len(output) > 1: # multi-class
        max = 0
        for i in range(len(output)):
            if output[i] > output[max]:
                max = i
    else:               # two-class, single output 0/1
        max = int(round(output[0]))
    return max

def Qk(cm, alpha):
    """ Compute the Q accuracy from a confusion matrix (see test method above) """
    Q = {}
    for a in alpha:
        i = alpha.index(a)
        Q[a] = (cm[i, i] / numpy.sum(cm[i])) * 100
    tp = 0; pos = 0
    for a in alpha:
        i = alpha.index(a)
        tp += cm[i, i]
        pos += sum(cm[i])
    return (float(tp) / float(pos)) * 100, Q

def readDenseDataFile(filename):
    """ Read data from file for training a neural network.
        The file follows the "dense" row-format:
        <i1> <i2> ... <im> | <o1> ... <on>
        where ix are m input values and ox are n output values """
    # first check format
    ninputs = None
    noutputs = None
    nexamples = 0
    f = open(filename)
    cnt = 0
    for row in f:
        cnt += 1
        inp, outp = row.split('|')
        indata = [ float(token) for token in inp.split() ]
        if ninputs:
            if len(indata) != ninputs:
                raise RuntimeError('Error reading file: Invalid input at row %d' % cnt)
        ninputs = len(indata)
        outdata = [ float(token) for token in outp.split() ]
        if noutputs:
            if len(outdata) != noutputs:
                raise RuntimeError('Error reading file: Invalid output at row %d' % cnt)
        noutputs = len(outdata)
    f.close()
    nexamples = cnt
    inm  = numpy.zeros((nexamples, ninputs))
    outm = numpy.zeros((nexamples, noutputs))
    f = open(filename)
    cnt = 0
    for row in f:
        inp, outp = row.split('|')
        inm[cnt]  = [ float(token) for token in inp.split()  ]
        outm[cnt] = [ float(token) for token in outp.split() ]
        cnt += 1
    f.close()
    return inm, outm


def fGaussian(x, mu = 0.0, sigma2 = 1.0):
    """ Gaussian PDF for numpy arrays """
    num = (x - mu) ** 2
    den = 2 * sigma2
    expon = numpy.exp(-num/den)
    return expon / numpy.sqrt(2.0 * numpy.pi * sigma2)

class KMeans():
    """
    K-means clustering is a special case of Expectation-Maximization (EM).
    In K-means clustering we consider samples
        x1,...,xn labeled with z1,...,zN with xt vectors in R^D and zt \in {1,...,K}.
    In other words, zt is a class label, or cluster label, for the data point xt.
    We can define a K-means probability model as follows where N(mu, I) denotes the
    D-dimensional Gaussian distribution with mean mu \in R^D and with the
    identity covariance matrix.
        theta = <mu_1,...,mu_K>, mu_k \in R^D
        P(x1, . . . , xn, z1, . . . zn) = PROD P(zt) P(xt|zt) = PROD 1/K N(mu_zt, I) (xt)
    We now consider the optimization problem defined for this model. For
    this model
        (mu_1,...,mu_K)* = argmin_mu min_z SUM || muzt - xt || ^2
    The optimization problem defines K-means clustering (under quadratic distortion).
    This problem is non-convex and in fact is NP-hard. The K-means algorithm is coordinate
    descent applied to this objective and is equivalent to EM under the above probability
    model. The K-means clustering algorithm can be written as follows where we specify a
    typical initialization step.
        1. Initialize mu_z to be equal to a randomly selected point xt.
        2. Repeat the following until (z1, . . . zn) stops changing.
            (a) zt   := argmin_z || mu_z - xt || ^2
            (b) Nz   := |{t: zt = z}|
            (c) mu_z := 1 / Nz SUM_t:zt=z xt
    In words, the K-means algorithm first assigns a class center mu_z for each class z.
    It then repeatedly classifies each point xt as belonging to the class whose center is
    nearest xt and then recomputes the class centers to be the mean of the point placed in that class.
    Because it is a coordinate descent algorithm, the sum of squares of the difference
    between each point and its class center is reduced by each update. This implies that the
    classification must eventually stabilize.
    The procedure terminates when the class labels stop changing.
    """
    def __init__(self, data):
        """ Construct a K-means classifier using the provided data.
            data: a two-dim numpy array, with one row corresponding to a data point.
            If training is not performed, the provided data is used as the "means". """
        assert len(data) > 0, "Data must be supplied"
        self.data  = data
        self.means = data.copy()
        self.samplelen = len(data[0])
        self.vars = numpy.empty((len(data), self.samplelen))

    def classify(self, sample):
        assert len(sample) == self.samplelen, "Sample vector has invalid length: " + str(len(sample))
        sqrdist = numpy.sum((self.means - sample) ** 2, 1)
        return sqrdist.argmin(0)

    def train(self, K):
        data = self.data
        N = len(data)
        clusters = numpy.zeros((N, 1))
        self.means = self.data[numpy.random.randint(N, size = K),:]  # pick K random samples
        while True:
            previous = clusters.copy()
            """ Compute cluster memberships GIVEN means """
            kdist = numpy.empty((len(data), K))
            for i in range(K):
                kdist[:,i] = numpy.sum((data - self.means[i]) ** 2, 1)
            clusters[:,0] = kdist.argmin(1)
            nsame = numpy.sum(previous[:,0] == clusters[:,0])
            if nsame == N:
                break
            """ Compute means GIVEN cluster memberships """
            for i in range(K):
                members = data[clusters[:,0] == i]
                self.means[i] = members.mean(0)  # mean over rows per column

def eucdist(v1, v2):
    diff = 0
    for i in range(len(v1)):
        diff += (v1[i] - v2[i])**2
    return math.sqrt(diff)




