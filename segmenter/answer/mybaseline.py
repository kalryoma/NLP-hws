import sys, codecs, optparse, os

# optparser : parse command-line
optparser = optparse.OptionParser()
optparser.add_option("-c", "--unigramcounts", dest='counts1w', default=os.path.join('data', 'count_1w.txt'), help="unigram counts")
optparser.add_option("-b", "--bigramcounts", dest='counts2w', default=os.path.join('data', 'count_2w.txt'), help="bigram counts")
optparser.add_option("-i", "--inputfile", dest="input", default=os.path.join('data', 'input'), help="input file to segment")
(opts, _) = optparser.parse_args()

# this is for baseline???
# Pdist is a python dictionary extract from counts1w or counts2w
class Pdist(dict):
    "A probability distribution estimated from counts in datafile."
    # read count1
    def __init__(self, filename, sep='\t', N=None, missingfn=None):
        self.maxlen = 0 
        for line in file(filename):
            # count1: key freq
            (key, freq) = line.split(sep)
            # convert key to utf8
            try:
                utf8key = unicode(key, 'utf-8')
            except:
                raise ValueError("Unexpected error %s" % (sys.exc_info()[0]))
            # add up numbers of current word
            # python2 dict.get(key, default) return x[key] or default
            self[utf8key] = self.get(utf8key, 0) + int(freq)
            # get the maxlen
            self.maxlen = max(len(utf8key), self.maxlen)
        # total number of freq
        # what does "or" do???
        self.N = float(N or sum(self.itervalues()))
        # for a single word, posibility = 1 / N
        self.missingfn = missingfn or (lambda k, N: 1./N)

    def __call__(self, key):
        if key in self: return float(self[key])/float(self.N)
        #else: return self.missingfn(key, self.N)
        elif len(key) == 1: return self.missingfn(key, self.N)
        else: return None

# the default segmenter does not use any probabilities, but you could ...
Pw  = Pdist(opts.counts1w)



old = sys.stdout
sys.stdout = codecs.lookup('utf-8')[-1](sys.stdout)
# ignoring the dictionary provided in opts.counts
with open(opts.input) as f:
    for line in f:
        utf8line = unicode(line.strip(), 'utf-8')
        output = [i for i in utf8line]  # segmentation is one word per character in the input
        print " ".join(output)
sys.stdout = old
