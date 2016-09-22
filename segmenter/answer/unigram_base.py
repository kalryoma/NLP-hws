import sys, codecs, optparse, os
import heapq, math, operator

# optparser : parse command-line
optparser = optparse.OptionParser()
optparser.add_option("-c", "--unigramcounts", dest='counts1w', default=os.path.join('data', 'count_1w.txt'), help="unigram counts")
optparser.add_option("-b", "--bigramcounts", dest='counts2w', default=os.path.join('data', 'count_2w.txt'), help="bigram counts")
optparser.add_option("-i", "--inputfile", dest="input", default=os.path.join('data', 'input'), help="input file to segment")
(opts, _) = optparser.parse_args()

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
        # what does "or" do???==>self.N=N, if N=None then N=sum(freq)
        self.N = float(N or sum(self.itervalues()))
        # for a single word, posibility = 1 / N
        self.missingfn = missingfn or (lambda k, N: 1./N)

    def __call__(self, key):
        if key in self: return float(self[key])/float(self.N)
        #else: return self.missingfn(key, self.N)
        elif len(key) == 1: return self.missingfn(key, self.N)
        else: return None

# define entry
class Entry(tuple):
    def __new__(self, minus, count, word, startposition, logprobability, backpointer):
        Entry.c = property(operator.itemgetter(1))
        Entry.w = property(operator.itemgetter(2))
        Entry.sp = property(operator.itemgetter(3))
        Entry.lp = property(operator.itemgetter(4))
        Entry.bp = property(operator.itemgetter(5))
        return tuple.__new__(Entry, (-logprobability, count, word, startposition, logprobability, backpointer))

def printsegment(index, processedline):
    if(index!=None):
        # find the highest index
        printsegment(chart[index].bp, processedline)
        processedline.append(chart[index].w)

def sameEntry(entry1, entry2):
    if(entry1.w == entry2.w and entry1.sp == entry2.sp):
        return True
    else:
        return False

def Weight(word):
    return len(word)

def createEntry(method, count, word, startposition, logprobability, backpointer):
    if(method == 1):
        return Entry(-logprobability, count, word, startposition, logprobability, backpointer)
    if(method == 2):
        return Entry(-logprobability, count, word, startposition, logprobability, backpointer)

old = sys.stdout
sys.stdout = codecs.lookup('utf-8')[-1](sys.stdout)

# the default segmenter does not use any probabilities, but you could ...
Pw  = Pdist(opts.counts1w)
Pw_bi = Pdist(opts.counts2w)

method = 2

#!



# start my own codes
# create an empty heap first
with open(opts.input) as f:
    num = 0
    for line in f:
        # if(num >= 2):
        #     break
        # print line
        # initialize the heap
        h = []
        utf8line = unicode(line.strip(), 'utf-8')
        # for each word that matches input at position 0
        find = False
        for word,freq in Pw.iteritems():
            if(utf8line.find(word) == 0):
                
                entry = createEntry(method, freq, word, 0, math.log10(Pw[word]/Pw.N), None)
                #entry = createEntry(method, freq, word, 0, math.log10(Pw[word]/Pw.N), None)
                heapq.heappush(h, entry)
                #print "push:" , entry.w, entry.lp
                find = True
        if not find:
            entry = createEntry(method, 0.5, utf8line[0], 0, math.log10(0.5/Pw.N), None)
            heapq.heappush(h, entry)
            #print "push:" , entry.w, entry.lp
        # iteratively fill in chart[i] for all i
        finalindex = len(utf8line) - 1
        endindex = -1
        chart = [None] * 2 * len(utf8line)
        while(len(h)!=0):
            # entry = top entry in the heap
            #print "heap:"
            #for item in h:
                #print item.w, item.lp
            entry = heapq.heappop(h)
            #print "pop:" , entry.w, entry.lp
            # get currtindex
            currtindex = entry.sp+len(entry.w)-1
            #print "index: " , entry.sp
            if currtindex > finalindex:
                break
            # if chart[currtindex] has a previous entry
            #print "current index: " , currtindex , endindex
            #for i in range(len(chart)-1):
                #if chart[i] != None:
                    #print i , chart[i]
            if chart[currtindex] != None:
                preventry = chart[currtindex]
                if(entry.lp > preventry.lp):
                    chart[currtindex] = entry
            else:
                chart[currtindex] = entry
            #if currtindex>endindex: 
            endindex = currtindex
            # for i in range(len(chart)):
            #     if (chart[i]!=None):
            #         print "chart[", i,"]:",chart[i]
            find = False
            for newword, freq in Pw.iteritems():
                if(utf8line.find(newword, endindex) == endindex + 1):
                    
                    newentry = createEntry(method, freq, newword, endindex + 1, entry.lp + math.log10(Pw[newword]/Pw.N), endindex)
                    
                    #print "newword:" , newword , math.log10(Pw[newword]/Pw.N)
                    checkexist = False
                    for ele in h:
                        if sameEntry(ele, newentry):
                            checkexist = True
                            break
                    if not checkexist:
                        heapq.heappush(h, newentry)
                        #print "push:" , newentry.w, newentry.lp
                    find = True
            if not find:
                if (endindex + 1) <= finalindex:
                    newentry = createEntry(method, 0.5, utf8line[endindex + 1], endindex + 1, entry.lp + math.log10(0.5/Pw.N), endindex)
                    checkexist = False
                    for ele in h:
                        if sameEntry(ele, newentry):
                            checkexist = True
                            break
                    if not checkexist:
                        heapq.heappush(h, newentry)
                        #print "push:" , newentry.w, newentry.lp
        #print "final chart:"
        #for i in range(len(chart)-1):
            #if chart[i] != None:
                #print i , chart[i] 
        processedline = []
        printsegment(finalindex, processedline)
        print " ".join(processedline)
        # num += 1              
    

# old = sys.stdout
# sys.stdout = codecs.lookup('utf-8')[-1](sys.stdout)
sys.stdout = old
