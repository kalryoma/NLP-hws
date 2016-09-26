# -*- coding:UTF-8 -*-
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
    def __new__(self, word, startposition, logprobability, backpointer):
        Entry.w = property(operator.itemgetter(1))
        Entry.sp = property(operator.itemgetter(2))
        Entry.lp = property(operator.itemgetter(3))
        Entry.bp = property(operator.itemgetter(4))
        return tuple.__new__(Entry, (-logprobability, word, startposition, logprobability, backpointer))

class _DIGIT:
    def __init__(self):
        self.value = [ u'０', u'１', u'２', u'３', u'４', u'５', u'６', u'７', u'８', u'９', u'·']
    def match(self, word):
        flag = True
        for i in range(len(word)):
            if(word[i] not in self.value):
                flag = False
                break
        return flag

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

def createEntry(method, word, startposition, logprobability, backpointer):
    if(method == 1):
        return Entry(word, startposition, logprobability, backpointer)
    if(method == 2):
        return Entry(word, startposition, logprobability/Weight(word), backpointer)

def mergeDigit(line, DIGIT):
    newline = []
    newword = []
    for i, word in enumerate(line):
        if(DIGIT.match(word) == False and len(newword)!= 0):
            newline.append("".join(newword))
            newline.append(word)
            newword = []
        elif(DIGIT.match(word) == False):
            newline.append(word)
        elif(DIGIT.match(word) == True):
            newword.append(word)
    return newline

old = sys.stdout
sys.stdout = codecs.lookup('utf-8')[-1](sys.stdout)

# the default segmenter does not use any probabilities, but you could ...
Pw  = Pdist(opts.counts1w)
DIGIT = _DIGIT()

method = 2

# start my own codes
# create an empty heap first
with open(opts.input) as f:
    num = 0
    for line in f:
        #if(num >= 8):
        #     break
        #if(num < 7):
        #    num = num + 1
        #    continue
        # initialize the heap
        h = []
        utf8line = unicode(line.strip(), 'utf-8')
        #print "%d:" %num + utf8line
        # for each word that matches input at position 0
        find = False
        for word,freq in Pw.iteritems():
            if(utf8line.find(word) == 0):
                entry = createEntry(method, word, 0, math.log10(Pw(word)), None)
                heapq.heappush(h, entry)
                # print "push:" , entry.w, entry.lp
                find = True
        if not find:
            entry = createEntry(method, utf8line[0], 0, 0, None)
            heapq.heappush(h, entry)
            # print "push:" , entry.w, entry.lp
        # iteratively fill in chart[i] for all i
        finalindex = len(utf8line) - 1
        endindex = -1
        chart = [None] * len(utf8line)
        while(len(h)!=0):
            # entry = top entry in the heap
            entry = heapq.heappop(h)
            # print "pop:" , entry.w, entry.lp
            # get currtindex
            currtindex = entry.sp+len(entry.w)-1
            if currtindex > finalindex:
                break
            # if chart[currtindex] has a previous entry
            if chart[currtindex] != None:
                preventry = chart[currtindex]
                if(entry.lp > preventry.lp):
                    chart[currtindex] = entry
            else:
                chart[currtindex] = entry
            if currtindex>endindex:
                endindex = currtindex
            # for i in range(len(chart)):
            #     if (chart[i]!=None):
            #         print "chart[", i,"]:",chart[i]
            find = False
            for newword, freq in Pw.iteritems():
                if(utf8line.find(newword, endindex) == endindex + 1):
                    newentry = createEntry(method, newword, endindex + 1, entry.lp + math.log10(Pw(newword)), endindex)
                    checkexist = False
                    for ele in h:
                        if sameEntry(ele, newentry):
                            checkexist = True
                            break
                    if not checkexist:
                        heapq.heappush(h, newentry)
                        # print "push:" , newentry.w, newentry.lp
                    find = True
            if not find:
                if (endindex + 1) <= finalindex:
                    newentry = createEntry(method, utf8line[endindex + 1], endindex + 1, entry.lp, endindex)
                    checkexist = False
                    for ele in h:
                        if sameEntry(ele, newentry):
                            checkexist = True
                            break
                    if not checkexist:
                        heapq.heappush(h, newentry)
                        # print "push:" , newentry.w, newentry.lp
                
        processedline = []
        printsegment(finalindex, processedline)
        #mergedline = mergeDigit(processedline, DIGIT)
        print " ".join(processedline)
        num += 1              
    

# old = sys.stdout
# sys.stdout = codecs.lookup('utf-8')[-1](sys.stdout)
sys.stdout = old
