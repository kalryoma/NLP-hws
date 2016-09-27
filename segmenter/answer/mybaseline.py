# -*- coding:UTF-8 -*-
import sys, codecs, optparse, os
import heapq, math, operator
import numpy as np

# optparser : parse command-line
optparser = optparse.OptionParser()
optparser.add_option("-c", "--unigramcounts", dest='counts1w', default=os.path.join('data', 'count_1w.txt'), help="unigram counts")
optparser.add_option("-b", "--bigramcounts", dest='counts2w', default=os.path.join('data', 'count_2w.txt'), help="bigram counts")
optparser.add_option("-i", "--inputfile", dest="input", default=os.path.join('data', 'input'), help="input file to segment")
optparser.add_option("-l", "--lambda", dest='ld', default=0.5, help="lambda")
(opts, _) = optparser.parse_args()


# Pdist is a python dictionary extract from counts1w or counts2w
class Pdist(dict):
    "A probability distribution estimated from counts in datafile."
    def __init__(self, filename, sep='\t', N=None, missingfn=None):
        self.maxlen = 0 
        for line in file(filename):
            (key, freq) = line.split(sep)
            try:
                utf8key = unicode(key, 'utf-8')
            except:
                raise ValueError("Unexpected error %s" % (sys.exc_info()[0]))
            self[utf8key] = self.get(utf8key, 0) + int(freq)
            self.maxlen = max(len(utf8key), self.maxlen)
        self.N = float(N or sum(self.itervalues()))
        self.missingfn = missingfn or (lambda k, N: 1./N)
        self.S0 = 0.5
        self.zero = self.S0 / float(self.N)

    def __call__(self, key):
        if key in self: return float(self[key])/float(self.N)
        #else: return self.missingfn(key, self.N)
        elif len(key) == 1: return self.missingfn(key, self.N)
        #else: return None
        else: 
            return self.zero
    def find(self, line, startindex):
        posibility = 0
        matchedword = []
        for i in range(self.maxlen + 1):
            word = line[startindex:startindex+i]
            if word in self:
                matchedword.append(word)
        return matchedword
    # aoSmooth, add one Smooth
    # def aoSmooth(self):
    # gtSmooth, Good Turing Smooth
    def gtSmooth(self):
        maxFreq = 0
        for word, freq in self.iteritems():
            if freq > maxFreq:
                maxFreq = freq
        Nr = [0] * (maxFreq + 2)
        for word, freq in self.iteritems():
            Nr[freq] = Nr[freq] + 1
        xp = []
        fp = []
        for i in range(len(Nr)):
            if Nr[i] != 0:
                xp.append(i)
                fp.append(Nr[i])
        for i in range(len(Nr)):
            Nr[i] = np.interp(i, xp, fp)
        Sr = [0] * (maxFreq + 1)
        for r in range(maxFreq + 1):
            if Nr[r]!=0:
                Sr[r] = 1.0 * (r + 1) * Nr[r + 1] / Nr[r]
        for word in self.keys():
            self[word] = Sr[self[word]]
        self.S0 = Sr[0]
        self.zero = self.S0/float(self.N)      
# define entry
class Entry(tuple):
    def __new__(self, word, startposition, logprobability, backpointer):
        Entry.w = property(operator.itemgetter(3))
        Entry.sp = property(operator.itemgetter(4))
        Entry.lp = property(operator.itemgetter(2))
        Entry.bp = property(operator.itemgetter(5))
        return tuple.__new__(Entry, (startposition, -logprobability, logprobability, word, startposition, backpointer))

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

def notPunctuation(word):
    punctuation = [u'０', u'１', u'２', u'３', u'４', u'５', u'６', u'７', u'８',     u'９', u'·', u'，', u'”', u'。', u'，', u'）', u'（', u'、']
    if word in punctuation:
        return False
    else:
        return True

def printsegment(index, processedline):
    if(index!=None):
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
        return Entry(word, startposition, logprobability / Weight(word), backpointer)

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

def mergeName(paragraph):
    possiblename = {}
    for line in paragraph:
        for i, word in enumerate(line):
            if len(word) == 1:
                possibleword = word
                j = i + 1
                while j < len(line) and len(line[j]) ==1 and notPunctuation(line[j]):
                    possibleword = possibleword + line[j]
                    if possibleword in possiblename:
                        possiblename[possibleword] +=1
                    else:
                        possiblename[possibleword] = 0
                    j +=1
    for key in possiblename.keys():
        if possiblename[key] <= 1 or len(key) <= 2:
            del possiblename[key]
    for key in possiblename.keys():
        flag = True
        for key2 in possiblename.keys():
            if key2!= key and key2.find(key)!=-1:
                del possiblename[key]
                break

    newparagraph = []
    newline = []
    for line in paragraph:
        i = 0
        while i < len(line):
            word = line[i]
            if len(word) == 1:
                possibleword = word
                j = i + 1
                flag = False
                while j<len(line) and len(line[j]) == 1:
                    possibleword = possibleword + line[j]
                    if possibleword in possiblename:
                        flag = True
                        break
                    j = j + 1
                if flag:
                    i = i + len(possibleword)
                    newline.append(possibleword)
                else:
                    i = i + 1
                    newline.append(word)
            else:
                i = i + 1
                newline.append(word)
        newparagraph.append(newline)
        newline = []
    return newparagraph
    
old = sys.stdout
sys.stdout = codecs.lookup('utf-8')[-1](sys.stdout)
DIGIT = _DIGIT()

# the default segmenter does not use any probabilities, but you could ...
Pw  = Pdist(opts.counts1w)
Pb = Pdist(opts.counts2w)

Pw.gtSmooth()
Pb.gtSmooth()

method = 1

ld = float(opts.ld)

# start my own codes
# create an empty heap first
with open(opts.input) as f:
    num = 0
    paragraph = []
    thisparagraph = []
    for line in f:
        #if(num >= 22):
        #    break
        # initialize the heap
        h = []
        utf8line = unicode(line.strip(), 'utf-8')
        # for each word that matches input at position 0
        matchedword = Pw.find(utf8line, 0)

        for word in matchedword:
            biword = "<s> "+word
            entry = createEntry(method, word, 0, math.log10(ld * Pb(biword) + (1-ld) * Pw(word)), None)
            heapq.heappush(h, entry)
        
        if len(matchedword) == 0:
            entry = createEntry(method, utf8line[0], 0, math.log10(ld * Pb.zero + (1 - ld) * Pw.zero), None)
            heapq.heappush(h, entry)
        
        # iteratively fill in chart[i] for all i
        finalindex = len(utf8line) - 1
        chart = [None] * len(utf8line)
        endindex = -1
        while(len(h)!=0):
            # entry = top entry in the heap
            entry = heapq.heappop(h)
            # get endindex
            currentindex = entry.sp + len(entry.w) - 1
            if currentindex > finalindex:
                break
            # if chart[endindex] has a previous entry
            if chart[currentindex] != None:
                preventry = chart[currentindex]
                if(entry.lp > preventry.lp):
                    chart[currentindex] = entry
            else:
                chart[currentindex] = entry
            endindex = currentindex
            newmatchedword = Pw.find(utf8line, endindex + 1)
            #print newword
            for newword in newmatchedword:
                biword = entry.w + " " + newword
                newentry = createEntry(method, newword, endindex + 1, entry.lp + math.log10(ld * Pb(biword) / Pw(entry.w) + (1 - ld) * Pw(newword)), endindex)
                checkexist = False
                for ele in h:
                    if sameEntry(ele, newentry):
                        checkexist = True
                        break
                if not checkexist:
                    heapq.heappush(h, newentry)
            if len(newmatchedword) == 0:
                if (endindex + 1) <= finalindex:
                    newentry = createEntry(method, utf8line[endindex + 1], endindex + 1, entry.lp + math.log10(ld * Pb.zero/Pw(entry.w) + (1 - ld) * Pw.zero), endindex)
                    checkexist = False
                    for ele in h:
                        if sameEntry(ele, newentry):
                            checkexist = True
                            break
                    if not checkexist:
                        heapq.heappush(h, newentry)
                
        processedline = []
        printsegment(finalindex, processedline)
        mergedline = mergeDigit(processedline, DIGIT)
        thisparagraph.append(mergedline)
        if(mergedline[1] == u"完"):
            paragraph.append(thisparagraph)
            thisparagraph = []
        #num += 1              
    for thisparagraph in paragraph:
        mergedparagraph = mergeName(thisparagraph)
        for mergedline in mergedparagraph:
            print " ".join(mergedline)
sys.stdout = old
