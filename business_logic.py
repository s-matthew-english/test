'''
Created on Sep 24, 2014

@author: matthias
'''
import sys
import re
from collections import defaultdict


class basic_tagger(object):    
    def __init__(self):

        #this holds the emission count for each word/tag pair
        self.emission_counts = defaultdict(int)
        
        #n grams are like 'maximum likelihood estimate' 
        self.n = 3
        #for unigram, bigram and trigram
        self.ngram_counts = [defaultdict(int) for i in range(self.n)]
        
        #this just holds all the tags, which is just 'O' and 'I-GENE'
        self.all_states = set()
        
        #holds the number of times the word has been observed
        #in the training data set
        self.word_counts = defaultdict(int)
        
        #number of times a tag is observed
        #used in determining the maximum likelihood
        #estimate
        self.ne_tag_counts = defaultdict(int)
        
        #this is useful for the filter in part 3
        self.rare_words_filter = {} 
                
    def read_counts(self, iterator, output):

        #responsible for populating the 
        #data structures
        for line in iterator:
            parts = line.strip().split(" ")
            #the number of times the word has been seen
            count = float(parts[0])
            #if it's a word under consideration
            if parts[1] == "WORDTAG":
                
                #data from the counts training file
                """
                19 WORDTAG I-GENE CBP
                6 WORDTAG O disappeared
                5 WORDTAG O modality
                """
                ne_tag = parts[2]   #this is the tag associated with the word
                                    #and which has generated that count

                word = parts[3] #grab the word
                
                #populate the data structures
                self.emission_counts[(word, ne_tag)] = count
                self.all_states.add(ne_tag)
                
                #you count all the times you've seen
                #a particular tag
                self.ne_tag_counts[ne_tag] += count
                
                #all the times
                #you've seen a particular word
                self.word_counts[word] += count
            
            elif parts[1].endswith("GRAM"):
                """  
                345128 1-GRAM O
                41072 1-GRAM I-GENE
                13 2-GRAM I-GENE STOP
                """
                #i.e. 1-gram, 2-gram, 3-gram
                n = int(parts[1].replace("-GRAM",""))
                
                #if it's a bi-gram, then
                #grab out the parts of it
                ngram = tuple(parts[2:])
                
                #how many times you've seen that n-gram
                self.ngram_counts[n-1][ngram] = count 
        
class simple_tagger(basic_tagger):
    
    #this part here is the crux of the simple tagger in fact
    #counting up the MAXIMUM LIKELIHOOD ESTIMATE for each word 
    def count_tag(self, word):
        count = 0
        c = 0
        
        
        #sets dont support specific grab operations
        #self.all_states is a set (maybe fix this)
        tag_list = list(self.all_states)
        
        
        #if it's rare for all the tags, return the rare value (happens to be I-GENE)
        if self.emission_counts[(word, tag_list[0])] < 5: 
            if self.emission_counts[(word, tag_list[1])] < 5:
                return self.count_tag("_RARE_")
        
        # for 'O' and 'I-GENE'
        for tag in self.all_states:
            
            #if it's appeared at all in the training data
            if self.emission_counts[(word, tag)] > 0:
                
                c = self.emission_counts[(word, tag)]/self.ne_tag_counts[tag]
                
            #if count < c:
            if c > count:

                count = c
                #assign the tag, the tag from which we've been iterating
                #the argmax starts at zero, the tag giver lives in this 
                #little castle, and at fist, just have to be bigger than zero, but next you have to be bigger 
                #than the biggest guy who came before
                ne_tag = tag

        #return the argmax tag
        return ne_tag

    def write_tags(self, iterator, output):

        for word in iterator:

            if word:
                #create the output file by
                #writing that word (that you just read in from gene.dev, and the 
                #tag by feeding it into the count_tag
                #method that you were just examining up top
                
                #because count_tag is returning the tag it gets from the count file
                output.write("%s %s\n" %(word, self.count_tag(word)))
            else:
                #and if it's not a word, just write a blank line
                #new sentence
                output.write("\n");
                          
                
class viterbi_tagger(basic_tagger):
    def set_rare_word_filter(self, filter):
        #populate the rare_words_filter data structure in the basic
        #tagger with a filter we feed in here.
        self.rare_words_filter = filter

    def tag_sentence(self, sentence, output):
        #execute the viterbi calculation
        viterbi = Viterbi(sentence, self)
        return viterbi.make_tag_sequence()

    def write_tags(self, iterator, output):
        k = 0
        for sentence in iterator:
            k += 1
            tag_sequence = self.tag_sentence(sentence, output)
            for i in range(len(sentence)):
                output.write('%s %s\n' % (sentence[i], tag_sequence[i]))
            output.write('\n');


             

class Viterbi(object):
    """
        Counts and stores coefficients along iteration path
    """
    def __init__(self, sentence, tagger):
        self.tagger = tagger
        self.sentence = sentence
        self.pi = defaultdict(int)
        self.bp = defaultdict(int)


    def get_factor(self, q_args, e_args):
        #refer to trigram function
        q = self.get_q(q_args)
        #refer to emission function
        e = self.get_e(e_args)
        
        
        return q*e

 
    def get_q(self, args):
        [w, u, v] = list(args)
        #the trigram probability
        #print("(",w,",", u,",", v,") / (",w,", ",u,") = ",self.tagger.ngram_counts[self.tagger.n-1][(w, u, v)]/self.tagger.ngram_counts[self.tagger.n-2][(w, u)],'\n') 
        #print("\n")

        return self.tagger.ngram_counts[self.tagger.n-1][(w, u, v)]/self.tagger.ngram_counts[self.tagger.n-2][(w, u)]


    def get_e(self, args):
        [word, v] = list(args)
        #the emission probability
        return self.tagger.emission_counts[(word, v)]/self.tagger.ne_tag_counts[v]


    def filter_rare_word(self, word):
        if self.tagger.word_counts[word] >= 5:
            return word
        # Use rare word classes
        #the tagger is the viterbi tagger if the regular expression
        #corresponds to a certain word, return the replacement tag
        #that is identified with the characteristics specified
        else:
            for [mark, regex] in self.tagger.rare_words_filter:
                if re.search(regex, word):
                    
                    return mark
            return '_RARE_'


    def count_step_coeff(self, step):
        word = self.filter_rare_word(self.sentence[step-1])
        if step == 1:
            print("HOT STEPPA", step)
            for v in self.tagger.all_states:
                self.pi[(step, '*', v)] = self.get_factor(('*', '*', v), (word, v))
        elif step == 2:
            #for 'O' and 'I-GENE' 
            for v in self.tagger.all_states:
                for u in self.tagger.all_states:
                    print("HOT STEPPA", step)
                    self.pi[(step, u, v)] = self.pi[(step - 1,'*', u)]*self.get_factor(('*', u, v), (word, v))
        else:
            #for 'O' and 'I-GENE' 
            for v in self.tagger.all_states:
                #for 'O' and 'I-GENE' 
                for u in self.tagger.all_states:
                    #for 'O' and 'I-GENE' 
                    for w in self.tagger.all_states:
                        #create a depth-first search and execute the following calculation
                        #for each contingent
                        #the emission parameter multiplied with the trigram probability
                        
                        #populate the default dict, filling it in at each step
                        print("step",step)
                        pi = self.pi[(step - 1, w, u)]*self.get_factor((w, u, v), (word, v))
                        print('\n',"*********************************")
                        print('\t',"self.pi[(",step - 1,",", w,",", u,")] =",self.pi[(step - 1, w, u)],'\n')
                        print('\t',"self.get_factor((",w,",", u,",", v,")",",", "(",word,",", v,")) =",'\n',self.get_factor((w, u, v), (word, v)))
                        print("*********************************",'\n')
                        
                        if pi > self.pi[(step, u, v)]:
                            #this is the most important part
                            #of the entire algorithm, where
                            #we populate the dictionary with the 
                            #backpointer argmax for ascribing tags
                            self.pi[(step, u, v)] = pi
                            self.bp[(step, u, v)] = w
                            print("self.bp[(step, u, v)]",self.bp[(step, u, v)])

    def count_coeff(self):
        n = len(self.sentence)
        #for each position in the sentence
        for i in range(n):
            self.count_step_coeff(i+1)


    def make_tag_sequence(self):
        self.count_coeff()
        n = len(self.sentence)
        max_val = 0
        #edge case
        if n == 1:
            max_val = 0
            for v in self.tagger.all_states:
                p = self.pi[(n, '*', v)]*self.get_q(('*', v, 'STOP'))
                if p > max_val:
                    max_val = p
        else:
            for u in self.tagger.all_states:
                for v in self.tagger.all_states:
                    #for all tags (twice) for each sentence
                    #the last position of the sentence times the STOP
                    p = self.pi[(n, u, v)]*self.get_q((u, v, 'STOP'))
                    if p > max_val:
                        max_val = p
                        #keep the argmax and return it
                        #we start from the back
                        tag_sequence = [u, v]
                        #edge case
                        if max_val == 0:
                            tag_sequence = [u, v]
                        if n == 2:
                            #edge case
                            return tag_sequence 
            #standard operating procedure
            for k in range(n-2, 0, -1):
                #start at n-2 because we begin by assigning the first 
                #two tags
                prev = self.bp[(k+2, tag_sequence[0], tag_sequence[1])]
                
                print ("prev",prev)
                print("k",k)
                #viola, add two to k because the first two already exist
                tag_sequence.insert(0, prev)
                #insert at zero because we began from the back
            return tag_sequence




"""
We need to predict emission probabilities for words in the test data
that do not occur in the traning data. One simple approach is to map
infrequent words in the training data to a common class and to treat
unseen words as members of this class. Replace infrequent words 
(Count(x) < 5) in the original training data file with a common 
symbol _RARE_.

HMM taggers can be improved by grouping words into informative word
classes rather than just into a single class of rare words. Here 
four rare word classes are implemented: Numeric, All Capital, Last
Capital and Rare. 
"""


class Replacer(object):
    """
    Find rare words in train corpus. 
    """

    def __init__(self):
        #create a container for all the words
        self.rare_word_list = defaultdict(int)
        self.rare_words_filter = {}

    def set_rare_word_filter(self, filter):
        self.rare_words_filter = filter

    def word_count(self, corpus_file, output):
        iterator = corpus_iterator(corpus_file)
        for word, ne_tag in iterator:
            #count the frequency with which
            #each word in the training data appears
            if word:
                #increment at each instance spotted
                self.rare_word_list[word] += 1
        #to all the words in the list
        for word in list(self.rare_word_list.keys()): 
            #if they're not rare
            if self.rare_word_list[word] >= 5:
                #remove them from the list
                del self.rare_word_list[word]

    def filter_rare_word(self, word):
        #if it's not a rare word
        if self.rare_word_list[word] == 0:
            #return it so you can
            #look for it's count
            return word
        else:
            #use rare word filter
            for [mark, regex] in self.rare_words_filter:
                if re.search(regex, word):
                    return mark
            #if no special characteristics
            #just return _RARE_
            return '_RARE_'
    
    def replace_rare(self, corpus_file, output):
        iterator = corpus_iterator(corpus_file)
        for word, ne_tag in iterator:
            if word is None:
                output.write('\n');
            else:
                output.write('%s %s\n' %(self.filter_rare_word(word), ne_tag))
   
   
   

def x_corpus_iterator(corpus_file):
    """
    Get an iterator object over the corpus file. The elements of the
    iterator contain sentence word. Blank lines, indicating
    sentence boundaries return None.
    """
    l = corpus_file.readline()
    while l:
        line = l.strip()
        if line: # Nonempty line
            yield line
        else: # Empty line
            yield None
        l = corpus_file.readline()


def corpus_iterator(corpus_file):
    """
    Get an iterator object over the corpus file. The elements of the
    iterator contain (word, ne_tag) tuple. Blank lines, indicating
    sentence boundaries return (None, None).
    """
    l = corpus_file.readline()
    while l:
        line = l.strip()
        if line: # Nonempty line
            # Extract information from line.
            # Each line has the format
            # word ne_tag
            fields = line.split(" ")
            ne_tag = fields[-1]
            word = " ".join(fields[:-1])
            yield (word, ne_tag)
        else: # Empty line
            yield (None, None)
        l = corpus_file.readline()


def sentence_iterator(corpus_iterator):
    """
    Return an iterator object that yields one sentence at a time.
    Sentences are represented as lists of words.
    """
    current_sentence = [] #Buffer for the current sentence
    for l in corpus_iterator:
            if l==None:
                if current_sentence:  #Reached the end of a sentence
                    yield current_sentence
                    current_sentence = [] #Reset buffer
                    
                else: # Got empty input stream
                    sys.stderr.write("WARNING: Got empty input file/stream.\n")
                    raise StopIteration
            else:
                current_sentence.append(l) #Add token to the buffer

    if current_sentence: # If the last line was blank, we're done
        yield current_sentence  #Otherwise when there is no more token
                                # in the stream return the last sentence.
