'''
Created on Sep 23, 2014

@author: matthias
'''
import sys
#from replace_rare_use_classes import Replacer
from business_logic import Replacer

from count_freqs import Hmm

#from viterbi_tagger import Tagger, sentence_iterator, v_corpus_iterator
from business_logic import viterbi_tagger, sentence_iterator, x_corpus_iterator

from eval_gene_tagger import Evaluator, corpus_iterator



def usage():
    """
    python replace_rare.py [input_file] > [output_file]
        Read in a gene tagged training input file and 
        output its content with rare words replaced with _RARE_.
    """

if __name__ == "__main__":

    
    usage()

    input1 = open('gene.train',"r")
    input2 = open('gene.train',"r")
    
    output1 = open('gene.replace.filter.train',"w")
    output2 = open('gene.replace.filter.train',"w")

    print("\n1. Generate training corpus with rare words replaced by (things other than) \'_RARE_\'."'\n')

    # Initialize a counter
    replacer = Replacer()
    # Use rare word classes
    replacer.set_rare_word_filter([['_NUMERIC_' , "[0-9]+"], ['_ALL_CAPITALS_' , "^[A-Z]+$"], ['_LAST_CAPITAL_' , "[A-Z]+$"]])
    # Collect counts
    replacer.word_count(input1, output1)
    # Replace rare words
    replacer.replace_rare(input2, output2)
    
    output1.flush()
    output2.flush()
    
    
###################################################

###################################################


    print("\n2. Generate word count file.\n")
    
    freqs_input = open('gene.replace.filter.train',"r")
    freqs_output = open('gene.counts', "w")

    
    # Initialize a trigram counter
    counter = Hmm(3)
    # Collect counts
    counter.train(freqs_input)
    # Write the counts
    counter.write_counts(freqs_output)
    
    freqs_output.flush()
    
    
    
###################################################

###################################################   

    
    print("\n3. Tag dev corpus with Viterbi tagger.\n")
    
    v_counts_file = open('gene.counts',"r")
    v_gene_file   = open('tea.dev',"r")
    
    v_output_file_1 = open('gene_dev.p3.out',"w")
    v_output_file_2 = open('gene_dev.p3.out',"w")


    # Initialize a simple gene tagger
    v_tagger = viterbi_tagger()
    # Use rare word classes
    v_tagger.set_rare_word_filter([['_NUMERIC_', "[0-9]+"], ['_ALL_CAPITALS_', "^[A-Z]+$"], ['_LAST_CAPITAL_', "[A-Z]+$"]])
    # Read counts, from the basic tagger
    v_tagger.read_counts(x_corpus_iterator(v_counts_file), v_output_file_1)
    # Tag , unique to the viterbi tagger, iterate over the sentences, iterate over the words within
    v_tagger.write_tags(sentence_iterator(x_corpus_iterator(v_gene_file)), v_output_file_2) 
    
    
    
    v_output_file_1.flush()
    v_output_file_2.flush()
    

#############################################################    
    

#############################################################    

    print("\n4. Evaluate results against \'golden\' standard.")

    gs_iterator = corpus_iterator(open('tea.key',"r"))
    pred_iterator = corpus_iterator(open('gene_dev.p3.out',"r"), with_logprob = False)
    evaluator = Evaluator()
    print('\n')
    evaluator.compare(gs_iterator, pred_iterator)
    evaluator.print_scores()
    print('\n')