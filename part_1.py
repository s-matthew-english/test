'''
Created on Sep 23, 2014

@author: matthias
'''
#from replace_rare import Replacer
from business_logic import Replacer

from count_freqs import Hmm, simple_conll_corpus_iterator

from business_logic import simple_tagger, x_corpus_iterator

from eval_gene_tagger import Evaluator, corpus_iterator


def usage():
    """
    python replace_rare.py [input_file] > [output_file]
        Read in a gene tagged training input file and 
        output its content with rare words replaced with _RARE_.
    """

if __name__ == "__main__":


    usage()
    

    print("\n1. Generate training corpus with rare words replaced by \'_RARE_\'.\n")
    
    input1 = open('gene.train',"r")
    #input1 = open('test.train',"r")
    input2 = open('gene.train',"r")
    #input2 = open('test.train',"r")

    output1 = open('gene.replace.train',"w")
    #output1 = open('test.replace.train',"w")
    output2 = open('gene.replace.train',"w")
    #output2 = open('test.replace.train',"w")

    # Initialize a counter
    replacer = Replacer()
    # Collect counts
    replacer.word_count(input1, output1)
    # Replace rare words
    replacer.replace_rare(input2, output2)
    
    output1.flush()
    output2.flush()
    
    
    
##############################################################    
    #commence p1a
    
    print("\n2. Generate word count file.\n")
    
    p1a_input = open('gene.replace.train',"r")
    #p1a_input = open('test.replace.train',"r")
    p1a_output = open('gene.counts', "w")
    #p1a_output = open('test.counts', "w")
    
    # Initialize a trigram counter
    counter = Hmm(3)
    # Collect counts
    counter.train(p1a_input)
    # Write the counts
    counter.write_counts(p1a_output)
    
    p1a_output.flush()
    
    
    
    
#############################################################    
    #commence p1b
    
    print("\n3. Tag dev corpus with simple tagger.\n")
    
    p1b_counts_file = open('gene.counts',"r")
    #p1b_counts_file = open('test.counts',"r")
    
    #p1b_gene_file   = open('gene.dev',"r")
    p1b_gene_file   = open('tea.dev',"r")
    
    #p1b_output_file_1 = open('bean_dev.p1.out',"w")
    p1b_output_file_1 = open('gene_dev.p1.out',"w")
    #p1b_output_file_2 = open('bean_dev.p1.out',"w")
    p1b_output_file_2 = open('gene_dev.p1.out',"w")


    # Initialize a simple gene tagger
    smp_tagger = simple_tagger()
    # Read counts
    smp_tagger.read_counts(x_corpus_iterator(p1b_counts_file), p1b_output_file_1)
    # Tag 
    smp_tagger.write_tags(x_corpus_iterator(p1b_gene_file), p1b_output_file_2)

    
    p1b_output_file_1.flush()
    p1b_output_file_2.flush()
    
    

#############################################################    
    #commence p1c

    print("\n4. Measure tagger score along with \'golden\' sample.\n")

    gs_iterator = corpus_iterator(open('tea.key',"r"))
    #gs_iterator = corpus_iterator(open('gene.key',"r"))
    #pred_iterator = corpus_iterator(open('bean_dev.p1.out',"r"), with_logprob = False)
    pred_iterator = corpus_iterator(open('gene_dev.p1.out',"r"), with_logprob = False)
    evaluator = Evaluator()
    evaluator.compare(gs_iterator, pred_iterator)
    evaluator.print_scores()
    print('\n')
