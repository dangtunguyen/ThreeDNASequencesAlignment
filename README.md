# Project Description

We will discuss a space-saving technique for pairwise sequence alignment
based on divide-and-conquer in the lectures. The technique can be easily
extended to multiple sequence alignment. Your task is to design an
O(m*n) space, O(m*n*k) time SP alignment algorithm for three DNA
sequences of lengths m, n, and k, respectively. Implement your algorithm
in C/C++/Java, and test it on real data.

You may find discussions on the space-saving technique for pairwise 
sequence comparison in the books by Jones and Pevzner, by Gusfield, 
and by Jiang, Xu, and Zhang. The following original papers provide more
detailed information:

D.S. Hirschberg. Algorithms for the longest common subsequence
problem. J.ACM, 24:664-75, 1977.

E. W. Myers and W. Miller. Optimal alignments in linear space.
Comp. Appl. Biosciences, 4:11-17, 1988.
At UCR Science Lib: call number is QH324.2 .C66 

Chao, Hardison and Miller. Recent developments in linear-space
alignment methods: a survey. Journal of Computational Biology.
Vol. 1-4, pp. 271-291. 1994. 

(i)  For simplicity, let's consider global alignment and the basic 
     SP score model where gaps are not specially treated. 

(ii) Your program should work for any score function on nucleotides.
     In other words, the user should be able to input a score function
in the form of a 5x5 matrix indexed by A, C, G, T, and space.
     The SP-score of a column of letters/spaces is the sum of the scores
     of each pair of letters/spaces in the column.

To test your program, use the Blast scores: Match = 5, Mismatch = -4,
and Indel = -8. The score between two spaces is 0.

(iii) Test your program on the following four sets of sequences:

     1.
     NM_013096.  Rattus norvegicus hemoglobin alpha, adult chain 2 (Hba-a2),
                 mRNA. 556 bps.
     NM_008218.  Mus musculus hemoglobin alpha, adult chain 1 (Hba-a1), 
                 mRNA. 569 bps.
     NM_000558.  Homo sapiens hemoglobin, alpha 1 (HBA1), mRNA. 627 bps.
     
     2.
     NM_010019.  Mus musculus death-associated protein kinase 2 (Dapk2),
                mRNA, 1792 bps.
     NM_001243563. Sus scrofa death-associated protein kinase 2 (DAPK2),
                mRNA, 1825 bps.
     NM_014326. Homo sapiens death-associated protein kinase 2 (DAPK2), 
                mRNA, 2628 bps.
                
     3.
     NM_000545. Homo sapiens HNF1 homeobox A (HNF1A), mRNA. 3417 bps
     NM_008261. Mus musculus hepatic nuclear factor 4 (Hnf4). 4391 bps
     NM_000457. Homo sapiens hepatocyte nuclear factor 4, alpha (HNF4A), 
                transcript variant 2, mRNA. 4737 bps
                
     4.
     NM_000492. Homo sapiens cystic fibrosis transmembrane conductance
                regulator (ATP-binding cassette sub-family C, member 7) 
                (CFTR), mRNA. 6132 bps
     NM_031506. Rattus norvegicus cystic fibrosis transmembrane 
                conductance regulator homolog (Cftr), mRNA. 6287 bps.
     NM_021050. Mus musculus cystic fibrosis transmembrane conductance
                regulator homolog (Cftr), mRNA. 6305 bps.

You may retrieve the sequences at http://www.ncbi.nlm.nih.gov/nucleotide/
using the accession numbers NM_000558.3, etc. Select "nucleotide" in 
the search option box. The sequence data is given at the bottom of 
the search result page. Note that the last dataset may take your algorithm 
quite some time to run, especially on an old computer.

For this question, submit a report (hard copy) along with HW3 with
                                               --------------
(a) a high-level description of your algorithm (i.e. high-level
    pseudo-code), and the data structures used,
(b) your source code, and
(c) the result of your program on each dataset, including the
    score of the optimal alignment obtained, the length of the alignment,
the number of columns with perfectly matched nucleotides in the
alignment, and the running time and memory (on a PC or laptop
    with specification of CPU and memory). Please do not include
the actual alignment in the report (because it is going to be quite
loooong :-). 
(d) Do you see any obvious conserved regions captured in your alignments?

Although this question is optional for HW2, it could be a very good idea
that you complete the dynamic programming routine for computing the 
optimal score of a 3-sequence alignment in O(m\*n) space, and test it on
some small data to make sure that it is correct. The recurrence relation is
given in Chapter 6.10 of the textbook. Note that this algorithm will need
the dynamic programming algorithm for pairwise sequence alignment to 
initialize its matrix.

Also, although your goal is to implement the O(mn)-space algorithm, it
will be useful for you to implement the standard O(mnk)-space dynamic 
programming algorithm for 3-sequence alignment and use it as a subroutine 
in the divide-and-conquer process when one of the sequences degenerates 
to one or zero letters.
