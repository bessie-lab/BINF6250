# Introduction
In this project, we implemented a Gibbs sampling algorithm to discover DNA motifs from a collection of nucleotide sequences. Gibbs sampling is a probabilistic method that iteratively updates candidate motif positions based on how well they match a Position Weight Matrix (PWM) constructed from the other sequences.The algorithm begins by randomly selecting a k-mer from each sequence. During each iteration, one sequence is temporarily excluded, and the remaining candidate motifs are used to construct a Position Frequency Matrix (PFM) and Position Weight Matrix (PWM). The PWM is then used to score every possible k-mer in the excluded sequence. These scores are converted into probabilities, allowing the algorithm to probabilistically select a new motif position. Repeating this process allows the motif predictions to gradually converge toward a shared sequence pattern.This project demonstrates the use of Python, NumPy, probability-based sampling, PFMs, PWMs, and sequence analysis to solve a biological pattern-discovery problem. The resulting motif can also be visualized using a sequence logo to examine the nucleotide conservation at each position.

# Pseudocode
def GibbsMotifFinder (seqs, k, seed=None):
    '''
    Function to find a pfm from a list of strings using a Gibbs sampler
    
    Args: 
        seqs (str list): a list of sequences, not necessarily in same lengths
        k (int): the length of motif to find
        seed (int, default=None): seed for np.random

    Returns:
        pfm (numpy array): dimensions are 4xlength
    '''
    # Use rng to make random samples/selections/numbers
    # Example: randint = rng.integer(1, 10)
    random.seed(seed)
    rng = np.random.default_rng(seed)

    pass

```
GibbsMotifFinder(seqs, k, seed)

START with the input sequences (seqs) and motif length k

SET the random seed

CREATE an empty list called motifs
    # stores one motif for each sequence


INITIALIZE motifs:

    LOOP through each sequence in seqs:

        RANDOMLY select a valid k-mer position

        RANDOMLY select a strand (+ or -)

        IF the reverse strand is selected:
            GET the reverse complement of the k-mer

        ADD the selected k-mer to motifs


START convergence loop:

    REPEAT until motifs stop changing OR 10,000 iterations are reached:


        RANDOMLY select one sequence index i

        REMOVE motifs[i] temporarily from motifs


        BUILD a PFM using all motifs except motifs[i]

            USE build_pfm() from motif_ops.py


        BUILD a PWM from the PFM

            USE build_pwm() from motif_ops.py


        CREATE an empty list called candidates

            # stores possible k-mers, scores, and strand information


        LOOP through every possible k-mer position in seqs[i]:


            GET the forward k-mer

            GET the reverse complement of the k-mer

                USE reverse_complement() from seq_ops.py


            SCORE the forward k-mer using the PWM

                USE score_kmer() from motif_ops.py


            SCORE the reverse complement using the PWM

                USE score_kmer() from motif_ops.py


            STORE the k-mer, score, position, and strand
            in candidates


        CONVERT the candidate scores into probabilities


        RANDOMLY SAMPLE one candidate using the probabilities

            # higher scoring candidates have a higher chance
            # but do not select the maximum score directly


        GET the sampled k-mer and strand information


        UPDATE motifs[i] with the sampled motif


        CHECK if motifs have converged


BUILD the final PFM using all motifs

    USE build_pfm() from motif_ops.py


RETURN the final PFM

```

# Successes
Description of the team's learning points

# Struggles
Description of the stumbling blocks the team experienced

# Personal Reflections
## Group Leader
Group leader's reflection on the project

## Other member
# Aamna
This project was definitely harder for me because there were more pieces that had to work together. I found it harder to tell where a problem was comign from, so breaking it into smaller pieces was especially important this time. Marcus reiterated to break big problems into smaller ones, so we tried to approached this project differently from the start. Instead of trying to get the whole sampler finished and then debugging it, we broke it into smaller pieces and wanted to get an understanding of each part before moving on. That was especially helpful with the reverse complement. It had us stuck for a while and we knew it was something we would eventually need to include, but we decided to leave it out for now so we could actually run the main part of the sampler and see what was working before adding this piece.

The scoring was another part that I liked working through. A group member first used a simple approach of adding 10 to the scores to make them positive for the weighted selection, which gave us something we could actually run and test. I thought that was a clever way to get the algorithm moving while we were still figuring out the scoring. From there, we changed it to keep the log2 scores and use np.exp2() to convert them into weights. Being able to run the updated version and see the motifs mostly have the shine-dalgarno motif was a good check that we were moving in the right direction.

I also feel a lot more comfortable working in notebooks and with gitHub now. I can move around the notebook pretty quickly, add or remove cells, and test one small change without feeling like I am going to break everything. I actually really like how interactive that makes the debugging process. The biggest challenge this time was probably the timing. Even though we had two weeks for the project, between the different time zones and finding a meeting time that worked for everyone, it still did not feel like a lot of time. I wish I could have had a few hours each day to work on it because I think having more time to test different things and talk through the algorithm would have helped a lot. Overall, I really enjoyed working on this and I feel like it tested me in so many different ways. I hope I can come back to this at some point and finish it off.

# Generative AI Appendix
As per the syllabus
