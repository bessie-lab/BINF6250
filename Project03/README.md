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

steps to be taken in writing the fuction
1.Randomly choose a motif from each sequence.
#make a list for all the motifs
motif = []
seq_length = length of the current sequence
k  = length of the motif
#now going to choose a starting position
for seq in seqs:
    seq_length = len(seq)
start_position = rng.integers(0, seq_length - k + 1)
# Extract the k-length motif
  motif = seq[start:start + k]

    # Store the motif
    motifs.append(motif)

return motifs# (was for practice only)

for iteration in range(1000):


2.Temporarily remove one motif.
3.Use the remaining motifs to construct a PFM.
4.Convert/use that PFM to obtain a PWM.
5. Use the PWM to score possible 10-mers in the removed sequence.
6.Use those scores to probabilistically choose a new motif.
7.Repeat.
8.At the end, create your final PFM.

```

# Successes
Description of the team's learning points

# Struggles
Description of the stumbling blocks the team experienced

# Personal Reflections
## Group Leader
Group leader's reflection on the project

## Other member
Other members' reflections on the project

# Generative AI Appendix
As per the syllabus
