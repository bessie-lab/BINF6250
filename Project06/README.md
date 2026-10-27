# BINF6250 - Project 06

# Introduction
This project builds a phylogenetic tree from a set of reverse transcriptase (RT) protein sequences using a distance-based approach. Pairwise similarity between sequences is measured with the Smith-Waterman local alignment algorithm, and those scores are converted into a distance matrix. The Neighbor-Joining (NJ)algorithm then uses that matrix to build an unrooted tree showing how thesequences relate to each other.Neighbor-Joining was chosen over UPGMA because it does not assume a constant rate of evolution across lineages, so branch lengths can differ between sequences. Smith-Waterman scoring finds the best matching region between twosequences, which suits sequences that vary in length or have diverged in parts.

# How it works
Read the data:sequences are loaded from a FASTA file(`lafayette_SARS_RT.fasta`).
Score pairs:every pair of sequences is aligned with Smith-Waterman.The best local score is normalized to a value   between 0 and 1, where 1 means a perfect match.
Build distances:each similarity is converted to a distance(`distance = 1 - similarity`), giving a symmetric distance matrix.
Join neighbors:Neighbor-Joining repeatedly merges the pair of nodes with the smallest Q-value and updates the distance matrix until two nodes remain.
Draw the tree:the result is written in Newick format and plotted with Biopython and Matplotlib.

# Files

- project06.ipynb - notebook with the implementation
- lafayette_SARS_RT.fasta -input sequences

#Requirements

Python 3, `numpy`, `matplotlib`, and `biopython`:

    pip install numpy matplotlib biopython

# Usage

Place the FASTA file in the same folder as the notebook, open it in VS Code or Jupyter, and choose **Run All**. The final cell reads the sequences, builds the distance matrix, runs Neighbor-Joining, and plots the tree.

# Pseudocode
def read_fasta(filename: str) -> Dict[str, str]:
    sequences = empty dictionary (name to sequence)
    current_name = nothing yet

    Open the file and read it line by line:
        strip spaces and the newline off the end of the line
        IF the line is empty:
            skip it
        IF the line starts with ">":
            current_name = last word of the line (i.e. "DM1")
            sequences[current_name] = empty string
        ELSE (it's part of a sequence):
            add the line to the end of sequences[current_name]

    Return sequences

Changes from Project 05 smith_waterman:
    1. cal_score is kept as-is; it still calculates each cell's score. The traceback arrow it records is not used here.
    2. Drop traceback, along with the traceback_matrix lines and the max_pos line, since nothing uses them anymore.
    3. Keep the original table-filling loop and the best_score line in a helper function, and remove both print statements.
    4. smith_waterman uses the Project 06 name and inputs, calls the helper, and returns the normalized score:
       normalized = raw(seq1, seq2) / smaller of [raw(seq1, seq1), raw(seq2, seq2)]

def cal_score(matrix, seq1, seq2, i, j, match, mismatch, gap):  #Original calc_score from Project05 untouched.  Needed to run Smith_Waterman function from Project05

def project5_smith_waterman_score(seq1, seq2, match, mismatch, gap):  #Original smith_waterman function from Project05, turned into a helper function for the new smith_waterman function in Project06

def smith_waterman(
    sequence1: str,
    sequence2: str,
    match_score: float = 2.0,
    mismatch_penalty: float = -1.0,
    gap_penalty: float = -1.0
) -> float:
    raw_score = project5_smith_waterman_score(sequence1 vs sequence2, using match_score, mismatch_penalty, gap_penalty)
    self_score1 = project5_smith_waterman_score(sequence1 vs sequence1, same scores)
    self_score2 = project5_smith_waterman_score(sequence2 vs sequence2, same scores)

    smaller_self_score = smaller of self_score1 and self_score2

    IF smaller_self_score is 0 (only happens if a sequence is empty):
        return 0.0 (avoids dividing by zero)
    
    Return (raw_score / smaller_self_score) Note: 1.0 = identical, closer to 0 = less similar

def build_distance_matrix(sequences: Dict[str, str]) -> Tuple[np.ndarray, List[str]]:
    ids = list of the sequence names, in file order (i.e. ["DM1", "DM2", ...])
    n = number of sequences
    distance_matrix = n x n grid of zeros (diagonal stays 0: each sequence vs itself)

    For each row i from 0 to n-1:
        FOR each column j from i+1 to n-1 (only pairs not yet done; skip the diagonal):
            similarity = smith_waterman (sequence for ids[i], sequence for ids[j])
            distance = 1 - similarity
            distance_matrix[i, j] = distance
            distance_matrix[j, i] = distance  (Table is symmetric, so fill both halves)

    Return distance_matrix, ids

def neighbor_joining(
    distance_matrix: np.ndarray,
    labels: List[str]
) -> str:
    PROPOSAL FOR DECISIONS IN PSEUDOCODE/CODE
        Final branch between the last two nodes is split evenly between them
        Branch lengths are written to 4 decimals places in the Newick text (calculations use full precision)
        All 20 sequences are used.  The docstring expects 19, note it in README
        Ties in the D* score: the first pair found is kept
        Negative branch lengths: set to 0, differences added to the other branch
        
    SETUP
    D = copy of distance_matrix (so original table isn't changed)
    current_labels = copy of labels (each starts as just a name, i.e. "DM1")

    WHILE more than 2 nodes remain:
        n = number of rows in D

        Step 1. Total Distance:  totals = sum of each for row of D

        Step 2. D* score for every pair, keep the lowest:  lowest_score = a very large number
        FOR each row i from 0 to n-1:
            FOR each column j from i+1 to n-1:
                score = (n-2) * D[i, j] - totals[i] - totals[j]
                IF score < lowest_score:
                    lowest_score = score
                    best_i = i
                    best_j = j

        Step 3. Branch lengths for best_i and best_j
        branch_i = (D[best_i, best_j] + (totals[best_i] - totals[best_j]) / (n-2)) / 2
        branch_j = D[best_i, best_j] - branch_i
        IF either branch is negative:
            set it to 0 and add the difference to the other branch
        
        Step 4. Distances from the new node to every other node k
        new_distances = empty list
        FOR each k that is not best_i or best_j:
            new_distances.append((D[best_i, k] + D[best_j, k] - D[best_i, best_j]) / 2)

        New_label for the joined pair:
        new_label = "(" + current_labels[best_i] + ":" + f"{branch_i:.4f}" + "," + current_labels[best_j] + ":" + f"{branch_j:.4f}" + ")"

        Step 5. Shrink the table
        keep = every index except best_i and best_j
        new_D = (n-1) x (n-1) grid of zeros
        copy D's row/columns for the kept indices into the top left of new_D
        put new_distances in the last row and last column of new_D (new node goes last)
        current_labels = labels at the kept indices, then new_label added at the end
        D = new_D

    Step 6. Exactly 2 nodes left: connect them  (The final branch is split evenly between two sides so the tree can be drawn)
    final_distance = D[0, 1] / 2
    tree = "(" + current_labels[0] + ":" + f"{final_distance:.4f}" + "," + current_labels[1] + ":" + f"{final_distance:.4f}" + ");"

    Return tree    

# Results Output



# Successes
The notebook reads the FASTA file, builds the distance matrix, runs Neighbor-Joining, and draws the tree. Running
Run All on the real data set produced a plotted phylogenetic tree with the sequence IDs as leaf labels.Instead of rewriting the alignment, we kept  `cal_score` function from project05 and wrapped it in the `smith_waterman` function that returned a single normalized score between 0 and 1.Dividing the best local score by the score of a perfect match of the shorter sequence meant identical sequences gave a similarity of exactly 1.0, which matched the docstring example and made the conversion to a distance (`1 - similarity`) simple.

# Struggles

# Reflection and Adaptive Learning

# Personal Reflections
## Group leader 

## Other member

## Other member

# Generative AI Appendix
