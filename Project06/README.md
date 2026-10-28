# BINF6250 - Project 06

# Introduction
This project implements a distance-based phylogenetic tree builder that applies the Neighbor-Joining algorithm to 20 HIV-1 reverse transcriptase (RT) protein sequences, each 229 amino acids long, from `lafayette_SARS_RT.fasta`. Every pair of sequences is aligned with the Smith-Waterman local alignment algorithm (match +2, mismatch −1, gap −1), and the best local score is normalized by the score of a perfect match of the shorter sequence, giving a similarity between 0 and 1. Each similarity is converted to a distance (`distance = 1 - similarity`) to build a symmetric 20 × 20 distance matrix. Neighbor-Joining then repeatedly joins the pair of nodes with the lowest D* (Q) score, calculates the two new branch lengths, and replaces the pair with a new node in a smaller matrix until two nodes remain. The result is an unrooted tree, written in Newick format and plotted with Biopython and Matplotlib. Neighbor-Joining was chosen over UPGMA because it does not assume a constant rate of evolution across lineages, so branch lengths can differ between sequences.

All required functionality described in the assignment instructions was completed. The `read_fasta()`, `smith_waterman()`, `build_distance_matrix()`, and `neighbor_joining()` functions were implemented, and the full pipeline runs in under a minute and produces a tree containing all 20 sequences. The Neighbor-Joining implementation was also tested on the 4-sequence example matrix from class, where it recovered the expected branch lengths of 11, 2, 6, 7, and 4.

Although the file name contains "SARS," every sequence header identifies the sequences as Human immunodeficiency virus 1 reverse transcriptase, which matches the assignment instructions. The file name was kept as provided so the notebook runs without changes.

# How it works
Read the data: sequences are loaded from a FASTA file(`lafayette_SARS_RT.fasta`).
Score pairs: every pair of sequences is aligned with Smith-Waterman.The best local score is normalized to a value   between 0 and 1, where 1 means a perfect match.
Build distances: each similarity is converted to a distance(`distance = 1 - similarity`), giving a symmetric distance matrix.
Join neighbors: Neighbor-Joining repeatedly merges the pair of nodes with the smallest Q-value and updates the distance matrix until two nodes remain.
Draw the tree: the result is written in Newick format and plotted with Biopython and Matplotlib.

# Files

- project06.ipynb - notebook with the implementation
- lafayette_SARS_RT.fasta -input sequences
- nj_tree.png - image of the final tree, shown in Results Output

# Requirements

Python 3, `numpy`, `matplotlib`, and `biopython`

pip install numpy matplotlib biopython

# Usage

Place the FASTA file in the same folder as the notebook, open it in VS Code or Jupyter, and choose **Run All**. The final cell reads the sequences, builds the distance matrix, runs Neighbor-Joining, and plots the tree.

# Pseudocode
## Final pseudocode that reflects submitted pseudocode  (See Appendix Starting Pseudocode at end of file for details on start and notes on changes)

```
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
    1. cal_score calculates each cell's score as in Project 05. The traceback move it records is not used.
    2. Traceback and the traceback matrix are dropped; only the highest score is needed.
    3. smith_waterman fills the score table itself (no helper function) and returns a normalized score:
       normalized = best local score / (match_score x length of the shorter sequence)

def cal_score(matrix, seq1, seq2, i, j, match, mismatch, gap):
    IF seq1[i-1] equals seq2[j-1]:
        match_score = match
    ELSE:
        match_score = mismatch
    diag_score = matrix[i-1, j-1] + match_score
    up_score   = matrix[i-1, j] + gap
    left_score = matrix[i, j-1] + gap
    score = biggest of 0, diag_score, up_score, left_score
    move = direction the score came from (0 END, 1 DIAG, 2 UP, 3 LEFT)
    Return score, move

def smith_waterman(sequence1, sequence2, match_score = 2.0, mismatch_penalty = -1.0, gap_penalty = -1.0) -> float:
    n = length of sequence1
    m = length of sequence2
    IF either sequence is empty:
        return 0.0

    matrix = (n+1) x (m+1) grid of zeros
    FOR each row i from 1 to n:
        FOR each column j from 1 to m:
            score = cal_score(matrix, sequence1, sequence2, i, j, match_score, mismatch_penalty, gap_penalty)  (move is ignored)
            matrix[i, j] = score

    Return highest number in matrix / (match_score x smaller of n and m)
        Note: 1.0 = identical, closer to 0 = less similar


def build_distance_matrix(sequences: Dict[str, str]) -> Tuple[np.ndarray, List[str]]:
    ids = list of the sequence names, in file order (i.e. ["DM1", "DM2", ...])
    n = number of sequences
    distance_matrix = n x n grid of zeros (diagonal stays 0: each sequence vs itself)

    FOR each row i from 0 to n-1:
        FOR each column j from i+1 to n-1 (only pairs not yet done; skip the diagonal):
            similarity = smith_waterman(sequence for ids[i], sequence for ids[j])
            distance = 1 - similarity
            distance_matrix[i, j] = distance
            distance_matrix[j, i] = distance  (table is symmetric, so fill both halves)

    Return distance_matrix, ids


def neighbor_joining(distance_matrix: np.ndarray, labels: List[str]) -> str:
    SETUP
    D = copy of distance_matrix (so the original table isn't changed)
    current_labels = copy of labels (each starts as just a name, i.e. "DM1")
    n = number of rows in D

    WHILE n > 2:
        Step 1. Total distance: totals = sum of each row of D

        Step 2. D* score for every pair, keep the lowest
        lowest_score = infinity
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
        IF branch_i is negative:
            branch_i = 0, branch_j = D[best_i, best_j]
        ELSE IF branch_j is negative:
            branch_j = 0, branch_i = D[best_i, best_j]

        Step 4. Distances from the new node to every other node k
        new_distances = empty list
        FOR each k from 0 to n-1:
            IF k is best_i or best_j: skip it
            new_distances.append((D[best_i, k] + D[best_j, k] - D[best_i, best_j]) / 2)

        new_label = "(" + label of best_i + ":" + branch_i (4 decimals) + ","
                        + label of best_j + ":" + branch_j (4 decimals) + ")"

        Step 5. Shrink the table
        keep = every index except best_i and best_j
        D = rows and columns of D at the kept indices, plus one extra row and column of zeros
        last row and last column of D = new_distances (new node goes last)
        current_labels = labels at the kept indices, then new_label added at the end
        n = number of rows in D

    Step 6. Exactly 2 nodes left: connect them (final branch split evenly)
    final_distance = D[0, 1] / 2
    tree = "(" + current_labels[0] + ":" + final_distance (4 decimals) + ","
               + current_labels[1] + ":" + final_distance (4 decimals) + ");"

    Return tree


def plot_tree(newick_tree: str) -> None:
    (provided by the template; one line added)
    save the figure as nj_tree.png before showing it, so it can be linked in the README
```

# Results Output

![Neighbor-Joining tree of 20 HIV-1 RT sequences](nj_tree.png)

Newick output from `neighbor_joining()`:

```
((P6:0.0278,((P5:0.0573,P12:0.0213):0.0034,(((P10:0.0003,P13:0.0128):0.0111,(P3:0.0021,P15:0.0241):0.0020):0.0050,(P8:0.0082,(DM1:0.0066,DM2:0.0000):0.0180):0.0008):0.0046):0.0022):0.0002,(P9:0.0245,((P11:0.0473,(DM3:0.0040,(JT1:0.0000,JT2:0.0000):0.0091):0.0116):0.0196,(P7:0.0493,(P2:0.0405,(P1:0.0767,(P4:0.0332,P14:0.0781):0.0073):0.0075):0.0170):0.0032):0.0046):0.0002);
```

The tree contains all 20 HIV-1 reverse transcriptase sequences. Overall, the sequences are highly similar: the largest pairwise distance in the matrix is 0.183, so most branches are short. JT1 and JT2 join with branches of length 0 because their sequences are identical, and DM1 and DM2 form the closest pair, with DM2's branch also at 0 (a tiny negative rounding error set to 0). DM3 does not group with DM1 and DM2; it clusters with JT1 and JT2 instead. The longest terminal branches belong to P14 (0.0781) and P1 (0.0767), meaning these sequences differ most from their nearest relatives in the tree. Because Neighbor-Joining produces an unrooted tree, the split at the far left of the plot (where the final branch was divided evenly) is not a true common ancestor; only the groupings and branch lengths are meaningful. 

# Successes
The notebook reads the FASTA file, builds the distance matrix, runs Neighbor-Joining, and draws the tree. Running the real data set produced a plotted phylogenetic tree with the sequence IDs as leaf labels.Instead of rewriting the alignment, we kept  `cal_score` function from project05 and wrapped it in the `smith_waterman` function that returned a single normalized score between 0 and 1. Dividing the best local score by the score of a perfect match of the shorter sequence meant identical sequences gave a similarity of exactly 1.0, which matched the docstring example and made the conversion to a distance (`1 - similarity`) simple.

# Struggles
Missing packages: The notebook failed at the import cell with `ModuleNotFoundError: No module named 'matplotlib'`. Installing `matplotlib` and `biopython` with `%pip install` inside the notebook, then restarting the kernel, fixed it. 

A near-zero negative branch: DM2's branch length came out as a tiny negative number (a floating-point rounding error) and was set to 0, so DM2 appears with no visible branch in the plot. With more time, we would separate the two possible causes of a negative branch. A value within a tiny tolerance of zero (for example, smaller than 1 × 10⁻¹⁰) is a rounding error and can safely be treated as exactly 0. A clearly negative value is different, it means the distances do not fit a tree perfectly, which is common with real sequence data. Our code handles both the same way, setting the branch to 0 and moving the difference to the other branch, which preserves the pair's total distance but does not record which case occurred. A future version could count and report real negative branches separately, since how often they occur could be a useful measure of how well the distances fit a tree. Using a distance correction like -log(similarity) instead of 1 - similarity could also make the distances fit a tree more closely and reduce negative branches.

Docstring examples that don't match the data: The template's docstrings expect 19 sequence IDs and 18 commas in the Newick string, but `lafayette_SARS_RT.fasta` contains 20 sequences, so our tree has 20 leaves and 19 commas. A binary tree with 20 leaves always has 19 commas, so the difference comes from the docstring, not the code. The docstring's example distance matrix (0.004, 0.012, 0.012 for the first three sequences) also could not be replicated: our values are 0.007, 0.046, and 0.039. The ordering matches (DM1 and DM2 are the closest pair), so the difference is most likely due to the template being made with different scoring or normalization. The docstrings were left as original.

# Personal Reflections
## Group leader
I could follow the Neighbor-Joining formulas, but I didn't fully understand why subtracting the row sums in the Q-matrix picks out true neighbors instead of just the two closest sequences. I also wasn't sure what it means when the algorithm produces negative branch lengths. I set them to zero so the tree would plot, but that hides the fact that my distances don't fit a perfect tree,and I don't know how much that affects the result. Finally, I converted Smith-Waterman scores into distances using `distance = 1 - normalized score`,and I'm not sure how well that simple formula represents real evolutionary distance. I also need to improve and get comfortable using the Numpy package it would make writing the neighbour-joining function easier and very efficient

## Other member: Justin Gubbens
I really enjoyed working on this project. Ildiko put together very in-depth pseudocode before our meeting over the weekend, which was very helpful in getting my bearings for the project as a whole. Although we started from the same pseudocode, I think we ended up with very different code in some areas, as I noticed a few areas where numpy functions could be used to make the code more efficient/concise and where some function calls could be hardcoded with consistent logic to decrease runtime. When I was reviewing Bessie's individually written code, I noticed that she used a lot of numpy functions that took some time for me to understand even with a background in computer science, so I recommended that we finalized the project using clearer (although maybe slightly less efficient) code using loops to ensure that everybody could understand it. Overall, I found this project to be a helpful exercise in practicing numpy and a great way to build off of our work from Project 5.

## Other member: Ildiko Polyak
The most valuable part of this project for me was comparing our three individual versions of the code. We all started from the same pseudocode, yet the three implementations ended up looking very different: different ways of normalizing the Smith-Waterman score, of shrinking the distance matrix, and of building the Newick string. Running my version side by side with the final group version and confirming they produced identical distance matrices and an identical tree showed me that very different code can implement the same algorithm correctly. Since I'm still a beginner at coding, reading how others solved the same problem taught me techniques I wouldn't have found on my own, like normalizing with a single alignment instead of three and using array operations to resize the matrix. It also showed me how much clear comments matter, code without them took me much longer to understand, even when it was well written.  Justin combined portions from everyone's code for the final program and it was interesting to see his decision process of what makes a program most efficient, yet understandable.  Code efficiency is still something I'm trying to learn more about and it feels like this vauge subjective topic that is developed over time and experience, rather than hardset rules that a beginner can memorize.

I'll admit that evolution and phylogeny have never been topics that interest me (BO-RING). What did catch my interest was applying the method to HIV reverse transcriptase sequences and comparing viral variants. It was a good lesson that the same algorithm can serve very different biological questions, from evolutionary history to tracking viral strains and drug resistance.  I often get so caught up in the coding that I forget about the biological aspect and this lesson was finally a little interesting and held my attention and refocused why we learn these hard algorithms in the first place.

# Appendix: Generative AI
Anthropic Claude Opus 5.5 prompt: "How do I add an phylogenetic tree image result output from my .ipynb file into my README.md file?"

Claude was used to adapt the Smith-Waterman function from project05.

Claude was used to create pseudocode from the final program (that had AI help) to compare to starting pseudocode that was human written.  The prompt used: "I've uploaded the starting pseudocode and the final code, now write a final version of pseudocode that reflects the final version of code, so that I can compare starting and final pseudocode."

Claude was asked about how to deal with near zero value branches in distance based phylogeny algorithms.

Claude was used to summarize all of the references provided at the end of the module and then asked the following prompt: "comparing all of the references to the information in the module, are there any parts that might disagree?"

# Appendix: Starting Pseudocode
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
## Comparison of Starting and Final Pseudocode
| Function | Starting pseudocode | Final code |
|---|---|---|
| `read_fasta` | As written | Unchanged |
| `cal_score` | As written | Same calculation, different variable names |
| `smith_waterman` | Helper function called 3 times; divide by the smaller self-alignment score | No helper; one alignment; divide by match_score × shorter length. Same result for this data. |
| `build_distance_matrix` | As written | Unchanged |
| `neighbor_joining` Steps 1–4, 6 | As written | Same logic; written more compactly |
| `neighbor_joining` Step 5 | New grid of zeros, copied in with loops | Kept rows and columns plus a padded extra row and column (NumPy) |
| `neighbor_joining` loop | `WHILE` more than 2 nodes, n set at the top of each round | `WHILE n > 2`, n updated at the end of each round |
| `plot_tree` | Template | Template + `savefig` line |