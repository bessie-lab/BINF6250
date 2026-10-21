# Project Description
The application of the **Smith-Waterman Algorithm** transformed Computational Biology and Bioinformatics. It established a mathematical framework to find the optimal `local sequence alignment` between two nucleotide or protein sequence, identifying pairs of segments across two long sequences that has maximum local similarity. The algorithm isolates highly conserved local regions or functional domains while ignoring non-conserved flanking regions.

The program has Key Features & Biological Significance:
* **Dynamic Programming Matrix**: Computes local alignment scores cell-by-cell by evaluating matches, mismatches, and insertion/deletion (indel) gap penalties.
* **Zero Threshold Option ("Free Start")**: Includes a minimum score threshold of `0` in the recurrence relation , preventing cumulative negative scores and allowing local alignments to start cleanly anywhere in the matrix.
* **Traceback Procedure**: Begins at the **global maximum score element** in the matrix (marking the end of the alignment) and steps backward along direction pointers until hitting a cell score of **`0`** (marking the start of the alignment).
* **Arbitrary Indel Handling**: Flexible scoring scheme accommodating match rewards, mismatch penalties, and gap penalties

## Pertinent documentation and usage write-up
### Requirements
- Python 3.11 or newer version 
- Jupyter Notebook (or an editor that runs notebook cells, like VS Code or PyCharm)
- Extra package Numpy needed. 


#Pseudocode
```
FUNCTION CAL_SCORE(matrix, seq1, seq2, i, j, match, mismatch, gap):
    # COMPARE SEQUENCE CHARACTERS
        IF seq1[i - 1] == seq2[j -1]
            match_score = match
        ELSE
            match_score =  mismatch
 
    # CALCULA PATH SCORES GROM THE PREDECESSORS
        diag_score = matrix[i - 1, j - 1] + match_score
        up_score   = matrix[i - 1, j] + gap
        left_score = matrix[i, j - 1] + gap
 
    # COMPUTE SCORE EITH THE ZERO THRESHOLD
        score = MAX(0, diag_score, up_score, left_score)
 
    # RECORD TRACEBACK DIRECTION
        IF score == 0
            move = (0) END
        ELSE IF score == diag_score
            move =  (1) DIAG  
        ELSE IF score == up_score
            move =  (2) UP    
        ELSE
            move = (3) LEFT  
 
        RETURN score, move
 
FUNCTION TRACEBACK(seq1, seq2, traceback_matrix, maximum_position):
 
    # Initialize lists:
    aligned_seq1 = empty list
    aligned_seq2 = empty list
 
    (current_row, current_col) = maximum_position
 
    WHILE current_row > 0 AND current_col > 0:
 
        current_move = traceback_matrix[current_row][current_col]
 
        IF current_move == END (0):
            BREAK
 
        ELSE IF current_move == DIAG (1)
            Append seq1[current_row - 1] TO aligned_seq1
            Append seq2[current_col - 1] TO aligned_seq2
 
            current_row = current_row - 1
            current_col = current_col - 1
 
        ELSE IF current_move == UP (2)
            Append seq1[current_row - 1] TO aligned_seq1
            Append "-" to aligned_seq2
 
            current_row = current_row - 1
 
        ELSE IF current_move == LEFT (3)
            Append "-" TO aligned_seq1
            Append seq2[current_col - 1] TO aligned_seq2
 
            current_col = current_col - 1
 
    # Reverse lists since traceback constructs the alignment backward
 
    aligned_seq1 = REVERSE(aligned_seq1) AND JOIN TO STRING
    aligned_seq2 = REVERSE(aligned_seq2) AND JOIN TO STRING
 
    RETURN aligned_seq1, aligned_seq2
 
FUNCTION SMITH_WATERMAN(seq1, seq2, match = +1, mismatch = -1, gap = -1):
    n = LENGTH(seq1)
    m = LENGTH(seq2)

    # Initialize scoring matrix and traceback matrix with zeros
    Initialize score_matrix of size (n + 1) × (m + 1) filled with 0s
    Initialize traceback_matrix of size (n + 1) × (m + 1) filled with 0s

    max_score = 0
    maximum_positions = []  // List to collect all cells sharing the maximum score

    # Forward Sweep: Fill score and traceback matrices cell-by-cell
    FOR i FROM 1 TO n
        FOR j FROM 1 TO m 
            (score, move) = CAL_SCORE(score_matrix, seq1, seq2, i, j, match, mismatch, gap)
            score_matrix[i][j] = score
            traceback_matrix[i][j] = move

            # Track global maximum cell positions (including ties)
            IF score > max_score 
                max_score = score
                maximum_positions = [(i, j)]        # Reset list with new highest position
            ELSE IF score == max_score AND score > 0 
                Append (i, j) to maximum_positions  # Add tie position

    # Backward Sweep: Traceback all optimal alignments
    optimal_alignment = []

    IF LENGTH(maximum_positions) == 1
        (aligned_seq1, aligned_seq2) = TRACEBACK(seq1, seq2, traceback_matrix, maximum_positions[0])
        RETURN aligned_seq1, aligned_seq2, score_matrix
    ELSE
        FOR EACH pos IN maximum_positions 
            Append TRACEBACK(seq1, seq2, traceback_matrix, pos) to optimal_alignment
        RETURN optimal_alignment, score_matrix
```
 
# Reflections
As we worked on this project, we developed a deeper understanding of the Smith-Waterman local alignment algorithm, particularly how its 4 option recurrence relation and zero reset mechanism help to identify regions of similarity between biological sequences. We learned how dynamic programming calculates alignment scores by evaluating matches, mismatches, and gaps while preventing negative scores from propagating through the scoring matrix.

While implementing the algorithm we were able to understand the importance of the traceback process on the reconstruction of the optimal local alignment from the highest score. Through this project, we strengthened our understanding of how dynamic programming transforms individual scoring decisions into meaningful biological sequence alignments.



 
# Personal Reflections
## Group Leader - Bessie
I really enjoyed working on this project . My group members were very supportive and very helpful when I did not understand the coding part which am still learning.Before this project, I understood “local” versus “global” alignment abstractly. Implementing it showed me why the traceback starts at the matrix maximum rather than the corner, and why poorly matching regions get trimmed away instead of dragging down the score.The hardest part for me was handling ties between the diagonal/up/left moves, getting the traceback logic right , handling the sequencing indices and determining which data type would work best for storing the optimal alignments initially had thought a tuple could work but it was not the best idea as shown by the challenges that I faced when trying to run the program. I resolved this by revisiting the pseudocode that we had worked on with my group members and trying to read more about the steps that we were taking in the algorithm and make sense of what each step means for the program .I want to keep on practicing how to build the logic necessary for all these algorithms and also improve my coding skills.

## Member
Ahmed Salman - For this assignment, we followed the professors plan to a T. We let our weakest programmer take the wheel while the other two supervised them. In our case, me and Graziano watched over Bessie as she coded the entire program from scratch. Bessie was really inquisitive and asked a lot of questions! But I have to say she should have more confidence in her answers, but also learn how to test her code. I think that if she improves her ability at googling what she is stuck with, she is going to become a much better programmer really fast. She already has the motivation to do it. Back to the assignment, thankfully I had good understanding of the algorithm from the start so I mainly had a problem with keeping track of index values in the matrix. I kept putting -1 where it doesn't belong or not adding one when I should. I then decided to re read my own code from a fresh perspective and "drew" the matrix every step of the way by hand (naturally it wasn't the full matrix just a short example so I can keep track of positions). That is how I fixed the index issue I had. There is still room for improvement, I would like to test my code on a larger sequence and see how it performs. I don't think this algorithm would work very well for aligning a small sequence to a whole genome as the `n` times `m` value sky rockets when using a whole genome as one strand. That is why the BLAST algorithm was created. I wonder how hard it would be to program a BLAST algorithm.
 
## Member - Graziano Peregrino  
### Issue
One of the challenges I faced during this project was understanding the small details and individual steps of the Smith-Waterman algorithm. Initially, I struggled to understand how the scoring calculation, recurrence relation, and traceback process would work together. As I was developing a better understanding of these functions, I was able to collaborate more effectively with my group. 

### What I Learned
My takeaway from this project it is how important to work together to reach the goal. As a group we were to help each other on our individual growth and development. We were able to work on the dynamic programming issue and how to get to the correct resolution of the problem. I also gained a better understanding of how breaking a complex problem into smaller one can be effective to find a solution. Working through the algorithm together helped both my technical reasoning and problem solving skills.

### Next Steps
As next steps I want to keep working on my reasoning skills and algorithmic thinking, specially when working with dynamic programming. I want to become more confident in my understanding of the methods so I can contribute more effectively to group and support my teammates. 

# Generative AI Appendix
 
Description of which generative AI was used and its version.
- NotebookLM 2.0 and Gemini 4 Argon
 
Updated the book and articles to notebookml and start asking about details of dynamic programming chapter.
 
The entire prompt that was used to generate the content.
- I need some help to understand dynamic programming and the walk of the nodes in the matrix.
- Can you show me a numerical example of filling a grid?
- I want to understand why the algorithm chooses the path with highest weight.
- Can you explain the diagonal reasoning?
- Can you explain dynamic programming in an arbitrary DAG?
- How can I calculate the weight of the diagonal?
- Can you explain the topological order on the dressing challenge problem?
- How would I modify the alignment graph to solve the gap penalty issue?
- How do affine gap penalties work?
- Can you provide some explanation on this The Four Option Recurrence Relation?
- How could we construct an alignment without storing any backtracking pointers?
 
An explanation of how it was used (e.g., to generate ideas).
- I used it as a tool to learn the concepts described on the chapter, asking for further clarification on the material and extra examples.
 
A justification for why generative AI was used.
- Sometimes just reading a concept and formula is not enough for me to have a grasp on the material, so I used it as an additional source to lern about dynamic programming.

# Reference
-Compeau, P., & Pevzner, P. (2015). Bioinformatics Algorithms: An Active Learning Approach (2nd ed., Vol. 1)
. Active Learning Publishers, LLC
