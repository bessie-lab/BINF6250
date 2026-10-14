# Introduction
This project implements a de Bruijn graph-based approach to genome assembly using Python.The implementation represents k-mers as directed edges between nodes formed by their prefixes and suffixes. By constructing a de Bruijn graph from sequencing reads, the program uses Eulerian walks to traverse the graph and assemble overlapping sequence fragments into longer contiguous sequences (contigs).The project explores key concepts in computational biology, including graph representation, k-mer analysis, graph traversal, and genome reconstruction. It provides a practical demonstration of how graph algorithms can be applied to solve problems in bioinformatics.

In this project we evaluate two datasets. The first is our toy dataset, with 9 reads of input which we are expected to build the deBrujin graph with k=6. The expected contig should be ATGGCGTACGTTACCATGGGCCTAAA. The second dataset involves reconstructing the yeast genome, with k-mer size in driver program as 31. It contains 2.4 million 150 bp reads with 30X coverage across the genome. With both datasets, we assume there is no genetic variation or sequencing errors involved. 

# Pseudocode

```
def add_edge(self, left, right):
  append right to graph[left]

def remove_edge(self, left, right):
  self.graph[left].remove(right)

def build_graph_from_reads(self, reads, k):
  for each read in reads:
            for i from 0 to len(read) − k:
                kmer  = read[i : i+k]
                left  = kmer[0 : k−1] # drop last base for prefix
                right = kmer[1 : k]   # drop first base for suffix
                add_edge(left, right)

def eulerian_walk(self, node, graph, seed=None, tour=None):
  tour = []
        while graph[node] still has edges:
            next_node = pick a random neighbor from graph[node]
            remove_edge(node, next_node) #consume edge                 
            tour.extend(eulerian_walk(next_node, graph)) #recurse 
        tour.append(node) # append on the way back out
        return tour # deepest node was added first, want to reverse it at the end

def assemble_contigs(self, seed=None):
  seed the RNG once
  graph = a copy of self.graph
  compute indegree of every node

  contigs = []
  #Start from ideal nodes which have more outgoing than incoming edges to start path
  for node in nodes where out > in:
      while graph[node] has edges:
          tour = reverse(eulerian_walk(node, graph))
          contigs.append(tour_to_sequence(tour))
  #Anything left over is a cycle, start anywhere that still has edges
  for node in all nodes:
      while graph[node] has edges:
          tour = reverse(eulerian_walk(node, graph))
          contigs.append(tour_to_sequence(tour))
  return contigs

def tour_to_sequence(self, tour):
  if tour is empty: return ""
  seq = [tour[0]]
  for node in tour[1:]:
      seq.append(node[-1])
  return "".join(seq)

def get_assembly_stats(self, contigs):
  lengths = sorted(len(c) for c in contigs, descending)
  num_contigs, total_length, longest = lengths[0], shortest = lengths[-1], mean
  n50: walk down the sorted lengths keeping a running sum;
      the first length where running sum >= total_length / 2 is the N50
  (handle empty contig case by returning zeros)


def write_fasta(self, contigs, filename):
  for i, contig in enumerate(contigs, 1):
    write ">contig_i length=N"
    write the sequence 


```

# Successes
The program successfully constructs a de Bruijn graph from DNA sequencing reads by breaking the sequences into k-mers and connecting their prefixes and suffixes to represent relationships between sequence fragments. It represents the graph using an adjacency list and supports graph manipulation by tracking connections between nodes. The program applies Eulerian walk concepts to traverse the graph, aiming to use each edge once to reconstruct sequences from the input data. Through this process, it can assemble overlapping sequence fragments into longer contigs. We developed an understanding of how de Bruijn graphs can be used to represent DNA sequence fragments and support genome assembly. Additionally, we were able to, through the code understand the usage of OOP, and modifying the deBrujin class and its attributes. For the toy dataset, we succesfully reconstructed the expected contig ATGGCGTACGTTACCATGGGCCTAAA , which gave us hope that the code could be implemented in the yeast dataset. The yeast dataset finished its run after 81 minutes, and produced 30 assembled contigs. Here are some of the key metrics we had built:

  • Total assembly: 22,958,782 bp
  • N50: 22,639,116 bp
  • Longest contig: 22,639,116 bp
  • Coverage: 30.0x

Throughout both projects, adding in `print` statements to indicate how many reads were processed throughout the run was very helpful for us to gauge that the program was successfully running, and where was taking the most time.

# Struggles
Our main struggles as a team were finding a time that worked for everyone given our schedules, and figuring out how to optimize the code for the yeast dataset. The toy dataset took a very short time to run and it outputted the contig as expected. The yeast dataset took 81 minutes to run compared to the expected 6-7 minutes, and finding ways to optimize the code was a struggle for us. One idea was to modify the `for` loop in the graph building step, although we did not get a chance to fully think that out. Also, conceptualizing recursion was a crucial barrier for us to overcome as it is an essential component of the Eulerian walk. Overall, it felt like there was still more to learn from this project that we maybe have not yet grasped, particularly in optimizing or modifying our code for it to successfully run in the yeast genome. 

# Personal Reflections
## Group Leader - Bessie
Through this project, I gained practical experience applying graph theory and Python programming to genome assembly in bioinformatics. I learned how to construct a de Bruijn graph by representing k-mers as directed edges between nodes and how Eulerian walks can be used to traverse the graph during sequence assembly. I strengthened my understanding of object-oriented programming, graph data structures, and algorithm design while developing methods to construct and manipulate graphs. Additionally, I improved my problem-solving and debugging skills by breaking the assembly process into smaller steps and addressing challenges associated with graph traversal and sequence reconstruction. Overall, this project deepened my understanding of how computational algorithms can be applied to biological problems, particularly the reconstruction of DNA sequences from smaller fragments.

## Yulia
This project was a challenge since we all had different schedules and it was difficult to work collaboratively as much. Bessie and Abigail had elected to meet first to brainstorm and start some work on code, and after I had finished the code based on their work. Bessie also wrote out her version of the code too, so we chatted on Teams to compare the codes and see how we could make the yeast dataset more efficient. It took my version 81 minutes for it to run. It took the longest to run at the building debrujin graph step. Assembling contigs took the second longest amount of time. I do suspect my memory could be a factor here, as I did not run it on an HPC. Marcus did mention it should take 6-7 minutes to run, but I wonder if he ran it on the HPC. I think the main bottleneck here was in the `for` loop. Also I wonder if forward and reverse strands are playing a role here. That's a thought that occured to me later in the project, so maybe I would try modifying the code to choose one of the strands and maybe that way the program would be more efficient.

Adding in `print` statements was helpful for me to visualize the progress of the datasets while they ran so it would not sit in silence. Through his project I refreshed my memory on OOP and working with class attributes. I had to constantly remind myself that they were not functions and therefore did not need to `return` anything. Once I wrapped my brain around that, it was easier. 

## Abby
This project was challenging to complete due to scheduling constraints. However, it was interesting to be able to apply graph theory in an algorithm as we just talked about very similar concepts in another class. It was helpful to be able to use these skills over multiple different problem sets and identify the nuances and differences between the code implementation of each. The eularian walk was the most difficult for me to conceptualize because I struggle to understand recursive theories.

## Other member
Other members' reflections on the project

# Generative AI Appendix
Claude was used to modify pseudocode logic, code syntax for eulerian_walk, assemble_contigs, tour_to_sequence, get_assembly_stats. It was also used to answer questions about how recursion works. 
