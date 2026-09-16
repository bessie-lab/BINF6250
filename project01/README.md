# Introduction
The program reads a variant call format file and creates a dictionary used to count the number of times a disease occurs

# Pseudocode

```
Parse line function 
1) The function takes a string as an argument
2) Clean and split the argument
3) grab the info column
2) Use conditionals to check if AF_EXAC is present
3) create a dictionary and convert it to a float and check if the < 0.001 condition is met
4) determine rarity of the disease using conditionals and return the list 

Read line function
1) Create an empty dictionary for counting disease occurrence
2) Open the vcf file
3) Read each line
4) Put each line as an argument into the parse line function
5) Check if the output of the parse line function is not an empty list
6) If so, take each disease in the list and check if it shows up in the dictionary
  a) If it does, increase its count by 1
  b) if it doesn't, add it to the dictionary and give it a count of 1
7) After doing this for all lines, return the dictionary so it can be printed at the end of the script
```

# Successes
- Learned how to use GitHub in a collaborative way (making changes and seeing it updated on a team member's github)
- Learned how to collaboratively write pseudocode
- Learned how to parse vcf files and extract important data

# Struggles
- Had trouble figuring out how to fork a project leader's GitHub and have it reflect all branches
- Had trouble seeing version history in a pull request 

# Personal Reflections
## Group Leader
I learned how to create a repo, create a branch, and work on that branch collaboratively with a teammate. This was really useful for learning how to isolate my work on a branch and I can see how I could do work independent of this collaboration on a different branch and why that would be useful. This was a really helpful exercise in becoming more comfortable with collaborative programming utilizing Github. 

## Other member
I learned how to fork a teammate's repository and how to commit changes and pull requests on Github. This was a great experience as it taught me on how we can work on the same program collaboratively . I also learnt how to open vcf files corrrectly and remove the tabs and splits in them in my code.

# Generative AI Appendix
Used Claude AI to debug code in the parse_file function
