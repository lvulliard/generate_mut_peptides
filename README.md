# generate_mut_peptides
## Generate mutated peptides in fasta format
Usage:
```
  generate_mut_peptides.py -n <nbSubstitutions> -p <refPeptidome> -m <substitutionMatrix> [-o <outputPrefix> -i]
  generate_mut_peptides.py --help
```
Options:
```
  -n --number=<nbSubstitutions>      Number of substitutions to perform
  -p --peptidome=<refPeptidome>      Reference peptidome to be mutated (fasta format) 
  -m --matrix=<substitutionMatrix>   Path to the substitution matrix pickle to use
  -o --output=<outputPrefix>         Prefix for the fasta output files
  -i --intervals-only                Output only the mutated regions, up to 10 flanking amino acids on each side of the mutations
  -h --help                          Show this screen.
```
