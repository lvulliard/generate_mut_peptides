#!/usr/bin/env python
# -*- coding: utf-8 -*-
##################################### Help ####################################
"""Generate mutated peptides in fasta format
Usage:
  generate_mut_peptides.py -n <nbSubstitutions> -p <refPeptidome> -m <substitutionMatrix>
  generate_mut_peptides.py --help

Options:
  -n --number=<nbSubstitutions>      Number of substitutions to perform
  -p --peptidome=<refPeptidome>      Reference peptidome to be mutated (fasta format) 
  -m --matrix=<substitutionMatrix>   Path to the substitution matrix pickle to use
  -h --help                          Show this screen.
"""

################################### Imports ###################################
from docopt import docopt
import numpy as np
import Bio
from Bio import SeqIO
from Bio.Seq import Seq, MutableSeq
from collections import defaultdict
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Alphabet import generic_protein
import pickle

################################## Functions ##################################

################################### Classes ###################################		
class peptide:
	"""Peptide to be mutated"""

	def __init__(self):
		self.nbSub = 0 # Number of substitutions to perform
		self.normPep = Seq("", generic_protein) # Normal peptide
		self.mutPep = MutableSeq("", generic_protein) # Reference peptide

	def addSub(self):
		self.nbSub += 1

	def setPep(self, pep):
		pep.seq.alphabet = generic_protein
		self.normPep = self.normPep + pep.seq
		self.mutPep = self.mutPep + pep.seq

	def __str__(self):
		return(str(self.nbSub)+"\n"+str(self.normPep)+"\n"+str(self.mutPep))

	def mutate(self, subMatrix):
		1

##################################### Main ####################################
arguments = docopt(__doc__)

aaList = 'ARNDCQEGHILKMFPSTWYVBZX'

# Number of amino acid substitutions to perform
nbSubstitutions = int(arguments["--number"])

# Number of peptides in the reference peptidome
nbPep = 0
with open(arguments["--peptidome"], "r") as handle:
	for i in SimpleFastaParser(handle):
		nbPep += 1

# Load substitution matrix
subMatrix = pickle.load(open(arguments["--matrix"], "r"))[:20,:20]

# Decide which mutation to perform
subProb = np.random.rand(nbSubstitutions)
subList = np.empty([nbSubstitutions, 2], dtype=int)
for k, p in zip(xrange(nbSubstitutions), subProb):
	cumProb = 0
	for i in xrange(20):
		for j in xrange(20):
			cumProb += subMatrix[i,j]
			if cumProb >= p:
				subList[k] =  [i,j]
				cumProb += 1
				break
		if cumProb >= 1:
			break

# Decide which peptide to mutate
pepDictIndices = np.random.randint(nbPep, size=nbSubstitutions)

# Store number of mutations
pepDict = defaultdict(peptide)
for i in pepDictIndices:
	pepDict[i].addSub()

# Read and store these peptides
subIndex = 0
for pep, i in zip(SeqIO.parse(arguments["--peptidome"], "fasta"), xrange(nbPep)):
	if i in pepDictIndices:
		# The residues to mutate need to be present in the selected sequence
		assert aaList[subList[subIndex,0]] in pep.seq
		subIndex += 1
		pepDict[i].setPep(pep)

# Perform mutations
for i in pepDictIndices:
	pepDict[i].mutate(subMatrix)
