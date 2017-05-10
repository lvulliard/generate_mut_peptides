#!/usr/bin/env python
# -*- coding: utf-8 -*-
##################################### Help ####################################
"""Generate mutated peptides in fasta format
Usage:
  generate_mut_peptides.py -n <nbSubstitutions> -p <refPeptidome> -m <substitutionMatrix> [-o <outputPrefix>]
  generate_mut_peptides.py --help

Options:
  -n --number=<nbSubstitutions>      Number of substitutions to perform
  -p --peptidome=<refPeptidome>      Reference peptidome to be mutated (fasta format) 
  -m --matrix=<substitutionMatrix>   Path to the substitution matrix pickle to use
  -o --output=<outputPrefix>         Prefix for the fasta output files
  -h --help                          Show this screen.
"""

################################### Imports ###################################
from __future__ import print_function
import sys
from docopt import docopt
import numpy as np
from collections import defaultdict
import pickle
import Bio
from Bio import SeqIO
from Bio.Seq import Seq, MutableSeq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Alphabet import generic_protein

################################## Functions ##################################
def eprint(*args, **kwargs):
	"""Print to stderr"""
	print(*args, file=sys.stderr, **kwargs)

################################### Classes ###################################		
class peptide:
	"""Peptide to be mutated"""

	def __init__(self):
		self.nbSub = 0 # Number of substitutions to perform
		self.sub = {} # Substitution realized
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

	def mutate(self, mutation):
		try:
			mutPos = np.random.choice([pos for pos, char in enumerate(self.mutPep) if char == aaList[mutation[0]] and pos not in self.sub.keys()])
			self.sub[mutPos] = mutation
			self.mutPep[mutPos] = aaList[mutation[1]]
		except ValueError:
			eprint("No "+aaList[mutation[0]]+" residues left in the mutated sequence.")
			return

##################################### Main ####################################
arguments = docopt(__doc__)
if not arguments['--output']:
	arguments['--output'] = ""

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
pepDictIndices.sort()

# Store number of mutations	
pepDict = defaultdict(peptide)
for i in pepDictIndices:
	pepDict[i].addSub()

# Read and store these peptides
subIndex = 0
for pep, i in zip(SeqIO.parse(arguments["--peptidome"], "fasta"), xrange(nbPep)):
	if i in pepDictIndices:
		# The residues to mutate need to be present in the selected sequence
		try:
			tempSubIndex = 0
			for j in pepDictIndices[pepDictIndices == i]:
				assert aaList[subList[subIndex+tempSubIndex,0]] in pep.seq
				tempSubIndex += 1
			subIndex += tempSubIndex
			pepDict[i].setPep(pep)
		except AssertionError:
			# Try the next peptide instead
			pepDictIndices[pepDictIndices == i] += 1
			# The current peptide will not have any mutation and the next one will have more
			pepDict[i+1].nbSub += pepDict[i].nbSub
			pepDict[i].nbSub = 0

# Perform mutations
for i, mut in zip(pepDictIndices, subList):
	pepDict[i].mutate(mut)

# Output
SeqIO.write([SeqRecord(s.mutPep) for s in pepDict.values()], arguments["--output"]+"mutatedRegion_mut.fa", "fasta")
SeqIO.write([SeqRecord(s.normPep) for s in pepDict.values()], arguments["--output"]+"mutatedRegion_norm.fa", "fasta")
