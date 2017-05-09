#!/usr/bin/env python
# -*- coding: utf-8 -*-
##################################### Help ####################################
"""Generate mutated peptides in fasta format
Usage:
  generate_mut_peptides.py <nbSubstitutions> [-p <refPeptidome>]
  generate_mut_peptides.py --help

Options:
  -p --peptidome     Reference peptidome to be mutated (fasta format) 
  -h --help          Show this screen.
"""

################################### Imports ###################################
from docopt import docopt
import numpy as np
import Bio
from Bio import SeqIO
from Bio.Seq import Seq, MutableSeq
from collections import defaultdict
from Bio.SeqIO.FastaIO import SimpleFastaParser

################################## Functions ##################################

################################### Classes ###################################		
class peptide:
	"""Peptide to be mutated"""

	def __init__(self):
		self.nbSub = 0 # Number of substitutions to perform
		self.normPep = Seq("") # Normal peptide
		self.mutPep = MutableSeq("") # Reference peptide

	def addSub(self):
		self.nbSub += 1

	def setPep(self, pep):
		self.normPep = self.normPep + pep
		self.mutPep = self.mutPep + pep

	def __str__(self):
		return(str(self.nbSub)+"\n"+str(self.normPep)+"\n"+str(self.mutPep))

##################################### Main ####################################
arguments = docopt(__doc__)

# Number of amino acid substitutions to perform
nbSubstitutions = int(arguments["<nbSubstitutions>"])

# Number of peptides in the reference peptidome
nbPep = 0
with open(arguments["<refPeptidome>"], "r") as handle:
	for i in SimpleFastaParser(handle):
		nbPep += 1

# Decide which peptide to mutate
pepDictIndices = np.random.randint(nbPep, size=nbSubstitutions)

# Store number of mutations
pepDict = defaultdict(peptide)
for i in pepDictIndices:
	pepDict[i].addSub()

