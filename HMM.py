#!/usr/bin/env python

#from __future__ import print_function

import os

if os.name == 'java':
    from Bio import MissingExternalDependencyError
    #This is a slight miss-use of MissingExternalDependencyError,
    #but it will do in the short term to skip this unit test on Jython
    raise MissingExternalDependencyError("This test can cause a fatal error "
        "on Jython with some versions of Java")

# standard modules
import random

from Bio import SeqIO

# biopython
from Bio import Alphabet
from Bio.Seq import MutableSeq
from Bio.Seq import Seq

# HMM stuff we are testing
from Bio.HMM import MarkovModel
from Bio.HMM import Trainer
from Bio.HMM import Utilities

# whether we should print everything out. Set this to zero for
# regression testing
VERBOSE = 0

Pfad ='D:\\UNI\\Semester5\\AlBi\\Praktikum\\Praktikum3\\'
gb_file = (Pfad + 'sequence.gb')
gene_list = open(Pfad + 'gen_list.txt', 'w')
seq_length = 0

seq_file = open(Pfad + 'sequence.fasta', 'r')

my_seq = ''		#the original sequence is saved as a my_seq string
for seq_record in SeqIO.parse(seq_file, 'fasta'):
	my_seq = seq_record.seq

#print my_seq[:0+5]
	

for record in SeqIO.parse(gb_file, 'genbank'):
	seq = str(record.seq)
	#print record.features
	#print seq[175], seq[176], seq[177], seq[178], seq[179], seq[180], seq[181], seq[182], seq[183], seq[184], seq[185]
	seq_length = len(seq)
	#print dir(record)

	for feature in record.features:
		#print dir(feature)
		if feature.type == 'CDS':
			if 'gene' in feature.qualifiers:
				if '(+)' in str(feature.location):
					geneName = feature.qualifiers['gene'][0]
					if 'locus_tag' in feature.qualifiers:
						#print ('%s\t%s\t%s\n' %(geneName, feature.location, feature.qualifiers['locus_tag']))
						gene_list.write('%s\t%s\t%s\n' %(geneName, feature.location, feature.qualifiers['locus_tag']))
					else:
						#print ('%s\t%s\n' %(geneName, feature.location))
						gene_list.write('%s\t%s\n' %(geneName, feature.location))
gene_list.close()


a = open(Pfad + 'gen_list.txt', 'r')
#code_bereich = open(Pfad + 'kodierenden_bereich.txt', 'w')
#separieren des kodierenden Bereichs des Genes und umschreiben der Bereiche in einer neuen Datei namens kodierenden_bereich.txt


list_anfang_codierende = []	#speichert start position aller kodierenden bereiche
list_end_codierende =  []	#speichert end position aller kodierenden bereiche
for line in a:
	cols = line.split('\t')
	
	if 'join' in cols[1]:
		bereich = cols[1].split(',')
		a_codierende = bereich[0].split(':')
		anfang_codierende = a_codierende[0].strip('join{[')
		e_codierende = bereich[len(bereich)-1].split(':')
		end_codierende = e_codierende[len(e_codierende)-1].strip('](+)}')
		
		list_anfang_codierende.append(anfang_codierende)
		list_end_codierende.append(end_codierende)
		#print '%s\t%s' %(anfang_codierende, end_codierende)
		#aa = int(end_codierende)-int(anfang_codierende)
		#code_bereich.write('%s\t%s\t%s\n' %(anfang_codierende, end_codierende, aa))
		
	else:
		bereich = cols[1].split(':')
		anfang_codierende = bereich[0].strip('[')
		end_codierende = bereich[1].strip('](+)')
		
		list_anfang_codierende.append(anfang_codierende)
		list_end_codierende.append(end_codierende)
		#print '%s\t%s' %(anfang_codierende, end_codierende)
		#aa = int(end_codierende)-int(anfang_codierende)
		#code_bereich.write('%s\t%s\t%s\n' %(anfang_codierende, end_codierende, aa))
	
#code_bereich.close()
#print list_anfang_codierende

#nicht_code_bereich = open(Pfad + 'nicht_kodierenden_bereich.txt', 'w')
list_anfang_nichtcodierende = []	#speichert start position aller nicht-kodierenden bereiche
list_end_nichtcodierende =  []	#speichert end position aller nicht-kodierenden bereiche

# In dem Fall dass der Sequenz mit einem nicht-kodierenden Bereich anfaengt
if list_anfang_codierende[0] != 1:
	anfang = 1
	ende = int(list_anfang_codierende[0]) - 1
	list_anfang_nichtcodierende.append(anfang)
	list_end_nichtcodierende.append(ende)
	#aa = int(ende) - int(anfang)
	#nicht_code_bereich.write('%s\t%s\t%s\n' %(anfang, ende, aa))

for i in range(0, len(list_anfang_codierende)-1):
	if (list_end_codierende[i] < list_anfang_codierende[i+1]):
		if int(list_end_codierende[i])+1 != int(list_anfang_codierende[i+1]):
			anfang = (int(list_end_codierende[i]))+1
			ende = (int(list_anfang_codierende[i+1]))-1
			
			list_anfang_nichtcodierende.append(anfang)
			list_end_nichtcodierende.append(ende)
			#aa = int(ende) - int(anfang)
			#nicht_code_bereich.write('%s\t%s\t%s\n' %(anfang, ende, aa))

#falls der Sequenz mit einem nicht-kodierenden Bereich endet
if int(list_end_codierende[-1]) < int(seq_length):
	list_anfang_nichtcodierende.append(int(list_end_codierende[-1])+1)
	list_end_nichtcodierende.append(int(seq_length))
	#aa = int(seq_length) - (int(list_end_codierende[-1])+1)
	#nicht_code_bereich.write('%s\t%s\t%s\n' %( int(list_end_codierende[-1])+1, int(seq_length), aa))

#nicht_code_bereich.close()

#print list_anfang_nichtcodierende
#print list_end_nichtcodierende
#print list_anfang_codierende

anzahlKodA = 0
anzahlKodT = 0
anzahlKodG = 0
anzahlKodC = 0
#findet die Anzahl jeder Nukleotide im kodierenden Bereich
for i in range(0, len(list_anfang_codierende)):
	#length = int(list_end_codierende[i]) - int(list_anfang_codierende[i])
	seq_now = my_seq[int(list_anfang_codierende[i]):int(list_end_codierende[i]) +1]
	anzahlKodA = anzahlKodA + seq_now.count('A')
	anzahlKodT = anzahlKodT + seq_now.count('T')
	anzahlKodG = anzahlKodG + seq_now.count('G')
	anzahlKodC = anzahlKodC + seq_now.count('C')

gesamtKodAnzahl = anzahlKodA + anzahlKodT + anzahlKodG + anzahlKodC


print '%s\t%s' %('anzahlKodA:', anzahlKodA)
print '%s\t%s' %('anzahlKodT:', anzahlKodT)
print '%s\t%s' %('anzahlKodG:', anzahlKodG)
print '%s\t%s' %('anzahlKodC:', anzahlKodC)
print '%s\t%s\n' %('gesamtKodAnzahl:', gesamtKodAnzahl)


anzahlNichtKodA = 0
anzahlNichtKodT = 0
anzahlNichtKodG = 0
anzahlNichtKodC = 0
#findet die Anzahl jeder Nukleotide im nichtkodierenden Bereich
#len(list_anfang_nichtcodierende)
for i in range(0,len(list_anfang_nichtcodierende)):
	#length = int(list_end_nichtcodierende[i]) - int(list_anfang_nichtcodierende[i])
	seq_now = my_seq[int(list_anfang_nichtcodierende[i]):int(list_end_nichtcodierende[i])+1]
	#print seq_now
	anzahlNichtKodA = anzahlNichtKodA + seq_now.count('A')
	anzahlNichtKodT = anzahlNichtKodT + seq_now.count('T')
	anzahlNichtKodG = anzahlNichtKodG + seq_now.count('G')
	anzahlNichtKodC = anzahlNichtKodC + seq_now.count('C')

gesamtNichtKodAnzahl = anzahlNichtKodA + anzahlNichtKodT + anzahlNichtKodG + anzahlNichtKodC


print '%s\t%s' %('anzahlNichtKodA:', anzahlNichtKodA)
print '%s\t%s' %('anzahlNichtKodT:', anzahlNichtKodT)
print '%s\t%s' %('anzahlNichtKodG:', anzahlNichtKodG)
print '%s\t%s' %('anzahlNichtKodC:', anzahlNichtKodC)
print '%s\t%s\n' %('gesamtNichtKodAnzahl:', gesamtNichtKodAnzahl)


#fuer Emission der Basen in einem Gen
emissionGenA = float(anzahlKodA)/float(gesamtKodAnzahl)
emissionGenT = float(anzahlKodT)/float(gesamtKodAnzahl)
emissionGenG = float(anzahlKodG)/float(gesamtKodAnzahl)
emissionGenC = float(anzahlKodC)/float(gesamtKodAnzahl)
print '%s\t%s' %('emissionGenA:', emissionGenA)
print '%s\t%s' %('emissionGenT:', emissionGenT)
print '%s\t%s' %('emissionGenG:', emissionGenG)
print '%s\t%s\n' %('emissionGenC:', emissionGenC)

#fuer Emission der Basen in einem Nicht-Gen
emissionNichtGenA = float(anzahlNichtKodA)/float(gesamtNichtKodAnzahl)
emissionNichtGenT = float(anzahlNichtKodT)/float(gesamtNichtKodAnzahl)
emissionNichtGenG = float(anzahlNichtKodG)/float(gesamtNichtKodAnzahl)
emissionNichtGenC = float(anzahlNichtKodC)/float(gesamtNichtKodAnzahl)
print '%s\t%s' %('emissionNichtGenA:', emissionNichtGenA)
print '%s\t%s' %('emissionNichtGenT:', emissionNichtGenT)
print '%s\t%s' %('emissionNichtGenG:', emissionNichtGenG)
print '%s\t%s' %('emissionNichtGenC:', emissionNichtGenC)


anzahlGenMet = 0
anzahlGenAla = 0
anzahlGenArg = 0
anzahlGenAsn = 0
anzahlGenAsp = 0
anzahlGenCys = 0
anzahlGenGlu = 0
anzahlGenGln = 0
anzahlGenGly = 0
anzahlGenHis = 0
anzahlGenIle = 0
anzahlGenLeu = 0
anzahlGenLys = 0
anzahlGenPhe = 0
anzahlGenPro = 0
anzahlGenSer = 0
anzahlGenThr = 0
anzahlGenTrp = 0
anzahlGenTyr = 0
anzahlGenVal = 0
anzahlGenStopp = 0
anzahlGenTriplett = 0

#findet die Anzahl jeder Aminosaeure im Gen

for i in range(0,len(list_anfang_codierende)):
	seq_now = my_seq[int(list_anfang_codierende[i]):int(list_end_codierende[i])]
	#print '%s\n' %(seq_now)
	for j in range(2, len(seq_now), 3):
		#print '%s%s%s' %(seq_now[j-2], seq_now[j-1], seq_now[j])
		if seq_now[j-2] == 'A':
			if seq_now[j-1] == 'T':
				if seq_now[j] == 'G':
					anzahlGenMet += 1
					anzahlGenTriplett += 1
					
		if seq_now[j-2] == 'G':
			if seq_now[j-1] == 'C':
				anzahlGenAla += 1
				anzahlGenTriplett += 1

		if seq_now[j-2] == 'C':
			if seq_now[j-1] == 'G':
				anzahlGenArg += 1
				anzahlGenTriplett += 1
		
		if seq_now[j-2] == 'A':
			if seq_now[j-1] == 'A':
				if (seq_now[j] == 'C'or seq_now[j] == 'T'):
					anzahlGenAsn += 1
					anzahlGenTriplett += 1
					
		if seq_now[j-2] == 'G':
			if seq_now[j-1] == 'A':
				if (seq_now[j] == 'C' or seq_now[j] == 'T'):
					anzahlGenAsp += 1
					anzahlGenTriplett += 1
		
		if seq_now[j-2] == 'T':
			if seq_now[j-1] == 'G':
				if (seq_now[j] == 'C' or seq_now[j] == 'T'):
					anzahlGenCys += 1
					anzahlGenTriplett += 1
		
		if seq_now[j-2] == 'G':
			if seq_now[j-1] == 'A':
				if (seq_now[j] == 'A' or seq_now[j] == 'G'):
					anzahlGenGlu += 1
					anzahlGenTriplett += 1
		
		if seq_now[j-2] == 'C':
			if seq_now[j-1] == 'A':
				if (seq_now[j] == 'A' or seq_now[j] == 'G'):
					anzahlGenGln += 1
					anzahlGenTriplett += 1
		
		if seq_now[j-2] == 'G':
			if seq_now[j-1] == 'G':
				anzahlGenGly += 1
				anzahlGenTriplett += 1
				
		if seq_now[j-2] == 'C':
			if seq_now[j-1] == 'A':
				if (seq_now[j] == 'C' or seq_now[j] == 'T'):
					anzahlGenHis += 1
					anzahlGenTriplett += 1
					
		if seq_now[j-2] == 'A':
			if seq_now[j-1] == 'T':
				if (seq_now[j] == 'C' or seq_now[j] == 'T' or seq_now[j] == 'A'):
					anzahlGenIle += 1
					anzahlGenTriplett += 1
		
		if seq_now[j-2] == 'T':
			if seq_now[j-1] == 'T':
				if (seq_now[j] == 'A' or seq_now[j] == 'G'):
					anzahlGenLeu += 1
					anzahlGenTriplett += 1
		
		if seq_now[j-2] == 'C':
			if seq_now[j-1] == 'T':
				anzahlGenLeu += 1
				anzahlGenTriplett += 1
		
		if seq_now[j-2] == 'A':
			if seq_now[j-1] == 'A':
				if (seq_now[j] == 'A' or seq_now[j] == 'G'):
					anzahlGenLys += 1
					anzahlGenTriplett += 1
		
		if seq_now[j-2] == 'T':
			if seq_now[j-1] == 'T':
				if (seq_now[j] == 'C' or seq_now[j] == 'T'):
					anzahlGenPhe += 1
					anzahlGenTriplett += 1

		if seq_now[j-2] == 'C':
			if seq_now[j-1] == 'C':
				anzahlGenPro += 1
				anzahlGenTriplett += 1
		
		if seq_now[j-2] == 'T':
			if seq_now[j-1] == 'C':
				anzahlGenSer += 1
				anzahlGenTriplett += 1

		if seq_now[j-2] == 'A':
			if seq_now[j-1] == 'G':
				if (seq_now[j] == 'C' or seq_now[j] == 'T'):
					anzahlGenSer += 1
					anzahlGenTriplett += 1
					
		if seq_now[j-2] == 'A':
			if seq_now[j-1] == 'C':
				anzahlGenThr += 1
				anzahlGenTriplett += 1
				
		if seq_now[j-2] == 'T':
			if seq_now[j-1] == 'G':
				if seq_now[j] == 'G':
					anzahlGenTrp += 1
					anzahlGenTriplett += 1
		
		if seq_now[j-2] == 'T':
			if seq_now[j-1] == 'A':
				if (seq_now[j] == 'C' or seq_now[j] == 'T'):
					anzahlGenTyr += 1
					anzahlGenTriplett += 1

		if seq_now[j-2] == 'G':
			if seq_now[j-1] == 'T':
				anzahlGenVal += 1
				anzahlGenTriplett += 1

		if seq_now[j-2] == 'T':
			if seq_now[j-1] == 'A':
				if (seq_now[j] == 'A' or seq_now[j] == 'G'):
					anzahlGenStopp += 1
					anzahlGenTriplett += 1
		if seq_now[j-2] == 'T':
			if seq_now[j-1] == 'G':
				if seq_now[j] == 'A':
					anzahlGenStopp += 1
					anzahlGenTriplett += 1
				

print '%s\t%s\n' %('anzahlGenTriplett:', anzahlGenTriplett)

emissionGenMet = float(anzahlGenMet)/float(anzahlGenTriplett)
emissionGenAla = float(anzahlGenAla)/float(anzahlGenTriplett)
emissionGenArg = float(anzahlGenArg)/float(anzahlGenTriplett)
emissionGenAsn = float(anzahlGenAsn)/float(anzahlGenTriplett)
emissionGenAsp = float(anzahlGenAsp)/float(anzahlGenTriplett)
emissionGenCys = float(anzahlGenCys)/float(anzahlGenTriplett)
emissionGenGlu = float(anzahlGenGlu)/float(anzahlGenTriplett)
emissionGenGln = float(anzahlGenGln)/float(anzahlGenTriplett)
emissionGenGly = float(anzahlGenGly)/float(anzahlGenTriplett)
emissionGenHis = float(anzahlGenHis)/float(anzahlGenTriplett)
emissionGenIle = float(anzahlGenIle)/float(anzahlGenTriplett)
emissionGenLeu = float(anzahlGenLeu)/float(anzahlGenTriplett)
emissionGenLys = float(anzahlGenLys)/float(anzahlGenTriplett)
emissionGenPhe = float(anzahlGenPhe)/float(anzahlGenTriplett)
emissionGenPro = float(anzahlGenPro)/float(anzahlGenTriplett)
emissionGenSer = float(anzahlGenSer)/float(anzahlGenTriplett)
emissionGenThr = float(anzahlGenThr)/float(anzahlGenTriplett)
emissionGenTrp = float(anzahlGenTrp)/float(anzahlGenTriplett)
emissionGenTyr = float(anzahlGenTyr)/float(anzahlGenTriplett)
emissionGenVal = float(anzahlGenVal)/float(anzahlGenTriplett)
emissionGenStopp = float(anzahlGenStopp)/float(anzahlGenTriplett)

print '%s\t%s' %('emissionGenMet:', emissionGenMet)
print '%s\t%s' %('emissionGenAla:', emissionGenAla)
print '%s\t%s' %('emissionGenArg:', emissionGenArg)
print '%s\t%s' %('emissionGenAsn:', emissionGenAsn)
print '%s\t%s' %('emissionGenAsp:', emissionGenAsp)
print '%s\t%s' %('emissionGenCys:', emissionGenCys)
print '%s\t%s' %('emissionGenGlu:', emissionGenGlu)
print '%s\t%s' %('emissionGenGln:', emissionGenGln)
print '%s\t%s' %('emissionGenGly:', emissionGenGly)
print '%s\t%s' %('emissionGenHis:', emissionGenHis)
print '%s\t%s' %('emissionGenIle:', emissionGenIle)
print '%s\t%s' %('emissionGenLeu:', emissionGenLeu)
print '%s\t%s' %('emissionGenLys:', emissionGenLys)
print '%s\t%s' %('emissionGenPhe:', emissionGenPhe)
print '%s\t%s' %('emissionGenPro:', emissionGenPro)
print '%s\t%s' %('emissionGenSer:', emissionGenSer)
print '%s\t%s' %('emissionGenThr:', emissionGenThr)
print '%s\t%s' %('emissionGenTrp:', emissionGenTrp)
print '%s\t%s' %('emissionGenTyr:', emissionGenTyr)
print '%s\t%s' %('emissionGenVal:', emissionGenVal)
print '%s\t%s\n' %('emissionGenStopp:', emissionGenStopp)




anzahlNichtGenMet = 0
anzahlNichtGenAla = 0
anzahlNichtGenArg = 0
anzahlNichtGenAsn = 0
anzahlNichtGenAsp = 0
anzahlNichtGenCys = 0
anzahlNichtGenGlu = 0
anzahlNichtGenGln = 0
anzahlNichtGenGly = 0
anzahlNichtGenHis = 0
anzahlNichtGenIle = 0
anzahlNichtGenLeu = 0
anzahlNichtGenLys = 0
anzahlNichtGenPhe = 0
anzahlNichtGenPro = 0
anzahlNichtGenSer = 0
anzahlNichtGenThr = 0
anzahlNichtGenTrp = 0
anzahlNichtGenTyr = 0
anzahlNichtGenVal = 0
anzahlNichtGenStopp = 0
anzahlNichtGenTriplett = 0
#findet die Anzahl jeder Aminosaeure in einem Nicht-Gen

for i in range(0,len(list_anfang_nichtcodierende)):
	seq_now = my_seq[int(list_anfang_nichtcodierende[i]):int(list_end_nichtcodierende[i])]
	#print '%s\n' %(seq_now)
	for j in range(2, len(seq_now), 3):
		#print '%s%s%s' %(seq_now[j-2], seq_now[j-1], seq_now[j])
		if seq_now[j-2] == 'A':
			if seq_now[j-1] == 'T':
				if seq_now[j] == 'G':
					anzahlNichtGenMet += 1
					anzahlNichtGenTriplett += 1
					
		if seq_now[j-2] == 'G':
			if seq_now[j-1] == 'C':
				anzahlNichtGenAla += 1
				anzahlNichtGenTriplett += 1

		if seq_now[j-2] == 'C':
			if seq_now[j-1] == 'G':
				anzahlNichtGenArg += 1
				anzahlNichtGenTriplett += 1
		
		if seq_now[j-2] == 'A':
			if seq_now[j-1] == 'A':
				if (seq_now[j] == 'C'or seq_now[j] == 'T'):
					anzahlNichtGenAsn += 1
					anzahlNichtGenTriplett += 1
					
		if seq_now[j-2] == 'G':
			if seq_now[j-1] == 'A':
				if (seq_now[j] == 'C' or seq_now[j] == 'T'):
					anzahlNichtGenAsp += 1
					anzahlNichtGenTriplett += 1
		
		if seq_now[j-2] == 'T':
			if seq_now[j-1] == 'G':
				if (seq_now[j] == 'C' or seq_now[j] == 'T'):
					anzahlNichtGenCys += 1
					anzahlNichtGenTriplett += 1
		
		if seq_now[j-2] == 'G':
			if seq_now[j-1] == 'A':
				if (seq_now[j] == 'A' or seq_now[j] == 'G'):
					anzahlNichtGenGlu += 1
					anzahlNichtGenTriplett += 1
		
		if seq_now[j-2] == 'C':
			if seq_now[j-1] == 'A':
				if (seq_now[j] == 'A' or seq_now[j] == 'G'):
					anzahlNichtGenGln += 1
					anzahlNichtGenTriplett += 1
		
		if seq_now[j-2] == 'G':
			if seq_now[j-1] == 'G':
				anzahlNichtGenGly += 1
				anzahlNichtGenTriplett += 1
				
		if seq_now[j-2] == 'C':
			if seq_now[j-1] == 'A':
				if (seq_now[j] == 'C' or seq_now[j] == 'T'):
					anzahlNichtGenHis += 1
					anzahlNichtGenTriplett += 1
					
		if seq_now[j-2] == 'A':
			if seq_now[j-1] == 'T':
				if (seq_now[j] == 'C' or seq_now[j] == 'T' or seq_now[j] == 'A'):
					anzahlNichtGenIle += 1
					anzahlNichtGenTriplett += 1
		
		if seq_now[j-2] == 'T':
			if seq_now[j-1] == 'T':
				if (seq_now[j] == 'A' or seq_now[j] == 'G'):
					anzahlNichtGenLeu += 1
					anzahlNichtGenTriplett += 1
		
		if seq_now[j-2] == 'C':
			if seq_now[j-1] == 'T':
				anzahlNichtGenLeu += 1
				anzahlNichtGenTriplett += 1
		
		if seq_now[j-2] == 'A':
			if seq_now[j-1] == 'A':
				if (seq_now[j] == 'A' or seq_now[j] == 'G'):
					anzahlNichtGenLys += 1
					anzahlNichtGenTriplett += 1
		
		if seq_now[j-2] == 'T':
			if seq_now[j-1] == 'T':
				if (seq_now[j] == 'C' or seq_now[j] == 'T'):
					anzahlNichtGenPhe += 1
					anzahlNichtGenTriplett += 1

		if seq_now[j-2] == 'C':
			if seq_now[j-1] == 'C':
				anzahlNichtGenPro += 1
				anzahlNichtGenTriplett += 1
		
		if seq_now[j-2] == 'T':
			if seq_now[j-1] == 'C':
				anzahlNichtGenSer += 1
				anzahlNichtGenTriplett += 1

		if seq_now[j-2] == 'A':
			if seq_now[j-1] == 'G':
				if (seq_now[j] == 'C' or seq_now[j] == 'T'):
					anzahlNichtGenSer += 1
					anzahlNichtGenTriplett += 1
					
		if seq_now[j-2] == 'A':
			if seq_now[j-1] == 'C':
				anzahlNichtGenThr += 1
				anzahlNichtGenTriplett += 1
				
		if seq_now[j-2] == 'T':
			if seq_now[j-1] == 'G':
				if seq_now[j] == 'G':
					anzahlNichtGenTrp += 1
					anzahlNichtGenTriplett += 1
		
		if seq_now[j-2] == 'T':
			if seq_now[j-1] == 'A':
				if (seq_now[j] == 'C' or seq_now[j] == 'T'):
					anzahlNichtGenTyr += 1
					anzahlNichtGenTriplett += 1

		if seq_now[j-2] == 'G':
			if seq_now[j-1] == 'T':
				anzahlNichtGenVal += 1
				anzahlNichtGenTriplett += 1

		if seq_now[j-2] == 'T':
			if seq_now[j-1] == 'A':
				if (seq_now[j] == 'A' or seq_now[j] == 'G'):
					anzahlNichtGenStopp += 1
					anzahlNichtGenTriplett += 1
		if seq_now[j-2] == 'T':
			if seq_now[j-1] == 'G':
				if seq_now[j] == 'A':
					anzahlNichtGenStopp += 1
					anzahlNichtGenTriplett += 1

print '%s\t%s\n' %('anzahlNichtGenTriplett:', anzahlNichtGenTriplett)

emissionNichtGenMet = float(anzahlNichtGenMet)/float(anzahlNichtGenTriplett)
emissionNichtGenAla = float(anzahlNichtGenAla)/float(anzahlNichtGenTriplett)
emissionNichtGenArg = float(anzahlNichtGenArg)/float(anzahlNichtGenTriplett)
emissionNichtGenAsn = float(anzahlNichtGenAsn)/float(anzahlNichtGenTriplett)
emissionNichtGenAsp = float(anzahlNichtGenAsp)/float(anzahlNichtGenTriplett)
emissionNichtGenCys = float(anzahlNichtGenCys)/float(anzahlNichtGenTriplett)
emissionNichtGenGlu = float(anzahlNichtGenGlu)/float(anzahlNichtGenTriplett)
emissionNichtGenGln = float(anzahlNichtGenGln)/float(anzahlNichtGenTriplett)
emissionNichtGenGly = float(anzahlNichtGenGly)/float(anzahlNichtGenTriplett)
emissionNichtGenHis = float(anzahlNichtGenHis)/float(anzahlNichtGenTriplett)
emissionNichtGenIle = float(anzahlNichtGenIle)/float(anzahlNichtGenTriplett)
emissionNichtGenLeu = float(anzahlNichtGenLeu)/float(anzahlNichtGenTriplett)
emissionNichtGenLys = float(anzahlNichtGenLys)/float(anzahlNichtGenTriplett)
emissionNichtGenPhe = float(anzahlNichtGenPhe)/float(anzahlNichtGenTriplett)
emissionNichtGenPro = float(anzahlNichtGenPro)/float(anzahlNichtGenTriplett)
emissionNichtGenSer = float(anzahlNichtGenSer)/float(anzahlNichtGenTriplett)
emissionNichtGenThr = float(anzahlNichtGenThr)/float(anzahlNichtGenTriplett)
emissionNichtGenTrp = float(anzahlNichtGenTrp)/float(anzahlNichtGenTriplett)
emissionNichtGenTyr = float(anzahlNichtGenTyr)/float(anzahlNichtGenTriplett)
emissionNichtGenVal = float(anzahlNichtGenVal)/float(anzahlNichtGenTriplett)
emissionNichtGenStopp = float(anzahlNichtGenStopp)/float(anzahlNichtGenTriplett)

print '%s\t%s' %('emissionNichtGenMet:', emissionNichtGenMet)
print '%s\t%s' %('emissionNichtGenAla:', emissionNichtGenAla)
print '%s\t%s' %('emissionNichtGenArg:', emissionNichtGenArg)
print '%s\t%s' %('emissionNichtGenAsn:', emissionNichtGenAsn)
print '%s\t%s' %('emissionNichtGenAsp:', emissionNichtGenAsp)
print '%s\t%s' %('emissionNichtGenCys:', emissionNichtGenCys)
print '%s\t%s' %('emissionNichtGenGlu:', emissionNichtGenGlu)
print '%s\t%s' %('emissionNichtGenGln:', emissionNichtGenGln)
print '%s\t%s' %('emissionNichtGenGly:', emissionNichtGenGly)
print '%s\t%s' %('emissionNichtGenHis:', emissionNichtGenHis)
print '%s\t%s' %('emissionNichtGenIle:', emissionNichtGenIle)
print '%s\t%s' %('emissionNichtGenLeu:', emissionNichtGenLeu)
print '%s\t%s' %('emissionNichtGenLys:', emissionNichtGenLys)
print '%s\t%s' %('emissionNichtGenPhe:', emissionNichtGenPhe)
print '%s\t%s' %('emissionNichtGenPro:', emissionNichtGenPro)
print '%s\t%s' %('emissionNichtGenSer:', emissionNichtGenSer)
print '%s\t%s' %('emissionNichtGenThr:', emissionNichtGenThr)
print '%s\t%s' %('emissionNichtGenTrp:', emissionNichtGenTrp)
print '%s\t%s' %('emissionNichtGenTyr:', emissionNichtGenTyr)
print '%s\t%s' %('emissionNichtGenVal:', emissionNichtGenVal)
print '%s\t%s\n' %('emissionNichtGenStopp:', emissionNichtGenStopp)


#merging list of start positions of coding and noncoding regions and sorting them in increasing order
all_starting_pos = (list_anfang_codierende) + (list_anfang_nichtcodierende)
all_start_pos = [int(x) for x in all_starting_pos]
all_start_pos = sorted(all_start_pos)

#print all_start_pos[:10]

#converting list of start positions into int values
anfang_codierende = [int(x) for x in list_anfang_codierende]
anfang_nichtcodierende = [int(x) for x in list_anfang_nichtcodierende]

#print anfang_codierende

genZuGen = 0
genZuNichtGen = 0
nichtGenZuGen = 0
nichtGenZuNichtGen = 0
startZuGen = 0
startZuNichtGen = 0

#Bestimmung der Anzahl der Transitionen
if anfang_codierende[0] == 1:
	startZuGen += 1
else:
	startZuNichtGen += 1
	
for i in range(1, len(all_start_pos)):
	if all_start_pos[i] in anfang_codierende:
		if all_start_pos[i-1] in anfang_codierende:
			genZuGen += 1
		if all_start_pos[i-1] in anfang_nichtcodierende:
			nichtGenZuGen += 1
	
	if all_start_pos[i] in anfang_nichtcodierende:
		if all_start_pos[i-1] in anfang_codierende:
			genZuNichtGen += 1
		if all_start_pos[i-1] in anfang_nichtcodierende:
			nichtGenZuNichtGen += 1

print '%s\t%s' %('startZuGen:', startZuGen)
print '%s\t%s' %('startZuNichtGen:', startZuNichtGen)
print '%s\t%s' %('genZuGen:', genZuGen)
print '%s\t%s' %('genZuNichtGen:', genZuNichtGen)
print '%s\t%s' %('nichtGenZuGen:', nichtGenZuGen)
print '%s\t%s\n' %('nichtGenZuNichtGen:', nichtGenZuNichtGen)

#Werte der Transition
trans_genZuGen = float(genZuGen)/float(genZuGen + genZuNichtGen)
trans_genZuNichtGen = float(genZuNichtGen)/float(genZuGen + genZuNichtGen)
trans_nichtGenZuGen = float(nichtGenZuGen)/float(nichtGenZuGen + nichtGenZuNichtGen)
trans_nichtGenZuNichtGen = float(nichtGenZuNichtGen)/float(nichtGenZuGen + nichtGenZuNichtGen)

print '%s\t%s' %('trans_genZuGen:', trans_genZuGen)
print '%s\t%s' %('trans_genZuNichtGen:', trans_genZuNichtGen)
print '%s\t%s' %('trans_nichtGenZuGen:', trans_nichtGenZuGen)
print '%s\t%s\n' %('trans_nichtGenZuNichtGen:', trans_nichtGenZuNichtGen)

class ObservationTypeAlphabet(Alphabet.Alphabet):
	letters = ['emissionMet', 'emissionAla', 'emissionArg', 'emissionAsn', 'emissionAsp', 'emissionCys', 'emissionGlu', 'emissionGln', 'emissionGly', 'emissionHis', 'emissionIle', 'emissionLeu', 'emissionLys', 'emissionPhe', 'emissionPro', 'emissionSer', 'emissionThr', 'emissionTrp', 'emissionTyr', 'emissionVal', 'emissionStopp']

class StateTypeAlphabet(Alphabet.Alphabet):
	letters =['G' , 'N']	#N--> NichtGen, G--> Gen
	
def _loaded_state(chance_observation, cur_state):
	if cur_state == 'G':
		if chance_observation == emissionGenMet:
			return 'emissionMet'
		elif chance_observation == emissionGenAla:
			return 'emissionAla'
		elif chance_observation == emissionGenArg:
			return 'emissionArg'
		elif chance_observation == emissionGenAsn:
			return 'emissionAsn'
		elif chance_observation == emissionGenAsp:
			return 'emissionAsp'
		elif chance_observation == emissionGenCys:
			return 'emissionCys'
		elif chance_observation == emissionGenGlu:
			return 'emissionGlu'
		elif chance_observation == emissionGenGln:
			return 'emissionGln'
		elif chance_observation == emissionGenGly:
			return 'emissionGly'
		elif chance_observation == emissionGenHis:
			return 'emissionHis'
		elif chance_observation == emissionGenIle:
			return 'emissionIle'
		elif chance_observation == emissionGenLeu:
			return 'emissionLeu'
		elif chance_observation == emissionGenLys:
			return 'emissionLys'
		elif chance_observation == emissionGenPhe:
			return 'emissionPhe'
		elif chance_observation == emissionGenPro:
			return 'emissionPro'
		elif chance_observation == emissionGenSer:
			return 'emissionSer'
		elif chance_observation == emissionGenThr:
			return 'emissionThr'
		elif chance_observation == emissionGenTrp:
			return 'emissionTrp'
		elif chance_observation == emissionGenTyr:
			return 'emissionTyr'
		elif chance_observation == emissionGenVal:
			return 'emissionVal'
		elif chance_observation == emissionStopp:
			return 'emissionStopp'
		else:
			return ValueError("Unexpected Gen_observation")
	
	elif cur_state == 'N':
		if chance_observation == emissionNichtGenMet:
			return 'emissionMet'
		elif chance_observation == emissionNichtGenAla:
			return 'emissionAla'
		elif chance_observation == emissionNichtGenArg:
			return 'emissionArg'
		elif chance_observation == emissionNichtGenAsn:
			return 'emissionAsn'
		elif chance_observation == emissionNichtGenAsp:
			return 'emissionAsp'
		elif chance_observation == emissionNichtGenCys:
			return 'emissionCys'
		elif chance_observation == emissionNichtGenGlu:
			return 'emissionGlu'
		elif chance_observation == emissionNichtGenGln:
			return 'emissionGln'
		elif chance_observation == emissionNichtGenGly:
			return 'emissionGly'
		elif chance_observation == emissionNichtGenHis:
			return 'emissionHis'
		elif chance_observation == emissionNichtGenIle:
			return 'emissionIle'
		elif chance_observation == emissionNichtGenLeu:
			return 'emissionLeu'
		elif chance_observation == emissionNichtGenLys:
			return 'emissionLys'
		elif chance_observation == emissionNichtGenPhe:
			return 'emissionPhe'
		elif chance_observation == emissionNichtGenPro:
			return 'emissionPro'
		elif chance_observation == emissionNichtGenSer:
			return 'emissionSer'
		elif chance_observation == emissionNichtGenThr:
			return 'emissionThr'
		elif chance_observation == emissionNichtGenTrp:
			return 'emissionTrp'
		elif chance_observation == emissionNichtGenTyr:
			return 'emissionTyr'
		elif chance_observation == emissionNichtGenVal:
			return 'emissionVal'
		elif chance_observation == emissionNichtGenStopp:
			return 'emissionStopp'
		else:
			return ValueError("Unexpected NichtGen_observation")
	
	else:
		return ValueError("Unexpected cur_state %s" % cur_state)



def generate_rolls(num_rolls):

    # start off in the Intron region
	
    cur_state = 'N'
    roll_seq = MutableSeq('', ObservationTypeAlphabet())
    state_seq = MutableSeq('', StateTypeAlphabet())
#    return roll_seq.toseq(), state_seq.toseq()

	
    # generate the sequence
    for roll in range(num_rolls):
        state_seq.append(cur_state)
        # generate a random number
        chance_num = random.random()
        #chance_num = chr(random.randrange(0,22))
        #numbers = random.randrange(0,22)
        #print numbers
    
        # add on a new roll to the sequence
        new_roll = _loaded_state(chance_num, cur_state)
        roll_seq.append(new_roll)

        # now give us a chance to switch to a new state
    
        chance_num = random.random()
		
        if cur_state == 'N':
			cur_state = 'G'
        elif cur_state == 'G' and chance_num > trans_genZuGen:
			cur_state = 'N'
        else:
            cur_state = 'G'
    
    return roll_seq.toseq(), state_seq.toseq()	


#build a HMM	
mm_builder = MarkovModel.MarkovModelBuilder(StateTypeAlphabet(), ObservationTypeAlphabet())

mm_builder.set_initial_probabilities({'G': startZuGen, 'N' : startZuNichtGen})
mm_builder.allow_all_transitions()

mm_builder.set_transition_score('G', 'G', trans_genZuGen)
mm_builder.set_transition_score('G', 'N', trans_genZuNichtGen)
mm_builder.set_transition_score('N', 'G', trans_nichtGenZuGen)
mm_builder.set_transition_score('N', 'N', trans_nichtGenZuNichtGen)

mm_builder.set_emission_score('G', 'emissionMet', emissionGenMet)
mm_builder.set_emission_score('G', 'emissionAla', emissionGenAla)
mm_builder.set_emission_score('G', 'emissionArg', emissionGenArg)
mm_builder.set_emission_score('G', 'emissionAsn', emissionGenAsn)
mm_builder.set_emission_score('G', 'emissionAsp', emissionGenAsp)
mm_builder.set_emission_score('G', 'emissionCys', emissionGenCys)
mm_builder.set_emission_score('G', 'emissionGlu', emissionGenGlu)
mm_builder.set_emission_score('G', 'emissionGln', emissionGenGln)
mm_builder.set_emission_score('G', 'emissionGly', emissionGenGly)
mm_builder.set_emission_score('G', 'emissionHis', emissionGenHis)
mm_builder.set_emission_score('G', 'emissionIle', emissionGenIle)
mm_builder.set_emission_score('G', 'emissionLeu', emissionGenLeu)
mm_builder.set_emission_score('G', 'emissionLys', emissionGenLys)
mm_builder.set_emission_score('G', 'emissionPhe', emissionGenPhe)
mm_builder.set_emission_score('G', 'emissionPro', emissionGenPro)
mm_builder.set_emission_score('G', 'emissionSer', emissionGenSer)
mm_builder.set_emission_score('G', 'emissionThr', emissionGenThr)
mm_builder.set_emission_score('G', 'emissionTrp', emissionGenTrp)
mm_builder.set_emission_score('G', 'emissionTyr', emissionGenTyr)
mm_builder.set_emission_score('G', 'emissionVal', emissionGenVal)
mm_builder.set_emission_score('G', 'emissionStopp', emissionGenStopp)

mm_builder.set_emission_score('N', 'emissionMet', emissionNichtGenMet)
mm_builder.set_emission_score('N', 'emissionAla', emissionNichtGenAla)
mm_builder.set_emission_score('N', 'emissionArg', emissionNichtGenArg)
mm_builder.set_emission_score('N', 'emissionAsn', emissionNichtGenAsn)
mm_builder.set_emission_score('N', 'emissionAsp', emissionNichtGenAsp)
mm_builder.set_emission_score('N', 'emissionCys', emissionNichtGenCys)
mm_builder.set_emission_score('N', 'emissionGlu', emissionNichtGenGlu)
mm_builder.set_emission_score('N', 'emissionGln', emissionNichtGenGln)
mm_builder.set_emission_score('N', 'emissionGly', emissionNichtGenGly)
mm_builder.set_emission_score('N', 'emissionHis', emissionNichtGenHis)
mm_builder.set_emission_score('N', 'emissionIle', emissionNichtGenIle)
mm_builder.set_emission_score('N', 'emissionLeu', emissionNichtGenLeu)
mm_builder.set_emission_score('N', 'emissionLys', emissionNichtGenLys)
mm_builder.set_emission_score('N', 'emissionPhe', emissionNichtGenPhe)
mm_builder.set_emission_score('N', 'emissionPro', emissionNichtGenPro)
mm_builder.set_emission_score('N', 'emissionSer', emissionNichtGenSer)
mm_builder.set_emission_score('N', 'emissionThr', emissionNichtGenThr)
mm_builder.set_emission_score('N', 'emissionTrp', emissionNichtGenTrp)
mm_builder.set_emission_score('N', 'emissionTyr', emissionNichtGenTyr)
mm_builder.set_emission_score('N', 'emissionVal', emissionNichtGenVal)
mm_builder.set_emission_score('N', 'emissionStopp', emissionNichtGenStopp)

#get_markov_model()

baum_welch_mm = mm_builder.get_markov_model()
standard_mm = mm_builder.get_markov_model()

rolls, states = generate_rolls(3000)

#predicted_states, prob = my_mm.viterbi(rolls, StateTypeAlphabet())
#print("prob: %f" % prob)