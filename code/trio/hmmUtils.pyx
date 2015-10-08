#!/usr/bin/env python
#Rick Dewey 8.14.12
#Last modified 8.14.12 Rick Dewey for Cython
#Python utils for HMM inheritance analysis of pedigrees
#usage: Cython internal toolset for HMM 

import sys
import os
import re
import fileinput
import getopt
import argparse
import numpy as np
cimport numpy as np
import math

DTYPEf = np.float
DTYPEi = np.int
ctypedef np.float_t DTYPEf_t
ctypedef np.int_t DTYPEi_t 

#loads transition matrix
def load_transition(path, selftrans, nonselftrans, sub_ped):  

     	#trio case
	if sub_ped == "t":
		tmat = np.zeros([3, 3], np.float)
		i = 0
		tfile = open(path, "r")  
		for line in tfile.readlines():
			linelist = line.replace("nonself", str(nonselftrans)).replace("self", str(selftrans)).replace("\n", "").split("\t")
			tmat[i][0] = math.log10(np.float(linelist[0]))
			tmat[i][1] = math.log10(np.float(linelist[1]))
			tmat[i][2] = math.log10(np.float(linelist[2]))
			i+=1
		return tmat
    
	#quartet case
	elif sub_ped == "q":
		tmat = np.zeros([6, 6], np.float)
		i = 0

		tfile = open(path, "r")
		for line in tfile.readlines():
			linelist = line.replace("nonself", str(nonselftrans)).replace("self", str(selftrans)).replace("\n", "").split("\t")
			tmat[i][0] = math.log10(np.float(linelist[0]))
			tmat[i][1] = math.log10(np.float(linelist[1]))
			tmat[i][2] = math.log10(np.float(linelist[2]))
			tmat[i][3] = math.log10(np.float(linelist[3]))
			tmat[i][4] = math.log10(np.float(linelist[4]))
			tmat[i][5] = math.log10(np.float(linelist[5]))
			i+=1
		return tmat
	
	else:
		print >> sys.stderr, "Error in hmmUtils.load_transition: Unknown sub-pedigree format: "+sub_ped
		exit(1)
        
#loads emission matrix for hmm
def load_emission(path, comp, mie, err, sub_ped):
    
	#trio case
	if sub_ped == "t":
		emat = np.zeros([3, 27], np.float)

		efile = open(path, "r")
		trash = efile.readline()

		#alleles are in father, mother, child order
		#there are 12 allele assortments consistent with MIEs and one uniformly heterozygous position indicating possible compressions
		erprobG = np.float(err)/np.float(12)
		goodprob = (1-np.float(err))/np.float(15)

		erprobC = (1-np.float(comp))/np.float(26)

		erprobM = (1-np.float(mie))/np.float(15)
		mieprob = np.float(mie)/np.float(12)
                        
		i = 0
		#0 is MIErich, 1 is compression, 2 is good data
		for line in efile.readlines():
			linelist = line.replace("erprobM", str(erprobM)).replace("erprobG", str(erprobG)).replace("goodprob", str(goodprob)).replace("comprob", str(comp)).replace("erprobC", str(erprobC)).replace("mieprob", str(mieprob)).split("\t")
			emat[0][i] = math.log10(np.float(linelist[2]))
			emat[1][i] = math.log10(np.float(linelist[3]))
			emat[2][i] = math.log10(np.float(linelist[4]))
			i+=1
		return emat

	#quartet case
	elif sub_ped == "q":
		efile = open(path, "r")
		trash = efile.readline()
    
		emat = np.zeros([6, 81], np.float)

		#alleles are in father, mother, child, child order
		#there are 52 allele assortments consistent with MIEs and one uniformly heterozygous position indicating possible compressions
		erprobHI = np.float(err)/np.float(65)
		erprobID = np.float(err)/np.float(66)
		goodprobID = (1-np.float(err))/np.float(16)
		goodprobHI = (1-np.float(err))/np.float(15)

		erprobC = (1-np.float(comp))/np.float(80)

		erprobM = (1-np.float(mie))/np.float(29)
        
		mieprob = np.float(mie)/np.float(52)
                        
		i = 0
		#0 is MIErich, 1 is compression, 2 is good data
		for line in efile.readlines():
			linelist = line.replace("erprobM", str(erprobM)).replace("erprobHI", str(erprobHI)).replace("erprobID", str(erprobID)).replace("goodprobHI", str(goodprobHI)).replace("goodprobID", str(goodprobID)).replace("comprob", str(comp)).replace("erprobC", str(erprobC)).replace("mieprob", str(mieprob)).split("\t")
			emat[0][i] = math.log10(np.float(linelist[2]))
			emat[1][i] = math.log10(np.float(linelist[3]))
			emat[2][i] = math.log10(np.float(linelist[4]))
			emat[3][i] = math.log10(np.float(linelist[5]))
			emat[4][i] = math.log10(np.float(linelist[6]))
			emat[5][i] = math.log10(np.float(linelist[7]))
			i+=1
		return emat

	else:
		print >> sys.stderr, "Error in hmmUtils.load_emission: Unknown sub-pedigree format: "+sub_ped
		exit(1)

#creates array of observations
def makeobs_set(numsubj):
	obs = []
	n_assort = 3**numsubj
	for i in range(0, n_assort):
		obs.append(i)
	return obs

#returns map of allele assortments to allele codes for hmm, allele order father, mother, child 
def input_map(path_efile, sub_ped):
	imap = {}

	#trio case
	if sub_ped == "t":
		tempin = open(path_efile+"EmissionTrioGen.txt", "r")
	elif sub_ped == "q":
		tempin = open(path_efile+"EmissionQuartetGen.txt", "r")
	else:
		print >> sys.stderr, "Error in hmmUtils.input_map: Unknown sub-pedigree format: "+sub_ped
		exit(1)
		
	trash = tempin.readline()
	for line in tempin.readlines():
		linelist = line.split("\t")
		imap[linelist[1]] = linelist[0]
	return imap
        
#creates and returns as a tuple emission, transition, start, state, and observation matrices
def make_hmm(mpath, sprob, comp, mprob, errprob, sub_ped):
	if sub_ped == "t":
		#distribute non-self transition probability to non-self states
		nprob = (1-np.float(sprob))/np.float(2)

		#create matrices for HMM, all probabilities in log10 space to avoid underflow from np.floating point
		transmat = load_transition(mpath+"TransitionTrioGen.txt", sprob, nprob, "t")
		emitmat = load_emission(mpath+"EmissionTrioGen.txt", comp, mprob, errprob, "t")
		startp = np.array([math.log10(np.float(1)/np.float(3)), math.log10(np.float(1)/np.float(3)), math.log10(np.float(1)/np.float(3))])
		states = np.array([0, 1, 2])
		inmap = input_map(mpath, "t")
        
	elif sub_ped == "q":
		#distribute non-self transition probability to non-self states
		nprob = (1-np.float(sprob))/np.float(5)

		#create matrices for HMM, all probabilities in log10 space to avoid underflow from np.floating point
		transmat = load_transition(mpath+"TransitionQuartetGen.txt", sprob, nprob, "q")
		emitmat = load_emission(mpath+"EmissionQuartetGen.txt", comp, mprob, errprob, "q")
		startp = np.array([math.log10(np.float(1)/np.float(6)), math.log10(np.float(1)/np.float(6)), math.log10(np.float(1)/np.float(6)), math.log10(np.float(1)/np.float(6)), math.log10(np.float(1)/np.float(6)), math.log10(np.float(1)/np.float(6))])
		states = np.array([0, 1, 2, 3, 4, 5])
		inmap = input_map(mpath, "q")
	else:
		print >> sys.stderr, "Error in hmmUtils.make_hmm: Unknown sub-pedigree format: "+sub_ped
		exit(1)
	return transmat, emitmat, startp, states, inmap

#uses sums of log-transformed probabilities to avoid np.floating point underflow problems
def viterbi(np.ndarray[DTYPEi_t] obs, np.ndarray[DTYPEi_t] states, np.ndarray[DTYPEf_t] start_p, np.ndarray[DTYPEf_t, ndim = 2] trans_p, np.ndarray[DTYPEf_t, ndim = 2] emit_p):
	cdef int y, y0, t, state
	cdef double prob
	cdef np.ndarray[DTYPEf_t, ndim = 2] V = np.zeros([len(obs), len(states)], np.float)
	cdef np.ndarray[DTYPEi_t, ndim = 2] path = np.zeros([len(states), len(obs)], np.int)
 	
	# Initialize base cases (t == 0)
	for y in states:
		V[0,y] = start_p[y] + emit_p[y,int(obs[0])]
		path[y,0] = y
 	
	# Run Viterbi for t > 0
	for t in range(1,len(obs)):
	
		for y in states:
			(prob, state) = max([(V[t-1,y0] + trans_p[y0,y] + emit_p[y,obs[t]], y0) for y0 in states])
			V[t,y] = prob
			path[y,0:t] = path[state,0:t]
			path[y,t] = y
 
	(prob, state) = max([(V[len(obs) - 1,y], y) for y in states])
	return (prob, path[state])              
