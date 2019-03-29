import sys
import math
import random
import scipy.stats
import numpy


PROPOSAL_VARIANCE = 0.01
SAMPLE_SIZE = 10000
MINIMUM_READ = 2
BURNIN = 1000
TRIM = 22
VERBOSE = False
PHI_START = 0.3
ALPHA = 1
BETA = 1
ALPHA_w = 600
BETA_w = 3400
w1 = 0.15
w0 = 1 - w1
mu_c = 0
sigma_c = PROPOSAL_VARIANCE
acctot = 0.0



Kn = 0.1
MEAN_DOMAIN_SPAN = 300
cache_nn = {}
def sigmoid(d,n):
	if d == 0: return(0.00)
	f = 1./(1.+math.exp(Kn*(MEAN_DOMAIN_SPAN-d)))
	if n in cache_nn: return f/cache_nn[n]
	tot = 0
	N = int(n+1)
	for i in range(1,N): tot += 1.0/(1.0+math.exp(Kn*(MEAN_DOMAIN_SPAN-i)))
	cache_nn[n] = tot
	return f/tot

def classify(n,r,p):
	q = 1-p; B = 1/math.log(1/p); u = math.log(n*q,1/p);
	pval = 1 - scipy.stats.gumbel_r.cdf(r,u,B)
	if pval < 0.05: return(1)
	else: return(0)


def regress(X,Y):
	N = len(X)
	xbar = numpy.average(X)
	ybar = numpy.average(Y)
	xybar = numpy.average([X[i]*Y[i] for i in range(N)])
	x2bar = numpy.average([X[i]*X[i] for i in range(N)])
	B = (xybar - xbar*ybar)/(x2bar - xbar*xbar)
	A0 = ybar - B*xbar

	yfit = [ A0 + B *X[i] for i in range(N)]
	yres = [Y[i] - (A0 + B *X[i]) for i in range(N)]
	var = sum([math.pow(yres[i],2) for i in range(N) ])/(N-2)
	std = math.sqrt(var)

	return(B, A0, std)

def getR1(n):
	"""Small Correction term. Defaults to 0.000016."""
	return(0.000016)

def getR2(n):
	"""Small Correction term. Defaults to 0.00006."""
	return(0.00006)

def getE1(n):
	"""Small Correction term. Defaults to 0.01."""
	return(0.01)
	
def getE2(n):
	"""Small Correction term. Defaults to 0.01."""
	return(0.01)

def getGamma():
	"""Euler-Mascheroni constant ~ 0.577215664901"""
	return(0.5772156649015328606)
	
def ExpectedRuns(n,p):
	"""ER_n =  log(1/p)(nq) + gamma/ln(1/p) -1/2 + r1(n) + E1(n) (Schilling, 1990)"""
	q = 1-p
	gamma = getGamma(); r1 = getR1(n); E1 = getE1(n)
	A = math.log(n*q,1.0/p)
	B = gamma/math.log(1.0/p)
	ER = A + B -0.5 + r1 + E1
	return( ER )

def good_orf(n, s):
	return (n >= 3 and s >= 150)

def read_IGV_file(PATH):
	#Chromosome  Start   End Feature reads   TAs
	orf_to_reads = {}
	for line in open(PATH,"r"):
		if line.startswith("#"): continue
		if line.upper().startswith("CHROMOSOME"): continue
		tmp = line.split()
		if tmp[3] in orf_to_reads:
			orf_to_reads[tmp[3]].append( ( int(tmp[1]) , int(tmp[4]) ) )
		else:
			orf_to_reads[tmp[3]] = [ ( int(tmp[1]) , int(tmp[4]) ) ]
	return(orf_to_reads)


def get_orf_data(reads, min_read, bad=set()):
	G = len([orf for orf in reads if orf not in bad]); g = 0;
	K = numpy.zeros(G); N = numpy.zeros(G); R = numpy.zeros(G);
	S = numpy.zeros(G); T = numpy.zeros(G); ORF = [];
	for orf in sorted(reads):
		if orf in bad: continue
		run = [0,0,0]; maxrun = [0,0,0]; k = 0; n = 0; s = 0;
		for site in reads[orf]:
			n += 1
			if site[1] < min_read:
				run[0] +=1
				if run[0] == 1:	run[1] = site[0]
				run[2] = site[0]
			else:
				k += 1
				maxrun = max(run,maxrun)
				run = [0,0,0]
		maxrun = max(run, maxrun)
		r = maxrun[0]
		if r > 0: s = maxrun[2] + 2  - maxrun[1]
		else: s = 0
		t = max(reads[orf])[0] + 2  - min(reads[orf])[0]
		K[g]=k; N[g]=n; R[g]=r; S[g]=s; T[g]=t;
		ORF.append(orf)
		g+=1
	return((ORF,K,N,R,S,T))

	
def F_non(p, N, R):
	q = 1.0 - p
	total = numpy.log(scipy.stats.beta.pdf(p,ALPHA,BETA))
	mu = numpy.log(N*q) / numpy.log(1/p)
	sigma = 1/math.log(1/p);
	total+= sum(numpy.log(scipy.stats.gumbel_r.pdf(R, mu, sigma)))
	return(total)
	

def sample_Z(p, w1, N, R, S, T, mu_s, sigma_s, mu_r, sigma_r):
	G = len(N); h1 = numpy.zeros(G)
	q = 1-p
	mu = numpy.log(N*q) / numpy.log(1/p)
	sigma = 1/math.log(1/p);
	h0 = ((scipy.stats.gumbel_r.pdf(R,mu,sigma)) * scipy.stats.norm.pdf(S, mu_s*R, sigma_s)  * (1-w1))
	for g in xrange(G):
		h1[g] = sigmoid(S[g], T[g]) * scipy.stats.norm.pdf(R[g], mu_r*S[g], sigma_r) * w1
	p_z1 = h1/(h0+h1)
	new_Z = scipy.stats.binom.rvs(1, p_z1, size=G)
	return(new_Z)	


def useMsg():
	print """
Usage: python MH.py -f <FILE> <Optional flags>"
	-f <text> :=\tIGV-Formatted File With Reads. Required.
	-s <int> :=\tSample Size. Optional. Default: -s 1000
	-b <int> :=\tBurn-in Size. Optional. Default: -b 100
	-t <int> :=\tTrim number. Optional. Default: -t 1
	-m <int> :=\tMinimum read count to consider as insertion. Optional. Default: -m 2
	-v :=\t\tVerbose. Prints sample of parameters. Optional. Default: NOT verbose 
	-p <float> :=\tPhi-0 starting value. Optional. Default: -p 0.5
	"""
	sys.exit()
