import sys
from MH_tools import *


try:
	if "-f" in sys.argv:
		PATH = sys.argv[sys.argv.index("-f")+1]
	else:
		raise Exception("Please Provide an Input File.")
	if "-l" in sys.argv:
		LANES = [int(x) for x in sys.argv[sys.argv.index("-l")+1].split(",")]
	if "-m" in sys.argv:
		MINIMUM_READ = int(sys.argv[sys.argv.index("-m")+1])
	if "-s" in sys.argv:
		SAMPLE_SIZE = int(sys.argv[sys.argv.index("-s")+1])
	if "-p" in sys.argv:
		PHI_START = float(sys.argv[sys.argv.index("-p")+1])	
	if "-b" in sys.argv:
		BURNIN = int(sys.argv[sys.argv.index("-b")+1])
	if "-t" in sys.argv:
		TRIM = int(sys.argv[sys.argv.index("-t")+1])
	if "-v" in sys.argv:
		VERBOSE = True
except:
	print ""
	e = sys.exc_info()[1]
	print 'Error Parsing Input Arguments (Msg: "%s")' % e,
	useMsg()
	sys.exit()

orf_to_reads = read_IGV_file(PATH)
(ORF_all, K_all, N_all, R_all, S_all, T_all) = get_orf_data(orf_to_reads, MINIMUM_READ)
bad_orf_set = set([ORF_all[g] for g in xrange(len(N_all)) if not good_orf(N_all[g], T_all[g])]);
bad_orf_set.add("Rvnr01");
(ORF, K, N, R, S, T) = get_orf_data(orf_to_reads, MINIMUM_READ, bad_orf_set)

mu_s, temp, sigma_s = regress(R,S) # Linear regression to estimate mu_s, sigma_s for span data
mu_r, temp, sigma_r = regress(S, R) # Linear regression to estimate mu_r, sigma_r for run data


N_GENES = len(N)
Z_sample = numpy.zeros((N_GENES, SAMPLE_SIZE))
Z = [classify(N[g], R[g], 0.5)   for g in xrange(N_GENES)]
Z_sample[:,0] = Z
N_ESS = numpy.sum(Z_sample[:,0] == 1)

phi_sample = numpy.zeros(SAMPLE_SIZE) #[]
phi_sample[0] = PHI_START
phi_old = PHI_START
phi_new = 0.00



i = 1; count = 0;
while i < SAMPLE_SIZE:

	# PHI
	acc = 1.0
	phi_new  = phi_old + random.gauss(mu_c, sigma_c)
	i0 = Z_sample[:,i-1] == 0
	if phi_new > 1 or phi_new <= 0 or (F_non(phi_new, N[i0], R[i0]) - F_non(phi_old, N[i0], R[i0])) < math.log(random.uniform(0,1)):
		phi_new = phi_old
		acc = 0.0
		flag = 0

	# Z
	Z = sample_Z(phi_new, w1, N, R, S, T, mu_s, sigma_s, mu_r, sigma_r)

	# w1
	N_ESS = sum(Z == 1)
	w1 = scipy.stats.beta.rvs(N_ESS + ALPHA_w, N_GENES - N_ESS + BETA_w)
	
	count +=1
	acctot+=acc

	if (count > BURNIN) and (count % TRIM == 0):
		phi_sample[i] = phi_new
		Z_sample[:,i] = Z
		i+=1

	
print "#Command used:\tpython %s" % (" ".join(sys.argv))
print "#MH Acceptance-Rate:\t%2.2f%%" % (100.0*acctot/count)
print "#Total Iterations Performed:\t%d" % count
print "#Sample Size:\t%d" % i
print "#phi estimate:\t%f" % numpy.average(phi_sample)
if VERBOSE: print "#Orf\tk\tn\tr\ts\tzbar\tSample"
else: print "#Orf\tk\tn\tr\ts\tzbar"
i = 0
for g in xrange(len(ORF_all)):
	k = K_all[g]; n = N_all[g];  r = R_all[g]; s = S_all[g]; orf=ORF_all[g];
	if VERBOSE:
		if ORF[g] not in bad_orf_set:			
			print "%s\t%d\t%d\t%d\t%d\t%f\t%s" % (orf, k, n, r, s, numpy.average(Z_sample[i,:]), ",".join(["%d" % x for x in Z_sample[i,:]]) )
			i+=1
		else:
			print "%s\t%d\t%d\t%d\t%d\t%f\t%d" % (orf, k, n, r, s, -1, -1)
	else:
		if orf not in bad_orf_set:
			print "%s\t%d\t%d\t%d\t%d\t%f" % (orf, k, n, r, s, numpy.average(Z_sample[i,:]))
			i+=1
		else:
			print "%s\t%d\t%d\t%d\t%d\t%f" % (orf, k, n, r, s, -1)
			
