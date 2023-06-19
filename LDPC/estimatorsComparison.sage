reset()
load_attach_path('./utils/')
load('list_sorting.sage')
load('MacKay_utils.sage')
load('stern_utils.sage')

q = 2
Fq = GF(q)
Fq_star = [Fq(i) for i in range(1,q)];
#inputFile = '96.3.963'
#inputFile = '252.252.3.252'
inputFile = '495.62.3.2915'
n,m,H = readFromFile("./input/" + inputFile)
k = n-m
w = 4
p = 2
ell = 6

# Compute the success probability of one iteration of Stern's ISD
succ_pr = binomial(floor((k+ell)/2), p)*binomial(ceil((k+ell)/2), p)*binomial(n-k-ell,w-2*p)/binomial(n,w)
t = 1000

## Variation: Set the confidence interval and the number of round
#confidence_interval = 10 #number of errors we admit (percentage of N_C(w))
#accuracy = 0.9 #accuracy of the estimate
#t= ceil(3 /  (confidence_interval^2 * succ_pr) * log(1/(1-accuracy)))
#print("Estimated round to reach the desired accuracy: ", t)


avg_hat_N = 0 #average number of solutions found
yes_Sol = 0 # number of positive answer from stern_isd() algorithm
p_star = 0
emp_Nw_old = 0
foundCodewords = Set([])
fileName = "./output/" + inputFile + "-" + str(w)

with open(fileName, "w") as f:
		f.write("Parameters" + "\n")
		f.write("Code ID: " + str(inputFile) + "\n")
		f.write("p: " + str(p) + " - k: " + str(k) + " - l: " + str(ell) + "\n")
f.close()

for num_attempt in range(1,t+1):
	with open(fileName, "a") as f:
			f.write(str(num_attempt) + "\n")
			X = stern_isd(H, Fq, Fq_star, n, k, p, ell, w) #Apply permutation and PGE

			# Previous estimator
			if (len(X)>0):
				print("Len(X)>0: ",len(X))
				yes_Sol+=1
			p_star = yes_Sol*1./num_attempt
			#emp_Nw_old = log(1-p_star) / (log(1-succ_pr*1.))
			emp_Nw_old = p_star / succ_pr
			f.write(" Estimate 0: " + str(emp_Nw_old*1.) + "\n")

			# Estimator 1
			avg_hat_N+=len(X)
			emp_Nw = avg_hat_N/(num_attempt*succ_pr)	
			f.write(" Estimate I: " + str(emp_Nw*1.) + "\n")

			# Estimator 2
			for u in X: 
				foundCodewords = foundCodewords.union(Set([u.str()]))
			p_star = (1 - (1-succ_pr)^num_attempt)
			emp_Nw = foundCodewords.cardinality() / p_star
			f.write(" Estimate II: " + str(emp_Nw*1.) + "\n")
	f.close()