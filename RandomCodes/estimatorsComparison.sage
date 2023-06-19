reset()
load_attach_path('../utils/')
load('list_sorting.sage')
load('MacKay_utils.sage')
load('stern_utils.sage')

n = 100 #code length
k = 50 #code dimension
q = 2 #finite field size
Fq = GF(q) #finite field
Fq_star = [x for x in Fq.list()[1:]]
H = random_matrix(Fq,n-k,n)

max_attempt = 20 # number of calls to ISD for each value of w
w_min = 10	#min value of w tested
w_max = 30 	#max value of w tested

p = 2 #p value in Stern's parameters
ell = 6 #l value in Stern's parameters

oldEstimatorList = [0 for i in range(n)]
estimatorOneList = [0 for i in range(n)]
theoreticalList = [0 for i in range(n)]

for w in range(w_min,w_max):
	succ_pr = binomial(floor((k+ell)/2), p)*binomial(ceil((k+ell)/2), p)*binomial(n-k-ell,w-2*p)/binomial(n,w)
	avg_hat_N = 0 #average number of solutions found
	yes_Sol = 0  #number of positive answer from stern_isd() algorithm
	p_star = 0
	emp_Nw_old = 0

	for num_attempt in range(1,max_attempt):
		X = stern_isd(H, Fq, Fq_star, n, k, p, ell, w) #Apply permutation and PGE
		
		#Old estimator
		if (len(X)>0):
			yes_Sol+=1
		p_star = yes_Sol*1./num_attempt
		emp_Nw_old = log(1-p_star) / (log(1-succ_pr*1.))

		#Estimator I
		avg_hat_N+=len(X)
		emp_Nw = avg_hat_N/(num_attempt*succ_pr)
	
	print("w: ", w)
	print(" Old Estimator: ", max(0,log(emp_Nw_old*1.)))
	print(" Estimator I: ",max(0,log(emp_Nw*1.)))
	print(" Expected:		 ", max(0,log(1.*binomial(n,w)*(q-1^w)/(q^(n-k)))))
	
	estimatorOneList[w] = max(0,log(emp_Nw*1.))
	oldEstimatorList[w] = max(0,log(emp_Nw_old*1.))
	theoreticalList[w] =  max(0,log(1.*binomial(n,w)*(q-1^w)/(q^(n-k))))

	fileName1 = "oldEstimatorList.dat"
	with open(fileName1, "w") as f:
		f.write("x y\n")
		for w in range(len(oldEstimatorList)):
			f.write(str(w) + " " + str(oldEstimatorList[w]) + "\n")
	f.close()

	fileName2 = "estimatorOneList.dat"
	with open(fileName2, "w") as f:
		f.write("x y\n")
		for w in range(len(estimatorOneList)):
			f.write(str(w) + " " + str(estimatorOneList[w]) + "\n")
	f.close()

	fileName3 = "theoreticalList.dat"
	with open(fileName3, "w") as f:
		f.write("x y\n")
		for w in range(len(theoreticalList)):
			f.write(str(w) + " " + str(theoreticalList[w]) + "\n")
	f.close()