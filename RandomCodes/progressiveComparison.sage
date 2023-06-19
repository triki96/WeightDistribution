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

w = 13
p = 2 #p value for Stern's parameters
ell = 6 #l value for Stern's parameters

max_attempt = 10000
oldEstimatorList = [0 for i in range(max_attempt)]
estimatorOneList = [0 for i in range(max_attempt)]
theoreticalList = [0 for i in range(max_attempt)]
estimatorTwoList = [0 for i in range(max_attempt)]

avg_hat_N = 0 #average number of solutions found
yes_Sol = 0  #number of positive answer from stern_isd() algorithm
foundCodewords = Set([]) #Different codewords found
p_star = 0
emp_Nw_old = 0

succ_pr = 1.*binomial(floor((k+ell)/2), p)*binomial(ceil((k+ell)/2), p)*binomial(n-k-ell,w-2*p)/binomial(n,w)

for num_attempt in range(1,max_attempt):	
	X = stern_isd(H, Fq, Fq_star, n, k, p, ell, w) #Apply permutation and PGE

	# Old estimator
	if (len(X)>0):
		yes_Sol+=1
	p_star = yes_Sol*1./num_attempt
	emp_Nw_old = log(1-p_star) / (log(1-succ_pr*1.))

	#Estimator 1
	avg_hat_N+=len(X)
	emp_Nw = avg_hat_N/(num_attempt*succ_pr)
	
	#Estimator 2
	for u in X: 
		foundCodewords = foundCodewords.union(Set([u.str()]))
	p_star2 = (1 - (1-succ_pr)^num_attempt)
	emp_Nw_stimatore2 = foundCodewords.cardinality() / p_star2

	print("n:",n,"  w:",w,"  q:",q,"  k:",k)
	print(str(num_attempt) + " of " + str(max_attempt))
	print(" Old algorithm estimate: ", max(0,log(emp_Nw_old*1.)))
	print(" Estimator I estimate:   ",max(0,log(emp_Nw*1.)))
	print(" Estimator 2 estimate:   ",max(0,log(emp_Nw_stimatore2*1.)))
	print(" Expected:		 ", max(0,log(1.*binomial(n,w)*(q-1^w)/(q^(n-k)))))
	print("----------------------------------------------------------")
	estimatorOneList[num_attempt] = max(0,log(emp_Nw*1.))
	oldEstimatorList[num_attempt] = max(0,log(emp_Nw_old*1.))
	estimatorTwoList[num_attempt] = max(0,log(emp_Nw_stimatore2*1.))
	theoreticalList[num_attempt] =  log(1.*binomial(n,w)*(q-1^w)/(q^(n-k)))


fileName1 = "oldEstimatorList.dat"
with open(fileName1, "w") as f:
	f.write("x y\n")
	for i in range(len(oldEstimatorList)):
		f.write(str(i) + " " + str(oldEstimatorList[i]) + "\n")
f.close()

fileName2 = "estimatorOneList.dat"
with open(fileName2, "w") as f:
	f.write("x y\n")
	for i in range(len(estimatorOneList)):
		f.write(str(i) + " " + str(estimatorOneList[i]) + "\n")
f.close()

fileName4 = "estimatorTwoList.dat"
with open(fileName4, "w") as f:
	f.write("x y\n")
	for i in range(len(estimatorTwoList)):
		f.write(str(i) + " " + str(estimatorTwoList[i]) + "\n")
f.close()

fileName3 = "theoreticalList.dat"
with open(fileName3, "w") as f:
	f.write("x y\n")
	for i in range(len(theoreticalList)):
		f.write(str(i) + " " + str(theoreticalList[i]) + "\n")
f.close()

