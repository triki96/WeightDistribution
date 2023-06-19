# coding=utf-8
#Vediamo come l'algoritmo che ha pensato Paolo migliora quello dei cinesi
reset()
load_attach_path('./utils/')
load('list_sorting.sage')
load('MacKay_utils.sage')
load('stern_utils.sage')

q = 2
Fq = GF(q)
Fq_star = [Fq(i) for i in range(1,q)];
inputFile = '495.62.3.2915'
#inputFile = '96.3.963'
n,m,H = readFromFile(inputFile)
k = n-m
w = 4
p = 2
ell = 20


# Set the confidence interval and the number of round
succ_pr = binomial(floor((k+ell)/2), p)*binomial(ceil((k+ell)/2), p)*binomial(n-k-ell,w-2*p)/binomial(n,w)
t = 10000
p_star = 0
emp_Nw_old = 0
fileName = "./outputEstimator2/" + inputFile + "-" + str(w)
with open(fileName, "w") as f:
		f.write("Parameters" + "\n")
		f.write("Code ID: " + str(inputFile) + "\n")
		f.write("p: " + str(p) + " - k: " + str(k) + " - l: " + str(ell) + "\n")
f.close()

foundCodewords = Set([])
for num_attempt in range(1,t+1):
	with open(fileName, "a") as f:
		f.write(str(num_attempt) + "\n")
		X = stern_isd(H, Fq, Fq_star, n, k, p, ell, w) #Apply permutation and PGE
		#print("Cardinalità soluzioni trovate: ",len(X))
		for u in X: 
			foundCodewords = foundCodewords.union(Set([u.str()]))
		#print("Cardinalità nuovo insieme: ", foundCodewords.cardinality())
		p_star = (1 - (1-succ_pr)^num_attempt)
		emp_Nw = foundCodewords.cardinality() / p_star
		f.write(" Estimate II: " + str(emp_Nw*1.) + "\n")
		#print(" Estimate II: " + str(emp_Nw*1.))
	f.close()