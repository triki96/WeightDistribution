#Controlliamo se l'algoritmo che ha pensato Paolo stima bene il numero di parole.
# Come codice di riferimento prendiamo un LDPC dal repertorio di McKay
reset()
load_attach_path('./utils/')
load('list_sorting.sage')
load('MacKay_utils.sage')
load('stern_utils.sage')

from sage.coding.reed_muller_code import QAryReedMullerCode

q = 2
Fq = GF(q)
Fq_star = [Fq(i) for i in range(1,q)];
inputFile = '495.62.3.2915'
n,m,H = readFromFile("./inputFiles/" + inputFile)
k = n-m
w = 4
p = 2
ell = 6

print(k)

# Set the confidence interval and the number of round
succ_pr = binomial(floor((k+ell)/2), p)*binomial(ceil((k+ell)/2), p)*binomial(n-k-ell,w-2*p)/binomial(n,w)
confidence_interval = 10 #number of errors we admit (percentage of N_C(w))
accuracy = 0.9 #accuracy of the estimate
t= ceil(3 /  (confidence_interval^2 * succ_pr) * log(1/(1-accuracy)))
print("t: ", t)


avg_hat_N = 0 #average number of solutions found
for num_attempt in range(t):
    print(num_attempt)
    X = stern_isd(H, Fq, Fq_star, n, k, p, ell, w) #Apply permutation and PGE
    avg_hat_N+=len(X)


Nw = 60 # Number of codeword of weight 4
emp_Nw = avg_hat_N/(t*succ_pr)
th_Nw = Nw
print("Emp. = ",emp_Nw*1.,", Th. = ", th_Nw*1.)