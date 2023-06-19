#Vediamo come l'algoritmo che ha pensato Paolo migliora quello dei cinesi
reset()
load_attach_path('./utils/')
load('list_sorting.sage')
load('MacKay_utils.sage')
load('stern_utils.sage')

q = 2
Fq = GF(q)
Fq_star = [Fq(i) for i in range(1,q)];
inputFile = '252.252.3.252'
n,m,H = readFromFile(inputFile)
k = n-m
w = 20
p = 2
ell = 6

# Set the confidence interval and the number of round
succ_pr = binomial(floor((k+ell)/2), p)*binomial(ceil((k+ell)/2), p)*binomial(n-k-ell,w-2*p)/binomial(n,w)
#confidence_interval = 10 #number of errors we admit (percentage of N_C(w))
#accuracy = 0.9 #accuracy of the estimate
#t= ceil(3 /  (confidence_interval^2 * succ_pr) * log(1/(1-accuracy)))
#print("Estimated round to reach the desired accuracy: ", t)
flag = True
i = 0
while flag:
    print(i)
    i+=1
    if succ_pr * i >= 1:
        flag = False
        print("m = ", i)

