library(pwr)
library(effectsize)

#
# checking a couple of examples from the paper
#
result = pwr.2p2n.test(h=4.705882353, n1=10, n2=10, sig.level=0.05, power=NULL)
result$power # one reported

result = pwr.2p2n.test(h=1.176470588, n1=10, n2=10, sig.level=0.00489, power=NULL)
result$power # 0.427193520923755

#
# our scenario, to apply on the CV
#
result = pwr.2p.test(h=NULL, n=4, sig.level=0.05, power=0.8)
result$h

#
# reproduce figure
#
# apparently h = effect size = (foldchange - 1) / CV
cv = 0.4
replicates = c() 
powersA = c()
powersB = c()
powersC = c()
powersD = c()

for (i in 2:12){
  replicates = append(replicates, i)
  
  fc = 1.5
  effect_size = (fc-1)/cv
  result = pwr.2p.test(h=effect_size, n=i, sig.level=0.05, power=NULL)
  powersA = append(powersA, result$power)
  
  fc = 2
  effect_size = (fc-1)/cv
  result = pwr.2p.test(h=effect_size, n=i, sig.level=0.05, power=NULL)
  powersB = append(powersB, result$power)
  
  fc = 3
  effect_size = (fc-1)/cv
  result = pwr.2p.test(h=effect_size, n=i, sig.level=0.05, power=NULL)
  powersC = append(powersC, result$power)
  
  fc = 4
  effect_size = (fc-1)/cv
  result = pwr.2p.test(h=effect_size, n=i, sig.level=0.05, power=NULL)
  powersD = append(powersD, result$power)
}
plot(replicates, powersA, type='b', pch=19, col='black', ylim=c(0,1))
lines(replicates, powersB, type='b', pch=19, col='red')
lines(replicates, powersC, type='b', pch=19, col='blue')
lines(replicates, powersD, type='b', pch=19, col='green')

#
# computation of appropriate d
#

# n1 = 2 and n2 = 4 for four hours
pwr.t2n.test(n1=2, n2=4, sig.level=0.05, power=0.80, d=NULL) # d = 3.257188

# n1 = 4 and n2 = 4, 3, 2 for twenty-four hours
pwr.t2n.test(n1=4, n2=4, sig.level=0.05, power=0.80, d=NULL) # d = 2.380757
pwr.t2n.test(n1=4, n2=3, sig.level=0.05, power=0.80, d=NULL) # d = 2.683785
pwr.t2n.test(n1=4, n2=2, sig.level=0.05, power=0.80, d=NULL) # d = 3.257188

#
# how to calculate effect size: Cohen's d
#

#define two samples
data1 <- c(6, 6, 7, 8, 8, 10, 11, 13, 15, 15, 16, 17, 19, 19, 21)
data2 <- c(10, 11, 13, 13, 15, 17, 17, 19, 20, 22, 24, 25, 27, 29, 29)

cohens_d(data1, data2)

a = (length(data1)-1) * var(data1)
b = (length(data2)-1) * var(data2)
num = a + b
denom = length(data1) + length(data2) - 2
pooled = sqrt(num / denom)
(mean(data2) - mean(data1)) / pooled

