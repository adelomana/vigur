library(pwr)

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
cv = 0.6
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

replicates
powersA
powersB

cv = 1.2
sigma = sqrt(log(cv**2 + 1))
sigma

a = 1; b=1.2
effect = abs(log(a) - log(b))
effect

#
# computation of appropriate d
#

# n1 = 2 and n2 = 4 for four hours
pwr.t2n.test(n1=2, n2=4, sig.level=0.05, power=0.80, d=NULL)
# d = 3.257188

# n1 = 4 and n2 = 4, 3, 2 for twenty-four hours
pwr.t2n.test(n1=4, n2=4, sig.level=0.05, power=0.80, d=NULL) # d = 2.380757
pwr.t2n.test(n1=4, n2=3, sig.level=0.05, power=0.80, d=NULL) # d = 2.683785
pwr.t2n.test(n1=4, n2=2, sig.level=0.05, power=0.80, d=NULL) # d = 3.257188





#poooled sd https://resources.wolframcloud.com/FormulaRepository/resources/Pooled-Standard-Deviation

