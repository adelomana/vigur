library(pwr)

# read-in table with data (see comment above for structure details).
data = read.table('/Users/adrian/scratch/pmic_201100033_sm_supplinfo/Supplemental_data_2_Example_data.txt', sep="\t",header=T);

# If you wish to change the significance level, change it here:
significance=0.05

power <- array(dim=c(length(data[,1]),6));
colNames=c("Peptide","Protein","Effect Size","Power","P value","Fold Change");
colnames(power)=colNames;
power[,5]=data[,3];
power[,6]=data[,4];
data[,4]=abs(data[,4])-1;

power[,1]=as.character(data[,1]);
power[,2]=as.character(data[,2]);
power[,3]=data[,4]/(data[,5]);

for (i in 1:length(data[,1])){
  print(as.numeric(power[i,3]))
  result <- pwr.2p2n.test(h=as.numeric(power[i,3]), sig.level=data[i,3], power=NULL, n1=data[i,6], n2=data[i,7]);
  power[i,4] <- result$power;
}
# The output table will be written to the path entered below.
# The columns are: 1 peptide name, 2 protein name, 3 effect size, 4 power.

write.table(power, "/Users/adrian/scratch/pmic_201100033_sm_supplinfo/new.txt", sep="\t");
