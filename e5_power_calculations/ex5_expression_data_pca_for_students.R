## Power calculations
## Example one: two independent groups (independent t-test)
# A proposed study wishes to investigate the effects of a hypertensive drug on patients with a genetic variant (experimental group) 
# compared to a patients without (control group). Previous studies show that the minimum clinically important difference is 15 mmHg 
# and a pilot study the pooled standard deviation (SD) is 20 mmHg and measurements were normally distributed. If we match the
# experimental group 1 to 1 with the control group, what sample size will we need for 80% power and alpha level p=0.05.

install.packages("pwr") # Useful website https://www.statmethods.net/stats/power.html
library(pwr)

eff = 15 / 20
sig = 0.05
pow = 0.8

print(pwr.t.test(type = "two.sample", d = eff, sig.level = sig, power = pow))

# After a further pilot study, we find the pooled standard deviation was 10mmHg


## Example two: a post-hoc power calculation for two paired groups (paired t-test)
# Is  there  evidence  that  clofibrate  changes  the  mean cholesterol  level?  Cholesterol  is  to  be  measured 
# before and after receiving clofibrate. From previous studies,  a  mean  difference  of  40mg/dl  is  deemed 
# clinically significant, with pooled standard deviation of 50 with measurements normally distributed. What  sample 
# size  of  is  required for  alpha level p = 0.05 and  power 80%? 

# How many samples would we need for alpha level p = 0.01

# As a result of the poor response/dropout rate, only 12 patients were recruited. A mean difference of 50 and pooled 
# SD of 60  mg/dl were found. What is the power for this study?




## For single cell data
# You wish to do a single cell project on rat brains, you wish to check how much power you have when using 12 samples, aiming
# for 10k cells per sample, with a read depth of 100k reads. The cell type you are most interested in is expected to be around 
# 20% of the the otherall cells. For this project you are looking at differential expression.
install.packages("devtools")
devtools::install_github("hadley/devtools")
devtools::install_github("heiniglab/scPower") # Useful website https://scpower.helmholtz-muenchen.de/
library(scPower)

# May need 'limma' as well
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("limma")

# Looking at Aim 1
a1_power<-power.general.withDoublets(nSamples=12,nCells=10000,readDepth=100000, ct.freq=0.2,type="de",
                                     #> Loading required package: pwr
                                     ref.study=scPower::de.ref.study, ref.study.name="Blueprint (CLL) iCLL-mCLL", samplesPerLane=1,
                                     read.umi.fit = scPower::read.umi.fit[
                                       read.umi.fit$type=="10X_PBMC_1",], gamma.mixed.fits = scPower::gamma.mixed.fits, ct="CD14+ Monocytes", disp.fun.param=scPower::disp.fun.param, mappingEfficiency = 0.8,
                                     min.UMI.counts = 3,
                                     perc.indiv.expr = 0.5,
                                     sign.threshold = 0.05,
                                     MTmethod="Bonferroni")
print(a1_power)

# Looking at Aim 2/3
# As a follow up, you are going to look at 12 samples again, but in this case enriched for you cell type of interest, so at 
# least 80%, though with fewer cells per sample: 2500 and lower read depth of 50k reads. What power do you have for differential
# expression?
