install.packages("devtools")
library(devtools)
install_github("komorowskilab/R.ROSETTA")
library(R.ROSETTA)

# 3 --------------------------------
autcon
attributes(autcon)
dim(autcon)

autcon$decision
table(autcon$decision)



# 4 --------------------------------
?rosetta
autcon_4 = rosetta(autcon)

rules_4 = autcon_4$main
rules_4
dim(rules_4)

quality = autcon_4$quality
quality

viewRules(rules_4)

dim(rules_4[which(rules_4$pValue < 0.05),])
dim(rules_4[which(rules_4$pValue < 0.005),])
dim(rules_4[which(rules_4$pValue < 0.001),])

viewRules(rules_4[which(rules_4$pValue < 0.05),])


rules_4$cuts
# 5 --------------------------------
autcon_5 = rosetta(autcon, clroc='autism', roc=TRUE)

rules_5 = autcon_5$main

plotMeanROC(autcon_5)

viewRules(rules_5)

# 6 --------------------------------
?getFeatures

gr = getFeatures(autcon, rules_4)

setdiff(gr$autism, gr$control)
setdiff(gr$control, gr$autism)

setdiff(attributes(autcon)$names, union(gr$autism, gr$control))
setdiff(attributes(autcon)$names, intersect(gr$autism, gr$control))
setdiff(attributes(autcon)$names, gr$autism)
setdiff(attributes(autcon)$names, gr$control)

# 7 --------------------------------
?recalculateRules
new_rules = recalculateRules(autcon, rules_4[which(rules_4$pValue < 0.05),])

autism = new_rules[which(new_rules$decision == 'autism'),]
most_sig = autism[1,]

most_sig_inx = which(new_rules$pValue == most_sig$pValue)


?plotRule
plotRule(autcon, autism, type="heatmap", ind=most_sig_inx)

#new_rules$supportSetRHS

plotRule(autcon, autism, type="boxmap", ind=most_sig_inx)

