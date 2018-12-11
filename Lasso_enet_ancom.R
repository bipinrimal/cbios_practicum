library(rhdf5)
library(phyloseq)
library(glmnet)
library(plyr)
source("merging.R")
source("read_hdf5.r")
source("ANCOM_updated_code.R")
##############################################################################
########################functions to read in phyloseq file####################
##############################################################################

###############################################################
#########################Analysis###############################
################################################################
# #Read dataset, extract otu table and the biome response
# data = readRDS("subsetted_physeq2.rds")
# OTU = otu_table(data)
# sample = sample_data(data)
# biome=sample$env_biome
# 
# #Pruned the OTUs that has abundance 0 in more than 99% of samples
# keep = (apply(OTU,1,function(x) sum(x==0)) <= (0.99*dim(OTU)[2]))
# sum(keep==TRUE)
# 
# pruned.data = prune_taxa(keep,data)
# pruned.OTU = otu_table(pruned.data)
# dim(pruned.OTU)
# 
# #Feed the newly pruned phyloseq data to the crawling algorithm
# ptm = proc.time()
# merged.data = merge_brothers(pruned.data)
# time = (proc.time()-ptm)[3]
# #take 14 hours to merge
# saveRDS(merged.data, "merged_data.rds")

merged.data = readRDS("merged_data.rds")
merged.OTU = otu_table(merged.data)
dim(merged.OTU)
#The merging with brother code produce error
# Error in merge_taxa.indices.internal(x, eqtaxa, archetype) : 
#   invalid OTU indices provided as eqtaxa or archetype.
# ptm = proc.time()
# merged.uncle = merge_brother_with_uncle(merged.data)
# time2 = (proc.time()-ptm)[3]

# ptm = proc.time()
# merge.correlated = filter_tree_correlation(merged.data,thres_corr=0.7)
# time2 = (proc.time()-ptm)[3]
# #takes 10 mins to merge
# dim(otu_table(merge.correlated))
# uncor.OTU = otu_table(merge.correlated)
# uncor.OTU[1,10:500]
# #negative OTU shows up here
# saveRDS(merge.correlated, "merged_correlated_data.rds")

#########For Lasso to work, we need to create a model matrix, i.e. the predictor matrix, columns are 
#########microbiomes, rows are observations/samples (i.e. human)
ptm=proc.time()
m.X = model.matrix(~.-1,t(merged.OTU))
time3 = (proc.time()-ptm)[3]

sample = sample_data(merged.data)
biome=sample$env_biome
Y = rep("0", 1300)
Y[biome == "urban biome"] = "1"

#########glmnet is a function to perform lasso, family = binomial means our response is binary
#########lambda is a penalization parameters, the larger the lambdas, the more coefficients/predictors
#########lasso forces to zero. This is how lasso obtain a sparse model. We are going to create a model
#########over a grid of lasso from small to large.
grid =10^ seq (10,-2, length =100)
lasso <- glmnet(x = m.X, y=Y, alpha=1, family="binomial",lambda = grid)
dim(coef(lasso))


########Divide data to training and testing set to evaluate the model ability
########to predict biome. That is, we train out model with training data, then use the 
########model to predict the test data. Here, we select testing size = 500
set.seed(1)
train = sample(1: nrow(m.X), 800)
test =(-train)
Y.test =Y[test]
########Select the best tuning/penalizing parameter lambda through 10-fold cross validation
########The function for this is cv.glmnet.
cv.out =cv.glmnet(m.X[train,], Y[train], alpha =1, family = "binomial")
plot(cv.out)
bestlam =cv.out$lambda.min
bestlam
########Now, let's find out which microbiome/predictor lasso keeps in the model
lasso.train = glmnet(m.X[train,], Y[train],alpha = 1, family = "binomial", lambda = grid)
lasso.coefs = predict(lasso.train, type = "coefficients",s= bestlam )
lasso.important.microbiomes = lasso.coefs@Dimnames[[1]][which(lasso.coefs!=0)]
lasso.important.microbiomes
length(lasso.important.microbiomes) 
#the model selects 73 "merged predictors" out of 1829. Note that there are many 
#nested microbiomes in 1 "merged predictor"

########Model's ability to predict biomes
lasso.train = glmnet(m.X[train,], Y[train],alpha = 1, family = "binomial", lambda = grid)
plot(lasso.train)
prob= predict(lasso.train,s=bestlam , newx=m.X[test,], type = "response")
pred = rep("1", 500)
pred[prob < 0.5] = "0"

sum(pred != Y[test])
lasso.result = table(pred, Y[test])
lasso.result      #the columns are test observation, the rows are model prediction
#######The table above tells us that our model wrongly predict 30 biomes out of 500, which is not so bad.
#######It underpredict village biome, which might be due to the fact that village biome makes up only 
#######7% of our data. It correcly predicts 437 urban biomes out of 449, and correctly predicts 33 biomes out of 51.
true.pos = lasso.result[1,1]/sum(lasso.result[,1])
true.neg = lasso.result[2,2]/sum(lasso.result[,2])
false.pos = lasso.result[2,1]/sum(lasso.result[,1])
false.neg = lasso.result[1,2]/sum(lasso.result[,2])
c(true.pos, true.neg, false.pos, false.neg)

############################################lasso is known to underperform in terms of prediction compared 
########to ridge and enet if multicollinearity is present. Let's try e-net
# ELASTIC NET WITH 0 < ALPHA < 1
#We will perform cross validation to select both alpha, and lambda within each alpha. 
#To do this we need to fix the folds across alpha
a <- seq(0.1, 0.9, 0.05)
foldid=sample(1:10,size=length(Y[train]),replace=TRUE)

ptm = proc.time()
CV = lapply(1:length(a), function(x){
  cv.glmnet(x= m.X[train,],y=Y[train],foldid=foldid,alpha=a[x], family = "binomial")
})
time4 = (proc.time()-ptm)[3]
#take about a minute to compute all CV results across alpha

par(mfrow=c(2,3))
min.mse = c()
for(i in 1:length(a)) {
  plot(CV[[i]])
  min.mse[i] = CV[[i]]$cvm[CV[[i]]$lambda == CV[[i]]$lambda.min]
}

# dev.off()
# plot(log(CV[[1]]$lambda.min), min.mse[1])
# for (i in 2:length(a)){
#   points(log(CV[[i]]$lambda.min), min.mse[i])
# }

plot(log(CV[[1]]$lambda),CV[[1]]$cvm,pch=19,col=1,xlab="log(Lambda)",ylab=CV[[1]]$name)
col=rgb(0,0,1:17,maxColorValue = 17)
for(i in 2:length(a)) points(log(CV[[i]]$lambda),CV[[i]]$cvm,pch=19,col=col[i])
legend("topleft",legend=c("alpha= 0.1","alpha= 0.9"),pch=19,col=c(1, col[i]))

#alpha = 0.1 has min.mse that dominates other alpha. Looking at this trend, we conclude that lasso (alpha = 1) 
#does not perform well. In fact, its min.mse = 0.29 on the same train data. We don't let alpha go all the way to 0,
#since alpha = 0 corresponds to ridge regression, and we want variable selection here
cv.out$cvm[cv.out$lambda == cv.out$lambda.min]

enet.train = glmnet(x= m.X[train,],y=Y[train],alpha=0.7, family = "binomial")
enet.coefs = predict(enet.train, type = "coefficients",s= CV[[which(a==0.5)]]$lambda.min)
enet.important.microbiomes = enet.coefs@Dimnames[[1]][which(enet.coefs!=0)]
enet.important.microbiomes
write(enet.important.microbiomes,"important_microbiomes.txt")

enet.prob= predict(enet.train,s=CV[[which(a==0.5)]]$lambda.min , newx=m.X[test,], type = "response")
enet.pred = rep("1", 500)
enet.pred[enet.prob < 0.5] = "0"

sum(enet.pred != Y[test])
enet.result = table(enet.pred, Y[test])
enet.result      
result

# install.packages("class")
# library(class)
# knn.pred=knn(m.X[train,], m.X[test,],Y[train],k=1)
# sum(Y[test]!= knn.pred )
# knn.result = table(knn.pred,Y[test])
# 
# which(lasso.coefs!=0)
# knn.pred2=knn(m.X[train,which(lasso.coefs!=0)], m.X[test,which(lasso.coefs!=0)],Y[train],k=1)
# sum(Y[test]!= knn.pred2 )
# knn.result2 = table(knn.pred2,Y[test])
# 
# knn.pred3=knn(m.X[train,which(enet.coefs!=0)], m.X[test,which(enet.coefs!=0)],Y[train],k=1)
# sum(Y[test]!= knn.pred3)
# knn.result3 = table(knn.pred3,Y[test])


#########################################################################
#########################################################################
#t test on raw abundances




########################################################################
########################################################################
#We have 1153 urban samples and only 147 village samples. See any problems?
#Comparing Microbial composition
#Use ANCOM method (analysis of composition of microbiomes)
Sample.ID = rownames(m.X)
OTUdat = data.frame(Sample.ID,m.X)
Vardat = data.frame(Sample.ID,Y)

ptm = proc.time()
comparison_test=ANCOM.main(OTUdat=OTUdat,
                             Vardat = Vardat,
                            adjusted=FALSE,
                             repeated=F,
                             main.var="Y",
                             adj.formula=NULL,
                             repeat.var=NULL,
                             longitudinal=FALSE,
                             random.formula=NULL,
                             multcorr=2,
                             sig=0.05,
                             prev.cut=0.99)
time5 = (proc.time()-ptm)[3]
#takes 24 hours if prev.cut = 0.99, takes 2.5 hhrs if prev.cut = 0.75
comparison_test$W.taxa
saveRDS(comparison_test$W.taxa, file = "ANCOM_results.rds")

#################################################################################################
#################################BALANCED SAMPLES#################################################
#################################################################################################
set.seed(1)
vil.idx = which(Y == 0)
urb.idx = which(Y == 1)
index = sample(1:7,1153, replace=T, prob = rep(1/7,7))
folds = list()
for (i in 1:7){
  folds[[i]] = c(vil.idx, urb.idx[which(index==i)])
}

#
ptm = proc.time()
results = lapply(1:7, function(i){
  keep = rep(FALSE,1300)
  #Keep only samples corresponding to index in folds[[i]]
  keep[folds[[i]]] = TRUE
  fold = prune_samples(keep,merged.data)
  
  fold.OTU = otu_table(fold)
  dim(fold.OTU)
  m.X = model.matrix(~.-1,t(fold.OTU))
  
  sample = sample_data(fold)
  biome=sample$env_biome
  Y.fold = rep("0", length(biome))
  Y.fold[biome == "urban biome"] = "1"
  
  n= length(biome)
  set.seed(1)
  train = sample(1: nrow(m.X), round(2*n/3))
  test =(-train)
  Y.fold.test =Y.fold[test]
  
  # grid =10^ seq (10,-2, length =100)
  # lasso <- glmnet(x = m.X, y=Y.fold, alpha=1, family="binomial",lambda = grid)
  # dim(coef(lasso))
  # 
  # cv.out =cv.glmnet(m.X[train,], Y[train], alpha =1, family = "binomial")
  # bestlam =cv.out$lambda.min
  # bestlam
  # ########Now, let's find out which microbiome/predictor lasso keeps in the model
  # lasso.train = glmnet(m.X[train,], Y.fold[train],alpha = 1, family = "binomial", lambda = grid)
  # lasso.coefs = predict(lasso.train, type = "coefficients",s= bestlam )
  # lasso.important.microbiomes = lasso.coefs@Dimnames[[1]][which(lasso.coefs!=0)]
  # lasso.important.microbiomes
  # length(lasso.important.microbiomes) 
  # #the model selects 73 "merged predictors" out of 1829. Note that there are many 
  # #nested microbiomes in 1 "merged predictor"
  # 
  # ########Model's ability to predict biomes
  # lasso.train = glmnet(m.X[train,], Y.fold[train],alpha = 1, family = "binomial", lambda = grid)
  # plot(lasso.train)
  # prob= predict(lasso.train,s=bestlam , newx=m.X[test,], type = "response")
  # pred = rep("1", n-round(2*n/3))
  # pred[prob < 0.5] = "0"
  # sum(pred != Y.fold[test])
  # lasso.result = table(pred, Y.fold[test])
  # lasso.result      #the columns are test observation, the rows are model prediction
  
  a <- seq(0.1, 0.9, 0.05)
  foldid=sample(1:10,size=length(Y.fold[train]),replace=TRUE)

  CV = lapply(1:length(a), function(x){
    cv.glmnet(x= m.X[train,],y=Y.fold[train],foldid=foldid,alpha=a[x], family = "binomial")
  })
  #take about half a minute to compute all CV results across alpha

  par(mfrow=c(2,4))
  # min.mse = c()
  # for(i in 1:length(a)) {
  #   plot(CV[[i]])
  #   min.mse[i] = CV[[i]]$cvm[CV[[i]]$lambda == CV[[i]]$lambda.min]
  # }


  # plot(log(CV[[1]]$lambda),CV[[1]]$cvm,pch=19,col=1,xlab="log(Lambda)",ylab=CV[[1]]$name)
  # col=rgb(0,0,1:17,maxColorValue = 17)
  # for(i in 2:length(a)) points(log(CV[[i]]$lambda),CV[[i]]$cvm,pch=19,col=col[i])
  # legend("topleft",legend=c("alpha= 0.1","alpha= 0.9"),pch=19,col=c(1, col[i]))
  # cv.out$cvm[cv.out$lambda == cv.out$lambda.min]

  enet.train = glmnet(x= m.X[train,],y=Y.fold[train],alpha=0.7, family = "binomial")
  enet.coefs = predict(enet.train, type = "coefficients",s= CV[[which(a==0.5)]]$lambda.min)
  enet.important.microbiomes = enet.coefs@Dimnames[[1]][which(enet.coefs!=0)]
  #return(enet.important.microbiomes)
  #write(enet.important.microbiomes,"important_microbiomes.txt")
  #
  enet.prob= predict(enet.train,s=CV[[which(a==0.5)]]$lambda.min , newx=m.X[test,], type = "response")
  enet.pred = rep("1", n-round(2*n/3))
  enet.pred[enet.prob < 0.5] = "0"

  sum(enet.pred != Y.fold[test])
  enet.result = table(enet.pred, Y.fold[test])
  return(enet.result)
})
enet.time = (proc.time()-ptm)[3]
enet.microbiomes = Reduce(intersect,enet.results)
write(enet.microbiomes,"enet_microbiomes.txt")
dev.off()
out = which(count(unlist(enet.results))$freq>4)
enet.microbiomes2 = count(unlist(enet.results))$x[out]
write(enet.microbiomes2,"enet_microbiomes_5andmore.txt")
