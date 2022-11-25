##------------------------------------------------------------
###Fit Small Aear Models using poll/election data
###Authors: Jiming Jiang, Yuanyuan Li, Peter X. K. Song

#This file generates all the results in the paper
##------------------------------------------------------------

################################Part 1. Support Functions####################
library(broom.mixed)
library(ggpubr)
library(MASS)
library(Matrix) 
library(lme4)
library(dplyr)
library(nlme)
setwd('..')
getwd()#workind directory should be the parent folder
source("R/sumca_funcs.R")

# Draw graphs using model output
plotall=function(lmm_result, data,labels="Abbrev"){
  #Capitalize variables
  data=rename(data,c("Candidate"="candidate", "Year"="year"))
  #labels includes the text to print in EBLUP figures
  result=list()
  ####Prediction interval
  data$fitted=lmm_result$EBLUP
  data$MSPE=lmm_result$MSPE
  result$fit=data.frame(EBLUP=data$fitted, MSPE=data$MSPE)
  #result$mspe_box=boxplot(data$MSPE, main="Distribution of MSPE")
  ct=qnorm(0.95)
  ##EBLUP
  result$eblup_pp=ggpaired(data, x = "Candidate", y = "fitted",id="id",
                           xlab = FALSE,point.size = 0,
                           ylab = "Predicted outcome",facet.by = "Year",
                           color = "white", line.color = "gray", line.size = 0.4,
                           short.panel.labs = F,label=labels,
                           font.label = list(size = 8),
                           label.select = list(criteria="(`y` >0.2  & `x` == 'R')"))+
    scale_color_manual(values = c("#0073C2FF", "#FC4E07"))+
    geom_hline(yintercept = 0, linetype = "dashed",size=0.3)+
    geom_point(aes(color = Candidate, shape=Candidate), position = position_dodge(0.3))+
    facet_grid(~Year)
  ##Plot MSPE
  mixed=data[c("Year","state","Candidate","fitted",
               "MSPE")]
  mixed=unique(mixed)
  mixed[,"state"]=as.numeric(mixed[,"state"])
  result$mspe_cycle=ggplot(mixed, aes(state, fitted)) +
    geom_errorbar(
      aes(ymin = fitted-ct*sqrt(MSPE), ymax = fitted+ct*sqrt(MSPE), color = Candidate,
          linetype=Candidate),
      position = position_dodge(0.3), width = 0.2)+
    geom_point(aes(color =Candidate, shape=Candidate), position = position_dodge(0.3)) +
    facet_wrap(~Year)+
    scale_color_manual(values = c("#0073C2FF", "#FC4E07"))+
    geom_hline(yintercept = 0, linetype = "dashed",size=0.3)+
    scale_x_discrete(labels = NULL, breaks = NULL)+
    ylab("Prediction interval")+xlab("State")+ theme(legend.position="top")
  return(result)
}

################################Part 2. SAE Model fitting####################
#### Read dataset####
model_data <- read.csv("data/model_data.csv")
model_data_5w <- read.csv("data/model_data_5w.csv")
election_2016 = read.csv("data/raw_data/election_2016.csv")
EV = election_2016$EV #get electoral college votes
states=unique(model_data$state) #all state names
#rename response variable to y
names(model_data)[names(model_data) == 'logp'] <- 'y' 
names(model_data_5w)[names(model_data_5w) == 'logp'] <- 'y'
# Convert integer columns to factor type
fac_coloumns=c("year","id","row_id","week")
model_data[,fac_coloumns]=lapply(model_data[,fac_coloumns],factor)
model_data_5w[,fac_coloumns]=lapply(model_data_5w[,fac_coloumns],factor)
# Create state-pollster label
model_data$pollandstate=paste(model_data$Abbrev,":", model_data$pollster)
model_data_5w$pollandstate=paste(model_data_5w$Abbrev,":", model_data_5w$pollster)

#### Fit SAE models####
#Model 1
raw_model1=lmer(y~candidate*year+(0+candidate|state),
                data=model_data)
# Summarize model output
model_output=function(model,data,state.label="Abbrev"){
lmm_result1=MSPEsumca(model,data = data)
model_plot=plotall(lmm_result1, data, state.label)
sum_stats=tidy(model)
#Bootstrap estimates for s.d. of random effects
sum_stats$std.error[sum_stats$effect=="ran_pars"]=lmm_result1$sd_sigma
sum_stats$statistic=sum_stats$estimate/sum_stats$std.error
return(list(eblup.fig=model_plot$eblup_pp,mspe.fig=model_plot$mspe_cycle,
            sum_table=sum_stats))
}
result1=model_output(raw_model1, model_data)
result1$eblup.fig
result1$mspe.fig
write.table(result1$sum_table,file="Results/Table1.csv",sep =
              ",", append=T,row.names = F)

#Model 2
raw_model2=lmer(y~candidate*year+(0+candidate|year:state),data=model_data)
#state effect: random slope of candidate; nested in year
result2=model_output(raw_model2, model_data)
result2$eblup.fig
result2$mspe.fig
write.table(result2$sum_table,file="Results/Table1.csv",sep =",",
            append=T,row.names = F)

#Model 3
raw_model3=lmer(y~candidate*year+(0+candidate|state)+(1|pollster),data=model_data)
result3=model_output(raw_model3, model_data, state.label="pollandstate")
result3$eblup.fig
result3$mspe.fig
write.table(result3$sum_table,file="Results/Table1.csv",sep =",",
            append=T,row.names = F)

#Model 4
raw_model4=lmer(y~candidate*year+(0+candidate|year:state)+(0+1|pollster), data=model_data)
result4=model_output(raw_model4, model_data, state.label="pollandstate")
result4$eblup.fig
result4$mspe.fig
write.table(result4$sum_table,file="Results/Table1.csv",sep =",",
            append=T,row.names = F)


#### Transfer learning########
# Predictor using mixed effects
ME_predictor=function(ME_train, polldata,fun){
  logp_bar= aggregate(log_poll~state+candidate, data = polldata, FUN=fun)
  logp_bar=merge(logp_bar,ME_train, by=c("state","candidate"))
  logp_bar$SAE= exp(logp_bar[,"log_poll"]-logp_bar[,"ME"])
  return(logp_bar[,c("state","candidate","SAE")])
}

# Transfer learning 
adjust_model=function(train_data, test_data,formula,weights=NA){
# get ME estimate
train_model = lmer(formula,
               data=train_data)
train_data$ME=predict(train_model)
ME_train=unique(train_data[c("state","candidate","ME")]) 
ME_train=aggregate(ME~state+candidate,data=ME_train,mean)
#Actual poll results by party
test_data$pollrate=test_data$poll_dem
test_data$pollrate[test_data$candidate=="R"]=test_data$poll_rep[test_data$candidate=="R"]
test_data$log_poll=log(test_data$pollrate)# log of poll rates
if(!any(is.na(weights))){
  # Weighted model-based support rate estimate (wSAE)
  test_data$log_poll=weights*test_data$log_poll
  SAE=ME_predictor(ME_train, test_data,fun=sum)
  # Weighted Poll of polls estimate(wPoP)
  test_data$pollrate=weights*test_data$pollrate
  PoP=aggregate(pollrate~state+candidate, data = test_data, sum)
}else{
SAE=ME_predictor(ME_train, test_data,fun=mean)
PoP=aggregate(pollrate~state+candidate, data = test_data, mean)
}
PoP=rename(PoP,c("PoP"="pollrate"))
#Compare real election results with estimates
test_data$realrate=test_data$real_dem
test_data$realrate[test_data$candidate=="R"]=test_data$real_rep[test_data$candidate=="R"]
real_test=unique(test_data[c("state","candidate","realrate")])
#Merge with SAE estimates
real_test=merge(real_test,SAE, by=c("state","candidate"))
# Merge with PoP estimates
real_test=merge(real_test,PoP, by=c("state","candidate"))
Init=factor(levels=c("D", "R"),rep("D",length(states)))
Est_eval=data.frame(state=states,Actual=Init,
                       SAE=Init, PoP=Init,
                       Votes=EV)
for (s in states){
  rs=real_test[real_test$state==s,]
  #actual winner
  Est_eval$Actual[Est_eval$state==s]=rs$candidate[which.max(rs$realrate)]
  #winner predicted by SAE
  Est_eval$SAE[Est_eval$state==s]=rs$candidate[which.max(rs$SAE)]
  #winner predicted by PoP
  Est_eval$PoP[Est_eval$state==s]=rs$candidate[which.max(rs$PoP)]
}
result=list()
# get all wrong predictions
result$wrong=Est_eval[(Est_eval$SAE!=Est_eval$Actual)|
              (Est_eval$PoP!=Est_eval$Actual),]
# national votes differences
result$Actual_votes_sum=aggregate(Votes~Actual,data=Est_eval,sum)
result$SAE_votes_sum=aggregate(Votes~SAE,data=Est_eval,sum)
result$PoP_votes_sum=aggregate(Votes~PoP,data=Est_eval,sum)
result$Est_eval=Est_eval
return(result)
}

formula_1 = as.formula(y~candidate+(0+candidate|state))#Model 1
formula_3 = as.formula(y~candidate+(0+candidate|state)+(1|pollster))#Model 3
data16=model_data[model_data$year==2016,]
data20=model_data[model_data$year==2020,]

# Predict 2020 winners using 2016 data (SAE using Model1)
model1_adjust1=adjust_model(data16,data20,formula_1)
model1_adjust1$wrong
write.table(model1_adjust1$wrong, sep=",",file="Results/Table2.csv",
            row.names = FALSE)
model1_adjust1$Actual_votes_sum
model1_adjust1$SAE_votes_sum
model1_adjust1$PoP_votes_sum
#model1_adjust1$Est_eval

# Predict 2016 winners using 2020 data (SAE using Model1)
model1_adjust2=adjust_model(data20,data16,formula_1)
model1_adjust2$wrong
write.table(model1_adjust2$wrong, sep=",",file="Results/Table2.csv",append = TRUE,
            row.names = FALSE)
model1_adjust2$Actual_votes_sum
model1_adjust2$SAE_votes_sum
model1_adjust2$PoP_votes_sum
#model1_adjust2$Est_eval

#Predict 2020 winners using 2016 data (SAE using Model3)
model3_adjust1=adjust_model(data16,data20,formula_3)
model3_adjust1$wrong
write.table(model3_adjust1$wrong, sep=",",file="Results/Table2.csv",append = TRUE,
            row.names = FALSE)
model3_adjust1$Actual_votes_sum
model3_adjust1$SAE_votes_sum
model3_adjust1$PoP_votes_sum

# Predict 2016 winners using 2020 data (SAE using Model3)
model3_adjust2=adjust_model(data20,data16,formula_3)
model3_adjust2$wrong
write.table(model3_adjust2$wrong, sep=",",file="Results/Table2.csv",append = TRUE,
            row.names = FALSE)
model3_adjust2$Actual_votes_sum
model3_adjust2$SAE_votes_sum
model3_adjust2$PoP_votes_sum

#### Ranking the pollsters ####
####1.Ranking pollster using random effects####
#Seperate
m3.16=lmer(formula_3,data=data16)
re1=ranef(m3.16)$pollster
re1 <- data.frame(pollster = row.names(re1), RE=re1[,1],RE_abs=abs(re1[,1]))
m3.20=lmer(formula_3,data=data20)
re2=ranef(m3.20)$pollster
re2 <- data.frame(pollster = row.names(re2), RE=re2[,1], RE_abs=abs(re2[,1]))
comb_re1=merge(re1,re2,by.x="pollster",by.y="pollster",
               all = T,suffixes = c(".2016",".2020"))
comb_re2=merge(re1,re2,by.x="pollster",by.y="pollster",
            all = F,suffixes = c(".2016",".2020"))

#Combined
m3=as.formula(y~candidate*year+(0+candidate|state)+(1|pollster))#Model 3
m3.all=lmer(m3,data=model_data)
re12=ranef(m3.all)$pollster
re12 <- data.frame(pollster = row.names(re12), RE=re12[,1], RE_abs=abs(re12[,1]))

t3=data.frame(comb=re12$pollster[order(re12$RE_abs)][1:10],
           rank2016=comb_re1$pollster[order(comb_re1$RE_abs.2016)][1:10],
           rank2020=comb_re1$pollster[order(comb_re1$RE_abs.2020)][1:10])
write.table(t3,file="Results/Table3.csv",sep =
              ",", row.names = F)
t3

t3_bottom=data.frame(comb=re12$pollster[order(re12$RE_abs,decreasing = T)][1:10],
              rank2016=comb_re1$pollster[order(comb_re1$RE_abs.2016,decreasing = T)][1:10],
              rank2020=comb_re1$pollster[order(comb_re1$RE_abs.2020,decreasing = T)][1:10])
write.table(t3_bottom,file="Results/Table3.csv",sep =
              ",", row.names = F,append = TRUE)
t3_bottom

t5=data.frame(rank2016=comb_re2$pollster[order(comb_re2$RE_abs.2016)][1:10],
              rank2020=comb_re2$pollster[order(comb_re2$RE_abs.2020)][1:10])
write.table(t5,file="Results/Table5.csv",sep =
              ",", row.names = F)
t5

t5_bottom=data.frame(rank2016=comb_re2$pollster[order(comb_re2$RE_abs.2016,decreasing = T)][1:10],
                     rank2020=comb_re2$pollster[order(comb_re2$RE_abs.2020,decreasing = T)][1:10])
write.table(t5_bottom,file="Results/Table5.csv",sep =
              ",", row.names = F,append = TRUE)
t5_bottom

####2.Ranking pollster using mixed effects####
pollster_analysis=function(formula,data){
model= lmer(formula, data=data)
data$ME=predict(model)
pollster_bias=aggregate(ME~pollster,data = data,mean)
pollster_bias$ME_abs=abs(pollster_bias$ME)
return(pollster_bias)
}

#seperate
p1=pollster_analysis(formula_3,data16)
p2=pollster_analysis(formula_3,data20)
comb_me1=merge(p1,p2,by.x="pollster",by.y="pollster",
           all = T,suffixes = c(".2016",".2020"))
comb_me2=merge(p1,p2,by.x="pollster",by.y="pollster",
               all = F,suffixes = c(".2016",".2020"))


###Combined
m3=as.formula(y~candidate*year+(0+candidate|state)+(1|pollster))#Model 3
me12=pollster_analysis(m3,model_data)
t4=data.frame(comb=me12$pollster[order(me12$ME_abs)][1:10],
              rank2016=comb_me1$pollster[order(comb_me1$ME_abs.2016)][1:10],
              rank2020=comb_me1$pollster[order(comb_me1$ME_abs.2020)][1:10])
write.table(t4,file="Results/Table4.csv",sep =
              ",", row.names = F)
t4

t4_bottom=data.frame(comb=me12$pollster[order(me12$ME_abs,decreasing = T)][1:10],
                     rank2016=comb_me1$pollster[order(comb_me1$ME_abs.2016,decreasing = T)][1:10],
                     rank2020=comb_me1$pollster[order(comb_me1$ME_abs.2020,decreasing = T)][1:10])
write.table(t4_bottom,file="Results/Table4.csv",sep =
              ",", row.names = F,append = TRUE)
t4_bottom

t6=data.frame(rank2016=comb_me2$pollster[order(comb_me2$ME_abs.2016)][1:10],
              rank2020=comb_me2$pollster[order(comb_me2$ME_abs.2020)][1:10])
write.table(t6,file="Results/Table6.csv",sep =
              ",", row.names = F)
t6

t6_bottom=data.frame(rank2016=comb_me2$pollster[order(comb_me2$ME_abs.2016,decreasing = T)][1:10],
                     rank2020=comb_me2$pollster[order(comb_me2$ME_abs.2020,decreasing = T)][1:10])
write.table(t6_bottom,file="Results/Table6.csv",sep =
              ",", row.names = F,append = TRUE)
t6_bottom



#### Weighted predictions based on historic performance####
# Use 1/EBLUP as weiights
weights_fun=function(data_EBLUP){
  EBLUP=data_EBLUP[,2]
  data_EBLUP$non_inform=1/mean(EBLUP)#non-informative prior
  data_EBLUP$inform=1/EBLUP
  return(data_EBLUP)
}
# Use 1/rank(EBLUP) as weiights
weights_fun_ranks=function(data_EBLUP){
  EBLUP=data_EBLUP[,2]
  data_EBLUP$inform=1/rank(EBLUP)
  data_EBLUP$non_inform=1/(mean(rank(EBLUP)))#non-informative prior
  return(data_EBLUP)
}
# Scale weiights
wts_scale=function(dat){
  total=aggregate(inform~state+candidate,FUN=sum,data=dat)
  total=rename(total,c("sum_w"="inform"))
  dat2=merge(dat,total,by=c("state","candidate"),all.x=TRUE)
  dat2$weights=dat2$inform/dat2$sum_w
  return(dat2)
}
# predict 2020 by 2016
weights_16=weights_fun(p1[c("pollster","ME_abs")])
weight_data=merge(data20,weights_16,by="pollster",all.x=TRUE)
weight_data$inform[is.na(weight_data$inform)]=weights_16$non_inform[1]
wt_df=wts_scale(weight_data)
weights=wt_df$weights[order(wt_df$row_id)]
model3_eblup_20=adjust_model(data16,data20,formula_3,weights=weights)
model3_eblup_20$wrong


weights_16=weights_fun_ranks(p1[c("pollster","ME_abs")])
weight_data=merge(data20,weights_16,by="pollster",all.x=TRUE)
weight_data$inform[is.na(weight_data$inform)]=weights_16$non_inform[1]
wt_df=wts_scale(weight_data)
weights=wt_df$weights[order(wt_df$row_id)]
model3_rank_20=adjust_model(data16,data20,formula_3,weights=weights)
model3_rank_20$wrong
write.table(model3_rank_20$wrong,file="Results/Table7.csv",sep =
              ",", row.names = F)

# predict 2016 by 2020
weights_20=weights_fun(p2[c("pollster","ME_abs")])
weight_data=merge(data16,weights_20,by="pollster",all.x=TRUE)
weight_data$inform[is.na(weight_data$inform)]=weights_20$non_inform[1]
wt_df=wts_scale(weight_data)
weights=wt_df$weights[order(wt_df$row_id)]
model3_eblup_16=adjust_model(data20,data16,formula_3,weights=weights)
model3_eblup_16$wrong

weights_20=weights_fun_ranks(p2[c("pollster","ME_abs")])
weight_data=merge(data16,weights_20,by="pollster",all.x=TRUE)
weight_data$inform[is.na(weight_data$inform)]=weights_20$non_inform[1]
wt_df=wts_scale(weight_data)
weights=wt_df$weights[order(wt_df$row_id)]
model3_rank_16=adjust_model(data20,data16,formula_3,weights=weights)
model3_rank_16$wrong
write.table(model3_rank_16$wrong,file="Results/Table7.csv",sep =
              ",", row.names = F,append = TRUE)

################################Part 3. Time series model fiitting#################
#### Exploration analysis ####
model_data_5w$pdiff=rep(0,nrow(model_data_5w))
for (i in 1:nrow(model_data_5w)){
  model_data_5w$pdiff[i]=c(model_data_5w[i,"poll_dem"]-model_data_5w[i,"real_dem"],
                        model_data_5w[i,"poll_rep"]-model_data_5w[i,"real_rep"])[model_data_5w[i,"candidate"]]/100
}
model_dat_avg=aggregate(y~year+state+candidate+pollster+week,data=model_data_5w,mean)
#add state abbreviation
model_dat_avg=merge(unique(model_data_5w[,c("state","Abbrev")]),model_dat_avg, by="state")
avg_actual=aggregate(pdiff~year+state+candidate+pollster+week,data=model_data_5w,mean)
model_dat_avg$pdiff=avg_actual$pdiff
df <- rbind(data.frame(x=model_dat_avg$pdiff, type='difference'),
            data.frame(x=model_dat_avg$y, type='log ratio'))
ggplot(df, aes(x, group=type, col=type)) + geom_density(position='dodge')+
  scale_x_continuous(limits=c(-0.5,0.5),breaks = seq(-0.5, 0.5, 0.1),
                     minor_breaks = seq(-0.5,0.5,0.05))
model_dat_avg$id=interaction(model_dat_avg[,c("year","state","pollster","week")])
model_dat_avg$id=as.numeric(model_dat_avg$id)

## Plot graphs on the same y axis
plotall_ylimit=plotall=function(lmm_result, data,labels="Abbrev"){
  #Capitalize variables
  data=rename(data,c("Candidate"="candidate", "Year"="year"))
  #labels includes the text to print in EBLUP figures
  result=list()
  ####Prediction interval
  data$fitted=lmm_result$EBLUP
  data$MSPE=lmm_result$MSPE
  result$fit=data.frame(EBLUP=data$fitted, MSPE=data$MSPE)
  #result$mspe_box=boxplot(data$MSPE, main="Distribution of MSPE")
  ct=qnorm(0.95)
  ##EBLUP
  result$eblup_pp=ggpaired(data, x = "Candidate", y = "fitted",id="id",
                           xlab = FALSE,point.size = 0,ylim=c(-0.3,0.75),
                           ylab = "Predicted outcome",facet.by = "Year",
                           color = "white", line.color = "gray", line.size = 0.4,
                           short.panel.labs = F,label=labels,
                           font.label = list(size = 8),
                           label.select = list(criteria="(`y` >0.2  & `x` == 'R')"))+
    scale_color_manual(values = c("#0073C2FF", "#FC4E07"))+
    geom_hline(yintercept = 0, linetype = "dashed",size=0.3)+
    geom_point(aes(color = Candidate, shape=Candidate), position = position_dodge(0.3))+
    facet_grid(~Year)
  ##Plot MSPE
  mixed=data[c("Year","state","Candidate","fitted",
               "MSPE")]
  mixed=unique(mixed)
  mixed[,"state"]=as.numeric(mixed[,"state"])
  result$mspe_cycle=ggplot(mixed, aes(state, fitted)) +
    geom_errorbar(
      aes(ymin = fitted-ct*sqrt(MSPE), ymax = fitted+ct*sqrt(MSPE), color = Candidate,
          linetype=Candidate),
      position = position_dodge(0.3), width = 0.2)+
    geom_point(aes(color =Candidate, shape=Candidate), position = position_dodge(0.3)) +
    facet_wrap(~Year)+
    scale_color_manual(values = c("#0073C2FF", "#FC4E07"))+
    geom_hline(yintercept = 0, linetype = "dashed",size=0.3)+ylim(-0.4,0.8)+
    scale_x_discrete(labels = NULL, breaks = NULL)+
    ylab("Prediction interval")+xlab("State")+ theme(legend.position="top")
  return(result)
}


# Model 1 without time correlation
raw_model1=lmer(y~candidate*year+(0+candidate|state),
                data=model_dat_avg)
lmm_result1=MSPEsumca(raw_model1,data = model_dat_avg)
model1=plotall_ylimit(lmm_result1, model_dat_avg)
model1$eblup_pp
model1$mspe_cycle

#  Model 1 with time correlation
model_dat_avg$week=as.numeric(model_dat_avg$week)
lme_fit1 <- lme(y~candidate*year,random=~0+candidate|state,
                data=model_dat_avg,correlation =
                  corAR1(form=~week|state/year/candidate/pollster))
lmm_result1.1=MSPEsumca(lme_fit1,data = model_dat_avg)
#summary(lmm_result1.1$EBLUP-predict(lme_fit1))
model1.1=plotall_ylimit(lmm_result1.1, model_dat_avg)
model1.1$eblup_pp
model1.1$mspe_cycle

save.image(file='Results/results.RData')

