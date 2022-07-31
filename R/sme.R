##------------------------------------------------------------
###Fit Small Aear Models using poll/election data
###Authors: Jiming Jiang, Yuanyuan Li, Peter X. K. Song
##------------------------------------------------------------

library(broom.mixed)
library(ggpubr)
library(MASS)
library(Matrix) 
library(lme4)
library(dplyr)
library(nlme)
#setwd("~/paper/Election/PEA")
source("R/sumca_funcs.R")

###read dataset
long_dat <- read.csv("data/long_dat.csv")
long_dat$row_id=1:nrow(long_dat)
election_2016 = read.csv("data/raw_data/election_2016.csv")
EV = election_2016$EV
long_dat$week=-long_dat$week
long_dat$candidate=factor(long_dat$candidate)
long_dat$year=factor(long_dat$year)
polls=levels(long_dat$pollster)
nlevels(long_dat$pollster)
states=unique(long_dat$state)
poll16=unique(long_dat$pollster[long_dat$year==2016])
poll20=unique(long_dat$pollster[long_dat$year==2020])
mean(poll16 %in% poll20)
mean(poll20 %in% poll16)
pdiff=rep(0,nrow(long_dat))
for (i in 1:nrow(long_dat)){
  pdiff[i]=c(long_dat[i,"poll_dem"]-long_dat[i,"real_dem"],
    long_dat[i,"poll_rep"]-long_dat[i,"real_rep"])[long_dat[i,"candidate"]]
}

long_dat$pdiff=pdiff/100
names(long_dat)[names(long_dat) == 'logp'] <- 'y'
long_dat_avg=aggregate(y~year+state+candidate+pollster+week,data=long_dat,mean)
avg2=aggregate(pdiff~year+state+candidate+pollster+week,data=long_dat,mean)
long_dat_avg$pdiff=avg2$pdiff
error_id=interaction(long_dat_avg[,c(1,2,3,4)])
long_dat_avg$error_id=as.numeric(error_id)

par(mfrow=c(1,1))
plot(density(long_dat_avg$pdiff),col="red",main="Density plots")
lines(density(long_dat_avg$y),col="blue")
legend("topleft",        # Add legend to plot
       c("differences", "log ratios"),
       col = c("red","blue"),
       lty = 1)
dat <- data.frame(dens = c(long_dat_avg$pdiff, long_dat_avg$y)
                  , lines = rep(c("a", "b")))

df <- rbind(data.frame(x=long_dat_avg$pdiff, type='difference'),
            data.frame(x=long_dat_avg$y, type='log ratio'))
ggplot(df, aes(x, group=type, col=type)) + geom_density(position='dodge')
ggplot(df, aes(x, group=type, col=type)) + geom_density(position='dodge')+
  scale_x_continuous(limits=c(-0.5,0.5),breaks = seq(-0.5, 0.5, 0.1),
                     minor_breaks = seq(-0.5,0.5,0.05))

vec1 <- data.frame(x=long_dat_avg$pdiff)
vec2 <- data.frame(x=long_dat_avg$y)

library(ggplot2)

ggplot() + geom_density(aes(x=x), colour="red", data=vec1) + 
  geom_density(aes(x=x), colour="blue", data=vec2)


#####################################Section 2. Fit SAE models
plotall=function(lmm_result, data,labels="state"){
  #labels includes the text to print in EBLUP figures
  result=list()
  ####Prediction interval
  data$fitted=lmm_result$EBLUP
  data$MSPE=lmm_result$MSPE
  result$fit=data.frame(EBLUP=data$fitted, MSPE=data$MSPE)
  result$mspe_box=boxplot(data$MSPE, main="Distribution of MSPE")
  ct=qnorm(0.95)
  ##EBLUP
  result$eblup_pp=ggpaired(data, x = "candidate", y = "fitted",id="id",
                           xlab = FALSE,point.size = 0,
                           ylab = "Predicted outcome",
                           color = "white", line.color = "gray", line.size = 0.4,
                           facet.by = "year", short.panel.labs = F,label=labels,
                           font.label = list(size = 8),
                           label.select = list(criteria="(`y` >0.2  & `x` == 'R')"))+
    scale_color_manual(values = c("#0073C2FF", "#FC4E07"))+geom_point(aes(color = candidate, shape=candidate), position = position_dodge(0.3))
  ##Plot by states
  mixed=data[c("year","state","candidate","fitted",
               "MSPE")]
  mixed=unique(mixed)
  result$mspe_states=ggplot(mixed, aes(year, fitted)) +
    geom_errorbar(
      aes(ymin = fitted-ct*sqrt(MSPE), ymax = fitted+ct*sqrt(MSPE), color = candidate),
      position = position_dodge(0.3), width = 0.2)+
    geom_point(aes(color = candidate), position = position_dodge(0.3)) +
    facet_wrap(~state)+
    scale_color_manual(values = c("#0073C2FF", "#FC4E07")) 
  mixed[,"state"]=as.numeric(mixed[,"state"])
  result$mspe_cycle=ggplot(mixed, aes(state, fitted)) +
    geom_errorbar(
      aes(ymin = fitted-ct*sqrt(MSPE), ymax = fitted+ct*sqrt(MSPE), color = candidate,
          linetype=candidate),
      position = position_dodge(0.3), width = 0.2)+
    geom_point(aes(color =candidate, shape=candidate), position = position_dodge(0.3)) +
    #facet_grid(rows = vars(year))+
    facet_wrap(~year)+
    scale_color_manual(values = c("#0073C2FF", "#FC4E07"))+
    ylab("Prediction Interval")
  return(result)
}


long_dat_avg=long_dat_avg[order(long_dat_avg[,1],long_dat_avg[,2],
                                long_dat_avg[,3],long_dat_avg[,4]),]
long_dat_avg$id=interaction(long_dat_avg[,c(1,2,4,5)])
long_dat_avg$id=as.numeric(long_dat_avg$id)

#Model 1
raw_model1=lmer(y~candidate*year+(0+candidate|state),
                data=long_dat)

lmm_result1=MSPEsumca(raw_model1,data = long_dat)
summary(lmm_result1$EBLUP-predict(raw_model1))
lmm_result1$sd_beta
lmm_result1$sd_sigma#Bootstrap estimators in Table 1
model1=plotall(lmm_result1, long_dat,"")
model1$eblup_pp#Figure 1
model1$mspe_cycle#Figure 5
write.table(tidy(raw_model1),file="Results/Table1.csv",sep =
              ",", append=T)

#Model 2
raw_model2=lmer(y~candidate*year+(0+candidate|year:state),data=long_dat_avg)
long_dat_avg$grp=interaction(long_dat_avg$year,long_dat_avg$state)
#state effect: random slope of candidate; nested in year
write.table(tidy(raw_model2),file="Results/Table1.csv",sep =
              ",", append=T)
lmm_result2=MSPEsumca(raw_model2,data = long_dat_avg)
summary(lmm_result1$EBLUP-predict(raw_model2))
lmm_result2$sd_beta
lmm_result2$sd_sigma#Bootstrap estimators in Table 1
model2=plotall(lmm_result2, long_dat_avg)
model2$eblup_pp#Figure 1
model2$mspe_cycle#Figure 5



#Model 3
raw_model3=lmer(y~candidate*year+(0+candidate|state)+(1|pollster),data=long_dat)
raw_model3

summary(raw_model3)
write.table(tidy(raw_model3),file="Results/Table1.csv",sep =
              ",", append=T,col.names = F)
lmm_result3=MSPEsumca(raw_model3,data = long_dat)
lmm_result3$sd_beta
lmm_result3$sd_sigma
long_dat$pollandstate=paste(long_dat$Abbrev,":", long_dat$pollster)
model3=plotall(lmm_result3, long_dat,"pollandstate")
model3$eblup_pp
model3$mspe_cycle

#Model 4
raw_model4=lmer(y~candidate*year+(0+candidate|year:state)+(0+1|pollster), data=long_dat)
summary(raw_model4)
write.table(tidy(raw_model4),file="Results/Table1.csv",sep =
              ",", append=T,col.names = F)
lmm_result4=MSPEsumca(raw_model4,data = long_dat)
lmm_result4$sd_beta
lmm_result4$sd_sigma
model4=plotall(lmm_result4, long_dat,"pollandstate")
model4$eblup_pp
model4$mspe_cycle


############################################Section 3.1: transfer learning
pred1=function(theta, polldata,fun){
  logp_bar= aggregate(y~state+candidate, data = polldata, FUN=fun)
  pi= exp(logp_bar[,3]-theta)
  return(pi)
}

adjust_model=function(train_data, test_data,formula,weights=NA){
  #weights need to be calculated as either 1/n_j or informative weights
train_model = lmer(formula,
               data=train_data)
train_data$ME=predict(train_model)
ME_train=unique(train_data[c("year","state","candidate","ME")]) 
ME_train=aggregate(ME~state+candidate,data=ME_train,mean)
test_data$pollrate=test_data$poll_dem
test_data$pollrate[test_data$candidate=="R"]=test_data$poll_rep[test_data$candidate=="R"]
theta=ME_train$ME
test_data$y=log(test_data$pollrate)
if(!any(is.na(weights))){
  test_data$y=weights*test_data$y
  p1=pred1(theta, test_data,fun=sum)#based on model
  test_data$weighted_pr=weights*test_data$pollrate
  p2=aggregate(weighted_pr~state+candidate, data = test_data, sum)[,3]#crude average
}else{
p1=pred1(theta, test_data,fun=mean)
p2=aggregate(pollrate~state+candidate, data = test_data, mean)[,3]
}
p3=aggregate(pollrate~state+candidate, data = test_data, mean)[,3]#crude average
test_data$realrate=test_data$real_dem
test_data$realrate[test_data$candidate=="R"]=test_data$real_rep[test_data$candidate=="R"]
real_test=unique(test_data[c("year","state","candidate","realrate")])
real_test$p1=p1
real_test$p2=p2
real_test$p3=p3
support=factor(levels=c("D", "R"),rep("D",length(states)))
callresult1=data.frame(state=states,call=support,
                       model_pred=support, adjusted_crude_pred=support,
                       crude_pred=support,
                       Votes=EV)
for (s in states){
  rs=real_test[real_test$state==s,]
  callresult1$call[callresult1$state==s]=rs$candidate[which.max(rs$realrate)]
  callresult1$model_pred[callresult1$state==s]=rs$candidate[which.max(rs$p1)]
  callresult1$adjusted_crude_pred[callresult1$state==s]=rs$candidate[which.max(rs$p2)]
  callresult1$crude_pred[callresult1$state==s]=rs$candidate[which.max(rs$p3)]
}
result=list()
result$full=real_test
result$correct_p1=sum(callresult1$model_pred==callresult1$call)
result$correct_p2=sum(callresult1$adjusted_crude==callresult1$call)
result$correct_p3=sum(callresult1$crude_pred==callresult1$call)
result$wrong=callresult1[(callresult1$model_pred!=callresult1$call)|
              (callresult1$adjusted_crude_pred!=callresult1$call)|
                (callresult1$crude_pred!=callresult1$call),]
result$mspe1=sum((p1-real_test$realrate)^2)
result$mspe2=sum((p2-real_test$realrate)^2)

result$sum0=aggregate(Votes~call,data=callresult1,sum)
result$sum1=aggregate(Votes~model_pred,data=callresult1,sum)
result$sum2=aggregate(Votes~crude_pred,data=callresult1,sum)
result$callresult=callresult1
return(result)
}


formula_1 = as.formula(y~candidate+(0+candidate|state))#Model 1
formula_3 = as.formula(y~candidate+(0+candidate|state)+(1|pollster))#Model 3
data16=long_dat[long_dat$year==2016,]
data20=long_dat[long_dat$year==2020,]

#Model 1
model1_adjust1=adjust_model(data16,data20,formula_1)
model1_adjust1$wrong

write.table(model1_adjust1$wrong, sep=",",file="Results/Table2.csv",
            row.names = FALSE)
model1_adjust1$sum0
model1_adjust1$sum1
model1_adjust1$sum2
model1_adjust1$callresult
dat1=model1_adjust1$full

model1_adjust2=adjust_model(data20,data16,formula_1)
model1_adjust2$wrong
write.table(model1_adjust2$wrong, sep=",",file="Results/Table2.csv",append = TRUE,
            row.names = FALSE)
model1_adjust2$sum0
model1_adjust2$sum1
model1_adjust2$sum2
model1_adjust2$callresult
dat2=model1_adjust2$full

#Model3
model3_adjust1=adjust_model(data16,data20,formula_3)
model3_adjust1$wrong
write.table(model3_adjust1$wrong, sep=",",file="Results/Table2.csv",append = TRUE,
            row.names = FALSE)
model3_adjust1$sum0
model3_adjust1$sum1
model3_adjust1$sum2
dat3=model3_adjust1$full

model3_adjust2=adjust_model(data20,data16,formula_3)
model3_adjust2$wrong
write.table(model3_adjust2$wrong, sep=",",file="Results/Table2.csv",append = TRUE,
            row.names = FALSE)
model3_adjust2$sum0
model3_adjust2$sum1
model3_adjust2$sum2
dat4=model3_adjust2$full

#########################################Section 3.2: Ranking the pollsters
##########1.Ranking pollster using random effects
####seperate
m3.16=lmer(formula_3,data=data16)
re1=ranef(m3.16)$pollster
re1 <- data.frame(pollster = row.names(re1), RE=re1[,1],RE_abs=abs(re1[,1]))
m3.20=lmer(formula_3,data=data20)
re2=ranef(m3.20)$pollster
re2 <- data.frame(pollster = row.names(re2), RE=re2[,1], RE_abs=abs(re2[,1]))
comb_re1=merge(re1,re2,by.x="pollster",by.y="pollster",
               all = T,suffixes = c(".2016",".2020"))
write.table(comb_re1,file="Results/re.csv",sep =
              ",", row.names = F)
comb_re2=merge(re1,re2,by.x="pollster",by.y="pollster",
            all = F,suffixes = c(".2016",".2020"))
write.table(comb_re2,file="Results/re_both.csv",sep =
              ",", row.names = F)
#####Combined
m3=as.formula(y~candidate*year+(0+candidate|state)+(1|pollster))#Model 3
m3.all=lmer(m3,data=long_dat)
re12=ranef(m3.all)$pollster
re12 <- data.frame(pollster = row.names(re12), RE=re12[,1], RE_abs=abs(re12[,1]))
write.table(re12,file="Results/comb_re.csv",sep =
              ",", row.names = F)

t3=data.frame(comb=re12$pollster[order(re12$RE_abs)][1:10],
           rank2016=comb_re1$pollster[order(comb_re1$RE_abs.2016)][1:10],
           rank2020=comb_re1$pollster[order(comb_re1$RE_abs.2020)][1:10])
write.table(t3,file="Results/Table3.csv",sep =
              ",", row.names = F)

t3_bottom=data.frame(comb=re12$pollster[order(re12$RE_abs,decreasing = T)][1:10],
              rank2016=comb_re1$pollster[order(comb_re1$RE_abs.2016,decreasing = T)][1:10],
              rank2020=comb_re1$pollster[order(comb_re1$RE_abs.2020,decreasing = T)][1:10])
write.table(t3_bottom,file="Results/Table3.csv",sep =
              ",", row.names = F,append = TRUE)

t5=data.frame(rank2016=comb_re2$pollster[order(comb_re2$RE_abs.2016)][1:10],
              rank2020=comb_re2$pollster[order(comb_re2$RE_abs.2020)][1:10])
write.table(t5,file="Results/Table5.csv",sep =
              ",", row.names = F)

t5_bottom=data.frame(rank2016=comb_re2$pollster[order(comb_re2$RE_abs.2016,decreasing = T)][1:10],
                     rank2020=comb_re2$pollster[order(comb_re2$RE_abs.2020,decreasing = T)][1:10])
write.table(t5_bottom,file="Results/Table5.csv",sep =
              ",", row.names = F,append = TRUE)

#############2.Ranking pollster using mixed effects
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
write.table(comb_me1,file="Results/me.csv",sep =
              ",", row.names = F)
head(comb_me1)
comb_me2=merge(p1,p2,by.x="pollster",by.y="pollster",
               all = F,suffixes = c(".2016",".2020"))
write.table(comb_me2,file="Results/me_both.csv",sep =
              ",", row.names = F)

###Combined
m3=as.formula(y~candidate*year+(0+candidate|state)+(1|pollster))#Model 3
me12=pollster_analysis(m3,long_dat)
write.table(me12,file="Results/comb_me.csv",sep =
              ",", row.names = F)

t4=data.frame(comb=me12$pollster[order(me12$ME_abs)][1:10],
              rank2016=comb_me1$pollster[order(comb_me1$ME_abs.2016)][1:10],
              rank2020=comb_me1$pollster[order(comb_me1$ME_abs.2020)][1:10])
write.table(t4,file="Results/Table4.csv",sep =
              ",", row.names = F)

t4_bottom=data.frame(comb=me12$pollster[order(me12$ME_abs,decreasing = T)][1:10],
                     rank2016=comb_me1$pollster[order(comb_me1$ME_abs.2016,decreasing = T)][1:10],
                     rank2020=comb_me1$pollster[order(comb_me1$ME_abs.2020,decreasing = T)][1:10])
write.table(t4_bottom,file="Results/Table4.csv",sep =
              ",", row.names = F,append = TRUE)

t6=data.frame(rank2016=comb_me2$pollster[order(comb_me2$ME_abs.2016)][1:10],
              rank2020=comb_me2$pollster[order(comb_me2$ME_abs.2020)][1:10])
write.table(t6,file="Results/Table6.csv",sep =
              ",", row.names = F)

t6_bottom=data.frame(rank2016=comb_me2$pollster[order(comb_me2$ME_abs.2016,decreasing = T)][1:10],
                     rank2020=comb_me2$pollster[order(comb_me2$ME_abs.2020,decreasing = T)][1:10])
write.table(t6_bottom,file="Results/Table6.csv",sep =
              ",", row.names = F,append = TRUE)



#############3.Use Ranks of pollster as weights to adjust predictions
weights_fun=function(data_EBLUP){
  #data_EBLUP includes pollsters+EBLUP
  EBLUP=data_EBLUP[,2]
  data_EBLUP$non_inform=1/mean(EBLUP)#non-informative prior
  data_EBLUP$inform=1/EBLUP
  return(data_EBLUP)
}

weights_fun_ranks=function(data_EBLUP){
  #data_EBLUP includes pollsters+EBLUP
  EBLUP=data_EBLUP[,2]
  data_EBLUP$inform=1/rank(EBLUP)
  data_EBLUP$non_inform=1/(mean(rank(EBLUP)))#non-informative prior
  return(data_EBLUP)
}

wts_scale=function(dat){
  total=aggregate(inform~state+candidate,FUN=sum,data=dat)
  total=rename(total,c("sum_w"="inform"))
  dat2=merge(dat,total,by=c("state","candidate"),all.x=TRUE)
  dat2$weights=dat2$inform/dat2$sum_w
  return(dat2)
}
#predict 2020 by 2016
weights_16=weights_fun(p1[c("pollster","ME_abs")])
weight_data=merge(data20,weights_16,by="pollster",all.x=TRUE)
weight_data$inform[is.na(weight_data$inform)]=weights_16$non_inform[1]
wt_df=wts_scale(weight_data)
weights=wt_df$weights[order(wt_df$row_id)]
model3_adjust1_up=adjust_model(data16,data20,formula_3,weights=weights)
model3_adjust1_up$wrong


weights_16=weights_fun_ranks(p1[c("pollster","ME_abs")])
weight_data=merge(data20,weights_16,by="pollster",all.x=TRUE)
weight_data$inform[is.na(weight_data$inform)]=weights_16$non_inform[1]
wt_df=wts_scale(weight_data)
weights=wt_df$weights[order(wt_df$row_id)]
model3_adjust1_up2=adjust_model(data16,data20,formula_3,weights=weights)
model3_adjust1_up2$wrong
write.table(model3_adjust1_up2$wrong,file="Results/Table7.csv",sep =
              ",", row.names = F)



#predict 2016 by 2020
weights_20=weights_fun(p2[c("pollster","ME_abs")])
weight_data=merge(data16,weights_20,by="pollster",all.x=TRUE)
weight_data$inform[is.na(weight_data$inform)]=weights_20$non_inform[1]
wt_df=wts_scale(weight_data)
weights=wt_df$weights[order(wt_df$row_id)]
model3_adjust1_up=adjust_model(data20,data16,formula_3,weights=weights)
model3_adjust1_up$wrong

weights_20=weights_fun_ranks(p2[c("pollster","ME_abs")])
weight_data=merge(data16,weights_20,by="pollster",all.x=TRUE)
weight_data$inform[is.na(weight_data$inform)]=weights_20$non_inform[1]
wt_df=wts_scale(weight_data)
weights=wt_df$weights[order(wt_df$row_id)]
model3_adjust1_up3=adjust_model(data20,data16,formula_3,weights=weights)
model3_adjust1_up3$wrong
write.table(model3_adjust1_up2$wrong,file="Results/Table7.csv",sep =
              ",", row.names = F,append = TRUE)



###############Time series analysis
##Repeat section 1 using long_dat_5w data preprocessed by preprocessing-week.
plotall_ylimit=function(lmm_result, data,labels="state"){
  #labels includes the text to print in EBLUP figures
  result=list()
  ####Prediction interval
  data$fitted=lmm_result$EBLUP
  data$MSPE=lmm_result$MSPE
  result$fit=data.frame(EBLUP=data$fitted, MSPE=data$MSPE)
  result$mspe_box=boxplot(data$MSPE, main="Distribution of MSPE")
  ct=qnorm(0.95)
  ##EBLUP
  result$eblup_pp=ggpaired(data, x = "candidate", y = "fitted",id="id",
                           xlab = FALSE,point.size = 0, ylim=c(-0.3,0.75),
                           ylab = "Predicted outcome",
                           color = "white", line.color = "gray", line.size = 0.4,
                           facet.by = "year", short.panel.labs = F,label=labels,
                           font.label = list(size = 8),
                           label.select = list(criteria="(`y` >0.2  & `x` == 'R')"))+
    scale_color_manual(values = c("#0073C2FF", "#FC4E07"))+geom_point(aes(color = candidate, shape=candidate), position = position_dodge(0.3))
  ##Plot by states
  mixed=data[c("year","state","candidate","fitted",
               "MSPE")]
  mixed=unique(mixed)
  result$mspe_states=ggplot(mixed, aes(year, fitted)) +
    geom_errorbar(
      aes(ymin = fitted-ct*sqrt(MSPE), ymax = fitted+ct*sqrt(MSPE), color = candidate),
      position = position_dodge(0.3), width = 0.2)+
    geom_point(aes(color = candidate), position = position_dodge(0.3)) +
    facet_wrap(~state)+
    scale_color_manual(values = c("#0073C2FF", "#FC4E07")) 
  mixed[,"state"]=as.numeric(mixed[,"state"])
  result$mspe_cycle=ggplot(mixed, aes(state, fitted)) +ylim(-0.4,0.8)+
    geom_errorbar(
      aes(ymin = fitted-ct*sqrt(MSPE), ymax = fitted+ct*sqrt(MSPE), color = candidate,
          linetype=candidate),
      position = position_dodge(0.3), width = 0.2)+
    geom_point(aes(color =candidate, shape=candidate), position = position_dodge(0.3)) +
    #facet_grid(rows = vars(year))+
    facet_wrap(~year)+
    scale_color_manual(values = c("#0073C2FF", "#FC4E07"))+
    ylab("Prediction Interval")
  return(result)
}
#Model 1

#write.table(tidy(raw_model1),file="data/Table1.csv",sep=",", append=T,col.names = F)
raw_model1=lmer(y~candidate*year+(0+candidate|state),
                data=long_dat_avg)

lmm_result1=MSPEsumca(raw_model1,data = long_dat_avg)
summary(lmm_result1$EBLUP-predict(raw_model1))
lmm_result1$sd_beta
lmm_result1$sd_sigma#Bootstrap estimators in Table 1
model1=plotall_ylimit(lmm_result1, long_dat_avg)
model1$eblup_pp#Figure 1
model1$mspe_cycle#Figure 5

lme_fit1 <- lme(y~candidate*year,random=~0+candidate|state,
                data=long_dat_avg,correlation =
                  corAR1(form=~week|state/year/candidate/pollster))

plot(ACF(lme_fit1),alpha=0.05)

lmm_result1.1=MSPEsumca(lme_fit1,data = long_dat_avg)
summary(lmm_result1.1$EBLUP-predict(lme_fit1))
lmm_result1.1$sd_beta
lmm_result1.1$sd_sigma#Bootstrap estimators in Table 1
model1.1=plotall_ylimit(lmm_result1.1, long_dat_avg)
model1.1$eblup_pp#Figure 1
model1.1$mspe_cycle#Figure 5

#Model 2
raw_model2=lmer(y~candidate*year+(0+candidate|year:state),data=long_dat_avg)
long_dat_avg$grp=interaction(long_dat_avg$year,long_dat_avg$state)
lme_fit2 <- lme(y~candidate*year,random=~0+candidate|grp,
                data=long_dat_avg,correlation =
                  corAR1(form=~week|grp/candidate/pollster))
lme_fit2
raw_model2
#state effect: random slope of candidate; nested in year
#write.table(tidy(raw_model2),file="data/Table1.csv",sep =
#              ",", append=T)
lmm_result2=MSPEsumca(raw_model2,data = long_dat_avg)
summary(lmm_result1$EBLUP-predict(raw_model2))
lmm_result2$sd_beta
lmm_result2$sd_sigma#Bootstrap estimators in Table 1
model2=plotall_ylimit(lmm_result2, long_dat_avg)
model2$eblup_pp#Figure 1
model2$mspe_cycle#Figure 5


lmm_result2.1=MSPEsumca(lme_fit2,data = long_dat_avg)
summary(lmm_result1$EBLUP-predict(lme_fit2))
lmm_result2.1$sd_beta
lmm_result2.1$sd_sigma#Bootstrap estimators in Table 1
model2.1=plotall_ylimit(lmm_result2.1, long_dat_avg)
model2.1$eblup_pp#Figure 1
model2.1$mspe_cycle#Figure 5
save.image(file='resultAR.RData')

