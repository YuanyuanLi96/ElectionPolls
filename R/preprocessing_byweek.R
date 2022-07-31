##------------------------------------------------
###Preprocessing the poll data and Election data
###Author: Yuanyuan Li
##----------------------------------------------------

######################Step 1: Download polls data from 538 website,
#find the polls conducted within 1 month ("start_date") before 
#2016/2020 election date, save as "polls_2016.csv" and "polls_2020.csv".

library(dplyr)
raw_polls_2016=read.csv("data/raw_data/raw_polls_2016.csv")
raw_polls_2016=raw_polls_2016[raw_polls_2016$type=="now-cast",]
raw_polls_2020=read.csv("data/raw_data/raw_polls_2020.csv")
polldate=function(data){
  data$start_date=as.Date(data$start_date, "%m/%d/%y")
  data$end_date=as.Date(data$end_date, "%m/%d/%y")
  data$poll_date=data$start_date + floor((data$end_date-data$start_date)/2)
  data$poll_date=as.Date(data$poll_date, "%m/%d/%y")
  return(data)
}
add_date_2016=polldate(raw_polls_2016)
add_date_2020=polldate(raw_polls_2020)

filter_poll=function(df,W=NULL){
  #W is the number of weeks before election
  df$election_date=as.Date(df$election_date, "%m/%d/%y")
  if(is.null(W)){
    cutoff=df$election_date[1]-31
    df=df[df$start_date>cutoff,]
  }else{
  cutoff=df$election_date[1]-7*W
  df=df[df$poll_date>cutoff,]
  }
  df=df[!df$state %in% c("U.S.","Maine CD-1","Maine CD-2","",
                         "Nebraska CD-1","Nebraska CD-2","Nebraska CD-3"),]
  df$state=droplevels(df$state)
  names(df)[names(df) == "cycle"] <- "year"
  return(df)
}

polls_2016=filter_poll(add_date_2016,5)
polls_2016=rename(polls_2016,
                c("poll_dem"="rawpoll_clinton", "poll_rep"="rawpoll_trump"))
polls_2016=polls_2016[,c("year","state","election_date","start_date","end_date",
                     "poll_date","pollster","sample_size","poll_dem","poll_rep")]
write.csv(polls_2016,file="data/polls_2016_5w.csv",row.names = F)

polls_2020=filter_poll(add_date_2020,5)
polls_2020=polls_2020[,c("question_id","year","state","election_date","start_date","end_date",
                     "poll_date","pollster","sample_size","candidate_party","pct")]
polls_2020=polls_2020[polls_2020$candidate_party%in%c("DEM","REP"),]
polls_2020=reshape(polls_2020, idvar=c("question_id","year","state","election_date","start_date","end_date","sample_size",
                                    "poll_date","pollster"), 
                 timevar = "candidate_party",
                 direction="wide")
polls_2020=rename(polls_2020,
                c("poll_dem"="pct.DEM", "poll_rep"="pct.REP"))
polls_2020=polls_2020[,c("year","state","election_date","start_date","end_date","sample_size",
                     "poll_date","pollster","poll_dem","poll_rep")]
write.csv(polls_2020,file="data/polls_2020_5w.csv",row.names = F)

##################Step 2: Find consistent names of same pollster
library("readxl")
polls_2016=read.csv("data/polls_2016_5w.csv")
polls_2020=read.csv("data/polls_2020_5w.csv")
election_2016=read.csv("data/raw_data/election_2016.csv")
election_2020=read.csv("data/raw_data/election_2020.csv")
pollster_rate=read.csv("data/raw_data/pollsters538.csv")#pollsters in 538 website
levels(polls_2016$state)
levels(polls_2020$state)

recom=function(data,newnames){
  polls.cha=as.character(unique(data$pollster[which(data$pollster %in% newnames==F)]))
  newpolls=lapply(polls.cha,function(x)strsplit(x, "\\, | |\\/")[[1]])
  pollsters=lapply(newnames, function(x)strsplit(x, "\\, | |\\/")[[1]])
  recommand=list()
  for (i in 1:length(newpolls)){
    rec=NULL
    for (j in 1:length(pollsters)){
      if (all(newpolls[[i]] %in% c(pollsters[[j]],"Inc.", "Co.","Company"))){
        rec=c(rec, newnames[j])
      }
    }
    name_length=sapply(rec,nchar)
    rec=rec[which.min(name_length)]
    if(!is.null(rec)){
    recommand[[length(recommand)+1]]=list(pollster=polls.cha[i], rec=rec)
    }
  }
  return(recommand)
}
levels(polls_2016$pollster)
newnames=unique(as.character(pollster_rate$Pollster))
recommand1=recom(polls_2016,newnames)
recommand2=recom(polls_2020,newnames)
repl=function(data,recommand){
 for(h in 1:length(recommand)){
   levels(data$pollster)[levels(data$pollster)==
                          recommand[[h]]$pollster]=recommand[[h]]$rec
  }
  return(data)
}
repl1=repl(polls_2016,recommand1)
repl2=repl(polls_2020,recommand2)

#Extract the nearest poll using poll_date
filter=function(data){
  filtered=NULL
  states=levels(data$state)
  for (s in states){
    data_s=data[data$state==s,]
    pollster_s=unique(data_s$pollster)
    for (p in pollster_s){
      data_sp=data_s[data_s$pollster==p,]
      filtered=rbind(filtered,data_sp[which.max(data_sp$poll_date),])
    }
  }
  return(filtered)
}

#repl1=filter(repl1)
#repl2=filter(repl2)

#cleaned polls data
raw_polls=rbind(repl1,repl2)
dim(raw_polls)
nlevels(raw_polls$pollster)


################################Step 3: Add real election results
#add real results
n1=length(which(raw_polls$year==2016))
raw_polls$real_dem=rep(0,nrow(raw_polls))
raw_polls$real_rep=rep(0,nrow(raw_polls))
ds1= raw_polls[1:n1,]
ds2=raw_polls[(n1+1):nrow(raw_polls),]
states=levels(election_2016$state)
for (s in states){
  ds1$real_dem[ds1$state==s]= election_2016$ratio_dem[election_2016$state==s]
  ds1$real_rep[ds1$state==s]= election_2016$ratio_rep[election_2016$state==s]
  ds2$real_dem[ds2$state==s]= election_2020$dem_percent[election_2020$state==s]
  ds2$real_rep[ds2$state==s]= election_2020$rep_percent[election_2020$state==s]
}

all_p= rbind(ds1,ds2)
nrow(all_p)
all_p$logp_dem=log(all_p$poll_dem/all_p$real_dem)
all_p$logp_rep=log(all_p$poll_rep/all_p$real_rep)
long_dat=reshape(all_p, varying =c("logp_dem","logp_rep"), 
                 v.names = "logp",timevar = "candidate",
                 times=c("D","R"),
                 direction="long")
head(long_dat)
dim(long_dat)#1686*15
#long_dat=long_dat[order(long_dat$id),]
long_dat$year=factor(long_dat$year)

#candidate: 1=democratic, 2= republican
#long_dat$candidate=factor(as.numeric(factor(long_dat$candidate)))
head(long_dat)
names=colnames(long_dat)
colnames(long_dat)=names
colnames(long_dat)
polls=levels(long_dat$pollster)
long_dat$poll_date=as.Date(long_dat$poll_date)
statenames=election_2016[,c("state","Abbrev")]
long_dat$row_id=1:nrow(long_dat)
long_dat=merge(long_dat,statenames, by="state")
long_dat=long_dat[order(long_dat$row_id),]
long_dat$week=ceiling(difftime(long_dat$election_date,long_dat$poll_date, 
                             units = "weeks"))  
hist(as.numeric(long_dat$week))
write.table(long_dat, file ="data/long_dat_5w.csv", sep=",",
            row.names=FALSE, col.names=TRUE)

