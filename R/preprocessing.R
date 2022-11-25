##----------------------------------------------------------------------------------
###Preprocessing the poll data and Election data
###Author: Yuanyuan Li

#This file consists two parts: 1. support functions; 2.preprocess data
#The output files are "model_data.csv", "model_data_5w.csv".
##----------------------------------------------------------------------------------

################################Part 1. Support Functions####################
# Transform date from character to "Date" class
polldate=function(data){
  data$start_date=as.Date(data$start_date, "%m/%d/%y")
  data$end_date=as.Date(data$end_date, "%m/%d/%y")
  data$poll_date=data$start_date + floor((data$end_date-data$start_date)/2)
  data$poll_date=as.Date(data$poll_date, "%m/%d/%y")
  return(data)
}

# Extract state-wide polls happened within [W] weeks before election date 
filter_poll=function(df,W=NULL){
  #W is the number of weeks before election
  df$election_date=as.Date(df$election_date, "%m/%d/%y")
  if(is.null(W)){
    #If W is null, use 1 month data
    cutoff=df$election_date[1]-31
    df=df[df$start_date>cutoff,]
  }else{
  cutoff=df$election_date[1]-7*W
  df=df[df$poll_date>cutoff,]
  }
  #remove non-statewide poll results
  df=df[!df$state %in% c("U.S.","Maine CD-1","Maine CD-2","",
                         "Nebraska CD-1","Nebraska CD-2","Nebraska CD-3"),]
  df$state=droplevels(df$state)
  names(df)[names(df) == "cycle"] <- "year"
  return(df)
}

# Replace the pollsters'names in the polls by the names in the pollster data. 
# This step is to avoid the inconsistencies of polllster names from different data sets. 
pollster_name_rec=function(data,pollster_names){
  polls.cha=as.character(unique(data$pollster[which(data$pollster %in% pollster_names==F)]))
  newpolls=lapply(polls.cha,function(x)strsplit(x, "\\, | |\\/")[[1]])
  pollsters=lapply(pollster_names, function(x)strsplit(x, "\\, | |\\/")[[1]])
  recommand=list()
  for (i in 1:length(newpolls)){
    rec=NULL
    for (j in 1:length(pollsters)){
      if (all(newpolls[[i]] %in% c(pollsters[[j]],"Inc.", "Co.","Company"))){
        rec=c(rec, pollster_names[j])
      }
    }
    name_length=sapply(rec,nchar)
    rec=rec[which.min(name_length)]
    if(!is.null(rec)){
      recommand[[length(recommand)+1]]=list(pollster=polls.cha[i], rec=rec)
    }
  }
  # replace the pollster names in poll data to the recommended names in pollster data
  for(h in 1:length(recommand)){
    levels(data$pollster)[levels(data$pollster)==
                            recommand[[h]]$pollster]=recommand[[h]]$rec
  }
  return(data)
}

# Extract the latest polls per state-pollster combination
state_pollster_latest=function(data){
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


# Add real election results to polls data
add_actual_results=function(polls,election){
  polls$real_dem=rep(0,nrow(polls))
  polls$real_rep=rep(0,nrow(polls))
  states=levels(election$state)
  for (s in states){
    polls$real_dem[polls$state==s]= election$ratio_dem[election$state==s]
    polls$real_rep[polls$state==s]= election$ratio_rep[election$state==s]
  }
  return(polls)
}


################################Part 2. Proprecess data####################
setwd('..')
getwd()#workind directory should be the parent folder
library(dplyr)
library("readxl")
###### Read in datasets #####
raw_polls_2016=read.csv("data/raw_data/raw_polls_2016.csv")
raw_polls_2020=read.csv("data/raw_data/raw_polls_2020.csv")
election_2016=read.csv("data/raw_data/election_2016.csv")
election_2020=read.csv("data/raw_data/election_2020.csv")
pollster_rate=read.csv("data/raw_data/pollsters538.csv")#pollsters in 538 website
pollster_names=unique(as.character(pollster_rate$Pollster))
state_names=election_2016[,c("state","Abbrev")]

###### Combine polls data preprocess steps #####
# Clean 2016 raw poll data
polls_clean_2016=function(raw_polls_2016,pollster_names,W=NULL,latest=TRUE){
  # 1. Transform date from character to "Date" class
  polls_2016 = polldate(raw_polls_2016)
  # 2. remove unused fields
  polls_2016=polls_2016[polls_2016$type=="now-cast",]
  polls_2016=rename(polls_2016,
                    c("poll_dem"="rawpoll_clinton", "poll_rep"="rawpoll_trump",
                      "year"="cycle"))
  polls_2016=polls_2016[,c("year","state","election_date","start_date","end_date",
                           "poll_date","pollster","sample_size","poll_dem","poll_rep")]
  
  # 3. Extract state-wide polls happened within [W] weeks before election date 
  polls_2016 = filter_poll(polls_2016, W)
  # 4. Replace the pollsters'names in the polls by the names in the pollster data
  polls_2016 = pollster_name_rec(polls_2016,pollster_names)
  # 5. Extract the latest polls per state-pollster combination
  if(latest){
  polls_2016 = state_pollster_latest(polls_2016)
  }
  return(polls_2016)
}

# Clean 2020 raw poll data
polls_clean_2020=function(raw_polls_2020,pollster_names,W=NULL,latest=TRUE){
  # 1. Transform date from character to "Date" class
  polls_2020=polldate(raw_polls_2020)
  # 2. remove unused fields
  polls_2020=rename(polls_2020, c("year"="cycle"))
  polls_2020=polls_2020[,c("question_id","year","state","election_date","start_date","end_date",
                               "poll_date","pollster","sample_size","candidate_party","pct")]
  polls_2020=polls_2020[polls_2020$candidate_party%in%c("DEM","REP"),]
  # 3. Extract state-wide polls happened within [W] weeks before election date 
  polls_2020=filter_poll(polls_2020, W)
  # 4. Reshape 2020 data to the same format with 2016 data
  polls_2020=reshape(polls_2020, idvar=c("question_id","year","state","election_date","start_date","end_date","sample_size",
                                         "poll_date","pollster"), 
                     timevar = "candidate_party",
                     direction="wide")
  polls_2020=rename(polls_2020,
                    c("poll_dem"="pct.DEM", "poll_rep"="pct.REP"))
  polls_2020=polls_2020[,c("year","state","election_date","start_date","end_date",
                           "poll_date","pollster","sample_size","poll_dem","poll_rep")]
  # 5. Replace the pollsters'names in the polls by the names in the pollster data
  polls_2020 = pollster_name_rec(polls_2020,pollster_names)
  # 6. Extract the latest polls per state-pollster combination
  if(latest){
  polls_2020 = state_pollster_latest(polls_2020)
  }
  return(polls_2020)
}

###### Get 1 month polls data and add actual election results #####
polls_2016=polls_clean_2016(raw_polls_2016, pollster_names)
polls_2020=polls_clean_2020(raw_polls_2020, pollster_names)

# Add real election results to polls data
ds1 = add_actual_results(polls_2016, election_2016)
ds2 = add_actual_results(polls_2020, election_2020)
all_data = rbind(ds1,ds2)

# Add response variable, weeks feature
model_data_func=function(all_data, state_names){
all_data$logp_dem=log(all_data$poll_dem/all_data$real_dem)
all_data$logp_rep=log(all_data$poll_rep/all_data$real_rep)
long_dat=reshape(all_data, varying =c("logp_dem","logp_rep"), 
                 v.names = "logp",timevar = "candidate",
                 times=c("D","R"),
                 direction="long")
long_dat=long_dat[order(long_dat$id),]
long_dat$week=ceiling(difftime(long_dat$poll_date,long_dat$election_date,
                               units = "weeks"))  
long_dat$row_id=1:nrow(long_dat)#use row_id to prevent order changing by merge
long_dat=merge(long_dat,state_names, by="state")
long_dat=long_dat[order(long_dat$row_id),]
return(long_dat)
}
model_data=model_data_func(all_data, state_names)
dim(model_data)#1634*18
write.table(model_data, file ="data/model_data.csv", sep=",",
            row.names=FALSE, col.names=TRUE)

###### Get 5-weeks polls data and add actual election results #####
polls_2016_5w=polls_clean_2016(raw_polls_2016, pollster_names, W=5,latest = FALSE)
polls_2020_5w=polls_clean_2020(raw_polls_2020, pollster_names, W=5,latest = FALSE)

# Add real election results to polls data
ds1 = add_actual_results(polls_2016_5w, election_2016)
ds2 = add_actual_results(polls_2020_5w, election_2020)
all_data_5w = rbind(ds1,ds2)
# Add response variable, weeks feature
model_data_5w=model_data_func(all_data_5w, state_names)
dim(model_data_5w)#8244*18
write.table(model_data_5w, file ="data/model_data_5w.csv", sep=",",
            row.names=FALSE, col.names=TRUE)

