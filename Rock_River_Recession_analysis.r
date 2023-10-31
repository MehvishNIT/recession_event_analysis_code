
rockproject<-read.csv("path to file")## import file to r 
rockproject$'y(mm/d)'<-rockproject$ Q.ft.3.*.0000989 # unit conversion

library(lubridate) # package for spliting date into sub columns
rockproject$year <- year(mdy(rockproject$DATE))
rockproject$month <- formatC(month(mdy(rockproject$DATE)), width=2, flag=0)# fancy operations to make each date 2 digit
rockproject$day <- formatC(day(mdy(rockproject$DATE)), digits=1, flag=0) 

rockproject$time<-paste(rockproject$month,rockproject$day,sep="_" ) #seperate date by - 
rockproject<-rockproject[rockproject$time >="04_15" & rockproject$time <="10_31",] #getting summer months  


Y<-rockproject$'y(mm/d)'# naming variable
mdy_dt<-array(data=NA, dim = c("row-length",1)) #loop for operation1 = -dy/dt
for (i in (2:(length(Y)-1)) ) {
  mdy_dt[i] <- (Y[i-1]-Y[i+1])/2  
}
rockproject<-cbind(rockproject,mdy_dt)

## getting annual lowest streamflows
mav<-function(x,n){filter(x, rep(1/n,n), sides=2)} #function for moving average
DMA7<-mav(rockproject$`y(mm/d)`,7)   # here we have 7day moving average
DMA7<-as.numeric(as.character(DMA7))
rockproject<-cbind(rockproject, DMA7)
yL7<-aggregate(rockproject$DMA7, by = list(rockproject$year),FUN = min,na.rm=T)

## code for B08 algorithm

q99 <- quantile(rockproject$`y(mm/d)`, probs = 0.99)
#(condition 1) identify all rising limb data  
cond1 <- which(rockprojects$mdy_dt <= 0)
#(condtion2)is two points before first rising discharge and(condition3)is 3 points after peak of hydrograph 
cond123 <- unique( unlist (lapply(cond1, FUN = function(x) seq(from=(x-2), to = (x+3), by=1))))
# identify above 99 percentile,major events
ind_major <- which(rockproject$`y(mm/d)` >= q99)
# condition 4, remove four points after major events
cond4 <- unique( unlist (lapply(ind_major, FUN = function(x) seq(from=(x), to = (x+4), by=1))))
rockproject$mdy_dt[1:(length(cond4)-2)] = "NA" # condition 4 added two extra days
#conditon 5, remove anomalous decreases in baselfow
cond5 <- which(rockproject$mdy_dt[1:(nrow(rockproject)-1)] <= 0.5*rockproject$mdy_dt[2:nrow(rockproject)])+1
cond6 <- (which(rockproject$mdy_dt[1:(nrow(rockproject)-1)] <= rockproject$mdy_dt[2:nrow(rockproject)] ))
all_conds <- unique(c(ind_major, cond123, cond4, cond5, cond6 ))## identify all above points uniquely
rockproject$mdy_dt[all_conds] = NA # dropping for further analysis. 

plot(rockproject$'Y(mm/day)'[!all_conds],rockproject$'-dy/dt'[!all_conds],xlab ='Y', ylab ='-dy/dt',log='xy') 
abline(rq(rockproject$'-dy/dt'[!all_conds]~ rockproject$'Y(mm/day)'[!all_conds],tau=.05))

## our agorithm
load("~exact.list")# list containing daily date,discharge,rainfall data
##to chk continue zero rain days(<2mm) for min 12 days length for study station,finding start and end of events 
new.df = function(i) {
  
  x <- rle(exact.list[[i]][[1]][["Rainfall"]])
  mycode = which(x[["values"]] <= 2 & x[["lengths"]] >= 12)
  ends = (cumsum(rle(exact.list[[i]][[1]][["Rainfall"]])$"lengths")[mycode])
  
  newindex = ifelse(mycode > 1, mycode - 1, 0)
  starts = cumsum(rle(exact.list[[i]][[1]][["Rainfall"]])$"lengths")[newindex] + 1
  if (0 %in% newindex){
    starts = c(1, starts)
  }
  m = matrix(data = NA,
             nrow = length(starts),
             ncol = 2)
  m = data.frame(starts, ends)
  
  return(list(m))
}

## write this in list format, list of sequence of indices following <2mm for min 12 days

final.list=list(){
  final.list[i] = new.df(i)
}

##extract all other variables from above indices and remove first 5 days to each event of <2mm rain for min 12 days 

dryevents = function(A) {
  tryflow = list()
  for (i in 1:nrow(final.list[[A]])) {
    tryflow[[i]] <-sapply(exact.list[[A]][[1]],"[",final.list[[A]]$starts[i]:final.list[[A]]$ends[i])
    tryflow[[i]] <- tryflow[[i]][-c(1:5), ]
  }
  return(list(tryflow))
}

## save values of recession events
DRY_EVT=list(){
  DRY_EVT[i]=dryevents(i)
}

combined_list <- do.call(c, DRY_EVT)# this makes a common list of all sub-events 

cond=list()
cond1=list()
cond[[i]] <- which(combined_list [[i]][,4][2:length(combined_list [[i]][,4])] <= combined_list [[i]][,4][2:nrow(combined_list [[i]])-1])
cond1[[i]]=(cond[[i]]+1)
reqd<-which(lengths(cond1)>=7)# finding length less than 7 obtained from cond1
filtered2_list<-combined_list [reqd]## has length min 7 ,apply criteria 2b on this list
cond2<-cond1[reqd]#

##extract decreasing flow events only from above extracted events  

outputs <- list()
for (k in 1:length(cond2)) {
  counter <- 1
  output.name <- names(cond2)[[k]]
  for (i in 2:length(cond2[[k]])) {
    if (cond2[[k]][i] == (cond2[[k]][i-1])+ 1) {
      counter <- counter + 1
    } else {
      if (counter >= 7) {
        outputs[[length(outputs) + 1]] <- cond2[[k]][(i-counter):(i-1)]
        names(outputs)[length(outputs)] <- output.name
      }
      counter <- 1
    }
  }  
  if (counter >= 7) {
    outputs[[length(outputs) + 1]] <- cond2[[k]][(length(cond2[[k]])-counter+1):length(cond2[[k]])]
    names(outputs)[length(outputs)] <- output.name
  }
}

names(filtered2_list) <- seq_along(filtered2_list)
names(cond2) <- seq_along(cond2)


# Add a predecessor to each element of the list bcz code drops first day of event 
outputs <- lapply(outputs, function(x) c(x[1]-1, x))

## to pick events that are decreasing and have length min of 7 days
filtered3=list()
Z=names(outputs)
Z2=names(filtered2_list)
for (i in 1:length(Z)) {
  if (Z[i] %in% Z2) {
    print(c(Z[i],outputs[[i]]))
    filtered3[[i]]<-filtered2_list[[Z[i]]][outputs[[i]],]
  }
}
names(filtered3)=Z

## get K values for seperated events
k_fun=function(i){
  if(all(filtered3[[i]][,4]>0)){
    logy <-log(filtered3[[i]][,4])
    logy = matrix(logy,nrow = length(logy),ncol = 1) %>% as.data.frame()
    mt = matrix(1:nrow(logy),ncol = 1)
    model <- lm(logy$V1 ~ mt)
    k <- -1/coef(model)[2]
    return(k)} else {0}
}

k.matrix=matrix(data=NA, nrow="nrow",ncol=1)
for (i in 1:"nrow"){
  k.matrix[i]=round(k_fun(i))
}
