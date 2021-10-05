require(graphics)
#some data
#d <- density(faithful$eruptions, bw = "sj")
files <- list.files(path="C:/Users/Halima/Desktop/phd_script/data1/", pattern="*.txt")
name<-c()
track<-list()
di<-list()
my.list=list()
nb=6 #data_dimension
pt<-c(1, 916,917, 1230,1231, 1361,1362, 2043,2044, 2418,2419, 3180)

par(mfrow=c(3,2))

#One-dimensional segmentation (working on each datatrack independently)
for (i in 1:nb) {

  d<- as.numeric(read.table(paste("C:/Users/Halima/Desktop/phd_script/data1/",files[i],sep = "")))
  ts_y<-ts(d)
  require(pastecs)
  tp<-turnpoints(ts_y)
  plot(d,type = "b")
  write(tp$tppos,"C:/Users/Halima/Desktop/phd_script/top_2.txt", ncolumns = 15 )
  difference<-diff(tp$tppos)
  num<- as.numeric(append(difference,-1,0))
  concatenate<-cbind(tp$tppos,d[tp$tppos],num)
  fram<-as.data.frame(x = concatenate)
  colnames(fram)<-c("indice","value","difference")
  
  fram<-fram[order(fram$difference,fram$indice,decreasing = TRUE),]
  h<-head(fram,19)
  f<-h[order(h$indice,decreasing = FALSE),]
  fram$mean<-as.integer(fram$indice-fram$difference/2)
  f$mean<-as.integer(f$indice-f$difference/2)
  inter<-matrix(pt, 20,byrow = TRUE)
  intervale <-as.data.frame(inter)[-1,]#
  
  
  #intervale2<-intervale[-1]
  #intervale <- matrix(intervale2, 9,byrow = TRUE)
  f$interval<-intervale#$V1
  points(f$mean,d[f$mean],col="green",type = "h")
  track[[i]]<-fram$mean
  di[[i]]<-fram$difference
  my.list[[i]]=fram
  name<-c(name,paste("f", i, sep = ""))
}

#find occurence in data, sorting data
freq<-sapply(di,table)
ss<-lapply(freq, function(i) data.frame(i[-1]))
ss2<-sapply(ss, function(i) i$logg<-log(i$Freq))

for(j in 1:nb){ss[[j]]$log<-log(ss[[j]]$Freq)}
for(j in 1:nb){ss[[j]]$Var1<-as.integer(as.vector(ss[[j]]$Var1))}
dii<-sapply(di,function(i) sort( i , decreasing = FALSE))

#P-Val to detect cut-points
pval<-sapply(dii,function(i) quantile( i , c(0.5,0.7,0.9, 0.95, 0.99,0.995, 0.996,0.997,0.998,0.999, 0.9999) ))
cut_pt<-lapply(ss, function(i) match(0,i$log))

#Multi-D segmentation aggregation
sapply(1:nb,function(i) plot(as.integer(as.vector(ss[[i]]$Var1)),ss[[i]]$Freq, type = "p",col=ifelse(as.integer(as.vector(ss[[i]]$Var1)) < pval[,i]['99.5%'], "blue", "red")) )
sapply(1:nb,function(i) {plot(ss[[i]]$Var1,log(ss[[i]]$Freq), type = "p",col=ifelse(ss[[i]]$Var1 < ss[[i]][cut_pt[[i]],]$Var1, "blue", "red"));abline(glm(log ~ Var1,data = ss[[i]][ss[[i]]$Var1 < ss[[i]][cut_pt[[i]],]$Var1,]), col="green");abline(glm(log ~ Var1,data = ss[[i]][ss[[i]]$Var1 > ss[[i]][cut_pt[[i]],]$Var1,]), col="yellow") })

selected_pt<-lapply(1:nb, function(i) ss[[i]][ss[[i]]$Var1 >= ss[[i]][cut_pt[[i]]-1,]$Var1,])
selected_indice<-lapply(1:nb,function(i) my.list[[i]]$mean[my.list[[i]]$difference>=ss[[i]][cut_pt[[i]]-1,]$Var1])
selected_indice<-lapply(1:nb, function(i) sort(selected_indice[[i]]))

#plotting multi-D segmentation
col <- rgb(runif(nb+1),runif(nb+1),runif(nb+1))
par(mfrow=c(1,1))
plot(selected_indice[[1]],rep(1,length(selected_indice[[1]])),xlab = "DNA position", ylab = "DataTracks",ylim=range(1:(nb+1)),col=col[1],xaxt="n")

sapply(2:nb,function(i) {lines(selected_indice[[i]],rep(i,length(selected_indice[[i]])),col=col[i],type = "p") })
cd<-pt
cd<-c(1,cd[c(FALSE,TRUE)])
length(cd)<-length(cd)-1
lines(cd,rep((nb+1),length(cd)),col=col[1],type = "h")

################################################
length(selected_indice[[which.max(lapply(selected_indice, function(i) length(i)))]])
occurence<-data.frame(ind=numeric(),occ=integer(),ens=list())
i<-0
count<- 0
nbb=0
stop<-FALSE

while (i < length(selected_indice[[which.max(lapply(selected_indice, function(i) length(i)))]])){
  #  for (i in 1:length(track[[1]])) {
  i=i+1
   
  count<-0
  for (j in 1:nb) {
    # print(paste0(j,"_________",i))
    if (!is.na(selected_indice[[j]][i])){
    if(length(occurence$ind[occurence$ind %in% c((selected_indice[[j]][i]-10):(selected_indice[[j]][i]+10))])!=0){
      print(paste0(selected_indice[[j]][i],"__jjjjjjj__",occurence$ind[occurence$ind %in% c((selected_indice[[j]][i]-10):(selected_indice[[j]][i]+10))]))
      #if(track[[j]][i] %in% occurence$ind){
      #print(paste0(j,"hg"))
      occurence$occ[occurence$ind==occurence$ind[occurence$ind %in% c((selected_indice[[j]][i]-10):(selected_indice[[j]][i]+10))]]<-occurence$occ[occurence$ind==occurence$ind[occurence$ind %in% c((selected_indice[[j]][i]-10):(selected_indice[[j]][i]+10))]]+1
      occurence$ens[occurence$ind==occurence$ind[occurence$ind %in% c((selected_indice[[j]][i]-10):(selected_indice[[j]][i]+10))]]<-list(c(occurence$ens[occurence$ind==occurence$ind[occurence$ind %in% c((selected_indice[[j]][i]-10):(selected_indice[[j]][i]+10))]],selected_indice[[j]][i]))
      count=count+1
      
    }
    else{
      
      occurence<-rbind(occurence,data.frame(ind=selected_indice[[j]][i],occ=1,ens=list(ens=c(selected_indice[[j]][i]))))
      #occurence$)
      
    }}
  }
  print(count)
}
seg<-occurence[occurence$occ >= (nb/2-1),]

bc<-sapply(1:length(seg$ind),function(i) {as.integer(mean(unlist(seg$ens[[i]]))) })
seg$res<-bc

bb<-sapply(1:length(seg$ind),function(i) {unique(seg$res[[i]])})
seg$must_occ<-bb
#sapply(freq,function(i) plot(i[-1], type = "b") )

my.list2<-sapply(1:length(my.list),function(i) my.list[[i]]$mean[my.list[[i]]$difference >= pval[,i]['99.5%']] )

prop.table(freq[[1]][-1])


library(vioplot)
sapply(dii,function(i) vioplot(i[-1]))
ccc<-sapply(freq,function(i) vioplot(i[-1]))
l<-lapply(freq, function(i) as.data.frame(i[-1]))

viop<- vioplot(freq[[2]][-1])

occurence<-data.frame(ind=numeric(),occ=integer(),ens=list())
i<-0
count<- 0
nb=0
stop<-FALSE

while (i < length(track[[1]]) && stop==FALSE){
  #  for (i in 1:length(track[[1]])) {
  i=i+1
  print(paste0("_________",i))
  
  count<-0
  for (j in 1:length(files)) {
    # print(paste0(j,"_________",i))
    if(length(occurence$ind[occurence$ind %in% c((track[[j]][i]-5):(track[[j]][i]+5))])!=0){
      # print(paste0(i,"__jjjjjjj__"))
      #if(track[[j]][i] %in% occurence$ind){
      #print(paste0(j,"hg"))
      occurence$occ[occurence$ind==occurence$ind[occurence$ind %in% c((track[[j]][i]-5):(track[[j]][i]+5))]]<-occurence$occ[occurence$ind==occurence$ind[occurence$ind %in% c((track[[j]][i]-5):(track[[j]][i]+5))]]+1
      occurence$ens[occurence$ind==occurence$ind[occurence$ind %in% c((track[[j]][i]-5):(track[[j]][i]+5))]]<-list(c(occurence$ens[occurence$ind==occurence$ind[occurence$ind %in% c((track[[j]][i]-5):(track[[j]][i]+5))]],track[[j]][i]))
      count=count+1
      # occurence<-rbind(occurence,data.frame(ind=32,occ=3))
      
      
    }
    else{
      
      occurence<-rbind(occurence,data.frame(ind=track[[j]][i],occ=1,ens=list(ens=c(track[[j]][i]))))
      #occurence$)
      
    }
  }
  if(count==0){
    nb=nb+1
  }
  
  if(nb>=2){
    print("halimaaaaaa")
    print(i)
    stop=TRUE
  }
  print(count)
}
seg<-occurence$ind[occurence$occ > (length(files)/2)]

sapply(freq,function(i) plot(i[-1], type = "b") )


