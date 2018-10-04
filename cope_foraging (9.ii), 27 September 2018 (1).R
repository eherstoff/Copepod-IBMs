#######################################################################
# Last updated: 2 October 2018, EMH and MU
#Individual based model of adult A. tonsa foraging in patches that differ in stoichiometric quality
#Parameterize using data from EMH's 2016 Helgoland experiments. Focus adults only, looking at movement/behavior by precondition type
#######################################################################
#
#FOCAL QUESTIONS:
#(1) how long copepods take to locate high-quality food within a patch from an overall matrix of low-quality food, and 
#(2) how well do copes stay there? 
# ... and then tie this back to length of time patches are known to persist in the real world so we can determine if the values
#the model outputs for patch size/length of time are at all relevant to what's possible
#
#
#######################################################################
#Setting up parameters for the simulation for the movement of individuals through the grid #####
nearq<-function(x){ #This function determines the nearest 1/16 for #s in calculations
  ans<-vector(length=length(x))
  for (i in 1:length(x)){
    a<-floor(x[i]/0.0625)
    if((x[i]-(0.0625*a))<0.03125){
      ans[i]<-0.0625*a
    }
    else{
      ans[i]<-0.0625*(a+1)
    }
  }
  return(ans)
}



decimalplaces <- function(x) { #This function determines decimal places to store/return
  if ((x %% 1) != 0) {
    nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed=TRUE)[[1]][[2]])
  } else {
    return(0)
  }
}

gcd<-function(x,y){ #This function determines the greatest common divisor?
  x1<-x*10^(max(decimalplaces(x),decimalplaces(y)))
  y1<-y*10^(max(decimalplaces(x),decimalplaces(y)))
  M<-max(x1,y1)
  N<-min(x1,y1)
  if((M-N*floor(M/N))>0){
    i<-0
    while(i==0){
      R<-M-N*floor(M/N)
      if(R>0){
        M<-N
        N<-R
      }
      else{
        return(N/(10^(max(decimalplaces(x),decimalplaces(y)))))
        i<-1
      }
    }
  }
}

#######################################################################
#
#Parameters about movement pulled from experiment, and standardized for movement in model #######
########################################################################
#Standardize the time units based on loop times in different quality foods
#Standardize the cell unit based on swimming speeds in different quality foods
#
#Input the time spent performing a loop in high quality (hq) or low-quality (lq) food. Times are in ms #######
#From "R code for movement model params, 17 Sept 2018

#This is taken from #GradSchool/Summer2016/Swimming analysis for EMH/Helix time budgeting, EMH experiments.xls
#This is the time spent swimming in a loop (in ms), divided by the frame rate I watched (every 10th frame, shot at 1/1250 fps)
helix.time.hq <- (c(0.08,0.032,0.032,0.112,0.096,0.104,0.048,0.056,0.184,0.192,0.032,0.048,0.048,0.048,0.032,0.032))/(10/1250)
helix.time.lq <- (c(0.056,0.032,0.072,0.032,0.032,0.032,0.064,0.032,0.04))/(10/1250)

loop.hq <- sample(x=helix.time.hq, size=1000, replace=T)
loop.lq <- sample(x=helix.time.lq, size=1000, replace=T)

#Find the standardized time unit in milliseconds
time.unit<-gcd(nearq(mean(loop.hq)),nearq(mean(loop.lq)) )

#How long (ms) do they stay in a cell while doing a loop #######
#This will be used to tell how many steps the simulation needs to keep a cope in a cell when helices happen
time.hq<-loop.hq/time.unit
time.lq<-loop.lq/time.unit 

#######################################################################
#Input the speeds under different treatments #####Speeds are in millimeters/ms  #######
#Take stepwise velocity data from this file:
dat.2 <- read.csv("/Users/emilypetchler/Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/Calc params/_Model input data for active swimmers, 18 Sept. 2018.csv")
f2.dat2 <- dat.2[c(1:1115),]     #all adults precon with f2
noN.dat2 <- dat.2[c(1116:1920),] #all adults precon with noN

speed.hq <- round(nearq(sample(x=f2.dat2$velocity, size=1000, replace = T)))
speed.lq <- round(nearq(sample(x=noN.dat2$velocity, size=1000, replace = T))) 

#Find the minimum cell size in mm  #######
cell.size<-gcd(nearq(mean(speed.hq)),nearq(mean(speed.lq)))

#How many cells do they travel in one time unit? #######
jump.hq<- time.unit* speed.hq/cell.size
jump.lq<- time.unit* speed.lq/cell.size

#Approximate the jump sizes again by rounding to the nearest whole number  #
jump.hq<-round(jump.hq) 
jump.lq<-round(jump.lq)
#######################################################################
#Frequency of swimming movement, taken from adult data  ###############
#This is the overall probability of helical swimming (but I have things broken down by behavior type too!)

#This is taken from #GradSchool/Summer2016/Swimming analysis for EMH/Helix time budgeting, EMH experiments.xls (see master active swimming tab)
#This is the ratio of total video time spent performing helical swimming
ratio.helix.hq <- c(0.294117647, 0.081632653,0.095238095,0.424242424,0.48,0.419354839,0.193548387,0.212121212,0.469387755,0.444444444,0.181818182,0.272727273,0.3,0.162162162,0.111111111,0.19047619)
ratio.helix.lq <- c(0.225806452, 0.129032258,0.147540984,0.173913043,0.070175439,0.105263158,0.380952381,0.06557377,0.131578947)

#This is the ratio of all *other* active swimming, excluding helical swimming behaviors.
ratio.otheractive.hq <- c(0.294117647,0.346938776,0.738095238,0.454545455,0.4,0.548387097,0.64516129,0.484848485,0.244897959,0.444444444,0.681818182,0.590909091,0.65,0.513513514,0.638888889,0.666666667)
ratio.otheractive.lq <- c(0.548387097,0.548387097,0.360655738,0.52173913,0.280701754,0.315789474,0.333333333,0.37704918,0.184210526)

#This is the ratio of time spent doing nothing/waiting
ratio.waiting.hq <- 1-(ratio.helix.hq+ratio.otheractive.hq)
ratio.waiting.lq <- 1-(ratio.helix.lq+ratio.otheractive.lq)

prob.hq <- nearq(sample(size=1000, x=(1-ratio.waiting.hq), replace=T)) #probability of doing any movement in HQ food
prob.lq <- nearq(sample(size=1000, x=(1-ratio.waiting.lq), replace=T)) #

#######################################################################
#Swimming angles ######################################################
#Angle probabilities: 
#at-90:-72, -72:-54, -54:-36, -36:-18, -18:0, _0_, 0:18, 17:36, 36:54, 54:72, 72:90 (in degrees)
#Angles at _0_ indicate small movement in Y (or Z) with no X movement.
#From excel sheet, "/Users/emilypetchler/Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/Calc params/adult angle calcs 9 July 2018/XY_angle.dat, or XZ_angle.dat tabs).
angles.HQ_XY <- 100*c(0.027855153,
                      0.080779944,
                      0.105849582,
                      0.105849582,
                      0.057103064,
                      0.249303621,
                      0.061281337,
                      0.100278552,
                      0.111420613,
                      0.077994429,
                      0.019498607) 
angles.LQ_XY <- 100*c(0.022964509,
                      0.077244259,
                      0.087682672,
                      0.106471816,
                      0.08559499,
                      0.217118998,
                      0.075156576,
                      0.114822547,
                      0.104384134,
                      0.073068894,
                      0.03131524)



angles.HQ_XZ <- 100*c(0.016713092,
                      0.055710306,
                      0.093314763,
                      0.079387187,
                      0.048746518,
                      0.406685237,
                      0.048746518,
                      0.077994429,
                      0.086350975,
                      0.058495822,
                      0.025069638)
angles.LQ_XZ <- 100*c(0.03131524,
                      0.054279749,
                      0.077244259,
                      0.077244259,
                      0.050104384,
                      0.438413361,
                      0.045929019,
                      0.06263048,
                      0.060542797,
                      0.039665971,
                      0.058455115)

#angles.null <- 100*c(.090909, 090909, 090909, 090909, 090909, 090909, 090909, 090909, 090909, 090909, 090909) #equal probability in all directions

#Average angles in radians (11 total angles possible)
angles<-c(-1.413716694, -1.099557429, -0.785398163, -0.471238898, -0.157079633,
          0,
          0.157079633, 0.471238898, 0.785398163, 1.099557429, 1.413716694) #Translated above/correct angles WRT pi
#avg. angles are c(-81,-63,-45,-27,-9, 0, 9, 27, 45, 63, 81) IN DEGREES
#######################################################################
#Running the IBM
# # # # # # # # # # # # # # # # # # # # 
#Form a grid space filled with LQ food  #####
grid.size <- 1000
grid2<-matrix(rep(0,grid.size^2),nrow=grid.size,ncol=grid.size) #the grid is created and filled with LQ food

#Add in HQ patches  #####

#grid2[100:250,100:250]<-1    # 9.1204 % has HQ food in 500 x 500 px
#grid2[100:400,100:300]<-1    # 32.6284 % has HQ food in 500 x 500 px
#grid2[75:425,75:425]<-1    #49.2804 % has HQ food in 500 x 500 px

grid2[200:550,200:550]<-1    # 12.3201 % has HQ food in 1000 x 1000 px




#Visualizing the HQ and LQ patches ####
library(SparseM)
image((grid2), col=c("steelblue","darkgoldenrod1")) # plot it; HQ food is gold :)

##############
#Setting up matrix######
timepoints <- 30000  #time steps in the model.
num.copes <- 500     #number of copes to simulate in the model
res.i<-matrix(ncol=num.copes, nrow=timepoints) #cols are replicates, rows are timepoints. **number of columns needs to match i in Individual based simulations below.

getdim<-function(x,m){
  rows<-nrow(m)
  cols<-ncol(m)
  ydim<-ceiling(x/rows)
  xdim<-x-(ydim-1)*rows
  return(cbind(xdim,ydim))
}

getind<-function(d,m){
  return((d[2]-1)*nrow(m)+d[1])
}
######
#Individual based simulations; NOTE that i in start 'for' must match ncol in res.i in 'Setting up matrix'  #####
ind<-sample(1:11,1) #sample from the distribution of angles in "swimming angles" above

for(i in 1:num.copes){  #number of replicates; needs to match ncol in 'Setting up matrix'
  refx<-c(0) #reference angle starting at = 0
  init<-c(sample(1:nrow(grid2),1),sample(1:ncol(grid2),1))  #Initiate from a random cell
  res.i[1,i]<-getind(init,grid2)
  j<-2
  a<-angles-refx[length(refx)]
  while(j<timepoints){ #time steps (must match time steps in 'Setting up matrix' above)  
    if(grid2[init[1],init[2]]==0){  #looping in LQ food, stay put for set number of timesteps
      if(runif(1)>mean(prob.lq)){ #draw from the distribution of prob.lq, sample with replacement
        sample.x <- sample(time.lq, size=1)
        res.i[j:min(((j+sample.x)-1),timepoints),i]<-getind(init,grid2)
        j<-j+sample.x
      }
       else{  #angles and movement for LQ
        sets<-rbind(c(cos(a[0]),sin(a[0])), c(cos(a[1]),sin(a[1])),c(cos(a[2]),sin(a[2])),c(cos(a[3]),sin(a[3])),c(cos(a[4]),sin(a[4])),c(cos(a[5]),sin(a[5])),
                    c(cos(a[6]),sin(a[6])), c(cos(a[7]),sin(a[7])), c(cos(a[8]),sin(a[8])), c(cos(a[9]),sin(a[9])), c(cos(a[10]),sin(a[10])), c(cos(a[11]),sin(a[11]))  )
        
        ind<-sample(1:11,1,prob=angles.LQ_XY)    #turning angle changes depending on food quality
        #ind<-sample(1:11,1,prob=angles.null)   #turning angle remains the same (equal in all directions), regardless of food quality
        
        x<-init[1]+floor(sample(jump.lq, size=1)*sets[ind,1])
        y<-init[2]+floor(sample(jump.lq, size=1)*sets[ind,2])
        if(x==0){x<-x-1}
        if(y==0){y<-y-1}
        if(x%in%1:nrow(grid2) & y%in%1:ncol(grid2)){init<-c(x,y)}
        else{init<-c(x%%nrow(grid2),y%%ncol(grid2))}
        res.i[j,i]<-getind(init,grid2)
        j<-j+1
        refx<-c(refx,angles[ind]-(pi*0.5-refx[length(refx)]))
      }
    }
    
    # Second part of the code #######
    if(grid2[init[1],init[2]] == 1){   #when they are looping in HQ, stay put for set number of timesteps
      if(runif(1)>mean(prob.hq)){
        sample.x <- sample(time.lq, size=1)
        res.i[j:min((j+sample.x-1),timepoints),i]<-getind(init,grid2)
        j<-j+sample.x
      }
      else{  #angles and movement for HQ
        sets<-rbind(c(cos(a[1]),sin(a[1])),c(cos(a[2]),sin(a[2])),c(cos(a[3]),sin(a[3])),c(cos(a[4]),sin(a[4])),c(cos(a[5]),sin(a[5])),
                    c(cos(a[6]),sin(a[6])), c(cos(a[7]),sin(a[7])), c(cos(a[8]),sin(a[8])), c(cos(a[9]),sin(a[9])), c(cos(a[10]),sin(a[10])), c(cos(a[11]),sin(a[11]))  )
        
        ind<-sample(1:11,1,prob=angles.HQ_XY)

        x<-init[1]+floor(sample(jump.hq, size=1)*sets[ind,1])
        y<-init[2]+floor(sample(jump.hq, size=1)*sets[ind,2])
        if(x==0){x<-x-1}
        if(y==0){y<-y-1}
        if(x%in%1:nrow(grid2) & y%in%1:ncol(grid2)){init<-c(x,y)}
        else{init<-c(x%%nrow(grid2),y%%ncol(grid2))}
        res.i[j,i]<-getind(init,grid2)
        j<-j+1
        refx<-c(refx,angles[ind]-(pi*0.5-refx[length(refx)]))
        
      }
    }
  }
}

##############
#After the IBM has run, pulling data together.
##############
#Plotting the mean of the X and Y coordinates in space that were visited by the cope runs ######
res.i<-res.i[-max(timepoints),] #take off the last row with NA's

meanx<-c()
varyx<-c()
meany<-c()
varyy<-c()
for(i in 1:nrow(res.i)){
  meanx<-c(meanx,mean(getdim(res.i[i,],grid2)[,1]))
  varyx<-c(varyx,var(getdim(res.i[i,],grid2)[,1]))
  
  meany<-c(meany,mean(getdim(res.i[i,],grid2)[,2]))
  varyy<-c(varyy,var(getdim(res.i[i,],grid2)[,2]))
}
par(mfrow=c(2,2))
plot(varyx)
plot(varyy)
plot(meany~meanx, type="l")
par(mfrow=c(1,1))

#Plot first trajectory #####
res.k<-matrix(nrow=nrow(res.i),ncol=2)
for(i in 1:nrow(res.i)){
  res.k[i,]<-getdim(res.i[i,6],grid2)
}
plot(res.k[,2]~res.k[,1],type="l",xlim=c(0,grid.size),ylim=c(0,grid.size)) 




#Combine all the replicates #####
nres<-c()
for(i in 1:ncol(res.i)){
  nres<-c(nres,res.i[,i])
}


#Frequency in each cell (that was visited atleast once) #####
freq<-as.data.frame(table(nres))

freq[,2]<-freq[,2]/length(nres)

#Find the time spent in HQ patches vs. LQ patches #####
patches<-which(grid2==1)
nonpatches<-which(grid2==0)

food.freq<-freq[which(freq$nres%in%patches),]
food.freq1<-as.numeric(food.freq[,1])
food.freq2<-as.numeric(food.freq[,2])
food.freq<-cbind(food.freq1,food.freq2)

nonfood.freq<-freq[which(freq$nres%in%nonpatches),]
nonfood.freq1<-as.numeric(nonfood.freq[,1])
nonfood.freq2<-as.numeric(nonfood.freq[,2])
nonfood.freq<-cbind(nonfood.freq1,nonfood.freq2)  

result.food<-matrix(nrow=nrow(food.freq),ncol=3)
for(i in 1:nrow(food.freq)){
  result.food[i,1:2]<-getdim(food.freq[i,1],grid2)
}
result.food[,3]<-food.freq[,2]
colnames(result.food)<-c("x","y","density")

result.nonfood<-matrix(nrow=nrow(nonfood.freq),ncol=3)
for(i in 1:nrow(nonfood.freq)){
  result.nonfood[i,1:2]<-getdim(nonfood.freq[i,1],grid2)
}
result.nonfood[,3]<-nonfood.freq[,2]
colnames(result.nonfood)<-c("x","y","density")


df.food<-as.data.frame(result.food)
colnames(df.food)<-c('x',"y","density")
df.nonfood<-as.data.frame(result.nonfood)
colnames(df.nonfood)<-c('x',"y","density")




#######################################################################

#TO DO: take means of 0 and 1's to show frequency of helix vs swim vs rest... 
#Need to show that what I put in to behavior/model is what comes out.
#Get net vs total displacement? Show/animate a trajectory...?
#(Increase run time and space if it all looks OK...)
#Heatmap of where seeded, and where final position is for copes in model.
#count HQ vs LQ cells in during loops... 


#######################################################################
#Run a t-test to see if there's differnces in copepod density in HQ vs LQ patches
t.test(x=df.food$density, y=df.nonfood$density)
#######################################################################
#Plotting data
#######################################################################
#Basic params for plotting #####
#Frequency in each cell (that was visited atleast once) ######
freq<-as.data.frame(table(nres))
freq[,2]<-100*(freq[,2]/length(nres)) #*100 for %time in each cell

#Find the time spent in food HQ patches vs. LQ patches #####
patches<-which(grid2==1)
nonpatches<-which(grid2==0)

food.freq<-freq[which(freq$nres%in%patches),]
nonfood.freq<-freq[which(freq$nres%in%nonpatches),]

#
#Plotting means and SEs #########
test.concat <-c(mean(food.freq$Freq), mean(nonfood.freq$Freq)) 
test.error <- c(sd(food.freq$Freq)/sqrt(length(food.freq$Freq)), sd(nonfood.freq$Freq)/sqrt(length(nonfood.freq$Freq)))
plotdat <- barplot(test.concat, ylim=c(0,(mean(nonfood.freq$Freq)+0.005)), col=c("darkgoldenrod1", "steelblue"), cex.axis=1.35, cex.lab=1.5,
                   xlab=c("Hi qual patches       Low qual patches"), ylab="Average frequency in patch type")
abline(h=0, lwd=2) #zero line
arrows(plotdat, test.concat-test.error, plotdat, test.concat, angle=90, code=1, length=0.05, lwd=2)
arrows(plotdat, test.concat+test.error, plotdat, test.concat, angle=90, code=1, length=0.05, lwd=2)

########################################################################
#
#Plotting frequency by raw data #########
par(mfrow=c(1,2))
stripchart(food.freq$Freq, method="jitter", jitter=0.05, vertical=T, main="frequency in cells with HQ food", ylim = c(0,0.15), col="steelblue")
stripchart(nonfood.freq$Freq, method="jitter", jitter=0.05, vertical=T, main="frequency in cells with LQ food", ylim = c(0,0.15), col="darkgoldenrod1")
par(mfrow=c(1,1))

########################################################################


########################################################################
#Time within patches
########################################################################
# Do copes find a patch, and how long do they stay in/around it?
#create 0-1 data, compare 0 vs 1? (do they find the patch)
#once 1, how many 0 vs 1's? (how long are they in it?)
########################################################################
##This is what % of patches within the grid had food (where hq=1)
patch.food = length(which(grid2==1)) / 
                 (length(which(grid2==1))+length(which(grid2==0)))  

########################################################################
# Looking at ability to find/stay within HQ vs LQ patches ##################################
########################################################################

#Residence time within HQ patches? ##################################
res.j<-matrix(nrow=nrow(res.i),ncol=ncol(res.i))

for(i in 1:ncol(res.i)){
  for(j in 1:nrow(res.i)){
   res.j[j,i]<-grid2[res.i[j,i]]
  }
}
#residence times for each trial
res.time<-apply(res.j,2,sum)/timepoints #divided by # steps run
hist(res.time, xlab="% time spent in HQ patches", col="steelblue", 
     main="Residence time within a HQ patch", breaks=10, xlim=c(0,1))
  #histogram of residence time within the patches
abline(v=mean(res.time, na.rm=T), lty=2, lwd=2) #This line is the mean of the time spent in patches
abline(v=(patch.food), col="chocolate1", lwd=2) #This line is what % of the area in the model has hq food available

mean(res.time, na.rm=T)
sd(res.time, na.rm=T)

#####
#determine when they first find a patch ##########
ones<-function(a){
  one<-sum(a!=1)
  if(one<length(a)){return(min(which(a==1)))}
  else(return(length(a)))
}

first.visit<-apply(res.j,2,ones)
hisfirst.visit[!is.finite(first.visit)] <- 0 #replace any INF values with 0's

hist(first.visit, xlab="Time steps", col="steelblue", xlim=c(0,timepoints),
     main="Number of time steps to first locate a patch", breaks=10) 
  #at which step do copes first enter a patch?
abline(v=mean(first.visit), lty=2, lwd=2) #This line is the mean of the time spent in patches
abline(v=(patch.food*timepoints), col="chocolate1", lwd=2) #This line is what % of the area in the model has hq food available
# NOT SURE THE Orange line above is correct! It's the % of the model filled with food, * the length of the run for each cope.

mean(first.visit)
sd(first.visit)

#####
#patch occupany at any one time point ##############
patchougue<-apply(res.j,1,sum)/1000

plot(patchougue, pch=20, col="steelblue", xlab="Steps through time", ylab="% of copepods within a patch (out of 1)") 
  # % of copes within a patch throughout time

hist(patchougue, col="steelblue", xlab="% of copepods within a patch (out of 1)", main="Patch occupancy throughout time", breaks=10, xlim=c(0,1)) 
  #frequency distribution of % copes in a patch throughout time
abline(v=mean(patchougue, na.rm=T), lty=2, lwd=2) #This line is the mean of the time spent in patches
abline(v=patch.food,  col="chocolate1", lwd=2) #putting a vertical line at % of patches within the grid
#if the # of copes in the patch is > % patches with food, get a +# and show copes stay longer in a patch than outside of it.


mean(patchougue, na.rm=T)
sd(patchougue, na.rm=T)
#####
#once they find a HQ patch, how long are they there? #####
afterstay<-c()
starts<-c()
for(i in 1:ncol(res.j)){
  one<-ones(res.j[,i])
  if(one<length(res.j[,i])){
  afterstay<-c(afterstay,length(which(res.j[((one+1):nrow(res.j)),i]==1)))
  starts<-c(starts,one)}
}

par(mfrow=c(2,1))
hist(afterstay/((timepoints-1)-starts), col="aliceblue")
summary(afterstay/((timepoints-1)-starts)) #ratio of time (out of 1) that copes spend in HQ food once they find it

hist(afterstay, col="aliceblue")
summary(afterstay) #this is the number of steps spend in HQ food across the 10,000 steps of the model
#this is the overall time in HQ food (not ness. one long visit.. )

par(mfrow=c(1,1))
########################################################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 




##