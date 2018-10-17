#######################################################################
# Last updated: 15 October 2018, EMH and MU
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
#
#WHAT FIGURES SHOULD BE PRODUCED FOR THE PAPER?
#Residence and location params --> injestion information.
#
#Plot size of patch vs. how long to find; relate to net benefit of being in patch (vs randomly feeding)
#Relate to conditions where you'd actually see this...
#Fraction of time in HQ vs LQ patch; how fast to find based on patch size...
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
#at 0, 0-36, 36-72, 72-108, 108-144, 144-180,216-252, 252-288, 288-324, 324-360, NA's (in degrees)
#Angles at _0_ indicate small movement in Y (or Z) with no X movement.
#From excel sheet, "/Users/emilypetchler/Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/Calc params/cope angles, 5 October 2018, emh and mu.xlsx
#     (look at the HQ.dat and LQ.dat tabs).

#Angle_XY	%s (out of 1) for HQ
angles.HQ_XY <- 100*c(0.060099133, 0.039653036, 0.042750929, 0.00929368, 0.038413879, 0.037174721, 
                      0.083023544, 0.041511772, 0.011771995, 0.044609665, 0.035315985) #NA's =  0.555142503%


#Angle_XY	%s (out of 1) for LQ
angles.LQ_XY <- 100*c(0.037037037, 0.03250189, 0.024943311, 0.003779289, 0.021164021, 0.037037037, 
                      0.077853364, 0.039304611, 0.015873016, 0.038548753, 0.03250189) #NA's = 0.637944067



#Angle_XZ	%s (out of 1) for HQ
angles.HQ_XZ <- 100*c(0.09417596, 0.04708798, 0.071251549, 0.12267658, 0.03283767, 0.010532838,
                      0.015489467, 0.041511772, 0.008054523, 0, 0) #NA's = 	0.555142503


#Angle_XZ	%s (out of 1) for LQ
angles.LQ_XZ <- 100*c(0.077853364, 0.040060469, 0.048374906, 0.108843537, 0.022675737, 0.015873016, 
                      0.012093726, 0.025699169, 0.009070295, 0, 0)  #NA's =	0.637944067

#angles.null <- 100*c(.090909, 090909, 090909, 090909, 090909, 090909, 090909, 090909, 090909, 090909, 090909) #equal probability in all directions

#Average angles in radians (11 total angles possible)
a<-c(0, 0.314159265, 0.942477796, 1.570796327, 2.199114858, 2.827433388,
     3.455751919,4.08407045, 4.71238898, 5.340707511, 5.969026042)  #Translated above/correct angles WRT pi
#avg. angles are c(0, 18, 54, 90, 126, 162, 198, 234, 270, 306, 342) IN DEGREES
#######################################################################
#Running the IBM
# # # # # # # # # # # # # # # # # # # # 
#Form a grid space filled with LQ food  #####
grid.size <- 1000
grid2<-matrix(rep(0,grid.size^2),nrow=grid.size,ncol=grid.size) #the grid is created and filled with LQ food

#Add in HQ patches  #####
#grid2[400:750,400:750]<-1    # 12.3201 % has HQ food in 1000 x 1000 px; HQ food patch in the middle of the model.
grid2[0:360,0:360]<-1    # 12.96 % has HQ food in 1000 x 1000 px; HQ food patch at the corner of the model.
#grid2[0:180,0:180]<-1    # 3.24 % has HQ food in 1000 x 1000 px; HQ food patch at the corner of the model.
#grid2[000:90,000:90]<-1    # 0.81 % has HQ food in 1000 x 1000 px; HQ food patch at the corner of the model.




##############
#Null model test: parameters for HQ and LQ are the same
#Cancel this bit when running the actual simulation
#angles.HQ_XY<-angles.LQ_XY
#prob.hq<-prob.lq
#jump.hq<-jump.lq
#loop.hq<-loop.lq

##############
#Setting up matrix######
timepoints <- 100000  #time steps in the model. Change back to 500,000?
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
#Individual based simulation  #####
ind<-sample(1:11,1) #sample from the distribution of angles in "swimming angles" above
sets<-rbind(c(cos(a[0]),sin(a[0])), c(cos(a[1]),sin(a[1])),c(cos(a[2]),sin(a[2])),c(cos(a[3]),sin(a[3])),c(cos(a[4]),sin(a[4])),c(cos(a[5]),sin(a[5])),
            c(cos(a[6]),sin(a[6])), c(cos(a[7]),sin(a[7])), c(cos(a[8]),sin(a[8])), c(cos(a[9]),sin(a[9])), c(cos(a[10]),sin(a[10])), c(cos(a[11]),sin(a[11]))  )

#Set the counter for behaviors. Columns are the copepods, rows are at each time step.
count<-matrix(nrow=timepoints,ncol=num.copes)    #System: 0 for straight runs and 1 for loops; NA= waiting during loops.
qual<-matrix(nrow=timepoints,ncol=num.copes)     #0 for LQ and 1 for HQ

for(i in 1:num.copes){  #number of replicates; needs to match ncol in 'Setting up matrix'
  init<-c(sample(1:nrow(grid2),1),sample(1:ncol(grid2),1))  #Initiate from a random cell
  res.i[1,i]<-getind(init,grid2)
  count[1,i]<-grid2[getind(init,grid2)]
  qual[1,i]<-1
  j<-2
  
  while(j<timepoints){ #time steps (must match time steps in 'Setting up matrix' above)  
    if(grid2[init[1],init[2]]==0){  #looping in LQ food, stay put for set number of timesteps
      if(runif(1)<mean(prob.lq)){ #draw from the distribution of prob.lq, sample with replacement
        count[j,i]<-1  #update the behavior counter
        sample.x <- floor(sample(time.lq, size=1))
        res.i[j:min(((j+sample.x)-1),timepoints),i]<-getind(init,grid2)
        qual[j:min(((j+sample.x)-1),timepoints),i]<-0
        j<-j+sample.x
        
      }
      else{  #angles and movement for LQ
        count[j,i]<-0  #Update the counter
        ind<-sample(1:11,1,prob=angles.LQ_XY)    #turning angle changes depending on food quality
        jump<-sample(jump.lq,size=1)
        x<-init[1]+floor(jump*sets[ind,1])
        y<-init[2]+floor(jump*sets[ind,2])
        if(x%in%1:nrow(grid2) & y%in%1:ncol(grid2)){init<-c(x,y)}
        else{init<-c(x%%nrow(grid2),y%%ncol(grid2))}
        if(init[1]==0){init[1]<-1}
        if(init[2]==0){init[2]<-1}
        res.i[j,i]<-getind(init,grid2)
        j<-j+1
        qual[j,i]<-0
      }
    }
    
    # Second part of the code #######
    if(grid2[init[1],init[2]] == 1){   #when they are looping in HQ, stay put for set number of timesteps
      if(runif(1)<mean(prob.hq)){
        count[j,i]<-1  #update the behavior counter
        sample.x <- floor(sample(time.lq, size=1))
        res.i[j:min((j+sample.x-1),timepoints),i]<-getind(init,grid2)
        qual[j:min(((j+sample.x)-1),timepoints),i]<-1
        j<-j+sample.x
        
      }
      else{  #angles and movement for HQ
        count[j,i]<-0  #Update the counter
        ind<-sample(1:11,1,prob=angles.HQ_XY)
        jump<-sample(jump.hq,size=1)
        x<-init[1]+floor(jump*sets[ind,1])
        y<-init[2]+floor(jump*sets[ind,2])
        if(x%in%1:nrow(grid2) & y%in%1:ncol(grid2)){init<-c(x,y)}
        else{init<-c(x%%nrow(grid2),y%%ncol(grid2))}
        if(init[1]==0){init[1]<-1}
        if(init[2]==0){init[2]<-1}
        res.i[j,i]<-getind(init,grid2)
        j<-j+1
        qual[j,i]<-1
      }
    }
  }
}
res.i<-res.i[-nrow(res.i),]
res.i<-res.i[-1,]
##############

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

#Residence time within HQ patches? (contains res.j, for % time in HQ vs LQ patchs.)   ##################################
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
hist(afterstay/((timepoints-1)-starts), col="aliceblue", breaks=20)
summary(afterstay/((timepoints-1)-starts)) #ratio of time (out of 1) that copes spend in HQ food once they find it

hist(afterstay, col="aliceblue", breaks=20)
summary(afterstay) #this is the number of steps spend in HQ food across the 10,000 steps of the model
#this is the overall time in HQ food (not ness. one long visit.. )

par(mfrow=c(1,1))
########################################################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 




#After the IBM has run, pulling data together. ##############
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

#######################################################################
#Plotting the mean of the X and Y coordinates in space that were visited by the cope runs ######
par(mfrow=c(2,2))
plot(varyx)
plot(varyy)
plot(meany~meanx, type="l", xlim=c(0,grid.size), ylim=c(0,grid.size)) #mean final locations
library(SparseM) #Visualizing the HQ and LQ patches #
image((grid2), col=c("steelblue","darkgoldenrod1")) # plot where HQ and LQ food is. HQ food is gold. LQ is blue.
par(mfrow=c(1,1))

#Plot the trajectory of the first copepod in the model #####
library(plot3D) #will plot in color

res.k<-matrix(nrow=nrow(res.i),ncol=2) #the mean of all of the copepod XY points.
for(i in 1:nrow(res.i)){
  res.k[i,]<-getdim(res.i[i,3],grid2)
}

#Plot the trajectory as B&W, and as color
par(mfrow=c(1,2))

plot(x=res.k[,2], y=res.k[,1],type="l",xlim=c(0,grid.size),ylim=c(0,grid.size), main="B&W movement plot") 

scatter2D(x=res.k[,2], y=res.k[,1], col = rainbow(n = length(res.k[,2])),
          xlim=c(0,grid.size),ylim=c(0,grid.size), type="l", main="Colorized movement", sub="Reds=start points; Blues=end points") 

par(mfrow=c(1,1))



#Create a numeric to hold which points were visited during the movement of all copes in the model run.
nres<-c()
for(i in 1:ncol(res.i)){
  nres<-c(nres,res.i[,i])
}


#######################################################################
#Plot the mean % of copepods in HQ vs. LQ over time #####
#use the data on 'count' from IBM, which says whether the animal is in HQ or LQ food
means<-c() 
for(i in 1:nrow(res.j)){ a
  means<-c(means,sum(res.j[i,])/length(res.j[1,]))
}
plot(means,type="l", ylim=c(0,1)) #average occupancy in patches over time.
abline(h=patch.food, col="red") #this line is the % of the model filled with HQ food

#with 12.96% HQ in corner, average of 44.95998 % copes in model are in HQ food through time.
#with 3.24% HQ in corner, average of  14.92647 % copes in model are in HQ food through time.

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
# concatinating data on count and qual to determine displacement and feeding rates ######
# Matrix of 'count' (0 for straight runs and 1 for loops; NA= waiting during loops); and 'qual (0 for LQ and 1 for HQ patch)
#Determine how many points copes spend in HQ vs LQ food #########
helix.LQ.dat<-matrix(nrow=num.copes, ncol=1)
helix.HQ.dat<-matrix(nrow=num.copes, ncol=1)
helix.LQ_na.dat<-matrix(nrow=num.copes, ncol=1)
helix.HQ_na.dat<-matrix(nrow=num.copes, ncol=1)
straight.LQ.dat <- matrix(nrow=num.copes, ncol=1)
straight.HQ.dat <- matrix(nrow=num.copes, ncol=1)

for (i in 1:num.copes){
  helix.LQ.dat[i,]= length(which(count[,i]==1&qual[,i]==0)) #number of helix in LQ food
  helix.HQ.dat[i,]= length(which(count[,i]==1&qual[,i]==1)) #number of helix in HQ food
  
  helix.LQ_na.dat[i,]= length(which(is.na(count[,i])&qual[,i]==0)) #number of helix in LQ food
  helix.HQ_na.dat[i,]= length(which(is.na(count[,i])&qual[,i]==1)) #number of helix in HQ food
  
  
  straight.LQ.dat[i,]= length(which(count[,i]==0&qual[,i]==0)) #number of straight runs in LQ food
  straight.HQ.dat[i,]= length(which(count[,i]==0&qual[,2]==1)) #number of straight runs in HQ food
  i=i+1
}

#Looking at hist of data on helical swimming by patch type ##########
par(mfrow=c(2,1))
hist((helix.LQ.dat+helix.LQ_na.dat)/timepoints, col="aliceblue") #number of timepoints helixing in LQ food
abline(v=(patch.food), col="chocolate1") #add a vertical line at % space with HQ food.
hist((helix.HQ.dat+helix.HQ_na.dat)/timepoints, col="aliceblue") #number of timepoints helixing in LQ food
abline(v=(patch.food), col="chocolate1")
par(mfrow=c(1,1))


#######################################################################
#Barplot of % time helical swimming in LQ vs HQ food. #######
hel.dat <- c(mean((helix.LQ.dat+helix.LQ_na.dat)/timepoints), mean(helix.HQ.dat+helix.HQ_na.dat)/timepoints)
hel.err <- c(sd((helix.LQ.dat+helix.LQ_na.dat)/timepoints)/sqrt(num.copes),
             sd((helix.HQ.dat+helix.HQ_na.dat)/timepoints)/sqrt(num.copes))

hel.plot <- barplot(hel.dat, ylim=c(0,1), xlab=c("LQ patch     HQ patch"), ylab="% time",
                    sub="% Time helical swimming", col=c("steelblue","goldenrod1"))
abline(h=0)
arrows(hel.plot, hel.dat, hel.plot, hel.dat+hel.err, length = 0, lwd=2)
arrows(hel.plot, hel.dat, hel.plot, hel.dat-hel.err, length = 0, lwd = 2)


#Run a t-test to see if there's differnces in behavior in HQ vs LQ patches
t.test(x=(helix.LQ.dat+helix.LQ_na.dat), y=(helix.HQ.dat+helix.HQ_na.dat))




#Barplot Num of helices in LQ vs HQ food. ####### 
hel.dat <- c(mean((helix.LQ.dat)), mean(helix.HQ.dat))
hel.err <- c(sd((helix.LQ.dat))/sqrt(num.copes),sd((helix.HQ.dat))/sqrt(num.copes))

hel.plot <- barplot(hel.dat, ylim=c(0, max(hel.dat+hel.err)), xlab=c("LQ patch     HQ patch"), ylab="Number of helices performed",
                    sub="Number of helices performed in each food type", col=c("steelblue","goldenrod1"))
abline(h=0)
arrows(hel.plot, hel.dat, hel.plot, hel.dat+hel.err, length = 0, lwd=2)
arrows(hel.plot, hel.dat, hel.plot, hel.dat-hel.err, length = 0, lwd = 2)

time.unit*hel.dat
cell.size*hel.dat

#Run a t-test to see if there's differnces in behavior in HQ vs LQ patches
t.test(x=(helix.LQ.dat), y=(helix.HQ.dat))


#######
#Barplot of % straight swimming in LQ vs HQ food. #######
str.dat <- c(mean((straight.LQ.dat)/timepoints), mean(straight.HQ.dat)/timepoints)
str.err <- c(sd((straight.LQ.dat)/timepoints)/sqrt(num.copes),
             sd((straight.HQ.dat)/timepoints)/sqrt(num.copes))

str.plot <- barplot(str.dat, ylim=c(0,.01), xlab=c("LQ patch     HQ patch"), ylab="% time",
                    sub="% Time straight swimming", col=c("steelblue","goldenrod1"))
abline(h=0)
arrows(str.plot, str.dat, str.plot, str.dat+str.err, length = 0, lwd=2)
arrows(str.plot, str.dat, str.plot, str.dat-str.err, length = 0, lwd=2)

#Num of straight jumps in LQ vs HQ food. ####### 
str.dat <- c(mean((straight.LQ.dat)), mean(straight.HQ.dat))
str.err <- c(sd((straight.LQ.dat))/sqrt(num.copes),sd((straight.HQ.dat))/sqrt(num.copes))

str.plot <- barplot(str.dat, ylim=c(0,max(str.dat + str.err)), xlab=c("LQ patch     HQ patch"), ylab="% time",
                    sub="% Time straight swimming", col=c("steelblue","goldenrod1"))
abline(h=0)
arrows(str.plot, str.dat, str.plot, str.dat+str.err, length = 0, lwd=2)
arrows(str.plot, str.dat, str.plot, str.dat-str.err, length = 0, lwd=2)

time.unit*str.dat
cell.size*str.dat

#Run a t-test to see if there's differnces in behavior in HQ vs LQ patches
t.test(x=(straight.LQ.dat), y=(straight.HQ.dat))

#######################################################################





















###########################################################################################################
###########################################################################################################
#Averages for text of manuscript ########
#####
#How many cells on average in the standardized cell.size?  #####

cell.dat <- c()
HQ.jump.dat <- c()
LQ.jump.dat <- c()
for(i in 1:1000) {
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
  
  # i-th element of `u1` squared into `i`-th position of `usq`
  cell.dat[i] <- cell.size
  HQ.jump.dat[i] <- mean(jump.hq)
  LQ.jump.dat[i] <- mean(jump.lq)
  
  i=i+1
}

max(cell.dat)
min(cell.dat)
mean(cell.dat)
sd(cell.dat)/sqrt(length(cell.dat))



mean(HQ.jump.dat)
sd(HQ.jump.dat)/sqrt(length(HQ.jump.dat))

mean(LQ.jump.dat)
sd(LQ.jump.dat)/sqrt(length(LQ.jump.dat))


###########################################################################################################
#How many cells on average in the standardized time.unit?  #####
time.dat <- c()
helix.HQ.dat <- c()
helix.LQ.dat <- c()
for(i in 1:1000) {
  helix.time.hq <- (c(0.08,0.032,0.032,0.112,0.096,0.104,0.048,0.056,0.184,0.192,0.032,0.048,0.048,0.048,0.032,0.032))/(10/1250)
  helix.time.lq <- (c(0.056,0.032,0.072,0.032,0.032,0.032,0.064,0.032,0.04))/(10/1250)
  
  loop.hq <- sample(x=helix.time.hq, size=1000, replace=T)
  loop.lq <- sample(x=helix.time.lq, size=1000, replace=T)
  
  #Find the standardized time unit in milliseconds
  time.unit<-gcd(nearq(mean(loop.hq)),nearq(mean(loop.lq)) )
  
  # i-th element of `u1` squared into `i`-th position of `usq`
  time.dat[i] <- time.unit
  helix.HQ.dat[i] <- mean(loop.hq)
  helix.LQ.dat[i] <- mean(loop.lq)
  i=i+1
}

max(time.dat)
min(time.dat)
mean(time.dat)
sd(time.dat)/sqrt(length(time.dat))


mean(helix.HQ.dat)
sd(helix.HQ.dat)/sqrt(length(helix.HQ.dat))


mean(helix.LQ.dat)
sd(helix.LQ.dat)/sqrt(length(helix.LQ.dat))



###########################################################################################################
#How much time moving vs waiting in HQ vs LQ foods?  ###############
HQ.move.dat <- c()
LQ.move.dat <- c()

for(i in 1:1000) {
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
  
  # i-th element of `u1` squared into `i`-th position of `usq`
  HQ.move.dat[i] <- prob.hq
  LQ.move.dat[i] <- prob.lq
  i=i+1
}

mean(HQ.move.dat)
sd(HQ.move.dat)/sqrt(length(HQ.move.dat))
min(HQ.move.dat)
max(HQ.move.dat)



mean(LQ.move.dat)
sd(LQ.move.dat)/sqrt(length(LQ.move.dat))
min(LQ.move.dat)
max(LQ.move.dat)


###########################################################################################################






###########################################################################################################
#######################################################################
# (Do I need the stuff below?)
#######################################################################

#######################################################################
#TO DO: take means of 0 and 1's to show frequency of helix vs swim vs rest... 
#Need to show that what I put in to behavior/model is what comes out.

#Get net vs total displacement? Show/animate a trajectory...?
#Heatmap of final position is for copes in model.
#count HQ vs LQ cells in during loops... 


#######################################################################
#Run a t-test to see if there's differnces in copepod density in HQ vs LQ patches
t.test(x=df.food$density, y=df.nonfood$density)
#######################################################################
#Barplots of copes in HQ vs LQ patches
#######################################################################
#Basic params for plotting #####
#Frequency in each cell (that was visited atleast once) ######
freq<-as.data.frame(table(nres))
freq<-cbind(freq,100*(freq[,2]/length(nres))) #*100 for %time in each cell
freq<-cbind(getdim(as.numeric(freq[,1]),grid2),freq) #putting X, Y, nres in each XY cell, Freq in each XY cell, and % time in each XY cell (=tiny fraction!)
#colnames(freq)
#[1] "xdim"                           "ydim"                           "nres"                          
#[4] "Freq"                           "100 * (freq[, 2]/length(nres))"

#Find the time spent in food HQ patches vs. LQ patches #####
patches<-which(grid2==1)
nonpatches<-which(grid2==0)

food.freq<-freq[which(freq$nres%in%patches),]
nonfood.freq<-freq[which(freq$nres%in%nonpatches),]

#Plotting means and SEs of frequency in each cell that was visited at least once; considering whether it's HQ or LQ #########
test.concat <-c(mean(food.freq[,5]), mean(nonfood.freq[,5])) 
test.error <- c(sd(food.freq[,5])/sqrt(length(food.freq[,5])), sd(nonfood.freq[,5])/sqrt(length(nonfood.freq[,5])))

plotdat <- barplot(test.concat, ylim=c(0,(max(test.concat)+max(test.error))), col=c("darkgoldenrod1", "steelblue"), cex.axis=1.35, cex.lab=1.5,
                   xlab=c("Hi qual patches       Low qual patches"), ylab="Average frequency of cope visity to particular cell of each type")
abline(h=0, lwd=2) #zero line
arrows(plotdat, test.concat-test.error, plotdat, test.concat, angle=90, code=1, length=0.05, lwd=2)
arrows(plotdat, test.concat+test.error, plotdat, test.concat, angle=90, code=1, length=0.05, lwd=2)

########################################################################
#plot first vs last half; how many in HQ vs LQ patch
#plot num copes in HQ vs LQ patches through time?




########################################################################
#Time within patches
########################################################################
# Do copes find a patch, and how long do they stay in/around it?
#create 0-1 data, compare 0 vs 1? (do they find the patch)
#once 1, how many 0 vs 1's? (how long are they in it?)

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
hist(afterstay/((timepoints-1)-starts), col="aliceblue", breaks=20)
summary(afterstay/((timepoints-1)-starts)) #ratio of time (out of 1) that copes spend in HQ food once they find it

hist(afterstay, col="aliceblue", breaks=20)
summary(afterstay) #this is the number of steps spend in HQ food across the 10,000 steps of the model
#this is the overall time in HQ food (not ness. one long visit.. )

par(mfrow=c(1,1))
########################################################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

########################################################################
# Looking at frequency in HQ and LQ patches ###
########################################################################
#combine all data on quality and counts from IBM. NOTES:
#Set the counter for behaviors. Columns are the copepods, rows are at each time step.
#in 'count' :  0 for straight runs and 1 for loops; NA= waiting during loops.p
#in 'qual'  :  0 for LQ and 1 for HQ
# These are n timesteps of rows, and n copepods wide/columns.
dim(count)
dim(qual)

time.dat <- cbind(count[c(1:25), c(1:5)], qual[c(1:25), c(1:5)])
na.contiguous(time.dat[,1])
which(count[,1]&&qual[,1]==1)


# change all colnames
colnames(time.dat) <- c("count","qual")

###########################################################################################################


###########################################################################################################


###########################################################################################################


###########################################################################################################


###########################################################################################################


###########################################################################################################
