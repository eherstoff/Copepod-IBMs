# Issue with the IBM loop not compiling; there's a problem with adding the helical counts (20s) to the counter loop without being out of bounds.



#######################################################################
# Last updated: 19 Feb 2019, EMH and MU
#Individual based model of adult A. tonsa foraging in patches that differ in stoichiometric quality
#Parameterize using data from EMH's 2016 Helgoland experiments. Focus adults only, looking at movement/behavior by precondition type
#######################################################################
#
#FOCAL QUESTIONS:
#(1) How big, how persistant do HQ patches need to be to find/use?
#   >> Halve HQ patch size 
#   >> Time to first locate HQ patch?
#   >> Residence time within HQ patch?
#(2) Do behaviors allow more possible feeding in HQ patches?
#   >> Cacluate average residency in HQ patch / patch size (= how more likely to be in HQ patch than 1:1?)
#   >> Calc. number of helicies in HQ vs. overall
#(3) How important is behavior?
#   >> Compare behavior/model with null (where LQ=HQ behaviors) 
#


#######################################################################
# Creating empty matricies so whole process can be looped
#######################################################################

cell.size.matrix2 <- matrix(nrow=3, ncol=1)
time.unit.matrix2 <- matrix(nrow=3, ncol=1)

#Residence time in HQ patch
mean.res.time_percent.matrix2 <- matrix(nrow=3, ncol=1)
sd.res.time_percent.matrix2   <- matrix(nrow=3, ncol=1)
mean.res.time_n.steps.matrix2 <- matrix(nrow=3, ncol=1)
sd.res.time_n.steps.matrix2   <- matrix(nrow=3, ncol=1)

#Steps until first found HQ patch
mean.first.visit.matrix2 <- matrix(nrow=3, ncol=1)
sd.first.visit.matrix2   <- matrix(nrow=3, ncol=1)
#Mean #copes found HQ patch through time
mean.means.matrix2 <- matrix(nrow=3, ncol=1)
sd.means.matrix2   <- matrix(nrow=3, ncol=1)

#Count helical swimming
mean.helix.LQ.dat.matrix2 <- matrix(nrow=3, ncol=1)
sd.helix.LQ.dat.matrix2   <- matrix(nrow=3, ncol=1)

mean.helix.HQ.dat.matrix2 <- matrix(nrow=3, ncol=1)
sd.helix.HQ.dat.matrix2   <- matrix(nrow=3, ncol=1)

#Time spent swimming in HQ and LQ food  
mean.helix.LQ_na.dat.matrix2 <- matrix(nrow=3, ncol=1)
sd.helix.LQ_na.dat.matrix2   <- matrix(nrow=3, ncol=1)

mean.helix.HQ_na.dat.matrix2 <- matrix(nrow=3, ncol=1)
sd.helix.HQ_na.dat.matrix2   <- matrix(nrow=3, ncol=1)

########
for (z in 1:length(cell.size.matrix2)) {
#######################################################################
#Setting up functions to calculate key params for the movement of individuals through the grid #####
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
#Parameters about movement pulled from experiment, and standardized for movement in model #######
########################################################################
#Standardize the time units based on loop times in different quality foods
#Standardize the cell unit based on swimming speeds in different quality foods
  #
#Input the time spent performing a loop in high quality (hq) or low-quality (lq) food. Times are in ms #######
#From "R code for movement model params, 17 Sept 2018
#This is taken from #GradSchool/Summer2016/Swimming analysis for EMH/Helix time budgeting, EMH experiments.xls
  
  #Data will be resampled 2000 times, and then the averages will be averaged... 
  
  #Create empty matricies to hold data
  time.unit <- matrix(nrow=2000, ncol=1)
  loop.hq   <- matrix(nrow=2000, ncol=1)
  loop.lq   <- matrix(nrow=2000, ncol=1)
  
  #For loop to resample study data
  for(i in 1:length(time.unit)) {
    #This is the time spent swimming in a loop (in ms), divided by the frame rate I watched (every 10th frame, shot at 1/1250 fps)
    helix.time.hq <- (c(0.08,0.032,0.032,0.112,0.096,0.104,0.048,0.056,0.184,0.192,0.032,0.048,0.048,0.048,0.032,0.032))/(10/1250)
    helix.time.lq <- (c(0.056,0.032,0.072,0.032,0.032,0.032,0.064,0.032,0.04))/(10/1250)
    
    #Resample with replacement from the data
    loop.hq[i,] <- mean(sample(x=helix.time.hq, size=2000, replace=T) )
    loop.lq[i,] <- mean(sample(x=helix.time.lq, size=2000, replace=T) )
    
    #Find the standardized time unit in milliseconds
    time.unit[i,] <-gcd(nearq(loop.hq[i,]),nearq(loop.lq[i,]) )
    
    i=i+1 }
  
  
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
  
  #Resample a bunch here, like for loop and time data      #DON'T DO THIS ANYMORE
  #speed.hq <- matrix(nrow=10000, ncol=1)
  #speed.lq <- matrix(nrow=10000, ncol=1)
  
  #for(i in 1:length(speed.hq)) {
  #speed.hq[i,] <- mean(round(nearq(sample(x=f2.dat2$velocity, size=1000, replace = T))))
  #speed.lq[i,] <- mean(round(nearq(sample(x=noN.dat2$velocity, size=1000, replace = T)))) 
  #i=i+1 }
  
  
#Find the model's cell size in mm  #######
  #cell.size<-gcd(nearq(mean(speed.hq)),nearq(mean(speed.lq)))
  cell.size<-gcd(nearq(mean(f2.dat2$velocity)), nearq(mean(noN.dat2$velocity)))  #Probablly should use this one instead! No need to resample things--plenty of points to work with (1115 f2, 805 noN)
  #with non-resampled data, cell.size = 0.3125
  
#How many cells do they travel in one time unit? #######
  jump.hq<- time.unit[1:1115]* (f2.dat2$velocity)/cell.size     #Need to take only n timeunits, as long as data in velocity
  jump.lq<- time.unit[1116:1920]* (noN.dat2$velocity)/cell.size 
  
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

# Probabilities of behaviors in HQ and LQ   ###############
  prob.wait.hq <- mean(ratio.waiting.hq) #probability of pause/waiting in HQ food = 0.2077741
  sd.wait.hq   <- sd(ratio.waiting.hq)
  
  prob.helix.hq <- mean(ratio.helix.hq) #probability of helical swimming in HQ food = 0.2707739
  sd.helix.hq   <- sd(ratio.helix.hq)
  
  prob.straight.hq <- mean(ratio.otheractive.hq) #probability of straight swimming in HQ food = 0.521452
  sd.straight.hq   <- sd(ratio.otheractive.hq)
  
  #So, if in HQ, numbers from   
  #                          0 - 21  = Wait @ 21 steps
  #                         22 - 49  = Helix @ 27 steps
  #                         50 - 100 = Straight @ 52 steps
  # # # #
  
  prob.wait.lq <- mean(ratio.waiting.lq) #probability of pause/waiting in lq food = 0.4555456
  sd.wait.lq   <- sd(ratio.waiting.lq)
  
  prob.helix.lq <- mean(ratio.helix.lq) #probability of helical swimming in lq food = 0.1588707
  sd.helix.lq   <- sd(ratio.helix.lq)
  
  prob.straight.lq <- mean(ratio.otheractive.lq) #probability of straight swimming in lq food = 0.3855837
  sd.straight.lq   <- sd(ratio.otheractive.lq)
  
  #So, if in LQ, numbers from   
  #                          0 - 45  = Wait @ 46 steps
  #                         46 - 62  = Helix @ 16 steps
  #                         63 - 100 = Straight @ 39 steps
  # # # #

#How long to spend waiting? ######
wait.times <-read.csv("Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/Calc params/*wait times for IBM, 19 Feb 2019.csv")
HQ.wait <- wait.times[1:275, 4]
LQ.wait <- wait.times[276:516, 4]
#Where XY movement was zero and they weren't performing helical swimming, calc. number of steps (and time in s) they wait.  

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
  grid.size <- 1000 #run at 500  1,000   5,000 px^2
  grid2<-matrix(rep(0,grid.size^2),nrow=grid.size,ncol=grid.size) #the grid is created and filled with LQ food
  
  #Add in HQ patches, where model grid.size is 500 px  #####
  #grid2[0:354,0:354]<-1    # ~50% % has HQ food; HQ food patch at the corner of the model.
  #grid2[0:250,0:250]<-1    # 25% has HQ food; HQ food patch at the corner of the model.
  #grid2[0:177,0:177]<-1    # ~12.5% has HQ food; HQ food patch at the corner of the model.
  #grid2[0:125,0:125]<-1  # 6.25% has HQ food; HQ food patch at the corner of the model.
  #grid2[0:88,0:88]<-1    # ~3.13% has HQ food; HQ food patch at the corner of the model.
  #grid2[0:62,0:62]<-1    # ~1.56% has HQ food; HQ food patch at the corner of the model.
  #grid2[0:44,0:44]<-1    #  ~0.77% has HQ food; HQ food patch at the corner of the model.
  #grid2[0:31,0:31]<-1    #  ~0.38% has HQ food; HQ food patch at the corner of the model.
  #grid2[0:22,0:22]<-1    #  ~0.1936% has HQ food; HQ food patch at the corner of the model.
  #grid2[0:16,0:16]<-1    #  ~0.0961% has HQ food; HQ food patch at the corner of the model.
  #grid2[0:11,0:11]<-1    #  ~0.04805% has HQ food; HQ food patch at the corner of the model.
  
  
  
  #Add in HQ patches, where model grid.size is 1,000 px  #####
  #grid2[0:866,0:866]<-1    # ~75% has HQ food in 1000 x 1000 px; HQ food patch at the corner of the model.
  #grid2[0:707,0:707]<-1    # ~50% % has HQ food in 1000 x 1000 px; HQ food patch at the corner of the model.
  #grid2[0:500,0:500]<-1    # 25% has HQ food in 1000 x 1000 px; HQ food patch at the corner of the model.
  #grid2[0:354,0:354]<-1    # ~12.5% has HQ food in 1000 x 1000 px; HQ food patch at the corner of the model.
  grid2[000:250,000:250]<-1  # 6.25% has HQ food in 1000 x 1000 px; HQ food patch at the corner of the model.
  #grid2[0:177,0:177]<-1    # ~3.13% has HQ food in 1000 x 1000 px; HQ food patch at the corner of the model.
  #grid2[0:125,0:125]<-1    # ~1.56% has HQ food in 1000 x 1000 px; HQ food patch at the corner of the model.
  #grid2[0:88,0:88]<-1    #  ~0.77% has HQ food in 1000 x 1000 px; HQ food patch at the corner of the model.
  #grid2[0:62,0:62]<-1    #  ~0.38% has HQ food in 1000 x 1000 px; HQ food patch at the corner of the model.
  #grid2[0:44,0:44]<-1    #  ~0.1936% has HQ food in 1000 x 1000 px; HQ food patch at the corner of the model.
  #grid2[0:31,0:31]<-1    #  ~0.0961% has HQ food in 1000 x 1000 px; HQ food patch at the corner of the model.
  #grid2[0:22,0:22]<-1    #  ~0.04805% has HQ food in 1000 x 1000 px; HQ food patch at the corner of the model.
  
  
  
  
  #Add in HQ patches, where model grid.size is 5,000 px  #####
  #grid2[0:3540,0:3540]<-1    # ~50% % has HQ food; HQ food patch at the corner of the model.
  #grid2[0:2500,0:2500]<-1    # 25% has HQ food; HQ food patch at the corner of the model.
  #grid2[0:1770,0:1770]<-1    # ~12.5% has HQ food; HQ food patch at the corner of the model.
  #grid2[0:1250,0:1250]<-1  # 6.25% has HQ food; HQ food patch at the corner of the model.
  #grid2[0:880,0:880]<-1    # ~3.13% has HQ food; HQ food patch at the corner of the model.
  #grid2[0:620,0:620]<-1    # ~1.56% has HQ food; HQ food patch at the corner of the model.
  #grid2[0:440,0:440]<-1    #  ~0.77% has HQ food; HQ food patch at the corner of the model.
  #grid2[0:310,0:310]<-1    #  ~0.38% has HQ food; HQ food patch at the corner of the model.
  #grid2[0:220,0:220]<-1    #  ~0.1936% has HQ food; HQ food patch at the corner of the model.
  #grid2[0:160,0:160]<-1    #  ~0.0961% has HQ food; HQ food patch at the corner of the model.
  #grid2[0:110,0:110]<-1    #  ~0.04805% has HQ food; HQ food patch at the corner of the model.
  
########################################################################
  ##This is what % of patches within the grid had food (where hq=1)   ##### 
  patch.food = length(which(grid2==1)) / 
    (length(which(grid2==1))+length(which(grid2==0)))  
########################################################################
  
##############
  #Null model test: parameters for HQ and LQ are the same  #####
  #Cancel this bit when running the actual simulation
  #angles.HQ_XY<-angles.LQ_XY
  #prob.hq<-prob.lq
  #jump.hq<-jump.lq
  #loop.hq<-loop.lq
  
##############
#Setting up matrix######
  timepoints <-25 #100000  #time steps in the model. 
  num.copes <- 5 #1000     #number of copes to simulate in the model
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
  
  rand.num <- runif(1) #Generate a number from 0-1 for determining behavior type.
  #So, if in HQ, numbers from   
  #                          0 - .21  = Wait @ 21 steps
  #                         .22 - .49  = Helix @ 27 steps
  #                         .50 - 1.00 = Straight @ 52 steps
  # # # #
  #So, if in LQ, numbers from   
  #                          0 - .45  = Wait @ 46 steps
  #                         .46 - .62  = Helix @ 16 steps
  #                         .63 - 1.00 = Straight @ 39 steps
  # # # #
  
  #Set the counter for behaviors. Columns are the copepods, rows are at each time step.
  count<-matrix(nrow=timepoints,ncol=num.copes)    #System: 0 for wait/stay,10 for keep waiting, 1 for straight runs and 2 for loops; 20= keep looping.
  qual<-matrix(nrow=timepoints,ncol=num.copes)     #0 for LQ and 1 for HQ
  
  for(i in 1:num.copes){  #number of replicates; needs to match ncol in 'Setting up matrix'
    init<-c(sample(1:nrow(grid2),1),sample(1:ncol(grid2),1))  #Initiate from a random cell
    res.i[1,i]<-getind(init,grid2)
    count[1,i]<-1
    qual[1,i]<-grid2[getind(init,grid2)]
    j<-2
    
    while(j<timepoints){ #time steps (must match time steps in 'Setting up matrix' above)  
      rand.num<-runif(1) #generate a random number to determine behavior
      
      if(grid2[init[1],init[2]]==0){  # if you're in LQ patch
        
        if(rand.num<0.62 && rand.num>0.46){ #If helical swimming in LQ; stay put for set number of timesteps
          count[j,i]<-2  #update the behavior counter
          sample.x <- floor(sample(time.lq, size=1))
          
          
          count[c((j+1):(j+sample.x)),i]<-20 #Error is here. "subscript out of bounds". j generating decimals sometimes??
          
          
          res.i[j:min(((j+sample.x)-1),timepoints),i]<-getind(init,grid2)
          qual[j:min(((j+sample.x)-1),timepoints),i]<-0
          j<-j+sample.x
          qual[j,i]<-0
        }
        
        else{  
          if(rand.num>=0.62){ #If Straight swim @ angle
            count[j,i]<-1  #Update the counter
            ind<-sample(1:11,1,prob=angles.LQ_XY)    #turning angle changes depending on food quality
            jump<-sample(jump.lq,size=1)
            x<-init[1]+floor(jump*sets[ind,1])
            y<-init[2]+floor(jump*sets[ind,2])
            if(x%in%1:nrow(grid2) & y%in%1:ncol(grid2)){init<-c(x,y)} #If they don't swim off the environment edge, record X,Y
            else{init<-c(x%%nrow(grid2),y%%ncol(grid2))} #(When they swim off the edge, wrap their position)
            if(init[1]==0){init[1]<-1} #Wrapper code for dealing with copes at 0 (boundary) => move over one pixel to 1
            if(init[2]==0){init[2]<-1}
            res.i[j,i]<-getind(init,grid2) 
            j<-j+1
            qual[j,i]<-0
          }
          
          else{  
            if(rand.num<=0.46){ #if wait/paused
              count[j,i]<-0
              sample.x<- sample(LQ.wait, size=1) #this is how long to wait in LQ
                count[(j+1):(j+sample.x-1),i]<-10  
              res.i[j:min(((j+sample.x)-1),timepoints),i]<-getind(init,grid2)
              qual[j:min(((j+sample.x)-1),timepoints),i]<-0
              j<-j+sample.x #pause one step
              qual[j,i]<-0
            }
          }
        }
      }#end of LQ patch statement
      
      # Second part of the code, for HQ food #######
      if(grid2[init[1],init[2]] == 1){   #when they are looping in HQ, stay put for set number of timesteps
        if(rand.num>0.22 && rand.num<0.5){  #If helical swimming
          count[j,i]<-2  #update the behavior counter
          sample.x <- floor(sample(time.hq, size=1))
          count[j+1:j+sample.x-1,i]<-20
          res.i[j:min((j+sample.x-1),timepoints),i]<-getind(init,grid2)
          qual[j:min(((j+sample.x)-1),timepoints),i]<-1
          j<-j+sample.x
          qual[j,i]<-1
        }
        else{  #angles and movement for HQ
          if(rand.num>=0.5){ #If Straight swim @ angle
            count[j,i]<-1  #Update the counter
            ind<-sample(1:11,1,prob=angles.HQ_XY)
            jump<-sample(jump.hq,size=1)
            x<-init[1]+floor(jump*sets[ind,1])
            y<-init[2]+floor(jump*sets[ind,2])
            if(x%in%1:nrow(grid2) & y%in%1:ncol(grid2)){init<-c(x,y)} #wrapper for tiled environment
            else{init<-c(x%%nrow(grid2),y%%ncol(grid2))}
            if(init[1]==0){init[1]<-1}
            if(init[2]==0){init[2]<-1}
            res.i[j,i]<-getind(init,grid2)
            j<-j+1
            qual[j,i]<-1
          }
          
          else{  
            if(rand.num>=0.22){ #if wait/paused
              count[j,i]<-0
              sample.x<-sample(HQ.wait, size=1) #How long to wait in HQ
                count[j+1:j+sample.x-1,i]<-10  
              res.i[j:min(((j+sample.x)-1),timepoints),i]<-getind(init,grid2)
              qual[j:min(((j+sample.x)-1),timepoints),i]<-1
              j<-j+sample.x #pause x step
              qual[j,i]<-1
            }  
          }
        }
      }
    }
  }


res.i<-res.i[-nrow(res.i),]
res.i<-res.i[-1,]
##############

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
hist(res.time, col="grey")

########################################################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#determine when they first find a patch ##########
ones<-function(a){
  one<-sum(a!=1)
  if(one < length(a) ){ return(min(which(a==1))) }
  else(return(length(a)))
}

first.visit<-apply(res.j,2,ones)

#hist(first.visit, col="aliceblue")
#abline(v=(patch.food*timepoints), col="chocolate1")



########################################################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#Plot the mean % of copepods in HQ vs. LQ over time #####
#use the data on 'count' from IBM, which says whether the animal is in HQ or LQ food
means<-c() 
for(i in 1:nrow(res.j)){ a
  means<-c(means,sum(res.j[i,])/length(res.j[1,]))
}

plot(means,type="l", ylim=c(0,1)) #average occupancy in patches over time.
abline(h=patch.food, col="red") #this line is the % of the model filled with HQ food

mean(means)
sd(means)/sqrt(length(means))


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



###########################################################################################################

###########################################################################################################
#Data needed for spreadsheet
###########################################################################################################
cell.size.matrix2[z,] <-mean(cell.size)
time.unit.matrix2[z,] <-mean(time.unit)

#Residence time in HQ patch
mean.res.time_percent.matrix2[z,] <-mean(res.time, na.rm=T)  #% time
sd.res.time_percent.matrix2[z,] <-sd(res.time, na.rm=T)
mean.res.time_n.steps.matrix2[z,] <-mean(apply(res.j,2,sum, na.rm=T)) #n steps
sd.res.time_n.steps.matrix2[z,] <-sd(apply(res.j,2,sum, na.rm=T))

#Steps until first found HQ patch
mean.first.visit.matrix2[z,] <-mean(first.visit) #number of steps in model
sd.first.visit.matrix2[z,] <-sd(first.visit)
#Mean #copes found HQ patch through time
mean.means.matrix2[z,] <-mean(means)
sd.means.matrix2[z,] <-sd(means)

#Count helical swimming
mean.helix.LQ.dat.matrix2[z,] <-mean(helix.LQ.dat)  #number helices in LQ
sd.helix.LQ.dat.matrix2[z,] <-sd(helix.LQ.dat)

mean.helix.HQ.dat.matrix2[z,] <-mean(helix.HQ.dat)  #number helices in HQ
sd.helix.HQ.dat.matrix2[z,] <-sd(helix.HQ.dat)

#Time spent swimming in HQ and LQ food  
mean.helix.LQ_na.dat.matrix2[z,] <-mean(helix.LQ_na.dat) #time spent helixing in LQ
sd.helix.LQ_na.dat.matrix2[z,] <-sd(helix.LQ_na.dat)

mean.helix.HQ_na.dat.matrix2[z,] <-mean(helix.HQ_na.dat) #time spent helixing in HQ
sd.helix.HQ_na.dat.matrix2[z,] <-sd(helix.HQ_na.dat)


# finish the grand for loop with a progress bar and adding to z  ######
#create a progress bar
print(z)
Sys.sleep(0.01)
flush.console()

z=z+1  }
#####



###########################################################################################################
#Outputs for spreadsheet
###########################################################################################################
final.dat <- cbind(cell.size.matrix2, time.unit.matrix2, 
                   #Residence time in HQ patch
                   mean.res.time_percent.matrix2,
                   sd.res.time_percent.matrix2,
                   mean.res.time_n.steps.matrix2,
                   sd.res.time_n.steps.matrix2,
                   #Steps until first found HQ patch
                   mean.first.visit.matrix2,
                   sd.first.visit.matrix2,
                   #Mean #copes found HQ patch through time
                   mean.means.matrix2,
                   sd.means.matrix2,
                   #Count helical swimming
                   mean.helix.LQ.dat.matrix2,
                   sd.helix.LQ.dat.matrix2,
                   mean.helix.HQ.dat.matrix2,
                   sd.helix.HQ.dat.matrix2,
                   #Time spent swimming in HQ and LQ food  
                   mean.helix.LQ_na.dat.matrix2,
                   sd.helix.LQ_na.dat.matrix2,
                   mean.helix.HQ_na.dat.matrix2,
                   sd.helix.HQ_na.dat.matrix2 )

#####
print(final.dat)
#####





###########################################################################################################
#Plotting the mean of the X and Y coordinates in space that were visited by the cope runs
#Run a simulation once, and plot the trajectory of the 1st animal in the model ######

library(plot3D) #will plot in color

#find which cope # had the minimum residence time
which(res.time==min(res.time))
min(res.time)                  #Here, #10 @ 0.01567
#find which cope # had the minimum residence time
which(res.time==max(res.time))  
max(res.time)                  #Here, #418 @ 0.28982

#average residence time
mean(res.time) #0.1193485
sd(res.time)   #0.04641549
#Scanning through the residence times, #34 has a res.time of 0.11985

res.k<-matrix(nrow=nrow(res.i),ncol=2) #the mean of all of the copepod XY points.
for(i in 1:nrow(res.i)){
  res.k[i,]<-getdim(res.i[i,34],grid2)
}


#Plot the trajectory as B&W, and as color
plot(x=res.k[,2], y=res.k[,1],type="l",xlim=c(0,grid.size),ylim=c(0,grid.size), main="B&W movement plot") 

scatter2D(x=res.k[,2], y=res.k[,1], col = rainbow(n = length(res.k[,2])),
          xlim=c(0,grid.size),ylim=c(0,grid.size), type="l", main="Colorized movement for res.time; red=start, blue=end",
          sub="~Avg res.time in model, cope #34 @ 0.11985 (avg=0.1193485)") 


min.cope <- cbind(res.j[,c(10)], res.k)  #min. cope(#10) helix/not (col 1), with X, Y position (col 2,3)
max.cope <- cbind(res.j[,c(418)], res.k) #max cope(#418) 
avg.cope <- cbind(res.j[,c(34)], res.k)  #cope 34 ~ avg.

library(xlsx)
write.csv(min.cope, "min_cope.xls")
write.csv(max.cope, "max_cope.xls")
write.csv(avg.cope, "avg_cope.xls")

###########################################################################################################
#plotting average cope (above) movement, without lines connecting from side-to-side.

#1-10 ####
set1<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 1.csv')
set2<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 2.csv')
set3<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 3.csv')
set4<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 4.csv')
set5<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 5.csv')
set6<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 6.csv')
set7<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 7.csv')
set8<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 8.csv')
set9<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 9.csv')
set10<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 10.csv')
#11-20 ####
set11<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 11.csv')
set12<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 12.csv')
set13<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 13.csv')
set14<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 14.csv')
set15<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 15.csv')
set16<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 16.csv')
set17<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 17.csv')
set18<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 18.csv')
set19<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 19.csv')
set20<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 20.csv')
#21-30 ####
set21<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 21.csv')
set22<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 22.csv')
set23<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 23.csv')
set24<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 24.csv')
set25<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 25.csv')
set26<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 26.csv')
set27<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 27.csv')
set28<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 28.csv')
set29<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 29.csv')
set30<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 30.csv')
#31-41 ####
set31<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 31.csv')
set32<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 32.csv')
set33<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 33.csv')
set34<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 34.csv')
set35<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 35.csv')
set36<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 36.csv')
set37<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 37.csv')
set38<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 38.csv')
set39<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 39.csv')
set40<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 40.csv')
#41-50 ####
set41<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 41.csv')
set42<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 42.csv')
set43<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 43.csv')
set44<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 44.csv')
set45<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 45.csv')
set46<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 46.csv')
set47<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 47.csv')
set48<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 48.csv')
set49<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 49.csv')
set50<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 50.csv')
#51-60 ####
set51<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 51.csv')
set52<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 52.csv')
set53<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 53.csv')
set54<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 54.csv')
set55<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 55.csv')
set56<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 56.csv')
set57<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 57.csv')
set58<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 58.csv')
set59<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 59.csv')
set60<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 60.csv')
#61-70 (with 61_ii, 61_iii splits too)####
set61<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 61.csv')
set61_ii <- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 61_ii.csv')
set61_iii <- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 61_ii.csv')
set62<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 62.csv')
set63<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 63.csv')
set63_ii<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 63_ii.csv')
set64<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 64.csv')
set65<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 65.csv')
set66<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 66.csv')
set67<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 67.csv')
set68<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 68.csv')
set69<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 69.csv')
set70<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 70.csv')
#71-80 ####
set71<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 71.csv')
set72<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 72.csv')
set73<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 73.csv')
set74<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 74.csv')
set75<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 75.csv')
set76<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 76.csv')
set77<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 77.csv')
set78<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 78.csv')
set79<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 79.csv')
set80<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 80.csv')
#81-90 ####
set81<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 81.csv')
set82<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 82.csv')
set83<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 83.csv')
set84<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 84.csv')
set85<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 85.csv')
set86<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 86.csv')
set87<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 87.csv')
set88<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 88.csv')
set89<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 89.csv')
set90<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 90.csv')
#91-100 ####
set91<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 91.csv')
set92<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 92.csv')
set93<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 93.csv')
set94<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 94.csv')
set95<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 95.csv')
set96<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 96.csv')
set97<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 97.csv')
set98<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 98.csv')
set99<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 99.csv')
set100<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 100.csv')
#101-103 ####
set101<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 101.csv')
set102<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 102.csv')
set103<- read.csv('Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/*Current model/Cope tracking/Avg cope data sets/Avg cope, set 103.csv')

######
#Color gradients: of 22 points, from red (start) to blue (end)...
#colfunc <- colorRampPalette(c("red", "blue"))
#colfunc(22)
#[1] "#FF0000" "#F2000C" "#E60018" "#DA0024" "#CE0030" "#C2003C" "#B60048" "#AA0055" "#9D0061" "#91006D" "#850079"
#[12] "#790085" "#6D0091" "#61009D" "#5500AA" "#4800B6" "#3C00C2" "#3000CE" "#2400DA" "#1800E6" "#0C00F2" "#0000FF"

#1-10:    "#FF0000" "#F2000C"
#11-20:   "#E60018" "#DA0024" 
#21-30:   "#CE0030" "#C2003C" 
#31-40:   "#B60048" "#AA0055" 
#41-50:   "#9D0061" "#91006D" 
#51-60:   "#850079" "#790085" 
#61-70:   "#6D0091" "#61009D" "#5500AA"
#71-80:   "#4800B6" "#3C00C2"
#81-90:   "#3000CE" "#2400DA"
#91-100:  "#1800E6" "#0C00F2"
#100-103:  "#0000FF"

#Create an empty plot space to put data from the model run####
plot(0~1, xlim=c(0,1000), ylim=c(0,1000), 
     type="l", pch=20, cex=0.75, col="white", xlab="  ", ylab="  ")
#Create a black box to show where HQ food is in the model space #####
HQfood <- rep(177, 177)
vertical <- c(1:177)
zeros <- rep(0,177)

points(HQfood ~ HQfood,  type="l", col="black", lwd=3.5) #top side
points(zeros ~ vertical, type="l", col="black", lwd=3.5) #bottom side

points(vertical ~ HQfood, type="l", col="black", lwd=3.5) #left side
points(vertical ~ zeros, type="l", col="black", lwd=3.5) #right side


#Plot 1-10 ####
points(set1$X.position~set1$Y.position, type="l", pch=20, cex=0.75, col="#FF0000")
points(set2$X.position~set2$Y.position, type="l", pch=20, cex=0.75, col="#FF0000")
points(set3$X.position~set3$Y.position, type="l", pch=20, cex=0.75, col="#FF0000")
points(set4$X.position~set4$Y.position, type="l", pch=20, cex=0.75, col="#FF0000")
points(set5$X.position~set5$Y.position, type="l", pch=20, cex=0.75, col="#FF0000")
points(set6$X.position~set6$Y.position, type="l", pch=20, cex=0.75, col="#F2000C")
points(set7$X.position~set7$Y.position, type="l", pch=20, cex=0.75, col="#F2000C")
points(set8$X.position~set8$Y.position, type="l", pch=20, cex=0.75, col="#F2000C")
points(set9$X.position~set9$Y.position, type="l", pch=20, cex=0.75, col="#F2000C")
points(set10$X.position~set10$Y.position, type="l", pch=20, cex=0.75, col="#F2000C")

#11-20 ####
points(set11$X.position~set11$Y.position, type="l", pch=20, cex=0.75, col="#E60018")
points(set12$X.position~set12$Y.position, type="l", pch=20, cex=0.75, col="#E60018")
points(set13$X.position~set13$Y.position, type="l", pch=20, cex=0.75, col="#E60018")
points(set14$X.position~set14$Y.position, type="l", pch=20, cex=0.75, col="#E60018")
points(set15$X.position~set15$Y.position, type="l", pch=20, cex=0.75, col="#E60018")
points(set16$X.position~set16$Y.position, type="l", pch=20, cex=0.75, col="#DA0024")
points(set17$X.position~set17$Y.position, type="l", pch=20, cex=0.75, col="#DA0024")
points(set18$X.position~set18$Y.position, type="l", pch=20, cex=0.75, col="#DA0024")
points(set19$X.position~set19$Y.position, type="l", pch=20, cex=0.75, col="#DA0024")
points(set20$X.position~set20$Y.position, type="l", pch=20, cex=0.75, col="#DA0024")

#21-30 ####
points(set21$X.position~set21$Y.position, type="l", pch=20, cex=0.75, col="#CE0030")
points(set22$X.position~set22$Y.position, type="l", pch=20, cex=0.75, col="#CE0030")
points(set23$X.position~set23$Y.position, type="l", pch=20, cex=0.75, col="#CE0030")
points(set24$X.position~set24$Y.position, type="l", pch=20, cex=0.75, col="#CE0030")
points(set25$X.position~set25$Y.position, type="l", pch=20, cex=0.75, col="#CE0030")
points(set26$X.position~set26$Y.position, type="l", pch=20, cex=0.75, col="#C2003C")
points(set27$X.position~set27$Y.position, type="l", pch=20, cex=0.75, col="#C2003C")
points(set28$X.position~set28$Y.position, type="l", pch=20, cex=0.75, col="#C2003C")
points(set29$X.position~set29$Y.position, type="l", pch=20, cex=0.75, col="#C2003C")
points(set30$X.position~set30$Y.position, type="l", pch=20, cex=0.75, col="coral")

#31-40 ####
points(set31$X.position~set31$Y.position, type="l", pch=20, cex=0.75, col="#B60048")
points(set32$X.position~set32$Y.position, type="l", pch=20, cex=0.75, col="#B60048")
points(set33$X.position~set33$Y.position, type="l", pch=20, cex=0.75, col="#B60048")
points(set34$X.position~set34$Y.position, type="l", pch=20, cex=0.75, col="#B60048")
points(set35$X.position~set35$Y.position, type="l", pch=20, cex=0.75, col="#B60048")
points(set36$X.position~set36$Y.position, type="l", pch=20, cex=0.75, col="#AA0055")
points(set37$X.position~set37$Y.position, type="l", pch=20, cex=0.75, col="#AA0055")
points(set38$X.position~set38$Y.position, type="l", pch=20, cex=0.75, col="#AA0055")
points(set39$X.position~set39$Y.position, type="l", pch=20, cex=0.75, col="#AA0055")
points(set40$X.position~set40$Y.position, type="l", pch=20, cex=0.75, col="#AA0055")
#41-50 ####
points(set41$X.position~set41$Y.position, type="l", pch=20, cex=0.75, col="#9D0061")
points(set42$X.position~set42$Y.position, type="l", pch=20, cex=0.75, col="#9D0061")
points(set43$X.position~set43$Y.position, type="l", pch=20, cex=0.75, col="#9D0061")
points(set44$X.position~set44$Y.position, type="l", pch=20, cex=0.75, col="#9D0061")
points(set45$X.position~set45$Y.position, type="l", pch=20, cex=0.75, col="#9D0061")
points(set46$X.position~set46$Y.position, type="l", pch=20, cex=0.75, col="#91006D")
points(set47$X.position~set47$Y.position, type="l", pch=20, cex=0.75, col="#91006D")
points(set48$X.position~set48$Y.position, type="l", pch=20, cex=0.75, col="#91006D")
points(set49$X.position~set49$Y.position, type="l", pch=20, cex=0.75, col="#91006D")
points(set50$X.position~set50$Y.position, type="l", pch=20, cex=0.75, col="#91006D")
#51-60 ####
points(set51$X.position~set51$Y.position, type="l", pch=20, cex=0.75, col="#850079")
points(set52$X.position~set52$Y.position, type="l", pch=20, cex=0.75, col="#850079")
points(set53$X.position~set53$Y.position, type="l", pch=20, cex=0.75, col="#850079")
points(set54$X.position~set54$Y.position, type="l", pch=20, cex=0.75, col="#850079")
points(set55$X.position~set55$Y.position, type="l", pch=20, cex=0.75, col="#850079")
points(set56$X.position~set56$Y.position, type="l", pch=20, cex=0.75, col="#790085")
points(set57$X.position~set57$Y.position, type="l", pch=20, cex=0.75, col="#790085")
points(set58$X.position~set58$Y.position, type="l", pch=20, cex=0.75, col="#790085")
points(set59$X.position~set59$Y.position, type="l", pch=20, cex=0.75, col="#790085")
points(set60$X.position~set60$Y.position, type="l", pch=20, cex=0.75, col="#790085")
#61-70 ####
points(set61$X.position~set61$Y.position, type="l", pch=20, cex=0.75, col="#6D0091")
points(set61_ii$X.position~set61_ii$Y.position, type="l", pch=20, cex=0.75, col="#6D0091")
points(set61_iii$X.position~set61_iii$Y.position, type="l", pch=20, cex=0.75, col="#6D0091")
points(set62$X.position~set62$Y.position, type="l", pch=20, cex=0.75, col="#6D0091")
points(set63$X.position~set63$Y.position, type="l", pch=20, cex=0.75, col="#6D0091")
points(set63_ii$X.position~set63_ii$Y.position, type="l", pch=20, cex=0.75, col="#61009D")
points(set64$X.position~set64$Y.position, type="l", pch=20, cex=0.75, col="#61009D")
points(set65$X.position~set65$Y.position, type="l", pch=20, cex=0.75, col="#61009D")
points(set66$X.position~set66$Y.position, type="l", pch=20, cex=0.75, col="#61009D")
points(set67$X.position~set67$Y.position, type="l", pch=20, cex=0.75, col="#5500AA")
points(set68$X.position~set68$Y.position, type="l", pch=20, cex=0.75, col="#5500AA")
points(set69$X.position~set69$Y.position, type="l", pch=20, cex=0.75, col="#5500AA")
points(set70$X.position~set70$Y.position, type="l", pch=20, cex=0.75, col="#5500AA")
#71-80 ####
points(set71$X.position~set71$Y.position, type="l", pch=20, cex=0.75, col="#4800B6")
points(set72$X.position~set72$Y.position, type="l", pch=20, cex=0.75, col="#4800B6")
points(set73$X.position~set73$Y.position, type="l", pch=20, cex=0.75, col="#4800B6")
points(set74$X.position~set74$Y.position, type="l", pch=20, cex=0.75, col="#4800B6")
points(set75$X.position~set75$Y.position, type="l", pch=20, cex=0.75, col="#4800B6")
points(set76$X.position~set76$Y.position, type="l", pch=20, cex=0.75, col="#3C00C2")
points(set77$X.position~set77$Y.position, type="l", pch=20, cex=0.75, col="#3C00C2")
points(set78$X.position~set78$Y.position, type="l", pch=20, cex=0.75, col="#3C00C2")
points(set79$X.position~set79$Y.position, type="l", pch=20, cex=0.75, col="#3C00C2")
points(set80$X.position~set80$Y.position, type="l", pch=20, cex=0.75, col="#3C00C2")
#81-90 ####
points(set81$X.position~set81$Y.position, type="l", pch=20, cex=0.75, col="#3000CE")
points(set82$X.position~set82$Y.position, type="l", pch=20, cex=0.75, col="#3000CE")
points(set83$X.position~set83$Y.position, type="l", pch=20, cex=0.75, col="#3000CE")
points(set84$X.position~set84$Y.position, type="l", pch=20, cex=0.75, col="#3000CE")
points(set85$X.position~set85$Y.position, type="l", pch=20, cex=0.75, col="#3000CE")
points(set86$X.position~set86$Y.position, type="l", pch=20, cex=0.75, col="#2400DA")
points(set87$X.position~set87$Y.position, type="l", pch=20, cex=0.75, col="#2400DA")
points(set88$X.position~set88$Y.position, type="l", pch=20, cex=0.75, col="#2400DA")
points(set89$X.position~set89$Y.position, type="l", pch=20, cex=0.75, col="#2400DA")
points(set90$X.position~set90$Y.position, type="l", pch=20, cex=0.75, col="#2400DA")
#91-100 ####
points(set91$X.position~set91$Y.position, type="l", pch=20, cex=0.75, col="#1800E6")
points(set92$X.position~set92$Y.position, type="l", pch=20, cex=0.75, col="#1800E6")
points(set93$X.position~set93$Y.position, type="l", pch=20, cex=0.75, col="#1800E6")
points(set94$X.position~set94$Y.position, type="l", pch=20, cex=0.75, col="#1800E6")
points(set95$X.position~set95$Y.position, type="l", pch=20, cex=0.75, col="#1800E6")
points(set96$X.position~set96$Y.position, type="l", pch=20, cex=0.75, col="#0C00F2")
points(set97$X.position~set97$Y.position, type="l", pch=20, cex=0.75, col="#0C00F2")
points(set98$X.position~set98$Y.position, type="l", pch=20, cex=0.75, col="#0C00F2")
points(set99$X.position~set99$Y.position, type="l", pch=20, cex=0.75, col="#0C00F2")
points(set100$X.position~set100$Y.position, type="l", pch=20, cex=0.75, col="#0C00F2")
#101-103 ####
points(set101$X.position~set101$Y.position, type="l", pch=20, cex=0.75, col="#0000FF")
points(set102$X.position~set102$Y.position, type="l", pch=20, cex=0.75, col="#0000FF")
points(set103$X.position~set103$Y.position, type="l", pch=20, cex=0.75, col="#0000FF")
########
# NOTE: Export as 868 x 868 px (square figure)
######

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
