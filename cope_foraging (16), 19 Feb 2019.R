#######################################################################
# Last updated: 20 Feb 2019, EMH and MU
#Individual based model of adult A. tonsa foraging in patches that differ in stoichiometric quality
#Parameterize using data from EMH's 2016 Helgoland experiments. Focus adults only, looking at movement/behavior by precondition type
#######################################################################
#
#######################################################################
# Creating empty matricies so whole process can be looped
#######################################################################

cell.size.matrix2 <- matrix(nrow=10, ncol=1)
time.unit.matrix2 <- matrix(nrow=10, ncol=1)

#Residence time in HQ patch
mean.res.time_percent.matrix2 <- matrix(nrow=10, ncol=1)
sd.res.time_percent.matrix2   <- matrix(nrow=10, ncol=1)
mean.res.time_n.steps.matrix2 <- matrix(nrow=10, ncol=1)
sd.res.time_n.steps.matrix2   <- matrix(nrow=10, ncol=1)

#Steps until first found HQ patch
mean.first.visit.matrix2 <- matrix(nrow=10, ncol=1)
sd.first.visit.matrix2   <- matrix(nrow=10, ncol=1)
#Mean #copes found HQ patch through time
mean.means.matrix2 <- matrix(nrow=10, ncol=1)
sd.means.matrix2   <- matrix(nrow=10, ncol=1)

#Count helices performed
mean.helix.LQ.dat.matrix2   <- matrix(nrow=10, ncol=1)
sd.helix.LQ.dat.matrix2   <- matrix(nrow=10, ncol=1)
mean.helix.HQ.dat.matrix2   <- matrix(nrow=10, ncol=1)
sd.helix.HQ.dat.matrix2   <- matrix(nrow=10, ncol=1)

#n steps spent helical swimming in HQ and LQ food  
helix.steps.LQ.dat.matrix2        <- matrix(nrow=10, ncol=1)
sd.helix.steps.LQ.dat.matrix2   <- matrix(nrow=10, ncol=1)
helix.steps.HQ.dat.matrix2        <- matrix(nrow=10, ncol=1)
sd.helix.steps.HQ.dat.matrix2   <- matrix(nrow=10, ncol=1)

#n steps straight swimming
straight.LQ.dat.matrix2        <- matrix(nrow=10, ncol=1)
sd.straight.LQ.dat.matrix2   <- matrix(nrow=10, ncol=1)
straight.HQ.dat.matrix2        <- matrix(nrow=10, ncol=1)
sd.straight.HQ.dat.matrix2   <- matrix(nrow=10, ncol=1)

#n steps paused/not moving
waiting.steps.LQ.dat.matrix2        <- matrix(nrow=10, ncol=1)
sd.waiting.steps.LQ.dat.matrix2   <- matrix(nrow=10, ncol=1)
waiting.steps.HQ.dat.matrix2        <- matrix(nrow=10, ncol=1)
sd.waiting.steps.HQ.dat.matrix2   <- matrix(nrow=10, ncol=1)  



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
  
  #time.unit <- mean(c(0.08,0.032,0.032,0.112,0.096,0.104,0.048,0.056,0.184,0.192,0.032,0.048,0.048,0.048,0.032,0.032, 0.056,0.032,0.072,0.032,0.032,0.032,0.064,0.032,0.04) / (10/1250))
  #this is the mean across all the helical swimming times in both HQ and lQ food. It will be used to set time.unit
  
#How long (ms) do they stay in a cell while doing a loop #######
  #This will be used to tell how many steps the simulation needs to keep a cope in a cell when helices happen
  time.hq<-helix.time.hq/time.unit
  time.lq<-helix.time.lq/time.unit 
  
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
  cell.size<-gcd(nearq(mean(f2.dat2$velocity)), nearq(mean(noN.dat2$velocity)))  #No need to resample things--plenty of points to work with (1115 f2, 805 noN)
  #with non-resampled data, cell.size = 0.3125
  
#How many cells do they travel in one time unit? #######
  jump.hq<-  time.unit[1:1115]* (f2.dat2$velocity)/cell.size     #Need to take only n timeunits, as long as data in velocity
  #time.unit * (f2.dat2$velocity)/cell.size  
  jump.lq<- time.unit[1116:1920]* (noN.dat2$velocity)/cell.size 
  #time.unit * (noN.dat2$velocity)/cell.size  
  
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
  
####################
#How long to spend waiting/not moving? ####
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
  grid.size <- 5000 #run at 500  1,000   5,000 px^2
  grid2<-matrix(rep(0,grid.size^2),nrow=grid.size,ncol=grid.size) #the grid is created and filled with LQ food
  
#Add in HQ patches, where model grid.size is 500 px  #####
  grid2[0:354,0:354]<-1    # ~50% % has HQ food; HQ food patch at the corner of the model.
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
  #grid2[0:707,0:707]<-1    # ~50% % has HQ food in 1000 x 1000 px; HQ food patch at the corner of the model.
  #grid2[0:500,0:500]<-1    # 25% has HQ food in 1000 x 1000 px; HQ food patch at the corner of the model.
  #grid2[0:354,0:354]<-1    # ~12.5% has HQ food in 1000 x 1000 px; HQ food patch at the corner of the model.
  #grid2[000:250,000:250]<-1  # 6.25% has HQ food in 1000 x 1000 px; HQ food patch at the corner of the model.
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
#Null model test: parameters for HQ and LQ are the same  #####
  #Cancel this bit when running the actual simulation
  #angles.HQ_XY<-angles.LQ_XY
  #jump.hq<-jump.lq
  #loop.hq<-loop.lq
  
##############
#Setting up matrix######
  timepoints <-100000 #100000 time steps in the model. 
  num.copes <- 100   #1000 number of copes to simulate in the model
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
  count<-matrix(nrow=timepoints,ncol=num.copes)    #System: 0 for wait/stay,10 for keep waiting, 1 for straight runs and 2 for loops; 20= keep looping.
  qual<-matrix(nrow=timepoints,ncol=num.copes)     #0 for LQ and 1 for HQ
  
  for(i in 1:num.copes){  #number of replicates; needs to match ncol in 'Setting up matrix'
    init<-c(sample(1:nrow(grid2),1),sample(1:ncol(grid2),1))  #Initiate from a random cell
    res.i[1,i]<-getind(init,grid2) 
    count[1,i]<-1 #Set initial count of behaviors to 1
    qual[1,i]<-grid2[getind(init,grid2)] #record quality of pixel
    j<-2
    
    while(j<timepoints){ #time steps (must match time steps in 'Setting up matrix' above)  
      rand.num<-runif(1) #generate a random number to determine behavior
      
      if(grid2[init[1],init[2]]==0){  # if you're in LQ patch
        
        if(rand.num<0.62 && rand.num>0.46){ #If helical swimming in LQ; stay put for set number of timesteps
          count[j,i]<-2  #update the behavior counter; 2 is the first step of helix performed
          sample.x <- floor(sample(time.lq, size=1)) #determine how long to helix
          count[(j+1):min(((j+sample.x)-1),timepoints),i]<-20 #add 20's for all subsequent steps performing a helix
          res.i[j:min(((j+sample.x)-1),timepoints),i]<-getind(init,grid2)
          qual[j:min(((j+sample.x)-1),timepoints),i]<-0
          qual[j:min(((j+sample.x)-1),timepoints),i]<-0
          j<-j+sample.x
        }
        
        if(rand.num>=0.62){ #If Straight swim @ angle
          count[j,i]<-1  #Update the counter
          ind<-sample(1:11,1,prob=angles.LQ_XY) #determine what angle to move at
          jump<-sample(jump.lq,size=1) #determine how far to move
          x<-init[1]+floor(jump*sets[ind,1])
          y<-init[2]+floor(jump*sets[ind,2])
          if(x%in%1:nrow(grid2) & y%in%1:ncol(grid2)){init<-c(x,y)} #If they don't swim off the environment edge, record X,Y
          else{init<-c(x%%nrow(grid2),y%%ncol(grid2))} #(When they swim off the edge, wrap their position)
          if(init[1]==0){init[1]<-1} #Wrapper code for dealing with copes at 0 (boundary) => move over one pixel to 1
          if(init[2]==0){init[2]<-1}
          res.i[j,i]<-getind(init,grid2) 
          qual[j,i]<-0
          j<-j+1
        }
        
        if(rand.num<=0.46){ #if wait/paused
          count[j,i]<-0 
          sample.x<- floor(sample(LQ.wait, size=1)) #determine how long to wait
          count[(j+1):min(((j+sample.x)-1),timepoints),i]<-10 #add 10's for all steps during wait 
          res.i[j:min(((j+sample.x)-1),timepoints),i]<-getind(init,grid2)
          qual[j:min(((j+sample.x)-1),timepoints),i]<-0
          j<-j+sample.x 
        }
      }#end of LQ patch statement
      
      # Second part of the code, for HQ food #######
      if(grid2[init[1],init[2]] == 1){   #if they are in HQ food patch...
        
        if(rand.num>0.22 && rand.num<0.5){  #helical swimming
          count[j,i]<-2  #update the behavior counter
          sample.x <- floor(sample(time.hq, size=1)) #determine how long to loop
          count[(j+1):min(((j+sample.x)-1),timepoints),i]<-20 #add 20's for all steps performing a helix
          res.i[j:min((j+sample.x-1),timepoints),i]<-getind(init,grid2)
          qual[j:min(((j+sample.x)-1),timepoints),i]<-1
          j<-j+sample.x
        }
        
        
        if(rand.num>=0.5){ #If Straight swim @ angle
          count[j,i]<-1  #Update the counter
          ind<-sample(1:11,1,prob=angles.HQ_XY) #what angle to swim
          jump<-sample(jump.hq,size=1) #how far to swim
          x<-init[1]+floor(jump*sets[ind,1])
          y<-init[2]+floor(jump*sets[ind,2])
          if(x%in%1:nrow(grid2) & y%in%1:ncol(grid2)){init<-c(x,y)} #wrapper for tiled environment
          else{init<-c(x%%nrow(grid2),y%%ncol(grid2))}
          if(init[1]==0){init[1]<-1}
          if(init[2]==0){init[2]<-1}
          res.i[j,i]<-getind(init,grid2)
          qual[j,i]<-1
          j<-j+1
        }
        
        if(rand.num<=0.22){ #if paused/waiting
          count[j,i]<-0
          sample.x<-floor(sample(HQ.wait, size=1)) #How long to wait in HQ
          count[(j+1):min(((j+sample.x)-1),timepoints),i]<-10 #add 10's for all steps performing a wait
          res.i[j:min(((j+sample.x)-1),timepoints),i]<-getind(init,grid2)
          qual[j:min(((j+sample.x)-1),timepoints),i]<-1
          j<-j+sample.x 
        }  
      }#end HQ wait
    }#end of while
  }#end of master for
  
  
  res.i<-res.i[-nrow(res.i),]
  res.i<-res.i[-1,]
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
# Matrix of 'count' (2 or 20 = helix; 1 straight; 0 or 10 for wait); and 'qual (0 for LQ and 1 for HQ patch)
#Determine how many timepoints copes spend in HQ vs LQ food #########
  helix.LQ.dat<-matrix(nrow=num.copes, ncol=1)
  helix.HQ.dat<-matrix(nrow=num.copes, ncol=1)
  
  helix.steps.LQ.dat<-matrix(nrow=num.copes, ncol=1)
  helix.steps.HQ.dat<-matrix(nrow=num.copes, ncol=1)
  
  straight.LQ.dat <- matrix(nrow=num.copes, ncol=1)
  straight.HQ.dat <- matrix(nrow=num.copes, ncol=1)
  
  waiting.steps.LQ.dat <- matrix(nrow=num.copes, ncol=1)
  waiting.steps.HQ.dat <- matrix(nrow=num.copes, ncol=1)
  
  
  for (i in 1:num.copes){
    helix.LQ.dat[i,]= length(which(count[,i]==2&qual[,i]==0)) #total number of helices in LQ food
    helix.HQ.dat[i,]= length(which(count[,i]==2&qual[,i]==1)) #... helix in HQ food
    
    helix.steps.LQ.dat[i,]= length(which(count[,i]==2|20&qual[,i]==0)) #n timesteps spent helical swimming of helix in LQ food
    helix.steps.HQ.dat[i,]= length(which(count[,i]==2|20&qual[,i]==1)) #... in HQ food
    
    straight.LQ.dat[i,]= length(which(count[,i]==1&qual[,i]==0)) #number of straight runs in LQ food
    straight.HQ.dat[i,]= length(which(count[,i]==1&qual[,i]==1)) #... in HQ food
    
    waiting.steps.LQ.dat[i,]= length(which(count[,i]==0|10&qual[,i]==0)) #n steps spent paused/waiting in LQ food
    waiting.steps.HQ.dat[i,]= length(which(count[,i]==0|10&qual[,i]==1)) #... in HQ food
    
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
  #Mean number of copes found HQ patch through time
  mean.means.matrix2[z,] <-mean(means)
  sd.means.matrix2[z,] <-sd(means)
  
  #Count behaviors
  #number of helices performed
  mean.helix.LQ.dat.matrix2[z,]   <- mean(helix.LQ.dat)
  sd.helix.LQ.dat.matrix2[z,]   <- sd(helix.LQ.dat)
  mean.helix.HQ.dat.matrix2[z,]   <- mean(helix.HQ.dat)
  sd.helix.HQ.dat.matrix2[z,]   <- sd(helix.HQ.dat)
  #n steps spent helical swimming in HQ and LQ food  
  helix.steps.LQ.dat.matrix2[z,]        <- mean(helix.steps.LQ.dat)
  sd.helix.steps.LQ.dat.matrix2[z,]   <- sd(helix.steps.LQ.dat)
  helix.steps.HQ.dat.matrix2[z,]        <- mean(helix.steps.HQ.dat)
  sd.helix.steps.HQ.dat.matrix2[z,]   <- sd(helix.steps.HQ.dat)
  #n steps straight swimming
  straight.LQ.dat.matrix2[z,]        <- mean(straight.LQ.dat)
  sd.straight.LQ.dat.matrix2[z,]   <- sd(straight.LQ.dat)
  straight.HQ.dat.matrix2[z,]        <- mean(straight.HQ.dat)
  sd.straight.HQ.dat.matrix2[z,]   <- sd(straight.HQ.dat)
  #n steps paused/not moving
  waiting.steps.LQ.dat.matrix2[z,]        <- mean(waiting.steps.LQ.dat)
  sd.waiting.steps.LQ.dat.matrix2[z,]   <- sd(waiting.steps.LQ.dat)
  waiting.steps.HQ.dat.matrix2[z,]        <- mean(waiting.steps.HQ.dat)
  sd.waiting.steps.HQ.dat.matrix2[z,]   <- sd(waiting.steps.HQ.dat) 
  
  
  
  # finish the grand for loop with a progress bar and adding to z  ######
  #create a progress bar
  print(z)
  Sys.sleep(0.01)
  flush.console()
  
  z=z+1  }
#####
