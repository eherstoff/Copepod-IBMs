########################################
#Using data from Helgoland MS to pull parameters for modeling copepod patch use.

#Code last updated @ 2 October, 2018 
########################################

################################################################################################
#Read in all data! #########
data <-read.csv("/Users/emilypetchler/Documents/Grad School Stuff/Summer 2016/Experiment data from Summer 2016!/Combo'd behaviors, EMH expt v3.csv") #Combo'd behaviors, EMH expt.csv
data<- data[-406,] #getting rid of random extra line.

#add on a line for blocks as numbers, for GLM. 
block.num <- rep(1:27, each=15)
data <- cbind(data,block.num)
#Convert Stage, Fed, and Precondition to numeric so summary doesn't freak out
fed.num     <- as.factor(data$Fed)
data <- cbind(data, fed.num)
stage.num   <- as.factor(data$Stage)
data <- cbind(data, stage.num)
precon.num  <- as.factor(data$Precondition)
data <- cbind(data, precon.num)


#Correcting net and total displacement by length of video.
net_displacement_time.corrected   <- data$net_displacement_not.time.corrected / (data$total.time*10) #need to x10 b/c every 10th frame read
total_displacement_time.corrected <- data$total_displacement_not.time.corrected / (data$total.time*10)
data <- cbind(data,net_displacement_time.corrected)
data <- cbind(data,total_displacement_time.corrected)
#Correcting speed data, by x10
avg_velocity_corrected <- 10*data$avg_velocity_from_stepwise
data <- cbind(data, avg_velocity_corrected)
#correcting length of time for video (*10, because every 10th frame)
correct.total.time <- data$total.time * 10
data <- cbind(data, correct.total.time)

colnames(data)
#1] "Video.ID"                      "Precondition"                  "Fed"                          
#[4] "Stage"                         "Block"                         "total.time"                   
#[7] "time.actively.swimming"        "percent.active.swimming.time"  "Xdisplacement"                
#[10] "Ydisplacement"                 "Zdisplacement"                 "Jump.sink"                    
#[13] "Unidirectional.jump"           "Helix.looping"                 "Swimming.drops"               
#[16] "Turning"                       "Net_X"                         "Net_Y"                        
#[19] "Net_Z"                         "V_overall"                     "Jump.time"                    
#[22] "Uni.time"                      "Helix.time"                    "Drop.time"                    
#[25] "Turn.time"                     "Helix.timing"                  "ratio_helix.timing.total.time"
#[28] "net.total"                     "V_overall.corrected"           net_displacement_not.time.corrected"  
#[31] "avg_velocity_from_netmove"     "total_displacement_not.time.corrected" "avg_velocity_from_stepwise"           
#[34] "block.num"                     "avg_velocity_corrected"                "net_displacement_time.corrected"      
#[37] "total_displacement_time.corrected"                        

dev.off() #resets graphical parameters in R plots 


#########
# Separate out the data by stage X food treatments
####  Adult data ############
######################################## 
#### f2 reared adults
data.a<-data[c(1:15),]     # f2 adults fed f2 Rho in trials
avg.a.Net_X       <- mean(data.a$Net_X)
sterr.a.Net_X     <- sd(data.a$Net_X)/(sqrt(15))
avg.a.Net_Y       <- mean(data.a$Net_Y)
sterr.a.Net_Y     <- sd(data.a$Net_Y)/(sqrt(15))
avg.a.Net_Z       <- mean(data.a$Net_Z)
sterr.a.Net_Z     <- sd(data.a$Net_Z)/(sqrt(15))
avg.a.V_overall   <- mean(data.a$V_overall)
sterr.a.V_overall <- sd(data.a$V_overall)/(sqrt(15))
avg.a.PAST        <- mean(data.a$percent.active.swimming.time)
sterr.a.PAST      <- sd(data.a$percent.active.swimming.time)/(sqrt(15)) 
avg.a.Jump.time       <- mean(data.a$Jump.time)
sterr.a.Jump.time     <- sd(data.a$Jump.time)/(sqrt(15))
avg.a.Uni.time       <- mean(data.a$Uni.time)
sterr.a.Uni.time     <- sd(data.a$Uni.time)/(sqrt(15))
avg.a.Helix.time       <- mean(data.a$Helix.time)
sterr.a.Helix.time     <- sd(data.a$Helix.time)/(sqrt(15))
avg.a.Drop.time       <- mean(data.a$Drop.time)
sterr.a.Drop.time     <- sd(data.a$Drop.time)/(sqrt(15))
avg.a.Turn.time       <- mean(data.a$Turn.time)
sterr.a.Turn.time     <- sd(data.a$Turn.time)/(sqrt(15))
mean.net.a            <- mean(data.a$net.total)
sterr.net.a           <- sd(data.a$net.total)/(sqrt(15))
mean.net.displace.a   <- mean(data.a$net_displacement_time.corrected)
sd.net.displace.a     <- sd(data.a$net_displacement_time.corrected)/(sqrt(15))
mean.vel.netmove.a    <- mean(data.a$avg_velocity_from_netmove)
sterr.vel.netmove.a   <- sd(data.a$avg_velocity_from_netmove)/(sqrt(15))
mean.total.displace.a <- mean(data.a$total_displacement_time.corrected)
sd.total.displace.a   <- sd(data.a$total_displacement_time.corrected)/(sqrt(15))
mean.vel.stepwise.a   <- mean(data.a$avg_velocity_from_stepwise)
sd.vel.stepwise.a     <- sd(data.a$avg_velocity_from_stepwise)/(sqrt(15))




data.b<-data[c(16:30),]    # f2 adults fed -N Rho in trials
avg.b.Net_X       <- mean(data.b$Net_X)
sterr.b.Net_X     <- sd(data.b$Net_X)/(sqrt(15))
avg.b.Net_Y       <- mean(data.b$Net_Y)
sterr.b.Net_Y     <- sd(data.b$Net_Y)/(sqrt(15))
avg.b.Net_Z       <- mean(data.b$Net_Z)
sterr.b.Net_Z     <- sd(data.b$Net_Z)/(sqrt(15))
avg.b.V_overall   <- mean(data.b$V_overall)
sterr.b.V_overall <- sd(data.b$V_overall)/(sqrt(15))
avg.b.PAST        <- mean(data.b$percent.active.swimming.time)
sterr.b.PAST      <- sd(data.b$percent.active.swimming.time)/(sqrt(15))
avg.b.Jump.time       <- mean(data.b$Jump.time)
sterr.b.Jump.time     <- sd(data.b$Jump.time)/(sqrt(15))
avg.b.Uni.time       <- mean(data.b$Uni.time)
sterr.b.Uni.time     <- sd(data.b$Uni.time)/(sqrt(15))
avg.b.Helix.time       <- mean(data.b$Helix.time)
sterr.b.Helix.time     <- sd(data.b$Helix.time)/(sqrt(15))
avg.b.Drop.time       <- mean(data.b$Drop.time)
sterr.b.Drop.time     <- sd(data.b$Drop.time)/(sqrt(15))
avg.b.Turn.time       <- mean(data.b$Turn.time)
sterr.b.Turn.time     <- sd(data.b$Turn.time)/(sqrt(15))
mean.net.b            <- mean(data.b$net.total)
sterr.net.b           <- sd(data.b$net.total)/(sqrt(15))
mean.net.displace.b   <- mean(data.b$net_displacement_time.corrected)
sd.net.displace.b     <- sd(data.b$net_displacement_time.corrected)/(sqrt(15))
mean.vel.netmove.b    <- mean(data.b$avg_velocity_from_netmove)
sterr.vel.netmove.b   <- sd(data.b$avg_velocity_from_netmove)/(sqrt(15))
mean.total.displace.b <- mean(data.b$total_displacement_time.corrected)
sd.total.displace.b   <- sd(data.b$total_displacement_time.corrected)/(sqrt(15))
mean.vel.stepwise.b   <- mean(data.b$avg_velocity_from_stepwise)
sd.vel.stepwise.b     <- sd(data.b$avg_velocity_from_stepwise)/(sqrt(15))


data.c<-data[c(31:45),]    # f2 adults fed -P Rho in trials
avg.c.Net_X       <- mean(data.c$Net_X)
sterr.c.Net_X     <- sd(data.c$Net_X)/(sqrt(15))
avg.c.Net_Y       <- mean(data.c$Net_Y)
sterr.c.Net_Y     <- sd(data.c$Net_Y)/(sqrt(15))
avg.c.Net_Z       <- mean(data.c$Net_Z)
sterr.c.Net_Z     <- sd(data.c$Net_Z)/(sqrt(15))
avg.c.V_overall   <- mean(data.c$V_overall)
sterr.c.V_overall <- sd(data.c$V_overall)/(sqrt(15))
avg.c.PAST        <- mean(data.c$percent.active.swimming.time)
sterr.c.PAST      <- sd(data.c$percent.active.swimming.time)/(sqrt(15))
avg.c.Jump.time       <- mean(data.c$Jump.time)
sterr.c.Jump.time     <- sd(data.c$Jump.time)/(sqrt(15))
avg.c.Uni.time       <- mean(data.c$Uni.time)
sterr.c.Uni.time     <- sd(data.c$Uni.time)/(sqrt(15))
avg.c.Helix.time       <- mean(data.c$Helix.time)
sterr.c.Helix.time     <- sd(data.c$Helix.time)/(sqrt(15))
avg.c.Drop.time       <- mean(data.c$Drop.time)
sterr.c.Drop.time     <- sd(data.c$Drop.time)/(sqrt(15))
avg.c.Turn.time       <- mean(data.c$Turn.time)
sterr.c.Turn.time     <- sd(data.c$Turn.time)/(sqrt(15))
mean.net.c            <- mean(data.c$net.total)
sterr.net.c           <- sd(data.c$net.total)/(sqrt(15))
mean.net.displace.c   <- mean(data.c$net_displacement_time.corrected)
sd.net.displace.c     <- sd(data.c$net_displacement_time.corrected)/(sqrt(15))
mean.vel.netmove.c    <- mean(data.c$avg_velocity_from_netmove)
sterr.vel.netmove.c   <- sd(data.c$avg_velocity_from_netmove)/(sqrt(15))
mean.total.displace.c <- mean(data.c$total_displacement_time.corrected)
sd.total.displace.c   <- sd(data.c$total_displacement_time.corrected)/(sqrt(15))
mean.vel.stepwise.c   <- mean(data.c$avg_velocity_from_stepwise)
sd.vel.stepwise.c     <- sd(data.c$avg_velocity_from_stepwise)/(sqrt(15))


#### -N reared adults
data.d<-data[c(46:60),]    # noN adults fed f2 Rho in trials
avg.d.Net_X       <- mean(data.d$Net_X)
sterr.d.Net_X     <- sd(data.d$Net_X)/(sqrt(15))
avg.d.Net_Y       <- mean(data.d$Net_Y)
sterr.d.Net_Y     <- sd(data.d$Net_Y)/(sqrt(15))
avg.d.Net_Z       <- mean(data.d$Net_Z)
sterr.d.Net_Z     <- sd(data.d$Net_Z)/(sqrt(15))
avg.d.V_overall   <- mean(data.d$V_overall)
sterr.d.V_overall <- sd(data.d$V_overall)/(sqrt(15))
avg.d.PAST        <- mean(data.d$percent.active.swimming.time)
sterr.d.PAST      <- sd(data.d$percent.active.swimming.time)/(sqrt(15))
avg.d.Jump.time       <- mean(data.d$Jump.time)
sterr.d.Jump.time     <- sd(data.d$Jump.time)/(sqrt(15))
avg.d.Uni.time       <- mean(data.d$Uni.time)
sterr.d.Uni.time     <- sd(data.d$Uni.time)/(sqrt(15))
avg.d.Helix.time       <- mean(data.d$Helix.time)
sterr.d.Helix.time     <- sd(data.d$Helix.time)/(sqrt(15))
avg.d.Drop.time       <- mean(data.d$Drop.time)
sterr.d.Drop.time     <- sd(data.d$Drop.time)/(sqrt(15))
avg.d.Turn.time       <- mean(data.d$Turn.time)
sterr.d.Turn.time     <- sd(data.d$Turn.time)/(sqrt(15))
mean.net.d            <- mean(data.d$net.total)
sterr.net.d           <- sd(data.d$net.total)/(sqrt(15))
mean.net.displace.d   <- mean(data.d$net_displacement_time.corrected)
sd.net.displace.d     <- sd(data.d$net_displacement_time.corrected)/(sqrt(15))
mean.vel.netmove.d    <- mean(data.d$avg_velocity_from_netmove)
sterr.vel.netmove.d   <- sd(data.d$avg_velocity_from_netmove)/(sqrt(15))
mean.total.displace.d <- mean(data.d$total_displacement_time.corrected)
sd.total.displace.d   <- sd(data.d$total_displacement_time.corrected)/(sqrt(15))
mean.vel.stepwise.d   <- mean(data.d$avg_velocity_from_stepwise)
sd.vel.stepwise.d     <- sd(data.d$avg_velocity_from_stepwise)/(sqrt(15))

data.e<-data[c(61:75),]    # noN adults fed -N Rho in trials
avg.e.Net_X       <- mean(data.e$Net_X)
sterr.e.Net_X     <- sd(data.e$Net_X)/(sqrt(15))
avg.e.Net_Y       <- mean(data.e$Net_Y)
sterr.e.Net_Y     <- sd(data.e$Net_Y)/(sqrt(15))
avg.e.Net_Z       <- mean(data.e$Net_Z)
sterr.e.Net_Z     <- sd(data.e$Net_Z)/(sqrt(15))
avg.e.V_overall   <- mean(data.e$V_overall)
sterr.e.V_overall <- sd(data.e$V_overall)/(sqrt(15))
avg.e.PAST        <- mean(data.e$percent.active.swimming.time)
sterr.e.PAST      <- sd(data.e$percent.active.swimming.time)/(sqrt(15))
avg.e.Jump.time       <- mean(data.e$Jump.time)
sterr.e.Jump.time     <- sd(data.e$Jump.time)/(sqrt(15))
avg.e.Uni.time       <- mean(data.e$Uni.time)
sterr.e.Uni.time     <- sd(data.e$Uni.time)/(sqrt(15))
avg.e.Helix.time       <- mean(data.e$Helix.time)
sterr.e.Helix.time     <- sd(data.e$Helix.time)/(sqrt(15))
avg.e.Drop.time       <- mean(data.e$Drop.time)
sterr.e.Drop.time     <- sd(data.e$Drop.time)/(sqrt(15))
avg.e.Turn.time       <- mean(data.e$Turn.time)
sterr.e.Turn.time     <- sd(data.e$Turn.time)/(sqrt(15))
mean.net.e            <- mean(data.e$net.total)
sterr.net.e           <- sd(data.e$net.total)/(sqrt(15))
mean.net.displace.e   <- mean(data.e$net_displacement_time.corrected)
sd.net.displace.e     <- sd(data.e$net_displacement_time.corrected)/(sqrt(15))
mean.vel.netmove.e    <- mean(data.e$avg_velocity_from_netmove)
sterr.vel.netmove.e   <- sd(data.e$avg_velocity_from_netmove)/(sqrt(15))
mean.total.displace.e <- mean(data.e$total_displacement_time.corrected)
sd.total.displace.e   <- sd(data.e$total_displacement_time.corrected)/(sqrt(15))
mean.vel.stepwise.e   <- mean(data.e$avg_velocity_from_stepwise)
sd.vel.stepwise.e     <- sd(data.e$avg_velocity_from_stepwise)/(sqrt(15))

data.f<-data[c(76:90),]    # noN adults fed -P Rho in trials
avg.f.Net_X       <- mean(data.f$Net_X)
sterr.f.Net_X     <- sd(data.f$Net_X)/(sqrt(15))
avg.f.Net_Y       <- mean(data.f$Net_Y)
sterr.f.Net_Y     <- sd(data.f$Net_Y)/(sqrt(15))
avg.f.Net_Z       <- mean(data.f$Net_Z)
sterr.f.Net_Z     <- sd(data.f$Net_Z)/(sqrt(15))
avg.f.V_overall   <- mean(data.f$V_overall)
sterr.f.V_overall <- sd(data.f$V_overall)/(sqrt(15))
avg.f.PAST        <- mean(data.f$percent.active.swimming.time)
sterr.f.PAST      <- sd(data.f$percent.active.swimming.time)/(sqrt(15))
avg.f.Jump.time       <- mean(data.f$Jump.time)
sterr.f.Jump.time     <- sd(data.f$Jump.time)/(sqrt(15))
avg.f.Uni.time       <- mean(data.f$Uni.time)
sterr.f.Uni.time     <- sd(data.f$Uni.time)/(sqrt(15))
avg.f.Helix.time       <- mean(data.f$Helix.time)
sterr.f.Helix.time     <- sd(data.f$Helix.time)/(sqrt(15))
avg.f.Drop.time       <- mean(data.f$Drop.time)
sterr.f.Drop.time     <- sd(data.f$Drop.time)/(sqrt(15))
avg.f.Turn.time       <- mean(data.f$Turn.time)
sterr.f.Turn.time     <- sd(data.f$Turn.time)/(sqrt(15))
mean.net.f            <- mean(data.f$net.total)
sterr.net.f           <- sd(data.f$net.total)/(sqrt(15))
mean.net.displace.f   <- mean(data.f$net_displacement_time.corrected)
sd.net.displace.f     <- sd(data.f$net_displacement_time.corrected)/(sqrt(15))
mean.vel.netmove.f    <- mean(data.f$avg_velocity_from_netmove)
sterr.vel.netmove.f   <- sd(data.f$avg_velocity_from_netmove)/(sqrt(15))
mean.total.displace.f <- mean(data.f$total_displacement_time.corrected)
sd.total.displace.f   <- sd(data.f$total_displacement_time.corrected)/(sqrt(15))
mean.vel.stepwise.f   <- mean(data.f$avg_velocity_from_stepwise)
sd.vel.stepwise.f     <- sd(data.f$avg_velocity_from_stepwise)/(sqrt(15))


#### -P reared adults
data.g<-data[c(91:105),]    # noP adults fed f2 Rho in trials
avg.g.Net_X       <- mean(data.g$Net_X)
sterr.g.Net_X     <- sd(data.g$Net_X)/(sqrt(15))
avg.g.Net_Y       <- mean(data.g$Net_Y)
sterr.g.Net_Y     <- sd(data.g$Net_Y)/(sqrt(15))
avg.g.Net_Z       <- mean(data.g$Net_Z)
sterr.g.Net_Z     <- sd(data.g$Net_Z)/(sqrt(15))
avg.g.V_overall   <- mean(data.g$V_overall)
sterr.g.V_overall <- sd(data.g$V_overall)/(sqrt(15))
avg.g.PAST        <- mean(data.g$percent.active.swimming.time)
sterr.g.PAST      <- sd(data.g$percent.active.swimming.time)/(sqrt(15))
avg.g.Jump.time       <- mean(data.g$Jump.time)
sterr.g.Jump.time     <- sd(data.g$Jump.time)/(sqrt(15))
avg.g.Uni.time       <- mean(data.g$Uni.time)
sterr.g.Uni.time     <- sd(data.g$Uni.time)/(sqrt(15))
avg.g.Helix.time       <- mean(data.g$Helix.time)
sterr.g.Helix.time     <- sd(data.g$Helix.time)/(sqrt(15))
avg.g.Drop.time       <- mean(data.g$Drop.time)
sterr.g.Drop.time     <- sd(data.g$Drop.time)/(sqrt(15))
avg.g.Turn.time       <- mean(data.g$Turn.time)
sterr.g.Turn.time     <- sd(data.g$Turn.time)/(sqrt(15))
mean.net.g            <- mean(data.g$net.total)
sterr.net.g           <- sd(data.g$net.total)/(sqrt(15))

data.h<-data[c(106:120),]   # noP adults hed -N Rho in trials
avg.h.Net_X       <- mean(data.h$Net_X)
sterr.h.Net_X     <- sd(data.h$Net_X)/(sqrt(15))
avg.h.Net_Y       <- mean(data.h$Net_Y)
sterr.h.Net_Y     <- sd(data.h$Net_Y)/(sqrt(15))
avg.h.Net_Z       <- mean(data.h$Net_Z)
sterr.h.Net_Z     <- sd(data.h$Net_Z)/(sqrt(15))
avg.h.V_overall   <- mean(data.h$V_overall)
sterr.h.V_overall <- sd(data.h$V_overall)/(sqrt(15))
avg.h.PAST        <- mean(data.h$percent.active.swimming.time)
sterr.h.PAST      <- sd(data.h$percent.active.swimming.time)/(sqrt(15))
avg.h.Jump.time       <- mean(data.h$Jump.time)
sterr.h.Jump.time     <- sd(data.h$Jump.time)/(sqrt(15))
avg.h.Uni.time       <- mean(data.h$Uni.time)
sterr.h.Uni.time     <- sd(data.h$Uni.time)/(sqrt(15))
avg.h.Helix.time       <- mean(data.h$Helix.time)
sterr.h.Helix.time     <- sd(data.h$Helix.time)/(sqrt(15))
avg.h.Drop.time       <- mean(data.h$Drop.time)
sterr.h.Drop.time     <- sd(data.h$Drop.time)/(sqrt(15))
avg.h.Turn.time       <- mean(data.h$Turn.time)
sterr.h.Turn.time     <- sd(data.h$Turn.time)/(sqrt(15))
mean.net.h            <- mean(data.h$net.total)
sterr.net.h           <- sd(data.h$net.total)/(sqrt(15))
mean.net.displace.h   <- mean(data.h$net_displacement_time.corrected)
sd.net.displace.h     <- sd(data.h$net_displacement_time.corrected)/(sqrt(15))
mean.vel.netmove.h    <- mean(data.h$avg_velocity_from_netmove)
sterr.vel.netmove.h   <- sd(data.h$avg_velocity_from_netmove)/(sqrt(15))
mean.total.displace.h <- mean(data.h$total_displacement_time.corrected)
sd.total.displace.h   <- sd(data.h$total_displacement_time.corrected)/(sqrt(15))
mean.vel.stepwise.h   <- mean(data.h$avg_velocity_from_stepwise)
sd.vel.stepwise.h     <- sd(data.h$avg_velocity_from_stepwise)/(sqrt(15))


data.i<-data[c(121:135),]   # noP adults fed -P Rho in trials
avg.i.Net_X       <- mean(data.i$Net_X)
sterr.i.Net_X     <- sd(data.i$Net_X)/(sqrt(15))
avg.i.Net_Y       <- mean(data.i$Net_Y)
sterr.i.Net_Y     <- sd(data.i$Net_Y)/(sqrt(15))
avg.i.Net_Z       <- mean(data.i$Net_Z)
sterr.i.Net_Z     <- sd(data.i$Net_Z)/(sqrt(15))
avg.i.V_overall   <- mean(data.i$V_overall)
sterr.i.V_overall <- sd(data.i$V_overall)/(sqrt(15))
avg.i.PAST        <- mean(data.i$percent.active.swimming.time)
sterr.i.PAST      <- sd(data.i$percent.active.swimming.time)/(sqrt(15))
avg.i.Jump.time       <- mean(data.i$Jump.time)
sterr.i.Jump.time     <- sd(data.i$Jump.time)/(sqrt(15))
avg.i.Uni.time       <- mean(data.i$Uni.time)
sterr.i.Uni.time     <- sd(data.i$Uni.time)/(sqrt(15))
avg.i.Helix.time       <- mean(data.i$Helix.time)
sterr.i.Helix.time     <- sd(data.i$Helix.time)/(sqrt(15))
avg.i.Drop.time       <- mean(data.i$Drop.time)
sterr.i.Drop.time     <- sd(data.i$Drop.time)/(sqrt(15))
avg.i.Turn.time       <- mean(data.i$Turn.time)
sterr.i.Turn.time     <- sd(data.i$Turn.time)/(sqrt(15))
mean.net.i            <- mean(data.i$net.total)
sterr.net.i           <- sd(data.i$net.total)/(sqrt(15))
mean.net.displace.i   <- mean(data.i$net_displacement_time.corrected)
sd.net.displace.i     <- sd(data.i$net_displacement_time.corrected)/(sqrt(15))
mean.vel.netmove.i    <- mean(data.i$avg_velocity_from_netmove)
sterr.vel.netmove.i   <- sd(data.i$avg_velocity_from_netmove)/(sqrt(15))
mean.total.displace.i <- mean(data.i$total_displacement_time.corrected)
sd.total.displace.i   <- sd(data.i$total_displacement_time.corrected)/(sqrt(15))
mean.vel.stepwise.i   <- mean(data.i$avg_velocity_from_stepwise)
sd.vel.stepwise.i     <- sd(data.i$avg_velocity_from_stepwise)/(sqrt(15))

###
#Concat all adult data
all.adult.data <- rbind.data.frame(data.a, data.b, data.c, data.d, data.e, data.f, data.g, data.h, data.i)

#Blocking by stage, fed, & precond. #########

adults.all <- rbind(data.a, data.b, data.c, data.d, data.e, data.f, data.g, data.h, data.i)
  adults.fed_f2 <- rbind(data.a, data.d, data.g)
  adults.fed_loN<- rbind(data.b, data.e, data.h)
  adults.fed_loP<- rbind(data.c, data.f, data.i)

########################################
################################################################################################


################################################################################################
########################################
# **ACROSS ALL FOOD TREATMENTS**, 
# how much time (seconds) was spent feeding, moving, or doing nothing, considering only by STAGE? 
# ADULT DATA  ############

#Doing nothing
behaviors <- (mean(all.adult.data$Helix.looping) + mean(all.adult.data$Jump.sink) +
                mean(all.adult.data$Unidirectional.jump) + mean(all.adult.data$Turning) + 
                mean(all.adult.data$Swimming.drops)) #avg. of 7.8 frames were spent doing various behaviors.
avg.length <- mean(all.adult.data$total.time/(.008)) #avg. 30.82963 total frames were watched.

(avg.length-behaviors)/avg.length  #0.7236905 time doing nothing


#Feeding behaviors
# mean raw count of helix (3 frames on avg) / avg. video length
hel <- (all.adult.data$Helix.looping)/(all.adult.data$total.time/(0.008))
  mean(hel) #0.009891753
saw <- (all.adult.data$Jump.sink)/(all.adult.data$total.time/(0.008))     
  mean(saw) #0.0332523
mean(hel)+mean(saw) #0.04314406
  
#movement behaviors
udr <- (all.adult.data$Unidirectional.jump)/(all.adult.data$total.time/(0.008))
  mean(udr) #0.1103752
tns <- (all.adult.data$Turning)/(all.adult.data$total.time/(0.008))
  mean(tns) #0.09224274
drp <- (all.adult.data$Swimming.drops)/(all.adult.data$total.time/(0.008))
  mean(drp) #0.0386318


sd.avg.length <- sd(all.adult.data$total.time/(.008)) #+/- 9.401128 frames
sd.behav <- sd((all.adult.data$Helix.looping) + (all.adult.data$Jump.sink) +
                   (all.adult.data$Unidirectional.jump) + (all.adult.data$Turning) + 
                   (all.adult.data$Swimming.drops)) #+/- 4.084649 frames
  
(sd.avg.length-sd.behav)/sd.avg.length #0.565515 SD doing nothing
sd.behav/sd.avg.length                 #0.434485 SD behaviors  

sd(hel+saw)     #SD Feed:   0.03241378
sd(udr+tns+drp) #SD Move:   0.09194174
  
#StDev in feeding 
sd.hel <- sd(hel) #  0.02161038
sd.saw <- sd(saw) #  0.04530166
#StDev in moving
sd.udr <- sd(udr) #  0.05779602
sd.tns <- sd(tns) #  0.06050043
sd.drp <- sd(drp) #  0.03722947

 
#FOR ADULTS, across all treatments:
  #0.7236905 time was doing nothing.
  #0.2843938 time moving
  #   ==> 0.04314405 time feeding
  #   ==> 0.2412497 time was moving

########################################


###############################################################################################  
# *Animals that were ONLY precon f2! **
    # how much time (seconds) was spent feeding, moving, or doing nothing, considering only by STAGE?   
    #(NOTE: code annotation is heaviest here; see this section for the most details on what I was doing/why)
########################################
# ADULT DATA  ############
adults.precon_f2<- c(data.a, data.b, data.c)
#Doing nothing
#First, determine how many video frames the animals spent doing these focal behaviors.
behaviors <- (mean(adults.precon_f2$Helix.looping) + mean(adults.precon_f2$Jump.sink) +
              mean(adults.precon_f2$Unidirectional.jump) + mean(adults.precon_f2$Turning) + 
              mean(adults.precon_f2$Swimming.drops)) #avg. of 11 frames were spent doing various behaviors.
#Determine how many frames I watched during analysis
avg.length <- mean(adults.precon_f2$total.time/(.008)) #avg. 35.6 total frames were watched.
#Determine how many frames had nothing going on.  
(avg.length-behaviors)/avg.length  #0.6910112 time doing nothing
  

# # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  
#Within the f/2 data, looking at the range of behavior times.
behav.nums <-((adults.precon_f2$Helix.looping) + (adults.precon_f2$Jump.sink) + (adults.precon_f2$Unidirectional.jump) +
                (adults.precon_f2$Turning) + (adults.precon_f2$Swimming.drops))
behav.times <- (adults.precon_f2$total.time/(.008)) 
  
#Determine how many frames had nothing going on.  
nothing.range <- (behav.times-behav.nums)/behav.times  
max(nothing.range)  #0.8518519
min(nothing.range)  #0.5
mean(nothing.range) #0.6776935
sd(nothing.range)   #0.09534529

# # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  
#Feeding behaviors
#Helical swimmins is the behavior most strongly associated with feeding.  
  # Count the number of frames performing a helix; divide by the length of that video.
  hel <- (adults.precon_f2$Helix.looping)/(adults.precon_f2$total.time/(0.008))
  mean(hel) #0.04781568
  max(hel) #0.12 s
  min(hel) #0 s
  sd(hel)  #0.0325157
#Sawtoothing is sometimes associated with feeding in copepods.
  #Count the number of frames where animals performed a jump-sink (i.e., sawtoothing), and divide by length of video
  saw <- (adults.precon_f2$Jump.sink)/(adults.precon_f2$total.time/(0.008))     
  mean(saw) #0.0597406
  sd(saw)   #0.04147533

# # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  
#movement behaviors
#Unidirectional jumps: basically, the animal is swimming in a straight line here.  
  udr <- (adults.precon_f2$Unidirectional.jump)/(adults.precon_f2$total.time/(0.008))
  mean(udr) #0.0829428
#These are the number of frames during which the animal turned.  
  tns <- (adults.precon_f2$Turning)/(adults.precon_f2$total.time/(0.008))
  mean(tns) #0.095618
#Swimming drops could be associated with anti-predator behavior; 
#this was more of a check for me to make sure I wasn't scaring them too much during filming.  
#HOWEVER, this is also a movement the copepod could use to try and quickly escape a patch.
  drp <- (adults.precon_f2$Swimming.drops)/(adults.precon_f2$total.time/(0.008))
  mean(drp) #0.03618942

  
mean(udr)+ mean(tns) +mean(drp) + mean(hel) + mean(saw)  #mean time moving across ALL behaviors was 0.3223065 seconds

# # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # 



 
########################################

################################################################################################
# *Animals that were ONLY precon loN food! **
# how much time (seconds) was spent feeding, moving, or doing nothing, considering only by STAGE?  
########################################
# ADULT DATA  ############
adults.precon_loN<- c(data.d, data.e, data.f)
#Doing nothing
behaviors <- (mean(adults.precon_loN$Helix.looping) + mean(adults.precon_loN$Jump.sink) +
                  mean(adults.precon_loN$Unidirectional.jump) + mean(adults.precon_loN$Turning) + 
                  mean(adults.precon_loN$Swimming.drops)) #avg. of  frames were spent doing various behaviors.
avg.length <- mean(adults.precon_loN$total.time/(.008)) #avg.  total frames were watched.
  
(avg.length-behaviors)/avg.length  # avg time doing nothing = 0.689243
  

# # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # 

#Within the f/2 data, looking at the range of behavior times.
behav.nums <-((adults.precon_loN$Helix.looping) + (adults.precon_loN$Jump.sink) + (adults.precon_loN$Unidirectional.jump) +
                (adults.precon_loN$Turning) + (adults.precon_loN$Swimming.drops))
behav.times <- (adults.precon_loN$total.time/(.008)) 

#Determine how many frames had nothing going on.  
nothing.range <- (behav.times-behav.nums)/behav.times  
max(nothing.range)  #0.8780488
min(nothing.range)  #0.4285714
mean(nothing.range) #0.6832327
sd(nothing.range)   #0.1447196

# # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # 


#Feeding behaviors
  # mean raw count of helix (3 frames on avg) / avg. video length
  hel <- (adults.precon_loN$Helix.looping)/(adults.precon_loN$total.time/(0.008))
  mean(hel) #0.015484
  sd(hel)   #0.0183577
  
  saw <- (adults.precon_loN$Jump.sink)/(adults.precon_loN$total.time/(0.008))     
  mean(saw) #0.03087836
  sd(saw)   #0.03831577
  
# # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # 

  #movement behaviors
  udr <- (adults.precon_loN$Unidirectional.jump)/(adults.precon_loN$total.time/(0.008))
  mean(udr) #0.1382186
  tns <- (adults.precon_loN$Turning)/(adults.precon_loN$total.time/(0.008))
  mean(tns) #0.09714026
  drp <- (adults.precon_loN$Swimming.drops)/(adults.precon_loN$total.time/(0.008))
  mean(drp) #0.270405
    mean(udr)+mean(tns)+mean(drp) #0.2205488

mean(udr)+ mean(tns) +mean(drp) + mean(hel) + mean(saw)  #mean time moving across ALL behaviors was 0.3167673 seconds


# # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # 

########################################


################################################################################################  
#Mean data for swimming speeds, from xxx$avg_velocity_from_stepwise speed calculations ... in mm/s. ######
mean(all.adult.data$avg_velocity_corrected)  #1.731182 mm/s

mean(adults.fed_f2$avg_velocity_corrected)   #1.807389
mean(adults.fed_loN$avg_velocity_corrected)  #1.659302

#StDev Speeds
sd(all.adult.data$avg_velocity_corrected)  #0.9695579 mm/s

sd(adults.fed_f2$avg_velocity_corrected)   #0.9744149
sd(adults.fed_loN$avg_velocity_corrected)  #1.054233

################################################################################################
#Mean displacement, from XXX total_displacement_time.corrected... in mm. ################
mean(all.adult.data$total_displacement_time.corrected)  #1.68484 mm

mean(adults.fed_f2$total_displacement_time.corrected)   #1.792933
mean(adults.fed_loN$total_displacement_time.corrected)  #1.602786

################################################################################################
#Mean data for helix vs total video length (=amount of feeding), from xxx$ratio_helix.timing.total.time  calculations ... ######
mean(all.adult.data$ratio_helix.timing.total.time)  #0.05240004 

mean(adults.fed_f2$ratio_helix.timing.total.time)   #0.1329213
mean(adults.fed_loN$ratio_helix.timing.total.time)  #0.01327578

################################################################################################
#Mean data for helix time (=amount of feeding), from xxx$Helix.timing  calculations ... ######
mean(all.adult.data$Helix.timing*10)  #0.1357037 seconds per helix

mean(adults.fed_f2$Helix.timing*10)   #0.3466667
mean(adults.fed_loN$Helix.timing*10)  #0.03555556

################################################################################################


################################################################################################
################################################################################################


################################################################################################
# Calculating angle and displacements for 3D stuff
################################################################################################
#Using data for only animals that are swimming (all animals that aren't moving were removed from these data)
dat.2 <- read.csv("/Users/emilypetchler/Documents/Grad School Stuff/2018-2019 Seventh Year/IBM based on Helgoland/Calc params/_Model input data for active swimmers, 18 Sept. 2018.csv")
f2.dat2 <- dat.2[c(1:1115),] #all adults precon with f2
noN.dat2 <- dat.2[c(1116:1920),]

mean(f2.dat2$distance)            #0.1701406
  sd(f2.dat2$distance)/sqrt(1613) # SE = 0.006460391
mean(noN.dat2$distance)            #0.2526409
  sd(noN.dat2$distance)/sqrt(1331) #SE = 0.01213405


mean(f2.dat2$velocity)            #21.26757
  sd(f2.dat2$velocity)/sqrt(1613) #SE = 0.8075489
  sd(f2.dat2$velocity)            #SD = 32.43292
  
mean(noN.dat2$velocity)            #31.58008
  sd(noN.dat2$velocity)/sqrt(1331) #SE = 1.516758
  sd(noN.dat2$velocity)            #SD = 55.33568

# Break angles into 1+ x10 equal chunks.
# 0 (swimming straight)
# 0.001 to 36 
#	36 to 72 
#	72 to 108 
#	108 to 144
# 144 to 180
# 180 to 216
# 216 to 252
# 252 to 288	
# 288 to 324	
# 324 to 360

a=length(which(f2.dat2$ANGLE ==0))/1115                             # 0.329148
a2=length(which(f2.dat2$ANGLE > 0.001 & f2.dat2$ANGLE < 36 ))/1115   #0.361435
b=length(which(f2.dat2$ANGLE > 36 & f2.dat2$ANGLE < 72 ))/1115   #0.04304933
c=length(which(f2.dat2$ANGLE > 72 & f2.dat2$ANGLE < 108 ))/1115  #0.05022422
d=length(which(f2.dat2$ANGLE > 108 & f2.dat2$ANGLE < 144 ))/1115 #0.05919283
e=length(which(f2.dat2$ANGLE > 144 & f2.dat2$ANGLE < 180 ))/1115 #0.05201794
f=length(which(f2.dat2$ANGLE > 180 & f2.dat2$ANGLE < 216 ))/1115 #0.0206278
g=length(which(f2.dat2$ANGLE > 216 & f2.dat2$ANGLE < 252 ))/1115 #0.01165919
h=length(which(f2.dat2$ANGLE > 252 & f2.dat2$ANGLE < 288 ))/1115 #0.02690583
i=(1+length(which(f2.dat2$ANGLE > 288 & f2.dat2$ANGLE < 324 )))/1115 #0.02421525
j=length(which(f2.dat2$ANGLE > 324 & f2.dat2$ANGLE < 360 ))/1115 #0.02152466


#When in LQ fod (loN Rho preconditioned animals)
a=length(which(noN.dat2$ANGLE < 0.001))/805                          #0.3875776
a2=length(which(noN.dat2$ANGLE > 0.001 & noN.dat2$ANGLE < 36 ))/805   #0.3614907
b=length(which(noN.dat2$ANGLE > 36 & noN.dat2$ANGLE < 72 ))/805   #0.06832298
c=length(which(noN.dat2$ANGLE > 72 & noN.dat2$ANGLE < 108 ))/805  #0.05714286
d=length(which(noN.dat2$ANGLE > 108 & noN.dat2$ANGLE < 144 ))/805 #0.02857143
e=length(which(noN.dat2$ANGLE > 144 & noN.dat2$ANGLE < 180 ))/805 #0.03478261
f=length(which(noN.dat2$ANGLE > 180 & noN.dat2$ANGLE < 216 ))/805 #0.02236025
g=length(which(noN.dat2$ANGLE > 216 & noN.dat2$ANGLE < 252 ))/805 #0.01614907
h=length(which(noN.dat2$ANGLE > 252 & noN.dat2$ANGLE < 288 ))/805 #0.01118012
i=length(which(noN.dat2$ANGLE > 288 & noN.dat2$ANGLE < 324 ))/805 #0.01118012
j=length(which(noN.dat2$ANGLE > 324 & noN.dat2$ANGLE < 360 ))/805 #0.001242236


#Looking for autocorrelation between points. Does point X determine X+1 position?
plot(dat.2$ANGLE~(dat.2$ANGLE+1))
#get rid of all zeros
new.angles <- which(dat.2$ANGLE==0)
new.angles <- dat.2[(-new.angles),] #getting rid of all zeros
temp<-new.angles[,6]
plot(temp[-1]~temp[1:length(temp)-1]) #1 time step lag
plot(temp[3:length(temp)]~temp[1:(length(temp)-2)]) #2 time step lag
plot(temp[4:length(temp)]~temp[1:(length(temp)-3)]) #3 time step lag
#Autocorrelation fcn
acf(temp, lag.max=50, type=c("partial"))





################################################################################################

#2d angle cals
################################################################################################
#see excel sheet