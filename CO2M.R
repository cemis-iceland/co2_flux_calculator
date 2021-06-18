#### inputDir = path to data
####first_m   = how many measurment to exclude at the start of the measurment
####d         = diameter chamber [cm]
####h         = height chamber [cm]
####vadd      = additional volume through collar etc [m^3]
CO2M <- function(inputDir,first_m,d,h,vadd,span){
  files <- list.files(inputDir, pattern = '.csv$', full.names = F)
  for (i in 1:length(files)) {
    
    #import data skip first row
    dataT <- read.csv(paste(inputDir,files[i],sep=""),header=FALSE,sep=",",skip=2)
    #delete empty columns
    dataT <- dataT[colSums(!is.na(dataT)) > 0]
    colnames(dataT) <- c("time","CO2","temp","hum")
    #dataT$time = gsub("\\:","",dataT$time )
    #dataT$time = gsub("\\:","",dataT$time )
    
        #t <- dataT[ ,1]
    #split character and save relevant part
    #only_t <- substring(dataT[ ,1],1,5)
    only_t <- as.numeric(gsub(":","",str_match(dataT$time, "[0-9]*[0-9]:[0-9][0-9].[0-9]*")))
    dataT$timeT <- c(0,cumsum(abs(diff(only_t))))
    filename=files[i]
    sample_name=gsub("\\.csv.*","",filename)
    ######
    loessMod10 <- loess(dataT$CO2 ~ dataT$timeT, data=dataT, span=span)
    dataT$smoothed10 <- predict(loessMod10)
    dir.create(file.path(inputDir,sample_name)) 
    outputDir <- file.path(inputDir,sample_name)
    png(filename = paste(outputDir,(gsub(".csv","_Smoothing.png",filename)),sep="/"),width = 780, height = 780, units = "px", pointsize = 12)#overall plot
    plot(dataT$CO2,x= dataT$timeT, type="l", main="Loess Smoothing and Prediction", xlab="Time", 
         ylab="CO2 level") +
      lines(dataT$smoothed10, x=dataT$timeT, col="red")
    dev.off()
    png(filename = paste(outputDir,(gsub(".csv","_local_min_max.png",filename)),sep="/"),width = 780, height = 780, units = "px", pointsize = 12)#overall plot
    lmm <- local.min.max(dataT$smoothed10, dev=mean, add.points=T,plot = T, 
                         main="Local Minima and Maxima")
    dev.off()
    ######dataT$smoothed10[1] 1st measurment, if you need to start at n measurment just mod to dataT$smoothed10[x]
    lmm_minima <- c(dataT$smoothed10[1],lmm$minima)
    split_list <- list()
    for (i in 1:length(lmm$maxima)) {
      split_list[[i]] <- dataT[(which(grepl(lmm_minima[i], dataT$smoothed10))+first_m):which(grepl(lmm$maxima[i], dataT$smoothed10)),]
    }
    new_names <- paste(sample_name,paste("Measurment",seq(1:length(lmm$maxima)),sep=""),sep="_")
    for (i in 1:length(lmm$maxima)) {
      names(split_list)[i] <- new_names[i]
    }
    for (i in 1:length(lmm$maxima)) {
      split_list[[i]]$timeTII <- c(0,cumsum(abs(diff(split_list[[i]]$timeT))))
    }
    png(filename = paste(outputDir,(gsub(".csv","_Regressions.png",filename)),sep="/"),width = 780, height = 780, units = "px", pointsize = 12)#overall plot
    par(mfrow = c(2, round(length(split_list)/2)))
    
    dat = data.frame(filename=rep(0, length(split_list)),
                     sample_name=rep(0, length(split_list)),
                     hum_mean=rep(0, length(split_list)),
                     T_mean=rep(0, length(split_list)),
                     lm_intercept=rep(0, length(split_list)),
                     lm_slope=rep(0, length(split_list)),
                     adjR2=rep(0, length(split_list)),
                     fluxmol=rep(0, length(split_list)),
                     fluxumol=rep(0, length(split_list)))
#For each list that we split we do a regression
    for (i in 1:length(split_list)) {
      mod1_lm<-lm(split_list[[i]]$CO2 ~ split_list[[i]]$timeTII, data=split_list[[i]]) #linear regression; "~" means "is a variable of"
      summary(mod1_lm)
      #get coefficient(slope) and R^2
      lm_intercept <- summary(mod1_lm)$coefficients[1]
      #lm_intercept <- unname(lm_intercept)
      lm_slope <- summary(mod1_lm)$coefficients[2]
      #lm_slope <- unname(lm_slope)
      adjR2 <- summary(mod1_lm)$adj.r.squared
      T_mean <- mean(split_list[[i]]$temp) #average temp. [C?]
      hum_mean <- mean(split_list[[i]]$hum) #average hum. 
      ###calculate flux
      #define variables 
      d  <- d                  #diameter chamber [cm]
      h  <- h                  #height chamber [cm]
      vc  <- (pi*(d/2)^2*h)/ 10^6    #Volume [m^3]: area[cm^2](pi*diameter/4^2)*height[cm]*10^-6
      v2 <- 22.4*10^-3          #molar volume for ideal gas at 273 K [m^3/mol]
      vadd <- vadd              #additional volume through collar etc [m^3]
      v <- vc+vadd              #total volume [m^3]
      T1 <- T_mean +273.15 #average air temperature in relevant time window [K]
      T2 <- 273.15              #standard temperature [K]
      A  <- pi*(d/2)^2* 10^-4       #footprint area chamber [m^2]
      m  <- lm_slope           #delta C/delta t [ppm/s]
      
      #calculation flux rate
      fluxmol=m*10^-6*(v/(v2*(T2/T1)))*1/A #soil CO2 flux [mol/(m^2*s)]
      fluxumol=fluxmol*10^6                #soil CO2 flux [umol/(m^2*s)]
      dat[i,] <- c(filename,names(split_list[i]),hum_mean,T_mean,lm_intercept,lm_slope,adjR2,fluxmol,fluxumol)
      write.csv(dat,paste(outputDir,(gsub(".csv","_output.csv",filename)),sep="/"))
      ###
      plot(split_list[[i]]$timeTII,split_list[[i]]$CO2,type="l",col="black",panel.first=grid(nx = length(seq(1, length(split_list[[i]]$timeTII), by = 1))), 
           main=paste("linear regression M",i,sep = ""),xlab = "Time (sec)",axes = F,ylab = "CO2 level") +
        box() 
      axis(2)
      abline(lm(split_list[[i]]$CO2~split_list[[i]]$timeTII),col = "red")
      legend("bottomright", bty="n", legend=c(paste("lm (red) R2=", format(summary(mod1_lm)$adj.r.squared, digits=5))))
      ## rounded coefficients for better output
      cf <- round(coef(mod1_lm), 4) 
      ## sign check to avoid having plus followed by minus for negative coefficients
      eq <- paste0("CO2 = ", cf[1],
                   ifelse(sign(cf[2])==1, " + ", " - "), abs(cf[2]), " x sec")
      ## printing of the equation
      mtext(eq, 3, line=-2,cex = 0.8)
      #model 2 compared to actual
      seq=1:split_list[[i]]$timeTII[length(split_list[[i]]$timeTII)]
      axis(side = 1, at=seq[c(TRUE, FALSE)])
    }
    dev.off()
  }}


