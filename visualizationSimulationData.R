library(Luminescence);library(viridis)
type = c('Mounding', 'Mixing')
func = c('Gradational', 'Exponential','Abrupt')

ma <- function(x, n = 7){filter(x, rep(1 / n, n), sides = 2)}

# Plot distributions ####
aggr.layers = 5
cols=viridis(3)
setwd('C:/Users/Marijn/sciebo/Research/Bioturbation/Bioturbation_code/')
png('Figures/ageDistributions_depthFunctions.png',width = 30, height = 15, res = 500, units = 'cm')
{
  par(mfrow = c(1,2))

  for(t in type)
  {
    plot(0,pch=NA, xlim = c(0,10000), ylim = c(1.5,0),xlab = 'Age [a]',ylab = ' Depth [m]')
    
    for (f in (1:length(func)-1))
    {
      leg.order = c(1,4,2,5,3,6)
      if(t==type[1])
      {
        legend('bottom', 
               lty = c(NA,NA,NA,1,1,2),
               pch = c(15,15,15,NA,NA,NA),
               legend = c(func, 'Modes','Bioturbated fraction','Distributions'),
               col = c(cols[c(1:3)],1,1,1),
               lwd = c(NA,NA,NA,3,1,1),
               ncol = 2,
               bty = 'n')
        legend('topright', legend = 'A', bty = 'n')
      }
      else
      {
        legend('topright', legend = 'B', bty = 'n')
      }
      
      axis(side = 3, at = c(0,10000),labels = c(0,1))
      mtext(side = 3,'Bioturbated fraction [-]', line = 1)
      title(t, line = 2.6)
      
      agedepth = matrix(nrow = 200, ncol = 4)
      soil = read.csv(paste0('C:/Simulations/Bioturbation/0_Final/',t,'_diffFunctions/',f,'_10000_out_allsoils.csv'))
      ages =  read.csv(paste0('C:/Simulations/Bioturbation/0_Final/',t,'_diffFunctions/',f,'_10000_out_OSL_ages.csv'))
      
      layers= unique(soil$nlayer)
      for(l in 0:(length(layers)/aggr.layers-1))
      {
        lay.sel = seq(l*aggr.layers, ((l+1)*aggr.layers-1))
        soil.sel = soil[soil$nlayer%in%lay.sel,]
        depth = soil.sel$cumth_m[nrow(soil.sel)] - sum(soil.sel$thick_m)/2
        age = ages[ages$layer%in%lay.sel,]
        fbio = 1-length(which(age$grain_age==10000))/nrow(age)
        count = nrow(age)
        age = age[age$grain_age<10000,]
        if(nrow(age)>1 & length(unique(age$grain_age))>3)
        {
          dens = density(age$grain_age,bw='SJ')
          mode = dens$x[dens$y==max(dens$y)][1]
          dens$y = dens$y/max(dens$y)
          dens$y[dens$y<0.05] = NA
          dens$y = (dens$y-0.05)/.95
          dens$y = -dens$y*0.05 + depth
          IQR=quantile(age$grain_age,.75)-quantile(age$grain_age,.25)
          # points(rep(depth, nrow(age))~age$grain_age, pch=16, col=adjustcolor('gray34', alpha.f = 0.025))
          if(f==0){lines(dens$y ~dens$x, col = cols[f+1],lty=2)}
        }
        else
        {
          mode = NA
          IQR=NA
        }
        agedepth[l+1,] = c(depth,mode,fbio,IQR)
      }
     
      lines(y = agedepth[,1], 
            x = (agedepth[,2]), col = cols[f+1],lwd = 3, lty = 1) # Modes
      
      lines(y=agedepth[,1], x=(agedepth[,3]*10000), col = cols[f+1], lty = 1) # fbio
      # lines(y=agedepth[,1], x=(agedepth[,4]), col = cols[f+1], lty = 2) # IQR
    }
  }
}; dev.off()


# Plot mix of mounding and mixing ####
# With linear depth function
perc = c(0
         ,0.25,0.5,0.75,.8,.9,.95,.99,1)
aggr.layers = 5
ageStatistics = array(dim = c(length(perc),200,10))
names_ageStatistics = c('depth','mean','median','sd','mode','rsd','IQR','skewness','fbio', 'OD')

for(i in 0:(length(perc)-1))
{
  soil = read.csv(paste0('C:/Simulations/Bioturbation/0_Final/ratioBetweenProcesses/',i,'_10000_out_allsoils.csv'))
  ages =  read.csv(paste0('C:/Simulations/Bioturbation/0_Final/ratioBetweenProcesses/',i,'_10000_out_OSL_ages.csv'))
  layers = unique(soil$nlayer)
  for(l in 0:(length(layers)/aggr.layers-1))
  {
    lay.sel = seq(l*aggr.layers, ((l+1)*aggr.layers-1))
    soil.sel = soil[soil$nlayer%in%lay.sel,]
    depth = soil.sel$cumth_m[nrow(soil.sel)] - sum(soil.sel$thick_m)/2
    age = ages[ages$layer%in%lay.sel,]

    age$grain_age = age$grain_age/max(age$grain_age)
    
    fbio = 1-length(which(age$grain_age==1))/nrow(age)
    count = nrow(age)
    age = age[age$grain_age<1,]
    if(nrow(age)>1 & l !=max(soil$nlayer) & length(unique(age$grain_age))>1)
    {
      # dens = KDE(age$grain_age)
      dens = density(age$grain_age,bw='SJ',na.rm=T)
      mode = dens$x[dens$y==max(dens$y)][1]
      mean = mean(age$grain_age)
      median = median(age$grain_age)
      sd = sd(age$grain_age)
      OD = calc_CentralDose(data.frame(age = age$grain_age, err = rep(0,nrow(age))),sigmab = 0,plot=F,verbose=F)$summary$rel_OD
      rsd = sd/mean
      IQR = quantile(age$grain_age,0.75) - quantile(age$grain_age,0.25)
      skewness = 3 * (mean(age$grain_age)-median(age$grain_age))/sd(age$grain_age) # alternative Pearson Mode skewness
      
      ageStatistics[i + 1, l+1,] = c(depth,mean,median,sd,mode,rsd,IQR,skewness,fbio,OD)
    }
  }
  # ageStatistics[i+1,,1] = ageStatistics[i+1,,1] / max(ageStatistics[i+1,,1],na.rm=T)
}
# ageStatistics[,,1] = ageStatistics[,,1]/max(ageStatistics[,,1], na.rm = T) # Normalize depth
ageStatistics[,,1][ageStatistics[,,1]>1] = NA
statisticsExperimentalData = read.csv('C:/Users/Marijn/sciebo/Research/Bioturbation/Data/distributionStatistics.csv')


stat_indices = c(5,7,9)
labels = c('A','B','C')
setwd("C:/Users/Marijn/sciebo/Research/Bioturbation/Bioturbation_code/")

png('Figures/statistics_ratioBetweenProcesses.png', height = 20, width = 20, units = 'cm', res = 500)
{
  par(mfrow = c(2,2))
  par(mar = c(4.1,4.1,3.1,1.1))
  par(oma = c(1.1,1.1,1.1,1.1))
  par(xpd = F)
  
  
  exp_datasets = c('De_termites', 'De_worms', 'De_ants')
  exp_plot = statisticsExperimentalData[statisticsExperimentalData$dataset%in%exp_datasets,]
  colors = grey.colors(length(perc),0.9,0.1)
  # for(i in 2:dim(ageStatistics)[3])
  ct=0
  for(i in stat_indices)
  {
    ct=ct+1
    ylab= "ẑ"
    xlab= "â"; if(i>=8){xlab="-"}
    plot(0,pch = NA, xlim = range(c(ageStatistics[,,i],exp_plot[,i]),na.rm = T),
         ylim = c(1,0), main = paste0(labels[ct],': ',names_ageStatistics[i]),
         ylab = ylab, xlab=xlab, xpd = T,
         cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
    for(j in 1:dim(ageStatistics)[1])
    {
      lines(ageStatistics[j,,1]~ageStatistics[j,,i], col = colors[j],lwd=2)
    }
    for(d in exp_datasets)
    {
      dd = exp_plot[exp_plot$dataset==d,]
      lines(dd[,1]~dd[,i], type = 'b', pch = which(exp_datasets%in% d),
            col = cols[which(exp_datasets%in% d)],
            lty = 2,cex=2,lwd=2)
    }
    
  }
  plot(0,pch=NA,axes=F,xlab="",ylab="")
  
  mtext(side = 3, line = -1, text = 'Experimental data',adj=0,cex=1.2)
  legend('topleft', legend = c(NA,substr(exp_datasets,4,nchar(exp_datasets))),
         pch = c(NA,1:length(exp_datasets)),
         col = c(NA,cols),bty = 'n', lty = 2,cex = 1.2)
  
  mtext(side = 1, line = -6.5, text = 'Simulated data\nFraction subsurface mixing',adj=0, cex = 1.2)
  legend('bottomleft', ncol = 3, legend = perc, lty = 1, col = colors, bty = 'n',cex=1.2)
};dev.off()








# Upheaval  ####
years = c(1:9,seq(10,90,10),seq(100,900,100),seq(1000,5000,1000))

NSF_modes = matrix(nrow = 4*length(years),ncol = 5)
par(mfrow = c(1,1))
count1=0
count_NSF = 0
mds=c(10,25,50,75)
for(md in mds)
{
  count1 = count1+1
  for(t in c(0:(length(years)-1)))
  {
    count_NSF = count_NSF + 1
    soil = read.csv(paste0('upheaval_',md,'cm/',t,'_10000_out_allsoils.csv'))
    ages = read.csv(paste0('upheaval_',md,'cm/',t,'_10000_out_OSL_ages.csv'))
    ages = ages[ages$layer<md,]

    fbio = 1-nrow(ages[ages$grain_age==10000,])/nrow(ages)
    ages = ages[ages$grain_age<10000 & ages$grain_age>0,]
    dens = density(ages$grain_age,bw='SJ')
    dens$y = dens$y/max(dens$y)
    mode = dens$x[dens$y==max(dens$y)]
    IQR = quantile(ages$grain_age,.75)-quantile(ages$grain_age,.25)
    NSF_modes[count_NSF,] = c(md, years[t+1], fbio, mode,IQR)
  }
}
NSF_modes = data.frame(NSF_modes)
names(NSF_modes) = c('md', 't', 'fbio', 'mode', 'IQR')

cols=viridis(4)
png('C:/Users/Marijn/sciebo/Research/Bioturbation/Bioturbation_code/Figures/resultsUpheaval.png',width = 15, height = 15, units = 'cm', res = 500)
{
  par(fig = c(0,1,0,1))
  par(mar = c(5.1,4.1,2.1,4.1))
  plot(0,pch=NA,xlab = 'Upheaval frequency [a]',ylab = 'Age [a]',
       xlim=c(0,1000),ylim=c(0,8000))
  axis(side=4,at = c(0,4000,8000),labels = c(0,0.5,1))
  mtext(side = 4, line = 2.6, text = 'Bioturbated fraction [-]')
  ct=0
  for(l in unique(NSF_modes$md))
  {
    ct=ct+1
    temp = NSF_modes[NSF_modes$md==l,]
    lines(IQR~t, temp,col = cols[ct])
    lines(x=temp$t, y=temp$fbio*8000,col = cols[ct],lty=2)
  }
  
  points(mode~t,NSF_modes[NSF_modes$t<=5000,],col = cols[ceiling(NSF_modes$md/25+.1)],
  pch = ceiling(NSF_modes$md/25+.1))
  
  legend('topright',col = c(cols,1,1),legend = c(paste('md:',mds,'cm'),'IQR', "fbio"),bty='n',
         lty=c(NA,NA,NA,NA,1,2),pch = c(1,2,3,4,NA,NA))

};dev.off()



setwd('C:/Simulations/Bioturbation/')
years = c(1:9,seq(10,90,10),seq(100,900,100),seq(1000,5000,1000))
years_ind = c(1,10,19,28)

NSF_modes = matrix(nrow = 4*length(years),ncol = 5)
par(mfrow = c(1,1))
count1=0
count_NSF = 0
mds=c(10,25,50,75)
setwd('C:/Simulations/Bioturbation/')
labs = c('A','B')
png('C:/Users/Marijn/sciebo/Research/Bioturbation/Bioturbation_code/Figures/resultsUpheaval_depthplots.png',width = 30, height = 15, units = 'cm', res = 500)
{
  par(mfrow = c(1,2))

  mds = c(25,75)
  cols = viridis(length(years_ind))
  for(md in mds)
  {
    plot(0,pch=NA,xlim=c(0,10000),ylim=c(1,0),   xlab='Age [a]', ylab='Depth [m]')
    axis(side = 3, at = c(0,10000),labels = c(0,1))
    mtext(side = 3,'Bioturbated fraction [-]', line = 1)
    title(paste0(labs[which(mds%in%md)],': mixing depth: ',md,' cm'), line = 2.6)
    
    count1 = count1+1
    # for(t in c(0:(length(years)-1)))
    for(t in years_ind-1)
    {
      # count_NSF = count_NSF + 1
      soil = read.csv(paste0('upheaval_',md,'cm/',t,'_10000_out_allsoils.csv'))
      ages = read.csv(paste0('upheaval_',md,'cm/',t,'_10000_out_OSL_ages.csv'))
      agedepth = matrix(nrow = 200, ncol = 5)
      
      layers= unique(soil$nlayer)
      for(l in 0:(length(layers)/aggr.layers-1))
      {
        lay.sel = seq(l*aggr.layers, ((l+1)*aggr.layers-1))
        soil.sel = soil[soil$nlayer%in%lay.sel,]
        depth = soil.sel$cumth_m[nrow(soil.sel)] - sum(soil.sel$thick_m)/2
        age = ages[ages$layer%in%lay.sel,]
        fbio = 1-length(which(age$grain_age==10000))/nrow(age)
        count = nrow(age)
        age = age[age$grain_age<10000,]
        if(nrow(age)>1 & length(unique(age$grain_age))>3)
        {
          dens = density(age$grain_age,bw='SJ')
          mode = dens$x[dens$y==max(dens$y)][1]
          dens$y = dens$y/max(dens$y)
          dens$y[dens$y<0.05] = NA
          dens$y = (dens$y-0.05)/.95
          dens$y = -dens$y*0.05 + depth
          IQR = quantile(age$grain_age,.75) - quantile(age$grain_age,.25)
          
          # points(rep(depth, nrow(age))~age$grain_age, pch=16, col=adjustcolor('gray34', alpha.f = 0.025))
          # lines(dens$y ~dens$x, col = cols[which(years_ind%in%(t+1))],lty=2)
        }
        else
        {
          mode = NA
          IQR=NA
        }
        agedepth[l+1,] = c(depth,mode,fbio,IQR,count)
      }
      agedepth = agedepth[!is.na(agedepth[,1]),]
      
      lines(y = agedepth[,1], 
            x = (agedepth[,2]), col = cols[which(years_ind%in%(t+1))],lwd = 3, lty = 1) # Modes
      
      lines(y=agedepth[,1], x=(agedepth[,3]*10000), col = cols[which(years_ind%in%(t+1))], lty = 1) # fbio
      # lines(y=agedepth[,1], x=(agedepth[,4]), col = cols[which(years_ind%in%(t+1))], lty = 2) # IQR
    }
    if(md==25)
    {
      legend('bottomright', 
             legend = c(paste(years[years_ind],'a'),'Modes','Bioturbated',' fraction'),
             col=c(cols,1,1),
             lty=c(rep(NA,length(years_ind)),1,1,NA),
             pch=c(rep(15,length(years_ind)),NA,NA,NA),
             lwd=c(rep(1,length(years_ind)),3,1), bty='n',ncol=2)
      text('Upheaval frequency',x=6000,y=1.15)
    }
  }

};dev.off()




# 
# 
# # age_distributions
# aggr.layers = 5
# png('Figures/statistics_ratioBetweenProcesses_allProfiles_SJ.png', height = 60, width = 60, units = 'cm', res = 500)
# par(mfrow = c(3,3))
# for(i in 0:(length(perc)-1))
# {
#   agedepth = matrix(nrow = 200, ncol = 4)
#   
#   soil = read.csv(paste0('C:/Simulations/Bioturbation/0_Final/ratioBetweenProcesses/',i,'_10000_out_allsoils.csv'))
#   ages =  read.csv(paste0('C:/Simulations/Bioturbation/0_Final/ratioBetweenProcesses/',i,'_10000_out_OSL_ages.csv'))
#   
#   plot(0,pch=NA, xlim = c(0,10000), ylim = c(1.5,0),xlab = 'Age [a]',ylab = ' Depth [m]', main = paste(perc[i+1]*100,"% mixing"))
#   
#   layers= unique(soil$nlayer)
#   for(l in 0:(length(layers)/aggr.layers-1))
#   {
#     lay.sel = seq(l*aggr.layers, ((l+1)*aggr.layers-1))
#     soil.sel = soil[soil$nlayer%in%lay.sel,]
#     depth = soil.sel$cumth_m[nrow(soil.sel)] - sum(soil.sel$thick_m)/2
#     age = ages[ages$layer%in%lay.sel,]
#     fbio = 1-length(which(age$grain_age==10000))/nrow(age)
#     count = nrow(age)
#     age = age[age$grain_age<10000,]
#     if(nrow(age)>1 & length(unique(age$grain_age))>3)
#     {
#       # dens = KDE(age$grain_age)
#       dens = density(age$grain_age, bw = 'SJ')
#       # dens = density(age$grain_age, bw = 'nrd0')
#       mode = dens$x[dens$y==max(dens$y)][1]
#       dens$y = dens$y/max(dens$y)
#       dens$y[dens$y<0.05] = NA
#       dens$y = -dens$y*0.05 + depth
#       
#       points(rep(depth, nrow(age))~age$grain_age, pch=16, col=adjustcolor('black', alpha.f = 0.01))
#       if(l%%3==0){lines(dens$y ~dens$x, col = 1,lty=1)}
#     }
#     else
#     {
#       mode = NA
#     }
#     agedepth[l+1,] = c(depth,mode,fbio,count)
#   }
#   lines(agedepth[,1]~agedepth[,2], lwd = 1.5, col = 2)
#   lines(y=agedepth[,1], x=(agedepth[,3]*10000), lwd = 1.5, col = 2, lty = 2)
# }
# dev.off()


# Mix all processes ####
uph_frequencies = c(10,100,1000,2000,5000,10000)
aggr.layers = 5
cols = viridis(n=6)

png('c:/Users/Marijn/sciebo/Research/Bioturbation/Bioturbation_code/Figures/combinationsWithUpheaval.png',
    height=15,width=30,res=500,units='cm')
par(mfrow = c(1,2))
mains = c('Mounding and upheaval', 'Mixing and upheaval')
wds = c('moundingUpheaval', 'mixingUpheaval')
labs=c('A','B')
for(ds in wds)
{
  setwd(paste0('C:/Simulations/Bioturbation/0_Final/',ds,'/'))
  plot(0,pch=NA,xlim=c(0,10000),ylim=c(1.5,0),main='',
       xlab='Age [a]', ylab='Depth [m]')
  
  axis(side = 3, at = c(0,10000),labels = c(0,1))
  mtext(side = 3,'Bioturbated fraction [-]', line = 1)
  title(mains[which(wds%in%ds)], line = 2.6)
  
  for (i in c(0,1,2,3,4,5))
  {
    agedepth = matrix(nrow = 200, ncol = 4)
    soil = read.csv(paste0(i,"_10000_out_allsoils.csv"))
    ages = read.csv(paste0(i,"_10000_out_OSL_ages.csv"))
    layers=unique(soil$nlayer)
    for(l in 0:(length(layers)/aggr.layers-1))
    {
      lay.sel = seq(l*aggr.layers, ((l+1)*aggr.layers-1))
      soil.sel = soil[soil$nlayer%in%lay.sel,]
      depth = soil.sel$cumth_m[nrow(soil.sel)] - sum(soil.sel$thick_m)/2
      age = ages[ages$layer%in%lay.sel,]
      
      fbio = 1-length(which(age$grain_age==10000))/nrow(age)
      count = nrow(age)
      age = age[age$grain_age<10000,]
      if(nrow(age)>1 & length(unique(age$grain_age))>3)
      {
        dens = density(age$grain_age,bw='SJ')
        mode = dens$x[dens$y==max(dens$y)][1]
        dens$y = dens$y/max(dens$y)
        dens$y[dens$y<0.05] = NA
        # dens$y = (dens$y-0.05)/.95
        dens$y = -dens$y*0.05 + depth
        IQR = quantile(age$grain_age,.75)-quantile(age$grain_age,.25)
        # points(rep(depth, nrow(age))~age$grain_age, pch=16, col=adjustcolor('gray34', alpha.f = 0.025))
        # lines(dens$y ~dens$x, col = i + 1,lty=2)
      }
      else
      {
        mode = NA
        IQR=NA
      }
      agedepth[l+1,] = c(depth,mode,fbio,IQR)
    }
    
    lines(y = agedepth[,1], 
          x = (agedepth[,2]),lwd = 3, lty = 1,col=cols[i+1]) # Modes
    lines(y=agedepth[,1], x=(agedepth[,3]*10000), lty = 1,col=cols[i+1]) # fbio
    # lines(y=agedepth[,1], x=(agedepth[,4]), lty = 2,col=cols[i+1]) # IQR
  }
  
  if(ds==wds[1])
  {
    legend('bottomright', 
           legend = c(paste(uph_frequencies[1:5],'a'),'none','Modes','fbio'),
           col=c(cols,1,1),lty=c(rep(NA,6),1,1),pch=c(rep(15,6),NA,NA),lwd=c(rep(1,6),3,1), bty='n',ncol=3)
    text('Upheaval frequency',x=5500,y=1.22)
  }
  legend('topright', labs[which(wds%in%ds)],bty='n')
}
dev.off()




# Diffrent rates ####
rates = c(1,2.5,5,7.5,10)
aggr.layers = 5
cols = viridis(n=length(rates))

png('c:/Users/Marijn/sciebo/Research/Bioturbation/Bioturbation_code/Figures/moundingDiffRates.png',
    height=15,width=15,res=500,units='cm')
{
setwd('C:/Simulations/Bioturbation/0_Final/mounding_diffRates/')

plot(0,pch=NA,xlim=c(0,10000),ylim=c(1.5,0),main='',
     xlab='Age [a]', ylab='Depth [m]')
axis(side = 3, at = c(0,10000),labels = c(0,1))
mtext(side = 3,'Bioturbated fraction [-]', line = 1)
title('Mounding with different rates', line = 2.6)
cols=viridis(length(rates))
for(i in (1:length(rates))-1)
{
  agedepth = matrix(nrow = 200, ncol = 4)
  soil = read.csv(paste0(i,"_10000_out_allsoils.csv"))
  ages = read.csv(paste0(i,"_10000_out_OSL_ages.csv"))
  layers=unique(soil$nlayer)
  for(l in 0:(length(layers)/aggr.layers-1))
  {
    lay.sel = seq(l*aggr.layers, ((l+1)*aggr.layers-1))
    soil.sel = soil[soil$nlayer%in%lay.sel,]
    depth = soil.sel$cumth_m[nrow(soil.sel)] - sum(soil.sel$thick_m)/2
    age = ages[ages$layer%in%lay.sel,]
    
    fbio = 1-length(which(age$grain_age==10000))/nrow(age)
    count = nrow(age)
    age = age[age$grain_age<10000,]
    if(nrow(age)>1 & length(unique(age$grain_age))>3)
    {
      dens = density(age$grain_age,bw='SJ')
      mode = dens$x[dens$y==max(dens$y)][1]
      dens$y = dens$y/max(dens$y)
      dens$y[dens$y<0.05] = NA
      # dens$y = (dens$y-0.05)/.95
      dens$y = -dens$y*0.05 + depth
      IQR = quantile(age$grain_age,.75)-quantile(age$grain_age,.25)
      # points(rep(depth, nrow(age))~age$grain_age, pch=16, col=adjustcolor('gray34', alpha.f = 0.025))
      # lines(dens$y ~dens$x, col = i + 1,lty=2)
    }
    else
    {
      mode = NA
      IQR=NA
    }
    agedepth[l+1,] = c(depth,mode,fbio,IQR)
  }
  
  lines(y = agedepth[,1], 
        x = (agedepth[,2]),lwd = 3, lty = 1,col=cols[i+1]) # Modes
  lines(y=agedepth[,1], x=(agedepth[,3]*10000), lty = 1,col=cols[i+1]) # fbio
  # lines(y=agedepth[,1], x=(agedepth[,4]), lty = 2,col=cols[i+1]) # IQR
}
text('Bioturbation rate',y=1.2, x=4000)
legend('bottomright', 
       legend = c(legend = paste(rates,'kg m-2'),NA,'Modes','Bioturbated',' fraction'),
       col=c(cols,NA,1,1),
       lty=c(rep(NA,length(rates)),NA,1,1),
       pch=c(rep(15,length(rates)),NA,NA,NA,NA),
       lwd=c(rep(1,length(rates)),NA,3,1), bty='n',ncol=3)

};dev.off()



