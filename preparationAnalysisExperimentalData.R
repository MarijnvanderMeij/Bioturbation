library(openxlsx)
library(Luminescence)

# PDF of a sample
PDF_mean_sd = function(means, sds,normalized = T)
{
  range = round(range(c(means-2*sds, means + 2*sds)))
  x = seq(range[1],range[2], by = 1)
  densities = matrix(nrow = length(means), ncol = length(x))
  for (i in 1:length(means))
  {
    mean = means[i]
    sd = sds[i]
    y = 1/(sd*sqrt(2*pi))*exp(-0.5*((x-mean)/sd)^2)
    densities[i,] = y
  }
  densities = colSums(densities)
  if(normalized) {  densities = densities/max(densities)}
  else{densities = densities/sum(densities)}
  return(data.frame(x = x, y = densities))
}

DR_properties = c('DR','U','Th','K','DR')

# Prepare datasets ####
# Tree throw
setwd('C:/Users/Marijn/sciebo/Research/Bioturbation/Data/Trees_Czech/')
trees_info = read.csv('DR_trees.csv')

De_trees = data.frame()
for(prof in c(1:4))
{
  files = list.files(paste0('Profile_',prof),full.names = T)
  for(depth in c(15,30,50,100,140))
  {
    sheetnames =  getSheetNames(files[grepl(paste0('-',depth,'_'),files)])
    for(method in c('pIRIRSL140', 'IRSL50'))
    {
      sn = sheetnames[grepl(method,sheetnames) & grepl('forR',sheetnames) & !grepl('DRT',sheetnames)]
      if(length(sn)>1)
      {
        sn = sn[grepl('_mod',sn)]
      }
      De = read.xlsx(xlsxFile = files[grepl(paste0('-',depth,'_'),files)],
                     sheet = sn)
      De$depth = depth/100
      De$profile = prof
      De$method = method
      DR.sub = trees_info[trees_info$Profile==prof & trees_info$Depth==depth/100,]
      De$NSF = DR.sub$rel_n_sat_grains
      for(par in DR_properties)
      {
        De[[par]] = DR.sub[[par]]
        De[[paste0(par,'_err')]] = DR.sub[[paste0(par,'_err')]]
      }
      De_trees = rbind(De_trees, De)
    }
  }
}
De_trees = De_trees[De_trees$profile==2 & De_trees$method=='pIRIRSL140',] # Representative profile

# Worms chernozem
setwd('C:/Users/Marijn/sciebo/Research/Bioturbation/Data/Worms_Chernozem/')
worms_info = read.csv('DR_Chernozem.csv')
files = list.files(getwd(),full.names = T)
files = files[grepl('New_DE', files) & !grepl('~', files)]

De_worms = data.frame()
for(f in files)
{
  De = read.xlsx(f)
  ID = as.numeric(substr(f,start = 80,nchar(f)-5))
  De = De[,c('Gray','Gray_Err')]
  names(De) = c('De', 'err')
  De$depth = worms_info[worms_info$ID==ID,'depth']
  De$profile = worms_info[worms_info$ID==ID,'profile']
  De$NSF = worms_info[worms_info$ID==ID,'rel_n_sat_grains']
  De$ID = ID
  De$method = 'pIRIR150'
  for(par in DR_properties)
  {
    De[[par]] = worms_info[worms_info$ID==ID,par]
    De[[paste0(par,'_err')]] = worms_info[worms_info$ID==ID,paste0(par,'_err')]
  }
  
  De_worms = rbind(De_worms,De)
}
De_worms = De_worms[De_worms$profile==2,] # representative profile

# Ants Spain
setwd('C:/Users/Marijn/sciebo/Research/Bioturbation/Data/Ants_Spain/')
ants_info = read.csv('DR_ants.csv')
De_ants = read.csv('De_ants_Spain.csv')

for(ID in unique(De_ants$Sample))
{
  if(ID%in%ants_info$ID)
  {
    De_ants[De_ants$Sample==ID,'profile'] = ants_info[which(ants_info$ID==ID),'profile']
    De_ants[De_ants$Sample==ID,'depth'] = ants_info[which(ants_info$ID==ID),'depth']
    De_ants[De_ants$Sample==ID,'NSF'] = ants_info[which(ants_info$ID==ID),'rel_n_sat_grains']
    for(par in DR_properties)
    {
      De_ants[De_ants$Sample==ID,par] = ants_info[which(ants_info$ID==ID),par]
      De_ants[De_ants$Sample==ID,paste0(par,'_err')] = ants_info[which(ants_info$ID==ID),paste0(par,'_err')]
    }
  }
}
De_ants = De_ants[!is.na(De_ants$profile),]
De_ants = De_ants[De_ants$profile=='SC-10',] # representative profile

# Termites Ghana
setwd('C:/Users/Marijn/sciebo/Research/Bioturbation/Data/Termites_Ghana/')
termites_info = read.csv('DR_termites.csv')
De_termites = read.csv('De_termites_Ghana.csv')

for(ID in unique(De_termites$Sample))
{
  if(ID%in%termites_info$ID)
  {
    De_termites[De_termites$Sample==ID,'profile'] = termites_info[which(termites_info$ID==ID),'profile']
    De_termites[De_termites$Sample==ID,'depth'] = termites_info[which(termites_info$ID==ID),'depth']
    De_termites[De_termites$Sample==ID,'NSF'] = termites_info[which(termites_info$ID==ID),'rel_n_sat_grains']
    for(par in DR_properties)
    {
      De_termites[De_termites$Sample==ID,par] = termites_info[which(termites_info$ID==ID),par]
      De_termites[De_termites$Sample==ID,paste0(par,'_err')] = termites_info[which(termites_info$ID==ID),paste0(par,'_err')]
    }
  }
}
De_termites = De_termites[!is.na(De_termites$profile),]
De_termites = De_termites[De_termites$profile==1 & De_termites$depth <1.03, ]
 
# write.csv(De_termites, "De_termites_Ghana_complete.csv")


# Plough Germany
setwd('C:/Users/Marijn/sciebo/Research/Bioturbation/Data/Plough_CarboZALF/')
De_plough = read.csv('De_plough.csv')
plough_info = read.csv('DR_plough.csv')
for(ID in unique(De_plough$Sample))
{
  if(ID%in%plough_info$ID)
  {
    De_plough[De_plough$Sample==ID,'profile'] = plough_info[which(plough_info$ID==ID),'profile']
    De_plough[De_plough$Sample==ID,'depth'] = plough_info[which(plough_info$ID==ID),'depth']
    De_plough[De_plough$Sample==ID,'NSF'] = plough_info[which(plough_info$ID==ID),'rel_n_sat_grains']
    for(par in DR_properties)
    {
      De_plough[De_plough$Sample==ID,par] = plough_info[which(plough_info$ID==ID),par]
      De_plough[De_plough$Sample==ID,paste0(par,'_err')] = plough_info[which(plough_info$ID==ID),paste0(par,'_err')]
    }
  }
}
De_plough = De_plough[!is.na(De_plough$profile),]

# Calculate distribution statistics ####
agestat = data.frame(matrix(nrow = 50, ncol = 11))
names(agestat) = c('depth','mean','median','sd','mode','rsd','IQR','skewness','fbio', 'OD', 'dataset')

count = 1
for(dataset in c('De_trees', 'De_worms', 'De_termites', 'De_ants','De_plough'))
{
  De = get(dataset)
  # De$De = De$De / quantile(De$De, 0.99)
  # De$depth = De$depth / max(De$depth)
  for(depth in unique(De$depth))
  {
    De.t = De[De$depth ==depth,]
    depth_norm = depth/max(De$depth)
    nsf = De.t$NSF[1]; if(is.null(nsf)){nsf=NA}
    
    # Normalization
    De.t$age = De.t$De/De.t$DR * 1000
    De.t$age_err = De.t$err/De.t$DR * 1000
    De.t$err = De.t$err/max(De.t$De)
    
    # Set thresholds for ages, or otherwise maximum age values
    if(dataset=='De_trees'){norm = max(De.t$age)}
    if(dataset=='De_worms'){norm = 13200}
    if(dataset=='De_termites'){norm = 4000}
    if(dataset=='De_ants'){norm = max(De.t$age)}
    if(dataset=='De_plough'){norm = 220}
    
    # Normalize ages
    De.t$age = De.t$age/norm
    De.t$age_err = De.t$age_err/norm
    
    # Add grains that are outside of the thresholds to the NSF and calculate fbio
    if(is.na(nsf)){nsf=0}
    nsf = nsf + length(which(De.t$age>1))/nrow(De.t)*(1-nsf)
    fbio=1-nsf
    
    data.temp=get(dataset)
    data.temp[data.temp$depth==depth,'fbio'] = fbio
    assign(dataset, data.temp)
    # Filter out the grains
    De.t = De.t[!(De.t$age>1),]    
    
    if(nrow(De.t)>3)
    {
      dens = density(De.t$age,bw='SJ',na.rm=T)
      mode = dens$x[dens$y==max(dens$y)][1]
      mean = mean(De.t$age)
      median = median(De.t$age)
      sd = sd(De.t$age)
      rsd = sd/mean
      IQR = quantile(De.t$age,0.75) - quantile(De.t$age,0.25)
      skewness = 3 * (mean-median)/sd # alternative Pearson Mode skewness
      dat = De.t[,c('De','err')]
      dat = dat[dat$De>0,]
      OD = calc_CentralDose(dat,sigmab = 0,plot=F,verbose=F)$summary$rel_OD
      agestat[count, ] = c(depth_norm,mean,median,sd,mode,rsd,IQR,skewness,fbio, OD, dataset)
    }
    agestat[count,'depth'] = depth_norm 
    agestat[count,'fbio'] = fbio
    agestat[count,'dataset'] = dataset
    
    count = count + 1
  }
}
agestat = agestat[!is.na(agestat$dataset),]
write.csv(agestat, 'C:/Users/Marijn/sciebo/Research/Bioturbation/Data/distributionStatistics.csv',row.names = F)




# De analysis ####
png('C:/Users/Marijn/sciebo/Research/Bioturbation/Figures/depthProfiles_allOrganisms.png',
    width = 18, height = 12, units = 'cm', res = 500)
# cairo_pdf('C:/Users/Marijn/sciebo/Research/Bioturbation/Figures/depthProfiles_allOrganisms.pdf',
    # width = 18/2.54, height = 12/2.54)
{
  bool_age = T
  labels = c('A','B','C','D','E')
  par(mfrow = c(2,3))
  ct=0
  for(dataset in c('termites', 'worms', 'ants','trees','plough'))
  {
    ct=ct+1
    De.t = get(paste0('De_',dataset))
    for(prof in unique(De.t$profile))
    {
      De = De.t[De.t$profile==prof,]
      if(bool_age)
      {
        De$age = De$De/De$DR
        xlim= range(De$age)
        
        if(dataset=='worms'){ xlim=c(0,25)}
        if(dataset=='termites'){xlim = c(0,10)}
        if(dataset=='plough'){xlim = c(0,8)}
        
        plot(0,pch=NA, ylim = c(max(De.t$depth,na.rm=T),0),
             xlim = xlim,
             xlab = 'Age [ka]', ylab = 'Depth [m]')
      } else
      {
        plot(0,pch=NA, xlim = quantile(De$De,na.rm=T,probs = c(0,.95)), 
             ylim = c(max(De.t$depth, na.rm=T), 0),
             xlab = 'De [Gy]', ylab = 'Depth [m]')
      }
      axis(side=3,at = c(0,xlim[2]*c(.5,1)),labels = c(0,0.5,1))
      mtext(side = 3, line = 2.6, text = 'Bioturbated fraction [-]',cex=0.66)
      
      y_corr = 0.8* min(diff(sort(unique(De$depth))))
      
      for(m in unique(De$method))
      {
        col = which(unique(De$method)%in%m)
        modes = data.frame()
        for (d in sort(unique(De$depth)))
        {
          fbio=De[De$depth==d,'fbio'][1]
          sub = De[De$depth==d & De$profile==prof &De$method==m,]
          sub_De = sub$De
          if(bool_age)
          {
            sub_De = sub_De/sub$DR[1]
          }
          dens = density(sub_De, bw = 'SJ')
          dens$y=dens$y/max(dens$y)
          dens$y[dens$y<0.05]=NA
          
          y=-dens$y*y_corr +d
          x = dens$x
          
          # abline(h=d, lty = 2, col = 'grey')
          lines(y~x, col = col)
          points(rep(d,length(sub_De))~sub_De,pch = 16,col = adjustcolor('black',alpha.f = 3/length(sub_De)))
          
          modes = rbind(modes,c(d, x[which(y==min(y,na.rm=T))],fbio))
        }
        points(modes[,1]~modes[,2], col = col, pch = 16)
        lines(modes[,1]~modes[,2], col = col, lwd = 2)
        lines(y=modes[,1],x=(modes[,3])*max(xlim), col = col, lwd = 1,lty=2)
        assign(paste0("modes_",dataset),modes)
        
        
        
      }
      if(dataset=='worms'){abline(v=13.4,col = 'red',lty=2)}
      if(dataset=='termites'){abline(v=4,col = 'red',lty=2)}
      if(dataset=='plough'){abline(v=.220,col = 'red',lty=2)}  
      # legend('topright', legend = paste0(labels[ct],": ",dataset),box.lty=0,bg='white')
      legend('topright', legend = paste0(labels[ct],": ",dataset),bty='n')
      box(which='plot')
    }
  }
  plot(0,pch=NA,xlab='',ylab='',axes=F)
  legend('center', bty = 'n', lty = c(1,NA,NA,2,2), pch = c(NA,16,16,NA,NA), 
         legend = c('Age distribution', 'Mode', 'Measurement', 'Bioturbated fraction','Bioturbation period'),
         col = c(1,1,adjustcolor('black',alpha.f=0.1),1,'red'))

};dev.off()

# entire profile
par(mfrow = c(2,3))
# for(dataset in c('De_trees', 'De_worms'))
for(dataset in c('De_trees', 'De_worms', 'De_ants', 'De_termites', 'De_plough'))
{
  De = get(dataset)
  dens = KDE(De$De/De$DR)
  # dens = PDF_mean_sd(means = De$De, De$err)
  dens$y=dens$y/max(dens$y)
  dens$x[dens$y<0.02] = NA
  plot(dens$y~dens$x, type = 'l', xlab = 'De [Gy]', ylab = 'density [-]', main = dataset)
  
}


# DR analysis ####
# for(dataset in c('De_trees', 'De_worms'))
for(dataset in c('De_trees', 'De_worms', 'De_termites', 'De_plough'))
{
  De = get(dataset)
  par(mfrow = c(2,2))
  for(par in c('U','Th','K','DR'))
  {
    plot(0,pch=NA,
         xlim = range(c(De[[par]]-De[[paste0(par,'_err')]],
                        De[[par]]+De[[paste0(par,'_err')]]), na.rm=T), 
         ylim = rev(range(De$depth)),
         ylab = 'Depth', 
         xlab = paste(par,'[?]'))
    if(par=='U'){  mtext(text = dataset, side = 3, outer = T, line = -2)}
    for(p in unique(De$prof))
    {
      DR.sub = De[De$prof==p,c('depth',par,paste0(par,'_err')),]
      DR.sub = DR.sub[!duplicated(DR.sub),]
      DR.sub = DR.sub[order(DR.sub$depth),]
      lines(DR.sub$depth~DR.sub[[par]], col = 1)
      points(DR.sub$depth~DR.sub[[par]], col = 1, pch = 16)
      arrows(x0=DR.sub[[par]]-DR.sub[[paste0(par,'_err')]],
             x1=DR.sub[[par]]+DR.sub[[paste0(par,'_err')]],
             y0=DR.sub$depth,
             y1=DR.sub$depth,
             col = 1, code = 3,length=0.1,
             angle = 90)
    }
    abline(v = mean(De[[par]]), lty=2, col = 'grey')
    legend('topright', lty=1, col = 1,
           legend = paste('Profile ',unique(De$prof)))
  }
  
}




