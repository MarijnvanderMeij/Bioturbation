library(viridis)
z = seq(0,3,0.001)
perc = 1 # percentage of mixing at depth of 1 meter
z_upp = 1 # reference for where perc has to be reached

# Gradational
dd_grd = (z_upp - sqrt(z_upp^2*(1-perc))) / z_upp^2 # find dd_grd for perc

-dd_grd/2*z_upp^2 + z_upp
perc*(-dd_grd/2*(1/dd_grd)^2+1/dd_grd)

BT_grd = 1-dd_grd*z
BT_grd[BT_grd < 0] = 0

# Exponential
dd_exp = -log(1-perc)/z_upp

dd_exp = 6

BT_exp = exp(-dd_exp*z)

# Abrupt
dd_abr = z_upp/perc
BT_all = c(1,1,0,0)
z_uni = c(0, dd_abr,dd_abr, max(z))

cols = viridis(3)
png('C:/Users/Marijn/sciebo/Research/Bioturbation/Figures/depthFunctions.png', height = 6, width = 6, units = 'in', res = 500)
{
  plot(0,pch=NA,
       xlim = c(0,1),
       ylim = c(z_upp*1.5,0),
       xlab = "Relative rate",
       ylab = "Depth [m]", xaxt = 'n')
  axis(side = 1, at = c(0,1))
  abline(h = z_upp, col = 'grey', lty = 2)
  
  lines(z~ BT_grd, col = cols[1], lty = 1, lwd = 3)
  lines(z~ BT_exp, col = cols[2], lty = 2, lwd = 3)
  lines(z_uni~ BT_all, col = cols[3], lty = 3, lwd = 3)
  
  legend('bottom', legend = c("Gradational", "Exponential", "Abrupt"), lty = c(1,2,3), 
         col = cols, ncol = 3,bty = 'n', lwd = 3)
};dev.off()




