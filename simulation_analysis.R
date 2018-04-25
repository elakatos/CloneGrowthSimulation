
dirList = c('18_04_24_2', '18_04_05_1')



vafdata <- read.table(paste0('~/CRCdata/Simulation_results/',dir,'/Vafdf.csv'), sep=',', header=T, row.names=1)

rsq <- scan(paste0('~/CRCdata/Simulation_results/',dir,'/Rsqs.txt'))
mutrs <- scan(paste0('~/CRCdata/Simulation_results/',dir,'/MutRatios.txt'))

vafSummary <- data.frame(invf = as.numeric(row.names(vafdata)), avg = apply(vafdata, 1, mean), sdev = apply(vafdata, 1, sd))
vafdata$invf <- as.numeric(row.names(vafdata))
vafInd <- melt(vafdata, id='invf')

p <- ggplot() +
  geom_ribbon(data=vafSummary, aes(x = invf,ymin = avg - 2*sdev, ymax = avg+ 2*sdev),fill='skyblue3', alpha=0.3) +
  geom_line(data=vafInd, aes(x=invf, y=value, group=variable),colour='grey20',alpha=0.1) +
  theme_bw() + geom_line(data=vafSummary, aes(x=invf, y=avg),size=1.2)

p2 <- ggplot(vafSummary, aes(x=invf, y=sdev/avg)) + geom_line() + theme_bw() +
  scale_x_continuous(limits=c(30, 200)) + scale_y_continuous(limits=c(0.04, 0.12))

p3 <- ggplot() + aes(x='', y=rsq) + geom_violin(fill='slategray3') + coord_flip() + theme_bw()

p4 <- ggplot() + geom_violin(aes(x='', y=mutrs),fill='firebrick3') + theme_bw()
