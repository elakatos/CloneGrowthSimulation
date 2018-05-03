library('ggplot2'); library('reshape2'); library(ggpubr)

dirList = c('18_04_26_3', '18_04_26_4', '18_04_27_3', '18_04_27_4')

p <- ggplot()
p2 <- ggplot() +scale_x_continuous(limits=c(10, 200)) + scale_y_continuous(limits=c(0.04, 1))

rsqDF <- data.frame(matrix(vector(), nrow=200))
mutrDF <- data.frame(matrix(vector(), nrow=200))

colList = c('skyblue3', 'darksalmon', 'burlywood3', 'darkseagreen3', 'mediumorchid3')
i=1

for (dir in dirList){

vafdata <- read.table(paste0('~/CRCdata/Simulation_results/',dir,'/Vafdf.csv'), sep=',', header=T, row.names=1)

rsq <- scan(paste0('~/CRCdata/Simulation_results/',dir,'/Rsqs.txt'))
mutrs <- scan(paste0('~/CRCdata/Simulation_results/',dir,'/MutRatios.txt'))

vafSummary <- data.frame(invf = as.numeric(row.names(vafdata)), avg = apply(vafdata, 1, mean), sdev = apply(vafdata, 1, sd))
vafdata$invf <- as.numeric(row.names(vafdata))
vafInd <- melt(vafdata, id='invf')


rsqDF[, dir] <- rsq
mutrDF[, dir] <- mutrs
#if (startsWith(dir, '18_04_25')){mutrDF[,dir] <- mutrDF[,dir]/2}

p <- p +
  geom_ribbon(data=vafSummary, aes(x = invf,ymin = avg - 2*sdev, ymax = avg+ 2*sdev),fill=colList[i], alpha=0.3) +
  #geom_line(data=vafInd, aes(x=invf, y=value, group=variable),colour='grey20',alpha=0.1) +
  theme_bw() + geom_line(data=vafSummary, aes(x=invf, y=avg),size=1.2)

p2 <- p2 +
  geom_line(data=vafSummary, aes(x=invf, y=sdev/avg), colour=colList[i]) + theme_bw()

i = i+1

}

mutrM <- melt(mutrDF)
rsqM <- melt(rsqDF)

p3 <- ggplot(rsqM, aes(x=variable, y=value, fill=variable)) + geom_violin() +
  coord_flip() + theme_bw() + scale_fill_manual(values=colList)

p4 <- ggplot(mutrM, aes(x=variable, y=value, fill=variable)) +
  geom_violin() + theme_bw() + scale_fill_manual(values=colList)
