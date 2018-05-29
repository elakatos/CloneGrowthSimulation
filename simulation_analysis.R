library('ggplot2'); library('reshape2'); library(ggpubr)
setwd('~/CRCdata/Simulation_results/')
colList = c('darksalmon', 'burlywood3', 'skyblue3', 'darkseagreen3', 'mediumorchid3')


#dirList = c('18_04_26_3', '18_04_26_4', '18_04_27_3', '18_04_27_4')
dirList = c('18_05_10_1', '18_05_10_2', '18_05_10_3', '18_05_10_4')
#dirList = c('18_05_15_1', '18_05_15_2', '18_05_15_3', '18_05_15_4')

N <- 200

p <- ggplot()
p2 <- ggplot() +scale_x_continuous(limits=c(2, 2000)) + scale_y_continuous(limits=c(0, 1))

rsqDF <- data.frame(matrix(vector(), nrow=N))
mutrDF <- data.frame(matrix(vector(), nrow=N))
rsqEpDF <- data.frame(matrix(vector(), nrow=N))
cisDF <- data.frame(matrix(vector(), nrow=N))
vafIndTotal <- data.frame(matrix(vector(), ncol=4))
names(vafIndTotal) <- c('invf', 'variable', 'value', 'dir')

i=1

for (dir in dirList){

vafdata <- read.table(paste0(dir,'/Vafdf.csv'), sep=',', header=T, row.names=1)

#vafsample <- sample(names(vafdata))

rsq <- scan(paste0(dir,'/Rsqs.txt'))
mutrs <- scan(paste0(dir,'/MutRatios.txt'))
rsqEp <- scan(paste0(dir,'/Rsqs_ep.txt'))

vafSummary <- data.frame(invf = as.numeric(row.names(vafdata)), avg = apply(vafdata, 1, mean), sdev = apply(vafdata, 1, sd))
vafdata$invf <- as.numeric(row.names(vafdata))
vafInd <- melt(vafdata, id='invf')

vafInd$dir <- dir
vafIndTotal <- rbind(vafIndTotal, vafInd)

rsqDF[, dir] <- rsq
mutrDF[, dir] <- mutrs
rsqEpDF[, dir] <- rsqEp

cellImm <- scan(file=paste0(dir,'/run_julia_simulations_batch.sh.o'), what=character(), sep='\n')  
cellImm <- cellImm[seq(1, length(cellImm), by=2)]
cellImmScore <- sapply(cellImm, function(x) (1e5-as.numeric(unlist(strsplit(x, ' '))[3]))/1e5 )

cisDF[,dir] <- cellImmScore

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
rsqEpM <- melt(rsqEpDF); rsqEpM <- subset(rsqEpM, value > 0)
cisM <- melt(cisDF)

p1 <- ggplot(vafIndTotal, aes(x=invf, y=value, group=variable, colour=dir)) +
  geom_line(alpha=0.1) + facet_grid(~dir) + scale_colour_manual(values=colList) +
  theme_bw() + labs(y = 'cumulative VAF')

p3 <- ggplot(rsqM, aes(x=variable, y=value, fill=variable)) + geom_violin() +
  coord_flip() + theme_bw() + scale_fill_manual(values=colList) +
  labs(y = 'R^2 value', x='')

p4 <- ggplot(rsqEpM, aes(x=variable, y=value, fill=variable)) + geom_violin() +
  coord_flip() + theme_bw() + scale_fill_manual(values=colList)

p5 <- ggplot(mutrM, aes(x=variable, y=value, fill=variable)) +
  geom_violin() + theme_bw() + scale_fill_manual(values=colList)
p6 <- ggplot(cisM, aes(x=value, fill=variable)) + facet_grid(~variable) +
  geom_histogram(alpha=0.8) + theme_bw() + scale_fill_manual(values=colList) +
  labs(x = '% immunogenic cells')

pdf('Neutral_negative_comparison_psmall.pdf', width=10, height=5)
p1
p6
p3
dev.off()



# Biomodality in epitope-containing samples


dir= '18_05_10_3'


  vafdata <- read.table(paste0(dir,'/Vafdf.csv'), sep=',', header=T, row.names=1)
  vafdata_ep <- read.table(paste0(dir,'/Vafdf_ep.csv'), sep=',', header=T, row.names=1)
  names(vafdata_ep) <- sapply(names(vafdata_ep), function(x) paste0(x, '_ep'))
  
  vafdata_ratio <- vafdata_ep/vafdata
  
  #special_ind <- which(vafdata[nrow(vafdata),]>3500)
  
  vafdata$invf <- as.numeric(row.names(vafdata))
  vafdata_ratio$invf <- as.numeric(row.names(vafdata_ratio))
  
  vafdata <- cbind(vafdata, vafdata_ep)
  vafInd <- melt(vafdata, id='invf')
  vafInd$type <- endsWith(as.character(vafInd$variable), "_ep")
  
  vafInd2 <- vafInd[vafInd$type,]
  vafInd1 <- vafInd[!vafInd$type,]
  
  vafRatioInd <- melt(vafdata_ratio, id='invf')
  
  pInd <- ggplot() +
    geom_line(data=vafInd1, aes(x=invf, y=value, group=variable, colour=type),alpha=0.1) +
    theme_bw() #+ geom_line(data=vafSummary, aes(x=invf, y=avg),size=1.2)
  
  
  pInd2 <- ggplot() +
    geom_line(data=vafInd2, aes(x=invf, y=value, group=variable, colour=type),alpha=0.1) +
    theme_bw() #+ geom_line(data=vafSummary, aes(x=invf, y=avg),size=1.2)
  pRat <- pRat + 
    geom_line(data=vafRatioInd, aes(x=invf, y=value, group=variable),color='red',alpha=0.1) +
    theme_bw() +  scale_x_continuous(limits=c(1, 10))
  



#Read in immunotherapy-sims

# cellImm <- scan(file='18_05_18_1/run_julia_simulations.sh.o', what=character(), sep='\n')  
# cellImm <- cellImm[seq(1, length(cellImm), by=2)]
# cellImmScore <- sapply(cellImm, function(x) (1e5-as.numeric(unlist(strsplit(x, ' '))[3]))/1e5 )

dirList <- c('18_05_18_5', '18_05_18_6', 
             '18_05_18_7', '18_05_18_8', '18_05_21_1', '18_05_21_2',
             '18_05_21_3', '18_05_21_4', '18_05_21_5', '18_05_21_6')

for (dir in dirList){

pIT <- ggplot() + scale_color_gradientn(colours=c('darkblue','skyblue4', 'grey80','darksalmon', 'darkred'), values=c(0, 0.2, 0.4, 0.7, 1)) + theme_bw()

for (i in 1:200){
Npost <- read.table(paste0(dir,'/postIT_sparse_',i,'.txt'), header=T, sep=',')
Npost$immScore <- 1-Npost$nonImm/Npost$N
pIT <- pIT + geom_line(data=Npost, aes(x=t, y=N, colour=immScore), alpha=0.5)
}

#pdf(paste0('PIT_N_', dir,'.pdf'), width=7, height=5)
print(pIT)
#dev.off()
}


