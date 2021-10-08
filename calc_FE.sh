#!/bin/bash

#cat << EOF > bw.R
##!/bin/Rscript
#library(stats)
#C <- read.table('BIAS')
#print(bw.nrd0(C\$V2))
#EOF
#
#bw=$(Rscript --vanilla bw.R | awk '{print $2}')
bw=0.02
echo $bw

cat << EOF > opes.dat
# Read COLVAR file
di:         READ FILE=BIAS2 IGNORE_TIME VALUES=di
opes:      READ FILE=BIAS2 IGNORE_TIME VALUES=opes.bias
uwall:      READ FILE=BIAS2 IGNORE_TIME VALUES=uwall.bias
w1:         READ FILE=BIAS2 IGNORE_TIME VALUES=LW.bias

# Define weights
weights: REWEIGHT_BIAS TEMP=313 ARG=opes.bias,uwall.bias,w1.bias #w2.bias,w3.bias,w4.bias #,w5.bias,w6.bias,w7.bias,w8.bias,w9.bias

HISTOGRAM ...
  ARG=di
  GRID_MIN=0.1
  GRID_MAX=1.2
  GRID_BIN=500
  BANDWIDTH=$bw
  LOGWEIGHTS=weights
  LABEL=hh
... HISTOGRAM

DUMPGRID GRID=hh  FILE=histo FMT=%24.16e
EOF

sed '2,10000d' Colvar.data > BIAS2
cat opes.dat | plumed driver --noatoms --plumed /dev/stdin --kt 2.603

cat << 'EOF' > plot.R
#!/bin/Rscript
pdf('FE.pdf')
H <- read.table('histo')
F <- -log(H$V2)+2*log(H$V1)
H2 <- read.table('/media/Data/work/gromacs/LiBNMB/MTD_final_all/LiTFSI_S/histo_wall2')
F2 <- -log(H2$V2)+2*log(H2$V1)
plot(H$V1,F-min(F),ylim=c(0,15),xlim=c(0,0.6))
lines(H2$V1,F2-min(F2),col="red")
dev.off()
EOF

cat << 'EOF' > plot.R
#!/bin/Rscript
pdf('FE.pdf')
H <- read.table('histo')
F <- -log(H$V2)+2*log(H$V1)
L <- which(H$V1>0.4)[1]
F_L <- F-min(F[1:L])
Min <- which(F_L[1:L]==min(F_L[1:L]))
Max <- Min + which(F_L[Min:L]==max(F_L[Min:L]))
Bar <- F[Max] - F[Min]
H2 <- read.table('/media/Data/work/gromacs/LiBNMB/MTD_final_all/LiTFSI_S/histo_wall2')
F2 <- -log(H2$V2)+2*log(H2$V1)

par(mar=c(5,5,2,2))
#plot(H$V1,F-min(F[1:L]),ylim=c(-5*Bar,3*Bar),xlim=c(0.15,0.65),type='l',col="red",ylab="FE (kT)",xlab="r (nm)",cex.lab=1.5,lwd=2,cex.axis=1.2,main="LiTFSI (EC)")
plot(H$V1,F-min(F[1:L]),ylim=c(-2,10),xlim=c(0.15,0.65),type='l',col="red",ylab="FE (kT)",xlab="r (nm)",cex.lab=1.5,lwd=2,cex.axis=1.2,main="LiTFSI (EC)")
lines(H2$V1,F2-min(F2[1:L]),lty=2,lwd=2)
text(0.65,1.5*Bar,"Upper \n Wall",pos=2,cex=2,col="blue")
abline(v=0.65,lwd=3,col="blue")
abline(h=0,lty=3)
legend("topright",c('Bulk','Cluster'),col=c('black','red'),bg="antiquewhite",lty=c(2,1),cex=1.8,lwd=c(3,3))
library(shape)
shape::Arrows(x0=H$V1[Min],y0=0,x1=H$V1[Min],y1=F_L[Max],code=2,arr.adj=1,lwd=2)
text(H$V1[Min],Bar,paste0(round(Bar,2)," kT"),pos=3,cex=1.5,col="green3")
write.table(x = Bar,row.names = FALSE,col.names = FALSE,file = 'barrier')
dev.off()
EOF

Rscript --vanilla plot.R

