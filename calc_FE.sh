#!/bin/bash
eval "$(conda shell.bash hook)"
conda activate CLIPS
#set -e
REF='./reference'
Ion1=$1
Ion2=$2
Solv=$3
nstepsmtd=$4
i=$5
Time=$(printf %.1f $(echo "0.000001*$nstepsmtd*$i/10" | bc -l))
HISTO="$REF/$Ion1$Ion2/histo_wall2"
rm -rf bck.*
#cat << EOF > bw.R
##!/bin/Rscript
#library(stats)
#C <- read.table('BIAS2')
#print(bw.nrd0(C\$V2))
#EOF

#sed '2,200d' Colvar.data > BIAS2
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
  GRID_MIN=0.001
  GRID_MAX=2
  GRID_BIN=500
  BANDWIDTH=$bw
  LOGWEIGHTS=weights
  LABEL=hh
... HISTOGRAM

DUMPGRID GRID=hh  FILE=histo FMT=%24.16e
EOF

L=$(sed "2,500d" Colvar.data | wc -l)
echo $L
NL=$(echo "$i*$L/10" | bc)
sed "2,500d" Colvar.data | head -n $NL | sed -n "1~10 p" > BIAS2
echo $(wc -l BIAS2)
cat opes.dat | plumed driver --noatoms --plumed /dev/stdin --kt 2.603

cat << 'EOF' > plot.R
#!/bin/Rscript

# Find Minima/Maxima
####################
local.min.max <- function(x, dev=mean, plot=TRUE, add.points=FALSE,  ...) {
  x <- stats::na.omit(x)
  r <- rle(x) 
  minima <- which(rep(x=diff(sign(diff(c(-Inf, r$values, -Inf)))) == 2, times=r$lengths))
  maxima <- which(rep(x=diff(sign(diff(c(-Inf, r$values, -Inf)))) == -2, times=r$lengths)) 
  if (plot == TRUE) {				 
    plot(x,type="l", ...)
    graphics::points(x[minima]~minima,pch=19,col="blue") 
    graphics::points(x[maxima]~maxima,pch=19,col="red")
    graphics::abline(h=dev(x, na.rm=TRUE), col="grey")
    if (add.points == TRUE) graphics::points(x, col="grey")
    graphics::legend("topleft", legend=c("Minima","Maxima"), pch=c(19,19), 
                     col=c("blue","red"), bg="white")
  }
  return( list(minima=x[minima], maxima=x[maxima],
               devmin=abs(dev(x) - x[minima]), 
               devmax=abs(dev(x) - x[maxima])) )
}
####################

pdf('FE.pdf')
H <- read.table('histo')
FE <- -log(H$V2)+2*log(H$V1)

if (file.exists('HISTO')){ 
H2 <- read.table('HISTO') 
FE2 <- -log(H2$V2)+2*log(H2$V1)
}

# Trim 
x <- H$V1[is.finite(FE)]
F <- FE[is.finite(FE)]

if (file.exists('HISTO')){ 
x2 <- H2$V1[is.finite(FE2)]
F2 <- FE2[is.finite(FE2)]
}

# Calc barrier
LMM <- local.min.max(F,plot = FALSE)
Min <- which(F==LMM$minima[1])
Max <- which(F==LMM$maxima[2])

if (file.exists('HISTO')){ 
LMM2 <- local.min.max(F2,plot = FALSE)
}

if(x[Max] < 0.65){ 
Bar <- LMM$maxima[2] - LMM$minima[1]
} else {
L <- which(x>0.65)[1]
Bar <- F[L] - LMM$minima[1]
}

# calc binding energy
BE <- NA
L <- which(x>0.65)[1]
if(length(LMM$minima) > 1){ 
  BE <- LMM$minima[2] - LMM$minima[1]
  Min2 <- which(F==LMM$minima[2])
} else{ 
    BE <- F[L] - LMM$minima[1]
    Min2 <- L
}

# Plot
library(latex2exp)
par(mar=c(5,5,2,2))
plot(x,F-LMM$minima[1],ylim=c(-10,20),xlim=c(0.15,0.75),type='l',col="red",ylab=TeX('$\\Delta G$ $(kT)$'),xlab="r (nm)",cex.main=2,cex.lab=2,lwd=2,cex.axis=1.5,main=TeX('PMF (Ion1$^+$ - Ion2$^-$ / Solv)'))

#plot(x,F-LMM$minima[1],ylim=c(-10,20),xlim=c(0.15,0.75),type='l',col="red",ylab="FE (kT)",xlab="r (nm)",cex.lab=1.5,lwd=2,cex.axis=1.2,main="PMF (Ion1 - Ion2 in Solv)")

if (file.exists('HISTO')){ 
lines(x2,F2 - LMM2$minima[1],lty=2,lwd=2)
}

text(0.65,15,"Upper\nWall",pos=2,cex=2,col="blue")
text(0.65,-2.5,"Time (ns)",pos=4,cex=1.5,col="red")
abline(v=0.65,lwd=3,col="blue")
abline(h=0,lty=3)
legend("bottomright",c('Bulk','Cluster'),col=c('black','red'),bg="antiquewhite",lty=c(2,1),cex=1.8,lwd=c(3,3))
library(shape)
shape::Arrows(x0=x[Min],y0=0,x1=x[Min],y1=Bar,code=2,arr.adj=1,lwd=2)
text(x[Min],Bar+1,paste0(round(Bar,2)," kT"),pos=3,cex=2,col="green3")

shape::Arrows(x0=x[Min2],y0=0,x1=x[Min2],y1=BE,code=2,arr.adj=1,lwd=2)
text(x[Min2],BE+1,paste0(round(BE,2)," kT"),pos=3,cex=2,col="green3")

write.table(x = Bar,row.names = FALSE,col.names = FALSE,file = 'barrier')
write.table(x = BE,row.names = FALSE,col.names = FALSE,file = 'bindE')
dev.off()
EOF
sed -i "s|HISTO|$HISTO|g" plot.R
sed -i "s|Ion1|$Ion1|g" plot.R
sed -i "s|Ion2|$Ion2|g" plot.R
sed -i "s|Solv|$Solv|g" plot.R
sed -i "s|Time |$Time |g" plot.R
Rscript --vanilla plot.R

