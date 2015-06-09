
DensityResponse <- function(Data,xlab="",ylab="",main="",alpha=0.2,
                            binwidth=NULL,histFill="white",histCol="black",
                            densityFill="#FF6666",TitleSize=15,TextSize=15,
                            XAxisSize=15,YAxisSize=15,AngleLab=30,LegendPosition="right",
                            TitleAxesSize=15,tmar=1,bmar=1,rmar=1,lmar=1){
  if (!is.vector(Data)) stop("Input data must be a numeric vector")
  Data <- data.frame(Values=Data)
  if (is.null(binwidth)) binwidth <- abs(range(Data)[1] - range(Data)[2]) / 20
  p <- ggplot(Data, aes(x=Values)) + theme_bw() + 
    geom_histogram(aes(y=..density..),binwidth=binwidth,colour=histCol, fill=histFill) + 
    geom_density(alpha=alpha, fill=densityFill)+ ylab(ylab) + xlab(xlab) + ggtitle(main) +
    theme(text = element_text(size=TextSize),axis.text.x = element_text(size=XAxisSize,angle = AngleLab, hjust = 1),
          axis.title.x=element_text(size=TitleAxesSize),axis.title.y=element_text(size=TitleAxesSize),axis.text.y=element_text(size=YAxisSize),
          legend.position=LegendPosition,plot.title=element_text(size=TitleSize),
          legend.key=element_blank(), plot.margin=unit(c(tmar,rmar,bmar,lmar), "cm"))
return(p)
}

#################################################################################
## PCA analysis 

PCA <- function (Data, RowNames = NULL,cor=TRUE, scale = TRUE, center = TRUE,...) {
isnot.null <- function(x) ! is.null(x)
isnot.vector <- function(x) ! is.vector(x)
  if (is.matrix(Data) || is.data.frame(Data)) {
      ana <- prcomp(t(Data), cor = cor, scale = scale, center = center,...)
      PC1 <- as.vector(ana$rotation[, 1])
      PC2 <- as.vector(ana$rotation[, 2])
      if (is.null(RowNames)) {
        Data <- data.frame(PC1, PC2)
      out <- list(Data=Data,PCs_all=ana$rot,Std=ana$sdev,Info=summary(ana))
      }
      else {
        Data <- data.frame(PC1, PC2, RowNames)
      out <- list(Data=Data,PCs_all=ana$rot,Std=ana$sdev,Info=summary(ana))
      }
      return(out)
    }
  else {
      stop("Input data must be a numeric matrix or data.frame")
    }
}




##############
PairwiseDist <- function(Data,method="jaccard",...){
  if (is.matrix(Data) || is.data.frame(Data)){
    Data <- unique(Data)
    methods <- c("manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard",
                 "gower", "altGower", "morisita", "horn", "mountford", "raup" , "binomial",
                 "chao", "cao")
    method <- match.arg(method,methods)
    pwdist <- vegdist(Data, method = method,...)
    pwdist <- data.frame(as.vector(pwdist))
    names(pwdist) <- "Distance"
    return(pwdist)
  } else {
    stop("Input data must be a numeric matrix or data.frame")
  }
}

##############
PairwiseDistPlot <- function(Data,xlab="",ylab="",main="",TextSize=15,TitleSize=15,XAxisSize=15,YAxisSize=15,
                             TitleAxesSize=15,tmar=1,bmar=1,rmar=1,lmar=1,AngleLab=30,
                             binwidth=NULL,fillCol="white",Colour="black",DensityFill="#FF6666",DensityAlpha=.2){
  if (is.matrix(Data) || is.data.frame(Data)){
    if (is.null(binwidth)) binwidth <- abs(range(Data)[1] - range(Data)[2]) / 20
    p <- ggplot(Data, aes(x=Distance)) + theme_bw() +
      geom_histogram(aes(y=..density..),binwidth=binwidth,colour=Colour, fill=fillCol) + 
      geom_density(alpha=DensityAlpha, fill=DensityFill)+ ylab(ylab) + xlab(xlab) + ggtitle(main) +
      theme(text = element_text(size=TextSize),axis.text.x = element_text(size=XAxisSize,angle = AngleLab, hjust = 1),
            axis.title.x=element_text(size=TitleAxesSize),axis.title.y=element_text(size=TitleAxesSize),
            axis.text.y=element_text(size=YAxisSize),plot.title=element_text(size=TitleSize),
            legend.key=element_blank(), plot.margin=unit(c(tmar,rmar,bmar,lmar),"cm")) 
    return(p)
  } else {
    stop("Input data must be a numeric matrix or data.frame")
  }
}

##############
##############
## Maximum Model Performance
#
#MaxPerf <- function(meanNoise=0,sdNoise,meanResp,sdResp,lenPred,iters=1000,
#                    filename=NULL,pdfW=10,pdfH=10,TextSize=15,TitleSize=15,
#                    XAxisSize=15,YAxisSize=15,TitleAxesSize=15,tmar=1,bmar=1,
#                    rmar=1,lmar=1,AngleLab=30,LegendPosition="right"){
#  vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
#  R2 <- c()
#  R02 <- c()
#  Q2 <- c()
#  rmsep <- c()
#  for (i in 1:iters){
#    x <- rnorm(lenPred,mean=meanResp,sd=sdResp)
#    noise <- rnorm(length(x),mean=meanNoise,sd=sdNoise)
#    y <- x+noise
#    R2[i] <- Rsquared(y,x)
#    Q2[i] <- Qsquared2(y,x)
#    R02[i] <- Rsquared0(y,x)
#    rmsep[i] <- RMSE(y,x)
#  }
#  
#  pp <- data.frame(R2=R2)
#  title <- expression(paste("R"^{2}))
#  p1 <- ggplot(pp, aes(x=R2)) + 
#    theme_bw() + 
#    geom_density(alpha=.2, fill="#FF6666") + ggtitle(title)+ aes(y = ..count..)+ 
#    theme(text = element_text(size=TextSize),axis.text.x = element_text(size=XAxisSize,angle = AngleLab, hjust = 1),
#          axis.title.x=element_text(size=TitleAxesSize),axis.title.y=element_text(size=TitleAxesSize),
#          axis.text.y=element_text(size=YAxisSize),legend.position=LegendPosition,plot.title=element_text(size=TitleSize),
#          legend.key=element_blank(), plot.margin=unit(c(tmar,rmar,bmar,lmar),"cm")) +
#    geom_vline(aes(xintercept=mean(R2, na.rm=T)), color="blue", linetype="dashed", size=1)
#  
#  pp <- data.frame(rmsep=rmsep)
#  p2 <- ggplot(pp, aes(x=rmsep)) + 
#    theme_bw() + 
#    geom_density(alpha=.2, fill="#FF6666") + ggtitle('RMSEP') +  aes(y = ..count..)+ 
#    theme(text = element_text(size=TextSize),axis.text.x = element_text(size=XAxisSize,angle = AngleLab, hjust = 1),
#          axis.title.x=element_text(size=TitleAxesSize),axis.title.y=element_text(size=TitleAxesSize),
#          axis.text.y=element_text(size=YAxisSize),legend.position=LegendPosition,plot.title=element_text(size=TitleSize),
#          legend.key=element_blank(), plot.margin=unit(c(tmar,rmar,bmar,lmar),"cm")) +
#    geom_vline(aes(xintercept=mean(rmsep, na.rm=T)), color="blue", linetype="dashed", size=1)
#  
#  pp <- data.frame(R02=R02)
#  title <- expression(paste("R"[0]^{2}))
#  
#  p3 <- ggplot(pp, aes(x=R02)) + 
#    theme_bw() + 
#    geom_density(alpha=.2, fill="#FF6666") + ggtitle(title) + aes(y = ..count..)+ 
#    theme(text = element_text(size=TextSize),axis.text.x = element_text(size=XAxisSize,angle = AngleLab, hjust = 1),
#          axis.title.x=element_text(size=TitleAxesSize),axis.title.y=element_text(size=TitleAxesSize),
#          axis.text.y=element_text(size=YAxisSize),legend.position=LegendPosition,plot.title=element_text(size=TitleSize),
#          legend.key=element_blank(), plot.margin=unit(c(tmar,rmar,bmar,lmar),"cm")) +
#    geom_vline(aes(xintercept=mean(R02, na.rm=T)), color="blue", linetype="dashed", size=1)
#  
#  title <- expression(paste("Q"^{2}))
#  pp <- data.frame(Q2=Q2)
#  p4 <- ggplot(pp, aes(x=Q2)) + 
#    theme_bw() + 
#    geom_density(alpha=.2, fill="#FF6666") + ggtitle(title) + aes(y = ..count..)+ 
#    theme(text = element_text(size=TextSize),axis.text.x = element_text(size=XAxisSize,angle = AngleLab, hjust = 1),
#          axis.title.x=element_text(size=TitleAxesSize),axis.title.y=element_text(size=TitleAxesSize),
#          axis.text.y=element_text(size=YAxisSize),legend.position=LegendPosition,plot.title=element_text(size=TitleSize),
#          legend.key=element_blank(), plot.margin=unit(c(tmar,rmar,bmar,lmar),"cm")) +
#    geom_vline(aes(xintercept=mean(Q2, na.rm=T)), color="blue", linetype="dashed", size=1)
#  
#  if (isnot.null(filename)){
#    pdfname=paste(filename,".pdf",sep="")
#    pdf(pdfname,width=pdfW,height=pdfH)
#    grid.newpage()
#    pushViewport(viewport(layout = grid.layout(2, 2)))
#    print(p1, vp = vplayout(1,1))
#    print(p2, vp = vplayout(1,2))
#    print(p3, vp = vplayout(2,1))
#    print(p4, vp = vplayout(2,2))
#    dev.off()
#  } #else {
#    #grid.newpage()
#    #pushViewport(viewport(layout = grid.layout(2, 2)))
#    #print(p1, vp = vplayout(1,1))
#    #print(p2, vp = vplayout(1,2))
#    #print(p3, vp = vplayout(2,1))
#    #print(p4, vp = vplayout(2,2))
#  #}
#  p <- list()
#  p$p1 <- p1
#  p$p2 <- p2
#  p$p3 <- p3
#  p$p4 <- p4
#  return(p)
#}
#
##########
GetLegend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
