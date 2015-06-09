#################################################################################
## Model Assessment and Results Visualization
#################################################################################

ErrorBarplot <- function(X, Y, err, fill = X,
          main = "", ylab = "", xlab = "",
            TextSize = 15, TitleSize = 15, XAxisSize = 15, YAxisSize = 15, 
          TitleAxesSize = 15, AngleLab = 35, barcol = "red", barSize = 1, 
            barWidth = 0.3, LegendName = "Legend", ColLegend = 1, RowLegend = NULL, 
                          LegendPosition = "right", tmar = 1, bmar = 1, rmar = 1, lmar = 1, 
            stat = "identity") {
  if (length(X) != length(Y) || length(X) != length(err) || length(X) != length(fill)){
     stop("X, Y, err and fill must have the same length")
    }
  if (length(X) != length(unique(X))) {stop("The X aesthetics cannot contain repeated values (groupes")}

   data <- data.frame(X=X,Y=Y,err=err,Group=fill)
     errbar <- aes(ymin = Y- err, ymax=Y+err)

   p <- ggplot(data, aes(x=X, y = Y,fill = Group)) + theme_bw() + 
       geom_bar(position = "dodge", stat = stat) + 
       geom_errorbar(mapping = errbar,position = position_dodge(0.9), 
                     width = barWidth, color = barcol,size = barSize) +
     theme(text = element_text(size = TextSize), 
             axis.text.x = element_text(size = XAxisSize, angle = AngleLab, 
              hjust = 1), axis.title.x = element_text(size = TitleAxesSize), 
              axis.title.y = element_text(size = TitleAxesSize), 
              axis.text.y = element_text(size = YAxisSize), legend.position = LegendPosition, 
              plot.title = element_text(size = TitleSize), legend.key = element_blank(), 
              plot.margin = unit(c(tmar, rmar, bmar, lmar), "cm")) +
     guides(colour = guide_legend(LegendName, 
              ncol = ColLegend, nrow = RowLegend), shape = guide_legend(LegendName, 
              ncol = ColLegend, nrow = RowLegend)) + ylab(ylab) + xlab(xlab)
   return(p)
}

##############
#plotGrid <- function(plots,NRows,NCols,HeightBlocks,MyLegend=NULL,LegendRight=NULL,filename=NULL,PDFheight=10,PDFwidth=10){ 
#	isnot.null <- function(x) ! is.null(x)
#	isnot.vector <- function(x) ! is.vector(x)
#  if(is.null(MyLegend) && length(HeightBlocks) != NRows){stop("The length of each column is given in HeightBlocks. Thus, the length of HeightBlocks should be equal to the number of columns")}
#  if(isnot.null(MyLegend) && length(HeightBlocks) != 2 && is.null(LegendRight)){stop("HeightBlocks defines the height of the plots and the legend. Therefore, its length has to be equal to 2")}
#  if(isnot.null(MyLegend) && isnot.null(LegendRight) && length(HeightBlocks) != NRows){stop("The length of each column is given in HeightBlocks. Thus, the length of HeightBlocks should be equal to the number of columns")}
#  
#  if (is.null(filename)){
#    t <- c("grid.arrange(arrangeGrob(")
#    for (i in 1:length(plots)){
#      t <- paste(t,plots[i]," + theme(legend.position='none'),",sep="")
#    }
#    if (isnot.null(MyLegend)){
#      if (length(HeightBlocks)>1){
#        jj <-  c(paste0(HeightBlocks[1:length(HeightBlocks)-1],","),HeightBlocks[length(HeightBlocks)])
#      } else {
#        jj <- HeightBlocks[1]
#      }
#      jj <- paste(jj,collapse=" ")
#      t <- paste(t,"nrow= ",NRows ,", ncol= ",NCols,"),arrangeGrob(MyLegend,nrow=1),heights=c(",jj,"))" )
#    } else {
#      if (length(HeightBlocks)>1){
#        jj <-  c(paste0(HeightBlocks[1:length(HeightBlocks)-1],","),HeightBlocks[length(HeightBlocks)])
#      } else {
#        jj <- HeightBlocks[1]
#      }
#      jj <- paste(jj,collapse=" ")
#      t <- paste(t,"nrow= ",NRows ,", ncol= ",NCols ,", heights=c(",jj,")))" )
#    }
#    if(isnot.null(LegendRight)){
#      substr(t,nchar(t),nchar(t)) = ""
#      t <- paste(t,",ncol=2)")
#    }
#    return(eval(parse(text=t)))
#  } else {
#    t <- c("grid.arrange(arrangeGrob(")
#    for (i in 1:length(plots)){
#      t <- paste(t,plots[i]," + theme(legend.position='none'),",sep="")
#    }
#    if (isnot.null(MyLegend)){
#      if (length(HeightBlocks)>1){
#        jj <-  c(paste0(HeightBlocks[1:length(HeightBlocks)-1],","),HeightBlocks[length(HeightBlocks)])
#      } else {
#        jj <- HeightBlocks[1]
#      }
#      jj <- paste(jj,collapse=" ")
#      t <- paste(t,"nrow= ",NRows ,", ncol= ",NCols,"),arrangeGrob(MyLegend,nrow=1),heights=c(",jj,"))" )
#    } else {
#      if (length(HeightBlocks)>1){
#        jj <-  c(paste0(HeightBlocks[1:length(HeightBlocks)-1],","),HeightBlocks[length(HeightBlocks)])
#      } else {
#        jj <- HeightBlocks[1]
#      }
#      jj <- paste(jj,collapse=" ")
#      t <- paste(t,"nrow= ",NRows ,", ncol= ",NCols ,", heights=c(",jj,")))" )
#    }
#    if(isnot.null(LegendRight)){
#      substr(t,nchar(t),nchar(t)) = ""
#      t <- paste(t,",ncol=2)")
#    }
#    pdfname <- paste(filename,".pdf",sep="")
#    pdf(file=pdfname,width=PDFwidth,height=PDFheight)
#    eval(parse(text=t))
#    dev.off()
#    return(1) 
#  }
#}

CorrelationPlot <- function (pred,obs,margin=NULL,main="",ylab="Predicted",xlab="Observed",
                   PointSize=4,ColMargin="blue",TextSize=15,TitleSize=15,
                   XAxisSize=15,YAxisSize=15,TitleAxesSize=15,tmar=1,bmar=1,
                   rmar=1,lmar=1,AngleLab=30,LegendPosition="right",
                   PointColor="black",PointAlpha=1,
                   PointShape=16,MarginWidth=1) 
{
	isnot.null <- function(x) ! is.null(x)
	isnot.vector <- function(x) ! is.vector(x)
  if (isnot.vector(obs) || isnot.vector(pred)){
    stop("The input data must be two vectors")
  } else if ( length(obs) != length(pred) ){
    stop("Both vectors have to be of equal length")
  } else if(isnot.null(margin)) {
    Data <- data.frame(Observed=obs,Predicted=pred)
    p <- ggplot(Data, aes(x=Observed, y=Predicted)) +
      geom_point(size=PointSize,colour=PointColor,shape=PointShape,alpha=PointAlpha) +
      geom_abline(slope=1,intercept=margin/2,colour=ColMargin,size=MarginWidth) + 
      geom_abline(slope=1,intercept=-(margin/2),colour=ColMargin,size=MarginWidth) + theme_bw() + 
      ggtitle(main) + ylab(ylab) + xlab(xlab)+
      ylim(c(min(c(obs,pred)), max(c(obs,pred)))) + xlim(c(min(c(obs,pred)), max(c(obs,pred)))) +
      theme(text = element_text(size=TextSize),axis.text.x = element_text(size=XAxisSize,angle = AngleLab, hjust = 1),
            axis.title.x=element_text(size=TitleAxesSize),axis.title.y=element_text(size=TitleAxesSize),
            axis.text.y=element_text(size=YAxisSize),legend.position=LegendPosition,plot.title=element_text(size=TitleSize),
            legend.key=element_blank(), plot.margin=unit(c(tmar,rmar,bmar,lmar),"cm")) 
  } else {
    Data <- data.frame(Observed=obs,Predicted=pred)
    p <- ggplot(Data, aes(x=Observed, y=Predicted)) + geom_point(size=PointSize,colour=PointColor,shape=PointShape,alpha=PointAlpha) + theme_bw() +
      ggtitle(main) + ylab(ylab) + xlab(xlab) +
      ylim(c(min(c(obs,pred)), max(c(obs,pred)))) + xlim(c(min(c(obs,pred)), max(c(obs,pred)))) +
      theme(text = element_text(size=TextSize),axis.text.x = element_text(size=XAxisSize,angle = AngleLab, hjust = 1),
            axis.title.x=element_text(size=TitleAxesSize),axis.title.y=element_text(size=TitleAxesSize),
            axis.text.y=element_text(size=YAxisSize),legend.position=LegendPosition,plot.title=element_text(size=TitleSize),
            legend.key=element_blank(), plot.margin=unit(c(tmar,rmar,bmar,lmar),"cm"))   
  }
  return(p) 
}

RMSE <- function(v1, v2) {
  i1 <- which(!is.na(v1))
  i2 <- which(!is.na(v2))
  is <- intersect(i1, i2)
  v1 <- v1[is]
  v2 <- v2[is]
  residuals <- abs(v1-v2)
  return(as.numeric(sqrt( (residuals%*%residuals)/length(v1) )))
}

RMSE_CV <- function(model, digits = 3,metric="RMSE") {
  signif(min(as.vector(na.omit(model$results[metric]))), digits=3)
}  

MAE <- function (v1, v2) {
  i1 <- which(!is.na(v1))
  i2 <- which(!is.na(v2))
  is <- intersect(i1, i2)
  v1 <- v1[is]
  v2 <- v2[is]
  residuals <- abs(v1 - v2)
  return(sum(residuals)/length(v1))
}

slope <- function(v1,v2){ # v1=z.test v2=y.test
  return(sum(v2*v1)/sum(v1*v1))
}

Rsquared0 <- function(v1,v2) { #v1=z.test (y), v2=y.test (x)
  if (is.vector(v1) && is.vector(v2) && length(v1)==length(v2)){
    y_obs_mean <- mean(v2)
    yr0 = v1 * slope(v1,v2)
    first_term = (v2 - yr0)*(v2 - yr0)
    second_term= (v2-y_obs_mean)*(v2-y_obs_mean)
    return(1-(sum(first_term)/sum(second_term)))
  }
  else {print("Wrong input: input arguments are not vector or have unequal length")}
}


Rsquared <- function(v1,v2) { # v1=z.test (y), v2=y.test (x)
  if (is.vector(v1) && is.vector(v2) && length(v1)==length(v2)){
    y_obs_mean <- mean(v2)
    y_pred_mean <- mean(v1)
    first_term <- sum((v2-y_obs_mean) * (v1 - y_pred_mean))
    second_term <- sqrt(sum((v2-y_obs_mean)*(v2-y_obs_mean)) * sum((v1 - y_pred_mean)*(v1 - y_pred_mean)))
    division <- first_term / second_term
    return(division * division)
  }
  else {print("Wrong input: input arguments are not vector or have unequal length")}
}

Rsquared_CV <- function(model, digits = 3,metric='RMSE') {
model$results$Rsquared[model$results[metric] == min(model$results[metric], na.rm=TRUE)]
#model$results$Rsquared[which(model$results[metric] %in% min(model$results[metric], na.rm=TRUE))]
} 

Qsquared1 <- function(v1, v2, resp_tr) {
  if (is.vector(v1) && is.vector(v2) && length(v1)==length(v2)){
    y_tr_mean <- mean(resp_tr)
    first_term <- abs(v1-v2)*abs(v1-v2)
    second_term <- abs(v2-y_tr_mean)*abs(v2-y_tr_mean)
    return(1-(sum(first_term)/sum(second_term)))
  }
  else {print("Wrong input: input arguments are not vector or have unequal length")}
}

Qsquared2 <- function(v1, v2) {
  if (is.vector(v1) && is.vector(v2) && length(v1)==length(v2)){
    y_obs_mean <- mean(v2)
    first_term <- abs(v1-v2)*abs(v1-v2)
    second_term <- abs(v2-y_obs_mean)*abs(v2-y_obs_mean)
    return(1-(sum(first_term)/sum(second_term)))
  }
  else {print("Wrong input: input arguments are not vector or have unequal length")}
}

Qsquared3 <- function(v1, v2,resp_tr) {
  if (is.vector(v1) && is.vector(v2) && length(v1)==length(v2)){
    y_obs_mean <- mean(v2)
    y_tr_mean <- mean(resp_tr)
    first_term <- abs(v1-v2)*abs(v1-v2)
    first_term <- sum(first_term)/length(v1)
    second_term <- abs(resp_tr-y_tr_mean)*abs(resp_tr-y_tr_mean)
    second_term <- sum(second_term)/length(resp_tr)
    return(1-(first_term/second_term))
  }
  else {print("Wrong input: input arguments are not vector or have unequal length")}
}

Validation <- function(pred, obs, resp_tr){
  if (is.vector(pred) && is.vector(obs) && length(pred)==length(obs)){
    metrics <- list(R2 = Rsquared(pred,obs), R02 = Rsquared0(pred,obs), 
   Q2_1 = Qsquared1(pred,obs,resp_tr),Q2_2 = Qsquared2(pred,obs),Q2_3 = Qsquared3(pred,obs,resp_tr),
   RMSE = RMSE(pred,obs), Slope=slope(pred,obs), MAE = MAE(pred, obs))
  } else {
    stop("Wrong input: input arguments are not vector or have unequal length")
  }
  return(metrics)
}


