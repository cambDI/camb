# Function to calculate the minimum model performance

MinPerf <- function (meanNoise = 0, sdNoise, resp, lenPred, stds=NULL,
                             iters = 1000, filename = NULL, pdfW = 10, pdfH = 10, TextSize = 15, 
                             TitleSize = 15, XAxisSize = 15, YAxisSize = 15, TitleAxesSize = 15, 
                             tmar = 1, bmar = 1, rmar = 1, lmar = 1, AngleLab = 30, LegendPosition = "right") 
{
  isnot.vector <- function(x) ! is.vector(x)
  isnot.null <- function(x) ! is.null(x)
  vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
  if (is.null(stds)){stds=rnorm(length(resp),mean=meanNoise,sd=sdNoise)}else{stds=stds*sample(c(1,-1),size=len(resp),replace=T)}
  R2 <- c()
  R02 <- c()
  Q2 <- c()
  rmsep <- c()
  for (i in 1:iters) {
    ##
    set.seed(i)
    idx <- sample(seq(1,length(resp)),lenPred, replace=FALSE)
    x <- resp[idx]
    noise <- stds[idx]
    y <- x + noise
    y <- sample(y,lenPred,replace=FALSE)
    ##
    R2[i] <- Rsquared(y, x)
    Q2[i] <- Qsquared2(y, x)
    R02[i] <- Rsquared0(y, x)
    rmsep[i] <- RMSE(y, x)
  }
  pp <- data.frame(R2 = R2)
  title <- expression(paste("R"^{
    2
  }))
  p1 <- ggplot(pp, aes(x = R2)) + theme_bw() + geom_density(alpha = 0.2, 
                                                            fill = "#FF6666") + ggtitle(title) + aes(y = ..count..) + 
    theme(text = element_text(size = TextSize), axis.text.x = element_text(size = XAxisSize, 
                                                                           angle = AngleLab, hjust = 1), axis.title.x = element_text(size = TitleAxesSize), 
          axis.title.y = element_text(size = TitleAxesSize), 
          axis.text.y = element_text(size = YAxisSize), legend.position = LegendPosition, 
          plot.title = element_text(size = TitleSize), legend.key = element_blank(), 
          plot.margin = unit(c(tmar, rmar, bmar, lmar), "cm")) + 
    geom_vline(aes(xintercept = mean(R2, na.rm = T)), color = "blue", 
               linetype = "dashed", size = 1)
  pp <- data.frame(rmsep = rmsep)
  p2 <- ggplot(pp, aes(x = rmsep)) + theme_bw() + geom_density(alpha = 0.2, 
                                                               fill = "#FF6666") + ggtitle("RMSEP") + aes(y = ..count..) + 
    theme(text = element_text(size = TextSize), axis.text.x = element_text(size = XAxisSize, 
                                                                           angle = AngleLab, hjust = 1), axis.title.x = element_text(size = TitleAxesSize), 
          axis.title.y = element_text(size = TitleAxesSize), 
          axis.text.y = element_text(size = YAxisSize), legend.position = LegendPosition, 
          plot.title = element_text(size = TitleSize), legend.key = element_blank(), 
          plot.margin = unit(c(tmar, rmar, bmar, lmar), "cm")) + 
    geom_vline(aes(xintercept = mean(rmsep, na.rm = T)), 
               color = "blue", linetype = "dashed", size = 1)
  pp <- data.frame(R02 = R02)
  title <- expression(paste("R"[0]^{
    2
  }))
  p3 <- ggplot(pp, aes(x = R02)) + theme_bw() + geom_density(alpha = 0.2, 
                                                             fill = "#FF6666") + ggtitle(title) + aes(y = ..count..) + 
    theme(text = element_text(size = TextSize), axis.text.x = element_text(size = XAxisSize, 
                                                                           angle = AngleLab, hjust = 1), axis.title.x = element_text(size = TitleAxesSize), 
          axis.title.y = element_text(size = TitleAxesSize), 
          axis.text.y = element_text(size = YAxisSize), legend.position = LegendPosition, 
          plot.title = element_text(size = TitleSize), legend.key = element_blank(), 
          plot.margin = unit(c(tmar, rmar, bmar, lmar), "cm")) + 
    geom_vline(aes(xintercept = mean(R02, na.rm = T)), color = "blue", 
               linetype = "dashed", size = 1)
  title <- expression(paste("Q"^{
    2
  }))
  pp <- data.frame(Q2 = Q2)
  p4 <- ggplot(pp, aes(x = Q2)) + theme_bw() + geom_density(alpha = 0.2, 
                                                            fill = "#FF6666") + ggtitle(title) + aes(y = ..count..) + 
    theme(text = element_text(size = TextSize), axis.text.x = element_text(size = XAxisSize, 
                                                                           angle = AngleLab, hjust = 1), axis.title.x = element_text(size = TitleAxesSize), 
          axis.title.y = element_text(size = TitleAxesSize), 
          axis.text.y = element_text(size = YAxisSize), legend.position = LegendPosition, 
          plot.title = element_text(size = TitleSize), legend.key = element_blank(), 
          plot.margin = unit(c(tmar, rmar, bmar, lmar), "cm")) + 
    geom_vline(aes(xintercept = mean(Q2, na.rm = T)), color = "blue", 
               linetype = "dashed", size = 1)
  if (isnot.null(filename)) {
    pdfname = paste(filename, ".pdf", sep = "")
    pdf(pdfname, width = pdfW, height = pdfH)
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(2, 2)))
    print(p1, vp = vplayout(1, 1))
    print(p2, vp = vplayout(1, 2))
    print(p3, vp = vplayout(2, 1))
    print(p4, vp = vplayout(2, 2))
    dev.off()
  }
  p <- list()
  p$p1 <- p1
  p$p2 <- p2
  p$p3 <- p3
  p$p4 <- p4
  return(p)
}
