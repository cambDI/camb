PCAPlot <- function (Data,labels = NULL, Colors=NULL, Shapes=NULL,
                     main = "", ylab = "PC2", xlab = "PC1",  
          PointSize = 4, LegendPosition = "right", LegendName = "", 
          ColLegend = 1, RowLegend = NULL, TitleSize = 15, TextSize = 15, 
          XAxisSize = 15, YAxisSize = 15, AngleLab = 0, TitleAxesSize = 15, 
          LegendTitleSize = 15, LegendTextSize = 15, tmar = 1, bmar = 1, 
          rmar = 1, lmar = 1) 
{
  isnot.null <- function(x) !is.null(x)
  isnot.vector <- function(x) !is.vector(x)
  if (length(names(Data)) < 2 || length(names(Data)) > 3) {
    stop("Two PCs required. The Data.frame provided has less than two columns (PCs) or more than 3")
  }
  else if (names(Data)[1] != "PC1" || names(Data)[2] != "PC2") {
    stop("Column names have to be {PC1, PC2}")
  }
  else if (length(names(Data)) == 2 && is.null(labels)) {
    print("No names provided")
    p <- ggplot(Data, aes(x = PC1, y = PC2)) + geom_point(size = PointSize) + 
      theme_bw() + ggtitle(main) + ylab(ylab) + xlab(xlab) + 
      theme(text = element_text(size = TextSize), axis.text.x = element_text(size = XAxisSize, 
                                                                             angle = AngleLab, hjust = 1), axis.title.x = element_text(size = TitleAxesSize), 
            axis.title.y = element_text(size = TitleAxesSize), 
            axis.text.y = element_text(size = YAxisSize), 
            legend.position = LegendPosition, plot.title = element_text(size = TitleSize), 
            legend.key = element_blank(), legend.text = element_text(size = LegendTextSize), 
            legend.title = element_text(size = LegendTitleSize), 
            plot.margin = unit(c(tmar, rmar, bmar, lmar), 
                               "cm")) + 
      guides(colour = guide_legend(LegendName,ncol = ColLegend, nrow = RowLegend), shape = guide_legend(LegendName, 
                                                                                                                               ncol = ColLegend, nrow = RowLegend))
  }
  else if (length(names(Data)) == 2 && isnot.null(labels)) {
    if (length(labels) != nrow(Data) || isnot.vector(labels)) {
      stop("Either the names are not in a vector, or its length is not equal to the number of datapoints (rows of the input data)")
    }
    else {
      if (is.null(Colors))
        Colors=labels
      if (is.null(Shapes))
        Shapes=labels
      print("Names provided")
      Data <- data.frame(Data, labels = labels,Colors=Colors,Shapes=Shapes)
      a <- length(unique(labels))
      p <- ggplot(Data, aes(x = PC1, y = PC2, color = Colors, 
                            shape = Shapes)) + geom_point(size = PointSize) + 
        scale_shape_manual(values = 1:a) + theme_bw() + 
        ggtitle(main) + ylab(ylab) + xlab(xlab) + theme(text = element_text(size = TextSize), 
                                                        axis.text.x = element_text(size = XAxisSize, 
                                                                                   angle = AngleLab, hjust = 1), axis.title.x = element_text(size = TitleAxesSize), 
                                                        axis.title.y = element_text(size = TitleAxesSize), 
                                                        axis.text.y = element_text(size = YAxisSize), 
                                                        legend.position = LegendPosition, legend.text = element_text(size = LegendTextSize), 
                                                        legend.title = element_text(size = LegendTitleSize), 
                                                        plot.title = element_text(size = TitleSize), 
                                                        legend.key = element_blank(), plot.margin = unit(c(tmar, 
                                                                                                           rmar, bmar, lmar), "cm")) + guides(colour = guide_legend(LegendName, 
                                                                                                                                                                    ncol = ColLegend, nrow = RowLegend), shape = guide_legend(LegendName, 
                                                                                                                                                                                                                              ncol = ColLegend, nrow = RowLegend))
    }
  }
  else {
    print("Names provided in the third column of the data.frame")
    col=NULL
    labels <- unlist(Data[names(Data)[3]])
    a <- length(unique(labels))
    shapes_default <- c(15,16,17,18)
    sha <- scale_shape_manual(values = rep(shapes_default,a)[1:a]) 
    namee <- names(Data)[3]
    if (!is.null(Colors)){
      if (length(Colors) != nrow(Data)){stop("If you provide Colors, the number of elements must equal the number of rows in Data")}
      col <- scale_colour_manual(values = Colors)}
    if (!is.null(Shapes)){
      if (length(Shapes) != nrow(Data)){stop("If you provide Shapes, the number of elements must equal the number of rows in Data")}
     sha <- scale_shape_manual(values = Shapes)}
      
    
    #Data <- data.frame(Data, Colors=Colors)#,Shapes=Shapes)
    print(head(Data))
    p <- ggplot(Data, aes_string(x = "PC1", y = "PC2", colour = namee, shape = namee))+ 
      geom_point(size = PointSize) + 
      theme_bw() + ggtitle(main) + 
      ylab(ylab) + xlab(xlab) + theme(text = element_text(size = TextSize), 
                                      axis.text.x = element_text(size = XAxisSize, angle = AngleLab, 
                                                                 hjust = 1), axis.title.x = element_text(size = TitleAxesSize), 
                                      axis.title.y = element_text(size = TitleAxesSize), 
                                      axis.text.y = element_text(size = YAxisSize), legend.position = LegendPosition, 
                                      plot.title = element_text(size = TitleSize), legend.text = element_text(size = LegendTextSize), 
                                      legend.title = element_text(size = LegendTitleSize), 
                                      legend.key = element_blank(), 
                                      plot.margin = unit(c(tmar,rmar, bmar, lmar), "cm")) + 
      guides(colour = guide_legend(LegendName,ncol = ColLegend, nrow = RowLegend), shape = guide_legend(LegendName, 
                        ncol = ColLegend, nrow = RowLegend))
   if(!is.null(sha))
      p <- p + sha
   if(!is.null(col))
      p<-p+col
    
  }
  return(p)
}
