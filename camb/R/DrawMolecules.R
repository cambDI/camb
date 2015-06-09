DrawMoleculeInSDF <- function(structures.file, structure.number, file.name, useNameAsTitle) {
  if (file.info(structures.file)$size  == 0) {stop("Input file is empty")}
  print(structure.number)
  .C("R_drawMoleculeInSDF", structures.file, as.integer(structure.number), file.name, useNameAsTitle)
}

# some testing code left here for now
#.C("R_drawMoleculeInSDF", system.file("test_structures", "structures_10.sdf", package = "camb"), 1, "output.png", FALSE)

PlotMolecules <- function(sdf.file, IDs,pdf.file=NULL,PDFMain=NULL,useNameAsTitle=TRUE) {
isnot.vector <- function(x) ! is.vector(x) 
isnot.null <- function(x) ! is.null(x)
  #suppressWarnings(require(pbapply,quietly=TRUE,warn.conflicts=FALSE)) || stop("Pacakge 'pbapply' is required")
  is.installed <- function(mypkg) is.element(mypkg, installed.packages()[,1])
  if (length(IDs) !=4) {stop("Only four compounds per plot supported at the moment..")}
  if (file.info(sdf.file)$size  == 0) {stop("Input file is empty")}
  #if (!is.installed('png')){stop("Pacakge 'png' is required")}
  #if (!is.installed('grid')){stop("Pacakge 'grid' is required")}

  temp.png <- tempfile("temp", fileext=".png")
  base <- theme(axis.line = element_blank(), 
                axis.text.x = element_blank(), 
                axis.text.y = element_blank(),
                axis.ticks = element_blank(), 
                axis.title.x = element_blank(), 
                axis.title.y = element_blank(),
                panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                plot.margin=unit(c(0,0,0,0), "cm"))
  i <- 1
  for(id in IDs) {
    DrawMoleculeInSDF(structures.file=sdf.file, structure.number=as.integer(id), temp.png,useNameAsTitle=TRUE)
    img <- readPNG(temp.png)
    g <- rasterGrob(img, interpolate=TRUE)
    p <- qplot(1, 1, geom="blank") + theme_bw() 
    assign(paste("p",i,sep=""), p +annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) + base)
    i <- i+1
  }
  if (isnot.null(pdf.file)){
  pdf(pdf.file)
  grid.arrange(p1,p2,p3,p4,nrow=2,main=PDFMain)
  dev.off()
  }
  return(list(p1,p2,p3,p4))
}


