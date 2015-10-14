#################################################################################
## Standardization and Descriptor Calculation of Molecules and Peptides / Sequences
#################################################################################

StandardiseMolecules <- function(structures.file, 
                                 standardised.file,
                                 removed.file = "",
                                 properties.file="standardisation_info.csv",
                                 remove.inorganic = FALSE, 
                                 fluorine.limit = -1,
                                 chlorine.limit = -1,
                                 bromine.limit = -1,
                                 iodine.limit = -1,
                                 min.mass.limit = -1,
                                 max.mass.limit = -1,                      
                                 number.processed = -1) {
  # handle non-existant file
  if (!file.exists(structures.file)) {
    stop("File does not exist")
  }
  if (file.info(structures.file)$size  == 0) {stop("Input file is empty")}
  
  # deal with sdf or smi difference
  split <- strsplit(structures.file, "\\.")[[1]]
  filetype <- split[length(split)]
  if(tolower(filetype) == "sdf") {
    print("Standardising Structures: Reading SDF (R)")
    sink(file="standardisation.log", append=FALSE, split=FALSE)
    .C("R_standardiseMolecules", 
       structures.file, 
       standardised.file, 
       removed.file,
       properties.file,
       as.integer(1), # process SDF
       as.integer(remove.inorganic), 
       as.integer(fluorine.limit),
       as.integer(chlorine.limit),
       as.integer(bromine.limit),
       as.integer(iodine.limit),
       as.integer(min.mass.limit),
       as.integer(max.mass.limit),
       as.integer(number.processed)) 
    sink()
  }
  else if(tolower(filetype) == "smi") {
    print("Standardising Structures: Reading SMILES (R)")
    sink(file="standardisation.log", append=FALSE, split=FALSE)
    .C("R_standardiseMolecules", 
       structures.file, 
       standardised.file, 
       removed.file,
       properties.file,
       as.integer(0), # process SMILES
       as.integer(remove.inorganic), 
       as.integer(fluorine.limit),
       as.integer(chlorine.limit),
       as.integer(bromine.limit),
       as.integer(iodine.limit),
       as.integer(min.mass.limit),
       as.integer(max.mass.limit),
       as.integer(number.processed)) 
    sink()
  }
  else {
    print("Unrecognised file type")
  }
  
  return(c(remove.inorganic, fluorine.limit, chlorine.limit, bromine.limit, iodine.limit, min.mass.limit, max.mass.limit))
}

GetPropertiesSDF <- function(structures.file,number_processed=-1){ 
if (!file.exists(structures.file)) {stop("File does not exist")}
if (file.info(structures.file)$size  == 0) {stop("Input file is empty")}
output <- tempfile("prop_temp",fileext=".csv",tempdir())
.C("R_GetPropertiesSDF",structures.file,as.integer(number_processed),output)
print("Reading..")
properties <- read.csv(output,sep="\t",header=TRUE,comment.char="",quote = "")
properties <- properties[, !apply(is.na(properties), 2, all)]
return(properties)
}

ShowPropertiesSDF <- function(structures.file,type=1){ ## 1 refers to not smiles
if (!file.exists(structures.file)) {stop("File does not exist")}
if (file.info(structures.file)$size  == 0) {stop("Input file is empty")}
output <- tempfile("props_temp",fileext=".csv")
.C("R_ShowPropertiesSDF",structures.file,output,as.integer(type))
if (file.info(output)$size  == 0) {stop("The molecules in the file provided do not contain any property")}
props <- as.vector(read.table(output,header=FALSE)$V1)
return(props)
}

GetPropertySDF <- function(structures.file, property="", number_processed=-1){
if (!file.exists(structures.file)) {stop("File does not exist")}
if (file.info(structures.file)$size == 0) {stop("Input file is empty")}
output <- tempfile("prop_temp",fileext=".csv")
.C("R_GetPropertySDF",structures.file,property,as.integer(number_processed),output)
if (file.info(output)$size  == 0) {stop("The molecules in the file provided do not contain the specified property")}
prop <- read.csv(output,comment.char="",quote ="",header=T)
names(prop) <- property
return(prop)
}


##############
install_and_load <- function (package1, ...)  {   
   # convert arguments to vector
      packages <- c(package1, ...)
     # start loop to determine if each package is installed
    for(package in packages){
       # if package is installed locally, load
          if(package %in% rownames(installed.packages()))
            do.call('library', list(package))
       # if package is not installed locally, download, then load
          else {
          cat("The package is not in the library..\nInstalling it..\n") 
          install.packages(package)
          do.call("library", list(package))
         }
    } 
}


## Compound Descriptors
RemoveStandardisedPrefix <- function(descriptors) {
  descriptors$Name <- sapply(descriptors$Name, function(x) {strsplit(as.character(x), "Standardised_")[[1]][2]})
  descriptors
}

GeneratePadelDescriptors <- function(standardised.file, types = c("2D"), threads = -1, limit = -1) {
  cat("Verifying that rJava is in the library. If not, it will be installed..\n")
  install_and_load("rJava")
  if (file.info(standardised.file)$size  == 0) {stop("Input file is empty")}
  descriptors.file <- tempfile("descriptors", fileext=".csv")
  GeneratePadelDescriptors.internal(standardised.file, descriptors.file, types, threads)
  read.csv(descriptors.file)
}

GeneratePadelDescriptorsFile <- function(standardised.file, descriptors.file, types = c("2D"), threads = -1, limit = -1) {
  if (file.info(standardised.file)$size  == 0) {stop("Input file is empty")}
  GeneratePadelDescriptors.internal(standardised.file, descriptors.file, types, threads)
}

GeneratePadelDescriptors.internal <- function(structures.file, descriptors.file, types, threads = -1) {
  print("Generating Descriptors")
  .jinit()
  .jcall("java/lang/System","S","setProperty","java.awt.headless","true")
  # add all JARs
  #.jaddClassPath(Sys.glob("lib/*.jar")) this was the old code but it stopped working for some reason
  .jaddClassPath(Sys.glob(paste(system.file(package = "camb"), "java", "*.jar", sep="/")))
  # call the main() method on the main class
  
  # replace the options in the padel_config.txt file
  padel_config_file <- system.file("extdata", "padel_config.txt", package="camb")
  readCon  <- file(padel_config_file, open = "r")
  conffile <- tempfile("padel_config", fileext=".txt")
  confCon <- file(conffile)
  lines <- c()
  while (length(line <- readLines(readCon, n = 1, warn = FALSE)) > 0) {
    vec <- unlist(strsplit(line, "="))
    if(vec[1]=="Compute2D") {
      if("2D" %in% types) {
        lines <- c(lines, "Compute2D=true")
      }
    }
    else if(vec[1]=="MaxThreads") {
      lines <- c(lines, paste("MaxThreads=", threads, sep=""))
    }
    else if(vec[1]=="DescriptorFile") {
      lines <- c(lines, paste("DescriptorFile=", descriptors.file, sep=""))
    } 
    else if(vec[1]=="Directory") {
      lines <- c(lines, paste("Directory=", structures.file, sep=""))
    }
    else {
      lines <- c(lines, line)
    }
  }
  writeLines(lines, confCon)
  close(readCon)
  close(confCon)
  
  # replace the options in the padel_types.xml file
  padel_types_file <- system.file("extdata", "padel_types.xml", package="camb")
  readCon  <- file(padel_types_file, open = "r")
  descfile <- tempfile("padel_types", fileext=".xml")
  descCon <- file(descfile)
  lines <- c()
  while (length(line <- readLines(readCon, n = 1, warn = FALSE)) > 0) {
    vec <- unlist(strsplit(line, "\""))
    if(vec[2] %in% types) {
      vec[4] <- "true"
      lines <- c(lines, paste(vec, collapse="\""))
    }
    else {
      lines <- c(lines, line)
    }
  }
  writeLines(lines, descCon)
  close(readCon)
  close(descCon)
  
  .jcall("padeldescriptor.PaDELDescriptorApp", , "main", c("-config", conffile, "-descriptortypes", descfile))
}

##############
# Check if an AA is natural
checkAA <- function(x) {
  AADict = c("A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", 
             "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V",
			 "-", ## for the gaps
             "ALA","CYS","ASP","GLU","PHE","GLY","HIS","ILE","LYS",
             "LEU","MET","ASN","PRO","GLN","ARG","SER","THR","VAL",
             "TRP","TYR")
  return(x %in% AADict)
}

##############
## Three letter to one letter AA code
convert31 <- function(AA) {  
  threeL <- c("ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN",
              "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE",
              "PRO", "SER", "THR", "TRP", "TYR", "VAL","-")  
  names(threeL) <- c("A", "R", "N", "D", "C", "E", "Q", "G", "H", "I",
                     "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V","-")  
  res <- names(threeL[which(threeL==toupper(AA))])
  return(res)
}

AADescs <- function(Data, type="Z5",...){
  if (!is.vector(Data) && !is.character(Data) && !is.data.frame(Data) && !is.matrix(Data)){
    stop("Input must be a character, vector, data frame or matrix")
  } else {    
    match_AA1 <- function(AA,sel.=sel){
      if (checkAA(toupper(AA))){
        AA <- toupper(AA)
        if (nchar(AA) == 1){
          d <- descs[descs$AA1 == AA,sel.+2]
        } else {
          AA <- convert31(AA)
          d <- descs[descs$AA1 == AA,sel.+2]
        }
        return(d)
      } else {
        stop("Non natural Amino acid provided")
      }
    }
    
    match_AA1_vec <-function(v,sel.=sel){
      if (is.vector(v)){
        des <- sapply(v,match_AA1,sel.=sel)
        namesDes <- row.names(des)
        des <- unlist(des)
        namesDes <- c(t(sapply(namesDes,paste0,"_",v,seq(1,length(v)))))
        names(des) <- namesDes
        return(des)
      } else {
        stop("Input is not a vector or character")
      }
    }
    
    match_AA1_df <- function(df,colNames,sel.){
      if (is.data.frame(df) || is.matrix(df)){
        des <- t(apply(df,1,match_AA1_vec))
        row.names(des) <- seq(1,nrow(df))
        tableDescs <- table(colNames)
        colNamesPos=c()
        for(i in 1:length(tableDescs)){
          now <- which(colNames %in% names(tableDescs)[i])
          colNamesPos <- append(colNamesPos,paste0(colNames[now],"_",seq(1,length(now))))
        }
        colnames(des) <- c(t(sapply(colNamesPos,paste0,"_Aa",seq(1,ncol(df)))))
        return(des)
      } else {
        stop("Input is neither a data.frame nor a matrix")
      }
    }
    

    if  (is.vector(Data)){
	  descs_path <- system.file("extdata", "aa_descs.rds", package="camb")
	  descs <- readRDS(descs_path)
      types <- c("ProtFP8","TScales","VHSE","STScales","BLOSUM","FASGAI","MSWHIM","Z5","Z3")
      type <- match.arg(type,types,several.ok=TRUE)
      root <- strsplit(names(descs)[3:ncol(descs)],"_")
      root <- unlist(root)[seq(1,length(unlist(root)),2)]
      sel <- which(root %in% type)
      return(match_AA1_vec(Data,sel))

    } else {
	  descs_path <- system.file("extdata", "aa_descs.rds", package="camb")
	  descs <- readRDS(descs_path)
	  types <- c("ProtFP8","TScales","Tscales","VHSE","STScales","BLOSUM","FASGAI","MSWHIM","Z5","Z3")
      type <- match.arg(type,types,several.ok=TRUE)
      root <- strsplit(names(descs)[3:ncol(descs)],"_")
      root <- unlist(root)[seq(1,length(unlist(root)),2)]
      typeExt <- root[which(root %in% type)]
      sel <- which(root %in% type)
      return(match_AA1_df(Data,colNames=typeExt,sel=sel))
    }
  }
}


SeqDescs <- function(data,UniProtID=TRUE,type="AAC",...){
  if (UniProtID){
    if (!is.vector(data)) stop("The input UNIPROT identifiers must be in a vector")
    d <- data.frame(unlist(lapply(data,getUniProt)))
    types <- c("AAC","DC","TC","MoreauBroto","Moran","Geary",
               "CTDC","CTDT","CTDD","CTriad","SOCN","QSO",
               "PACC","APAAC")
    type <- match.arg(type,types,several.ok=TRUE)
    type <- paste0("extract",type,"(x)")
    type <- paste("c(",paste(type,collapse=","),")",sep="")
    t <- paste("t(apply(d,1,FUN=function(x)",type,"))",sep="") 
    des=eval(parse(text=t))
    row.names(des) <- data
    return(des) 
  } else {
    if (!is.data.frame(data) && !is.matrix(data)) stop("The input sequences must be a dataframe or a matrix")
    types <- c("AAC","DC","TC","MoreauBroto","Moran","Geary",
               "CTDC","CTDT","CTDD","CTriad","SOCN","QSO",
               "PACC","APAAC")
    type <- match.arg(type,types,several.ok=TRUE)
    type <- paste0("extract",type,"(x)")
    type <- paste("c(",paste(type,collapse=","),")",sep="")
    t <- paste("t(apply(data,1,FUN=function(x)",type,"))",sep="")
    des=eval(parse(text=t))
    row.names(des) <- rownames(data)
    return(des) 
  }
}



##############
## Circular Morgan Fingerprints as specified in RDkit
############## 

MorganFPs <- function (bits=512,radii=c(0,1,2),mols,output,keep="hashed_binary",images=FALSE,unhashed=FALSE,
                       RDkitPath="/usr/local/share/RDKit",PythonPath="/usr/local/lib/python2.7/site-packages",
					   extMols=FALSE,unhashedExt=FALSE,logFile=FALSE) {
  if (is.null(RDkitPath)) RDkitPath <- "/usr/local/share/RDKit"
  if (is.null(PythonPath)) RDkitPath <- "/usr/local/lib/python2.7/site-packages"
  
  pyfpath <- system.file("extdata", "FingerprintCalculator.py", package="camb")
  t <- paste(pyfpath," -bits",bits,"-radii",radii,#"--f",type,
             "-mols",mols,"-output",output,"-RDkitPath",RDkitPath,sep=" ")
  if (images) t <- paste(t,"-image",sep=" ")
  if (unhashed) t <- paste(t,"-unhashed",sep=" ")
  if (extMols!=FALSE) t <- paste(t,"-molsEXT",extMols,sep=" ")
  if (unhashedExt!=FALSE) t <- paste(t, "-unhashedEXT",sep=" ")
  if (logFile) t <- paste(t, " > ",output,".log",sep="")
  system(t) 
  types_keep <- c("hashed_binary","hashed_counts","unhashed_binary","unhashed_counts",
                  "hashed_binaryEXT","hashed_countsEXT","unhashed_binaryEXT","unhashed_countsEXT")
  keep <- match.arg(keep,types_keep,several.ok=FALSE)
  if (is.na(keep)){
    stop("Cannot load the fingerprints file. Possible types are:
         hashed_binary, hashed_counts, unhashed_binary, unhashed_counts, hashed_binaryEXT, hashed_countsEXT, unhashed_binaryEXT and unhashed_countsEXT")
  }  
  loadFile=paste(output,"_",keep,".csv",sep="")
  p=read.table(loadFile,sep=",")
  return(p)
}

##############


#loadMorganFPs <- function (type="hashed_binary",output){
#  types_keep <- c("hashed_binary","hashed_counts","unhashed_binary","unhashed_counts",
#                  "hashed_binaryEXT","hashed_countsEXT","unhashed_binaryEXT","unhashed_countsEXT")
#  type <- match.arg(type,types_keep,several.ok=FALSE)
#  if (is.na(type)){
#    stop("Cannot load the fingerprints file. Possible types are:
#         hashed_binary, hashed_counts, unhashed_binary, unhashed_counts, hashed_binaryEXT, hashed_countsEXT, unhashed_binaryEXT and unhashed_countsEXT")
#  }  
#  loadFile=paste(output,type,".csv",sep="_")
#  p=read.table(loadFile,sep=",")
#  return(p)
#}


mergeData <- function(a1,a2,a3,a4=NULL,a5=NULL,a6=NULL){
  return(cbind(a1,a2,a3,a4,a5,a6,deparse.level=0))
}

