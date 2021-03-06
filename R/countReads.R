## main function that takes a gtf file and bam files to count read for continuous windows
### The function requires for each sample, input and IP named for the same prefix (filename). input and IP each have postfix of input.bam/m6A.bam
#### If multiple IP sample used a shared input library, sharedInput can be specified.
#' @title countReads
#' @description This is the very first function in MeRIP-seq data analysis that initianize a `MeRIP` object. This function takes BAM files of Input/IP library of each samples as input and use given GTF file as gene annotation to divide genes into consecutive bins of user defined size.  
#' @param samplenames The names of each sample (prefix for bam files).
#' @param gtf The gtf format gene annotation file
#' @param fragmentLength The RNA fragment length (insert size of the library).
#' @param modification The modification used to name the BAM files.
#' @param bamFolder Path to the folder where bam file locates
#' @param binSize The size of consecutive bins to slice the transcripts
#' @param threads The number of threads to use for hyperthreading
#' @param strandToKeep According to library preparation protocol, choose which strand to count. Stranded RNA library usually seq the "ooposite" strand. Small RNA library seq the "same" strand.
#' @param outputDir The directory to save output files
#' @param saveOutput Logical option indicating whether to save output as an RDS file.
#' @param paired Logical indicating whether the input bam files are from paired end sequencing. Default is FALSE. If using paired end data, the read length will be estimated from the data and only good mate are counted.
#' @export
countReads<-function(
  samplenames,# file name of samples
  gtf, # gtf file used for peak calling
  fragmentLength = 150,
  bamFolder,
  outputDir=NA,
  modification = "m6A",
  binSize = 50,
  strandToKeep = "opposite",
  paired = FALSE,
  threads = 1,
  saveOutput = FALSE
){

  ##read bam files
  bamPath.input = paste(bamFolder,"/",samplenames,".input.bam",sep="")
  bamPath.IP = paste(bamFolder,"/",samplenames,".",modification,".bam",sep="")
  no.samples = length(samplenames)

  ## Check for missing files and index bam files
  if( !all(file.exists(bamPath.input)) ) stop( "input bam file missing!!!" )
  if( !all(file.exists(bamPath.IP)) ) stop( "IP bam file missing!!!" )
  num_bam_files <- length(bamPath.input)
  for (ibam in 1:num_bam_files) {
    inputfile = bamPath.input[ibam]
    IPfile = bamPath.IP[ibam]
    if (! file.exists(paste(inputfile,'.bai',sep=""))) {
      print(paste("Stage: index bam file", inputfile))
      indexBam(inputfile)
    }
    if (! file.exists(paste(IPfile,'.bai',sep=""))) {
      print(paste("Stage: index bam file", IPfile))
      indexBam(IPfile)
    }
  }


  
  ## This step removes ambiguous annotations and returns gene model
  cat("Reading gtf file to obtain gene model\nFilter out ambiguous model...\n")
  geneGRList = gtfToGeneModel(gtf) #get the gene model in GRList with only single chromosome and strand.
  cat("Gene model obtained from gtf file...\n")

  ## Check BAM headers and remove chr in geneModel that is not in BAM file. 
  bamHeader <- scanBamHeader(bamPath.input, what=c("targets") )
  seqLevels <- unique( unlist( lapply( bamHeader, function(x) names( x$targets) ) ) )
  geneGRList <- geneGRList[ unlist( runValue( seqnames( geneGRList ) ) ) %in% seqLevels ]
  
  
  no.genes=length(geneGRList)## define number of genes

  cat("counting reads for each genes, this step may takes a few hours....\n")
  start_time <- Sys.time()
  registerDoParallel( cores = threads)
  cat(paste("Hyper-thread registered:",getDoParRegistered(),"\n"))
  cat(paste("Using",getDoParWorkers(),"thread(s) to count reads in continuous bins...\n"))
  reads <- foreach(i = 1:no.genes, .combine = rbind) %dopar%{

    geneName = names(geneGRList)[i]
    geneModel =reduce( geneGRList[geneName][[1]] )## merge overlapping exons

    # DNA location to gene location conversion
    df.geneModel= as.data.frame(geneModel) ##data frame of gene model
    dna.range = as.data.frame(range(geneModel))
    df.geneModel$end = df.geneModel$end - dna.range$start + 1
    df.geneModel$start = df.geneModel$start - dna.range$start + 1
    DNA2RNA = rep(0,dna.range$end - dna.range$start +1)
    no.exon = dim(df.geneModel)[1]
    for (j in 1:no.exon){DNA2RNA[df.geneModel$start[j]:df.geneModel$end[j]]=1}
    exon.length = sum(DNA2RNA)
    DNA2RNA=cumsum(DNA2RNA)*DNA2RNA

    ## skip any gene with smaller than 200bp transcript
    if(exon.length < 200) {return(NULL)}

    #creat a corresponding map from RNA to DNA
    #RNA2DNA = 1:exon.length
    #pointer = 1
    #for (j in 1:no.exon){
    #  RNA2DNA[pointer:(pointer+df.geneModel$width[j]-1) ]= RNA2DNA[pointer:(pointer+df.geneModel$width[j]-1)] + df.geneModel$start[j] -pointer
    #  pointer = pointer + df.geneModel$width[j]
    #}
    #RNA2DNA = RNA2DNA + dna.range$start -1 #back to chromosome coordinates

    ## switch strand because stranded RNA library protocol sequence reverse strand
    if(strandToKeep == "opposite"){
      reads.strand = character()
      if(dna.range$strand == "+"){reads.strand = "-"}else{reads.strand = "+"} ## switch strand on RNA reads for Truseq protocol
    }else if(strandToKeep == "same"){
      reads.strand = as.character(dna.range$strand)
    }else{
      cat("Currently m6Amonter only support strand specific RNA-seq data.\nCounting reads at opposite strand by defalt...\n")
      reads.strand = character()
      if(dna.range$strand == "+"){reads.strand = "-"}else{reads.strand = "+"}
    }

    #creat center points of continuous window
    if(exon.length <= binSize){
      slidingStart= exon.length/2
      #mapping = data.frame(start = RNA2DNA[slidingStart-exon.length/2+1], end = RNA2DNA[slidingStart + exon.length/2]  )
    }else{
      slidingStart= round(seq(from = binSize/2, to = (exon.length - binSize/2), length.out = ceiling(exon.length/binSize) ) )
      #mapping = data.frame(start = RNA2DNA[slidingStart - binSize/2 +1], end = RNA2DNA[slidingStart + binSize/2 ]  )
    }

    #mapping$chr = as.character(dna.range$seqnames)
    #mapping$strand = as.character(dna.range$strand)
    #rownames(mapping) = paste(geneName,slidingStart,sep = ",")
    #geneRNA2DNA= rbind(geneRNA2DNA,mapping[c("chr","start","end","strand")])

    #count reads in all samples
    ba.IP = sapply(bamPath.IP,.countReadFromBam,which = range(geneModel),reads.strand = reads.strand,DNA2RNA = DNA2RNA,fragmentLength=fragmentLength,left=dna.range$start,sliding = slidingStart, binSize = binSize, paired = paired)
    ba.input = sapply(bamPath.input,.countReadFromBam,which = range(geneModel),reads.strand = reads.strand,DNA2RNA = DNA2RNA,fragmentLength=fragmentLength,left=dna.range$start,sliding = slidingStart, binSize = binSize, paired = paired)

    if(is.vector(ba.IP) ){# if there is only one window for this gene, make it a matrix to avoid bug
      ba.IP = matrix(ba.IP, nrow = 1)
      ba.input = matrix( ba.input, nrow = 1 )
    }
    ba.counts <- cbind(ba.input,ba.IP)
    rownames(ba.counts) <-  paste(geneName,slidingStart,sep = ",")

    ba.counts
  }
  rm(list=ls(name=foreach:::.foreachGlobals), pos=foreach:::.foreachGlobals)
  end_time <- Sys.time()
  cat(paste("Time used to count reads:",difftime(end_time, start_time, units = "mins"),"mins... \n"))

  colnames(reads) <- c(paste(samplenames,"input",sep = "-"),paste(samplenames,"IP",sep = "-"))


  data.out <- MeRIP(reads = reads, binSize = binSize, gtf = gtf, geneModel = geneGRList, bamPath.input = bamPath.input, bamPath.ip = bamPath.IP, samplenames = samplenames)
  if(saveOutput){
    ## create output directory
    dir.create(outputDir, showWarnings = FALSE, recursive = TRUE)
    saveRDS(data.out,paste0(outputDir,"/MeRIP_readCounts.RDS"))
  }


  return(data.out)
}
