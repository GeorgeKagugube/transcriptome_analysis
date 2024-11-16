staroutput_preprocessing <- function(path_to_star_files, one_star_file){
  # Function takes in two inputs
  # 1.path to directory with the star out put files
  # 2. On of these files to read and continue to create/build the data matrix
  # I would define these outside of the function as variables that I supply to
  # the function
  ## Set the working directory/project directory here
  setwd(path_to_star_files)
  
  ## Loading data into R
  raw_data <- read.table(file = one_star_file,
                         sep = "\t",
                         header = F,
                         as.is = T)
  
  ## Gather all the data files in the directory of interest here
  data_files <- dir(pattern = "*ReadsPerGene.out.tab")
  
  # Create an empty vetor here
  countmatrix <- c()
  
  ## Loop through the files and extract the counts from star
  for (sample in seq_along(data_files)) {
    input_data <- read.table(file = data_files[sample], 
                             sep = "\t", 
                             header = F, 
                             as.is = T)
    
    ## Combine teh samples into a count table here 
    countmatrix <- cbind(countmatrix, input_data[,2])
  }
  
  ## Convert the countmatrix to a data frame here
  CountMatrix <- as.data.frame(countmatrix)
  
  ## Add rw names (ENSEMBL Gene IdS)
  rownames(CountMatrix) <- raw_data[,1]
  
  ## Rename the columns with the sample names 
  colnames(CountMatrix) <- sub("ReadsPerGene.out.tab","",data_files)
  
  CountMatrix <- CountMatrix[-c(1:4),]
  ## Return the count matrix here
  return(CountMatrix)
}

