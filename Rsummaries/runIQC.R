args = commandArgs(trailingOnly=T)

source('plotInstrumentQC.R')



print(args)


files<-list.files(args[1],pattern = '*.hdf5')

for (myfile in files){
  if (!grepl('result', myfile) & grepl('.hdf5', myfile)){
    pdfout<-gsub(".hdf5", ".pdf", myfile)
    runall(paste(args[1],myfile,sep='/'),paste(args[1], pdfout,sep='/'))


  }
}
