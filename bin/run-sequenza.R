library(sequenza)
library(stringr)

library(sequenza)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)

sample_name<-as.character(args[2])

job_name<-str_replace_all(sample_name,"-","_")
print(job_name)

chromosome.list.male = paste0("chr",c(1:22,"X","Y"))
chromosome.list.female = paste0("chr",c(1:22,"X"))


if(as.character(args[3])=="male")
{
	job_name<-sequenza.extract(args[1],normalization.method="median",chromosome.list=chromosome.list.male,kmin=500,gamma=100)
	print("Enters Male")
} else {
	job_name<-sequenza.extract(args[1],normalization.method="median",chromosome.list=chromosome.list.female,kmin=500,gamma=100)
}

if(as.character(args[3])=="male")
{
	job_name.CP<-sequenza.fit(sequenza.extract=job_name,ratio.priority = FALSE,female=FALSE,chromosome.list=chromosome.list.male,XY = c(X = "chrX", Y = "chrY"))
} else {
	job_name.CP<-sequenza.fit(sequenza.extract=job_name,ratio.priority = FALSE,female=TRUE,chromosome.list=chromosome.list.female,XY = c(X = "chrX", Y = "chrY"))
}

out_dir<-paste(sample_name,"_OUTPUT",sep="")
print(out_dir)

if(as.character(args[3])=="male")
{
	sequenza.results(sequenza.extract = job_name, cp.table = job_name.CP,sample.id =sample_name, out.dir=out_dir, ratio.priority=FALSE,female=FALSE,chromosome.list=chromosome.list.male,XY = c(X = "chrX", Y = "chrY"))
} else {
	sequenza.results(sequenza.extract = job_name, cp.table = job_name.CP,sample.id =sample_name, out.dir=out_dir, ratio.priority=FALSE,female=TRUE,chromosome.list=chromosome.list.female,XY = c(X = "chrX", Y = "chrY"))
}

warnings()
