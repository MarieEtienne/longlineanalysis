args <- commandArgs(TRUE)

	file.in <- args[1]
	empty <- type.convert(args[2])
	
  print(file.in)
setwd("C:/Users/Marie/Rockfish/Longline")
source("SvnRCode/CLonglineData.R")
source("SvnRCode/DrawSamples.R")
source("SvnRCode/Utils.R")
source("SvnRCode/SEM.R")
source("SvnRCode/MEM.R")
source("SvnRCode/CPUE.R")
source("SvnRCode/SimuStudy.R")
source("SvnRCode/Output.R")

SimStudy(file.in, forget.empty=empty)