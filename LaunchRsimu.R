args <- commandArgs(TRUE)

	file.in <- args[1]
	empty <- type.convert(args[2])
	
  print(file.in)
setwd("C:/Users/Marie/Rockfish/Longline")
source("RCode/CLonglineData.R")
source("RCode/DrawSamples.R")
source("RCode/Utils.R")
#source("RCode/LongLineJags.R")
source("RCode/SEM.R")
source("RCode/MEM.R")
source("RCode/CPUE.R")
source("RCode/SimuStudy.R")
source("RCode/SweptArea.R")
source("RCode/Output.R")

SimStudy(file.in, forget.empty=empty)