#script to combine sighting data from inshore and offshore files
#sets effort = the greatest amount of effort in either the inshore or the offshore
#sums number of mom-calf sightings from inshore + offshore

rm(list=ls())

library(tidyr)
library(dplyr)
library(tidyverse)

source("ExcelFileParams.R")

data.ext <- "v3" # "v2"  # or v2
data.path <- paste0("data\\Formatted Annual Data ", data.ext, "\\")

out.dir <- paste0("data\\Formatted Annual Data Combined ", data.ext, "\\")  # for revising outputs
if (!exists(out.dir)) dir.create(out.dir)

FILES <- list.files(path = data.path, ".csv")
years <- c(1994:2019, 2021:last.yr)

for (i in 1:length(years)) {
  year<-years[i]
  flist<-list.files(path = data.path, as.character(year))
  if (length(flist) > 1){
    off.dat<-read.csv(paste0(data.path, flist[grep("offshore", flist, 
            fixed = TRUE)]))
    in.dat<-read.csv(paste0(data.path, flist[grep("inshore", flist, 
                                                  fixed = TRUE)]))
    combined.dat<-in.dat %>%
      full_join(off.dat, by = c("Date", "Shift")) %>%
      mutate(Effort = pmax(Effort.x, Effort.y, na.rm = TRUE), 
             Sightings = rowSums(across(c(Sightings.x, Sightings.y)), na.rm = TRUE)) %>%
      select(-c(Effort.x, Sightings.x, Effort.y, Sightings.y, Week.y)) %>%
      rename(Week = Week.x) #use weeks from inshore file. #issues with 1995 file
    write.csv(combined.dat, paste0(out.dir, year, " ", 
              "Formatted_combined_", data.ext, ".csv"
                                   ), row.names = FALSE)
  } else {
    in.dat<-read.csv(paste0(data.path, flist[grep("inshore", flist, 
                                                  fixed = TRUE)]))
    write.csv(in.dat, paste0(out.dir, year, " ", 
                                   "Formatted_combined_", data.ext, ".csv"
    ), row.names = FALSE)
  }
  
}



