# create script

ReactomePathwaysRelation <- read.table("data-raw/ReactomePathwaysRelation.txt", sep="\t", header=F, quote="\"", stringsAsFactors = F, check.names = F)
colnames(df) <- c("parent", "child")

usethis::use_data(ReactomePathwaysRelation)
