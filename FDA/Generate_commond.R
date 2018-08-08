
unlink('command_list')
for (i in 1:10){
	for(Type in c("M")){
		for(p in c(500)){
				for(n in c(150)){
							write(paste("R CMD BATCH --no-save --no-restore '--args case.id=", 
											i, " Type=\"",  Type, "\"", " p=", p, " n=", n,  "' FDA_study.R temp/GEV.output._", Type, 
											"_",  "p", p, "_n",n, "_id", i, 
											".txt  &", sep=""), file = "command_list", append = TRUE)
						}
					}
				}
}
