
unlink('command_list')
for (i in 11:30){
	for(Type in c("I")){
		for(p in c(500)){
				for(n in c(2000)){
							write(paste("R CMD BATCH --no-save --no-restore '--args case.id=", 
											i, " Type=\"",  Type, "\"", " p=", p, " n=", n,  "' GEV_study.R temp/GEV.output._", Type, 
											"_",  "p", p, "_n",n, "_id", i, 
											".txt  &", sep=""), file = "command_list", append = TRUE)
						}
					}
				}
}
