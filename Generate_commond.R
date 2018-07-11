
unlink('command_list')
for (i in 1:3){
	for(Type in c("Toeplitz")){
		for(p in c(50)){
				for(n in c(300)){
							write(paste("bsub -q day -o  temp/temp", i 
											," -M 4 R CMD BATCH --no-save --no-restore '--args case.id=", 
											i, " Type=\"",  Type, " p=", p, " n=", n,  "' GEV_study.R GEV.output._", Type, 
											"_",  "p", p, "_n",n, "_id", i, 
											".txt", sep=""), file = "command_list", append = TRUE)
						}
					}
				}
			}
