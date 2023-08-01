sig=read.csv("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/Validation_Selection/res.sig.validation.csv")
sig2=read.csv("/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/Validation_Selection/res.sig.validation.latest.csv")

sig2=sig2[which(sig2$Hit + sig2$New.Hit >= 1) , ]
table(sig2$Hit, sig2$New.Hit)
sig2=sig2[,-1]

sig2$X=sig$Validate[match(sig2$Pair, sig$Pair)]
write.csv(sig2, "/Volumes/share/mnt/Data0/PROJECTS/CROPSeq/FullScale/Results/Validation_Selection/res.sig.validation.combined.csv")
