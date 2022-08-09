
##TO DO:
# change the selection of high and low individuals by looping through the famlies within each site*provenance combination, and picking one highest and one lowest individual.

#CODES
# CD	crown density# CL	crown color# NS	north snc# NR	NORTH RHAB# SS	SOUTH SNC# SR	SOUTH RHAB# SRT	NEEDLE RETENTION# TCM 	TREE COMMENT


#Crown density was rated on a scale ranging from 1-4; a 1 corresponding to a sickly sparse crown lacking in needle retention, and a 4 corresponding to a full healthy crown. Crown color was rated on a scale from 1-3; a one corresponding to a green healthy crown color and a 3 corresponding to highly chlorotic, sickly crown color. Needle retention was rated on a secondary lateral branch located on the fourth whorl and on the south side of the tree. Needle retention was estimated as the proportion of needles retained in each year of growth e.g. 1.5 refers to a tree holding 100% of its current year needles, 50% of its second year needles and no third year needles.

#SNC severity was rated on both the north and south side of each tree on a secondary lateral branch on the fourth whorl from the top of the tree. Ratings ranged from 0 corresponding to no pseudothecia present on the underside of the needle, to 3 corresponding to greater that 66% of stomata occluded by pseudothecia. 

require(gdata)

fam <- read.xls ("~/Dropbox/desktop/adaptree/activity2_planning_gwas/Tag IDS study trees only 10-17-2017.xlsx")
fam2 <- fam[,1:2]

options(stringsAsFactors = F)
setwd("~/Dropbox/desktop/adaptree/activity2_planning_gwas/")
data <- read.table ("~/Dropbox/desktop/adaptree/activity2_planning_gwas/SNC_SSMT_9-14-15.csv",sep = ",", header = T)
data <- data[,1:12]
data <- data[data[,1] < 9434,] # avoid trees with TAGs less than 9434 because higher than that have no info.

data <- merge (data,fam2,by.x = "TAG",by.y = "TAG")
data <- data[is.na(data$FAMILY) == F,]


exclude1 <- c ("Evans","Stone","Slice")
data <- data[which(!(data$Site_name %in% exclude1)),]

exclude2 <- "CASIERRA"
data <- data[which(!(data$REGION %in% exclude2)),]



#sort first by crown density, and then by (colour + retention)
#susceptible + intolerant: Give precedence to high SNC and low crown density, then identify subs based on (high colour + low retention) 
#for susceptible + tolerant:  Give precedence to high SNC and high crown density, then identify subs based on (low colour + high retention)




data$SNC_BOTH <- (data$SS + data$NS)/2
data$R_BOTH <- (data$NR + data$SR) / 2

data <- data[which(data$R_BOTH <= 1),] 

data$index <- (4 - data$CL) + data$SRT
data$index2 <- (data$CD * 2) + (4 - data$CL) + data$SRT 

tapply (data$CD,list(data$Site_name,data$REGION),mean,na.rm = T)

data <- data[which (is.na(data$index2) == F),]


#pdf ("gwas_planning_index2_all_sites.pdf")
sites <- unique (data$Site_name)
regions <- unique (data$REGION)

results <- array (NA,c((length(sites)*length(regions)),19))
count <- 0


to_sample <- NULL
to_sample_resist <- NULL

bottoms_fam <- NULL
tops_fam <- NULL

for (ii in 1:length (sites)){
	
	subsite <- data[which(data$Site_name == sites[ii]),]
	
	
	for (jj in 1:length (regions)){
		count <- count + 1
		subregion <- subsite[which(subsite$REGION == regions[jj]),]
		
		results[count,1] <- sites[ii]
		results[count,2] <- regions[jj]
		results[count,3] <- sum (subregion$SNC_BOTH == 0,na.rm = T)
		results[count,4] <- sum (subregion$SNC_BOTH == 1 | subregion$SNC_BOTH == 0.5,na.rm = T)

	

		sus <- subregion[which (subregion$SNC_BOTH >= quantile(subregion$SNC_BOTH,0.25,na.rm = T)),]

		remain1 <- sus

		bottom_rank_fam <- NULL
		top_rank_fam <- NULL

		fam_rank_test <- F
		num_include <- 0
		
		while (fam_rank_test == F){
			
			the_biggest <- tapply (remain1$index2,remain1$FAMILY,max)
			the_smallest <- tapply (remain1$index2,remain1$FAMILY,min)
			
			the_range <- the_biggest - the_smallest
			
			the_nam <- names (the_range[the_range == max(the_range)])
			
			subfam <- remain1[remain1$FAMILY == the_nam[1],] #pick the first one, arbitrarily
			
			max_ind <- which (subfam$index2 == max (subfam$index2))
			min_ind <- which (subfam$index2 == min (subfam$index2))
			this_max <- subfam[max_ind[1],] #pick the first one, arbitrarily
			this_min <- subfam[min_ind[1],]
			
			this_val <- this_max$index2 - this_min$index2
			
			if (this_val > 2 && num_include < 10){ ##stop whenever there are a total of 10 entries for this site*provenance or when this_val falls is equal or below 2
				bottom_rank_fam <- rbind (bottom_rank_fam,this_min)
				top_rank_fam <- rbind (top_rank_fam,this_max)
			
				taken <- c(bottom_rank_fam$TAG,top_rank_fam$TAG)
			
				remain1 <- remain1[!(remain1$TAG %in% taken),]
				num_include <- num_include + 1
			} else {
				fam_rank_test <- T
			}
			

			
		}
			
		bottom_rank_fam$fam_diff <- top_rank_fam$index2 - bottom_rank_fam$index2
		top_rank_fam$fam_diff <- top_rank_fam$index2 - bottom_rank_fam$index2
		fam_diff <- top_rank_fam$index2 - bottom_rank_fam$index2
		
		bottom_rank_fam_rank <- bottom_rank_fam[order(fam_diff,decreasing = T),]
		top_rank_fam_rank <- top_rank_fam[order(fam_diff,decreasing = T),]
		
		bottom_rank_fam_rank$priority <- 1:nrow (bottom_rank_fam_rank)
		top_rank_fam_rank$priority <- 1:nrow (top_rank_fam_rank)
		
		bottoms_fam <- rbind (bottoms_fam, bottom_rank_fam_rank)
		tops_fam <- rbind (tops_fam, top_rank_fam_rank)


		# results[count,5] <- nrow (top_extra)
		# results[count,6] <- nrow (bottom_extra)
		# results[count,7] <- mean (top_rank$index2)
		# results[count,8] <- mean (bottom_rank$index2)
		# results[count,9] <- mean (top_rank$index2) - mean (bottom_rank$index2)
		# results[count,10] <- nrow(subregion)
		# results[count,11] <- mean (bottom_rank$SNC_BOTH)
		# results[count,12] <- mean (top_rank$SNC_BOTH)
		# results[count,13] <- quantile(subregion$SNC_BOTH,0.25,na.rm = T)
		# results[count,14] <- length (unique(top_3extra$FAMILY))
		# results[count,15] <- length (unique(bottom_3extra$FAMILY))
		# results[count,16] <- sum (bottom_3extra$FAMILY %in% top_3extra$FAMILY)
		# results[count,17] <- length (unique (sus_ord$FAMILY))

			
		

# # 		results[count,18] <- mean(top_rank_fam_rank$index2[1:8])
		# results[count,19] <- mean(bottom_rank_fam_rank$index2[1:8])
			

# # 		
		# bottoms <- cbind (bottom_3extra$Site_name, bottom_3extra$REGION, bottom_3extra$TAG,bottom_3extra$index2,array("intolerant",11),1:11,paste (bottom_3extra$REGION, bottom_3extra$Site_name,"I",sep = "_"))
		# tops <- cbind (top_3extra$Site_name, top_3extra$REGION, top_3extra$TAG,top_3extra$index2,array("tolerant",11),11:1,paste (top_3extra$REGION, top_3extra$Site_name,"T",sep = "_"))
		
		# to_sample <- rbind(to_sample,bottoms,tops)
				
		
		
		
		###RESISTANCE:
		
		r0 <- subregion[subregion$SNC_BOTH == 0,]
		r1 <- subregion[subregion$SNC_BOTH == 1,]
		
		row0 <- nrow (r0)
		row1 <- nrow (r1)

		if (nrow (r0) >= 8){
							
			if ((row0) >= 10){
				num_samp <- 10
			} else {
				num_samp <- row0 
			}
			
			samp1 <- r0[1:num_samp,]
			
			to_add <- cbind (samp1$Site_name,samp1$REGION,samp1$TAG,samp1$SNC_BOTH,1:num_samp,paste (samp1$REGION, samp1$Site_name,"R",sep = "_"),samp1$FAMILY)
			
			to_sample_resist <- rbind (to_sample_resist,to_add)
			
		} else if ((row0 + row1) >= 8 & row0 > 6){
						
			if ((row0 + row1) >= 10){
				num_samp <- 10
			} else {
				num_samp <- row0 + row1
			}

			samp1 <- rbind (r0, r1[1:(num_samp-row0),])
			
			to_add <- cbind (samp1$Site_name,samp1$REGION,samp1$TAG,samp1$SNC_BOTH,1:num_samp,paste (samp1$REGION, samp1$Site_name,"R",sep = "_"),samp1$FAMILY)
			
			to_sample_resist <- rbind (to_sample_resist,to_add)

		} 
		
		
		
		
		
	}
}
#dev.off()

bottoms_to_sample <- cbind (bottoms_fam$Site_name, bottoms_fam$REGION, bottoms_fam$TAG,bottoms_fam$FAMILY, bottoms_fam$index2,array("intolerant",nrow (bottoms_fam)),paste (bottoms_fam$REGION, bottoms_fam$Site_name,"I",sep = "_"),bottoms_fam$priority)
tops_to_sample <- cbind (tops_fam$Site_name, tops_fam$REGION, tops_fam$TAG, tops_fam$FAMILY, tops_fam$index2,array("tolerant",nrow (tops_fam)),paste (tops_fam$REGION, tops_fam$Site_name,"T",sep = "_"),tops_fam$priority)

to_sample <- rbind (bottoms_to_sample,tops_to_sample)


##PLOTTING:
##PLOTTING:
##PLOTTING:
plot (0,xlim = c (0.5,2.5),ylim = c (0,14),col = "white",xaxt = "n",xlab = "",ylab = "case-control index")
axis (1,c(1,2),c("ranked","family_ranked"))

for (i in 1:nrow(results)){
	points (c(1,2),c(as.numeric (results[i,7]),as.numeric (results[i,18])), type = "o", col = "red")
	points (c(1,2),c(as.numeric (results[i,8]),as.numeric (results[i,19])), type = "o", col = "blue")
	
}


bottoms_fam$code1 <- paste (bottoms_fam$REGION,bottoms_fam$Site_name,sep = "__")

par (mar = c (10,5,5,5))
stripchart (bottoms_fam$fam_diff ~ bottoms_fam$code1, vertical = T, ylab = "difference between top and bottom ranked individual in family", pch = 1)
par (las = 3)
#axis (1,1:length (unique(bottoms_fam$code1)),labels = unique (bottoms_fam$code1))
arrows (-1000,2,1000,2,col = "red")
arrows (-1000,8,1000,8,col = "blue")

##DONE PLOTTING
##DONE PLOTTING
##DONE PLOTTING



colnames (to_sample) <- c("Site_name","REGION","TAG","FAMILY","index","pool_type","pool_name","priority")

write.table (to_sample, "to_sample_tolerant_vs_intolerant_withFAMILY.csv", sep = ",",col.names = T, row.names = F, quote = F)




##SNC resistance
##SNC resistance
##SNC resistance
##SNC resistance
##SNC resistance
##SNC resistance
##SNC resistance

colnames (to_sample_resist) <- c("Site_name","REGION","TAG","SNC_score","priority","pool_name","family")


#check out SNC resistance (0's and 1's)
dd <- tapply ((as.numeric (results[,3]) + as.numeric (results[,4])),list(results[,1],results[,2]),sum)
ee <- tapply ((as.numeric (results[,10])),list(results[,1],results[,2]),sum)


#check out SNC resistance (0's only)
dd0 <- tapply ((as.numeric (results[,3])),list(results[,1],results[,2]),sum)
dd1 <- tapply ((as.numeric (results[,4])),list(results[,1],results[,2]),sum)
ee <- tapply ((as.numeric (results[,10])),list(results[,1],results[,2]),sum)


dd0[,!(colnames (dd0) %in% to_sample_resist[,2])]
dd1[,!(colnames (dd1) %in% to_sample_resist[,2])]


samp1 <- data[data$SNC_BOTH == 0 & data$REGION == "ORCASL",]
to_add <- cbind (samp1$Site_name,samp1$REGION,samp1$TAG,samp1$SNC_BOTH,1:nrow(samp1),paste (samp1$REGION,"combined","R",sep = "_"),samp1$FAMILY)			
to_sample_resist <- rbind (to_sample_resist,to_add)


samp1 <- data[data$SNC_BOTH == 0 & data$REGION == "ORSISL",]
to_add <- cbind (samp1$Site_name,samp1$REGION,samp1$TAG,samp1$SNC_BOTH,1:nrow(samp1),paste (samp1$REGION,"combined","R",sep = "_"),samp1$FAMILY)			
to_sample_resist <- rbind (to_sample_resist,to_add)


int1 <- data[data$SNC_BOTH == 0 & data$REGION == "WACASH",]
int2 <- data[data$SNC_BOTH == 1 & data$REGION == "WACASH",]
samp1 <- rbind (int1,int2[c(1,3),])
to_add <- cbind (samp1$Site_name,samp1$REGION,samp1$TAG,samp1$SNC_BOTH,1:nrow(samp1),paste (samp1$REGION,"combined","R",sep = "_"),samp1$FAMILY)			
to_sample_resist <- rbind (to_sample_resist,to_add)


##check family representation:
ind111 <- paste (to_sample_resist[,6],to_sample_resist[,7],sep = "__")
tab1 <- table (ind111)  #how many families are over-represented?
tab1[tab1 >= 3]

#manually remove samples from tab1 to bring down to 2:
ind1 <- which (ind111 == "ORSISH_Nortons_R__1171")
ind2 <- which (ind111 == "ORSISL_combined_R__3094")

to_sample_resist_good <- to_sample_resist[-c(15,86,87),]


##re-check family representation:
ind111 <- paste (to_sample_resist_good[,6], to_sample_resist_good[,7],sep = "__")
tab1 <- table (ind111)  #how many families are over-represented?
tab1[tab1 >= 3]

write.table (to_sample_resist_good, "to_sample_resistant_withFAMILY.csv", sep = ",",col.names = T, row.names = F, quote = F)




##Rhabdocline planning
##Rhabdocline planning
##Rhabdocline planning
##Rhabdocline planning
##Rhabdocline planning
##Rhabdocline planning
##Rhabdocline planning
##Rhabdocline planning
##Rhabdocline planning


options(stringsAsFactors = F)
setwd("~/Dropbox/desktop/adaptree/activity2_planning_gwas/")
data <- read.table ("~/Dropbox/desktop/adaptree/activity2_planning_gwas/SNC_SSMT_9-14-15.csv",sep = ",", header = T)
data <- data[,1:12]
data <- data[data[,1] < 9434,] # avoid trees with TAGs less than 9434 because higher than that have no info.

data <- merge (data,fam2,by.x = "TAG",by.y = "TAG")
data <- data[is.na(data$FAMILY) == F,]

#only focus on CASIERRA, and look at all sites except Nortons, Evans, and Stone, where there is very little rhabdocline.

exclude1 <- c ("Evans","Stone","Nortons")
data <- data[which(!(data$Site_name %in% exclude1)),]


include1 <- c ("CASIERRA")
data <- data[which((data$REGION %in% include1)),]


data$SNC_BOTH <- (data$SS + data$NS)/2
data$R_BOTH <- (data$NR + data$SR) / 2

data <- data[is.na(data$R_BOTH) == F,]


to_sample_rhab <- NULL

#pdf ("gwas_planning_rhabodcline.pdf")
sites <- unique (data$Site_name)

count <- 0

for (ii in 1:length (sites)){
	
	subsite <- data[which(data$Site_name == sites[ii]),]
	
	fams <- unique (subsite$FAMILY)
	
	for (jj in 1:length (fams)){
		
		subfam <- subsite[subsite$FAMILY == fams[jj],]
		
		
		if (min (subfam$R_BOTH) <= 1){
			
			subfam_ord <- subfam[order(subfam$R_BOTH),]
			
		
			samp1 <- subfam_ord[1,]		
			to_add <- cbind (samp1$Site_name,samp1$REGION,samp1$TAG,samp1$R_BOTH,array("resistant",nrow(samp1)),1:nrow(samp1),paste (samp1$REGION, samp1$Site_name,"rhabR",sep = "_"),samp1$FAMILY)
			to_sample_rhab <- rbind (to_sample_rhab,to_add)
	
	
	
			samp1 <- subfam_ord[nrow(subfam_ord),]
			to_add <- cbind (samp1$Site_name,samp1$REGION,samp1$TAG,samp1$R_BOTH,array("susceptible",nrow(samp1)),1:nrow(samp1),paste (samp1$REGION, samp1$Site_name,"rhabS",sep = "_"),samp1$FAMILY)	
			to_sample_rhab <- rbind (to_sample_rhab,to_add)

		
			
			
		}
		
		
	}
	

	
}




colnames (to_sample_rhab) <- c("Site_name","REGION","TAG","Rhab_infection","pool_type","priority","pool_name","family")

write.table (to_sample_rhab, "to_sample_rhabdocline_susceptible_resistant_withFAMILY.csv", sep = ",",col.names = T, row.names = F, quote = F)

int1 <- which (to_sample_rhab[,4] == "0")
to_samp_zero_resist <- to_sample_rhab[int1,]
to_samp_zero_sus <- to_sample_rhab[(int1+1),]

to_samp_zero <- rbind (to_samp_zero_resist, to_samp_zero_sus)
to_samp_zero <- to_samp_zero[order(to_samp_zero[,8],to_samp_zero[,1],to_samp_zero[,2]),]

write.table (to_samp_zero, "to_sample_rhabdocline_susceptible_resistant_zeros_withFAMILY.csv", sep = ",",col.names = T, row.names = F, quote = F)











