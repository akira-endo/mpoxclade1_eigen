## Install dependencies if needed ##############################################
if(!require (factoextra)) install.packages("factoextra")
if(!require (countrycode)) install.packages("countrycode")
if(!require (grid)) install.packages("grid")

## Standardisation #############################################################

load("../data/pca_indicators/dataset_final.RData")

means_df = data.frame(array(dim=c(nrow(dataset_final),ncol(dataset_final)-3)))
sd_df= data.frame(array(dim=c(nrow(dataset_final),ncol(dataset_final)-3)))


for (i in (1:nrow(dataset_final))){
  means_df[i,] <- colMeans(dataset_final[,4:ncol(dataset_final)])
  sd_df[i,] <- sapply(dataset_final[,4:ncol(dataset_final)],sd)
}

all_std_wide = dataset_final
all_std_wide[,4:ncol(dataset_final)] = (dataset_final[,4:ncol(dataset_final)]-means_df)/sd_df

## Conduct PCA #################################################################

library(factoextra)

pca <- prcomp(all_std_wide[,4:ncol(dataset_final)], scale = TRUE)
eigenvalue <- round(get_eigenvalue(pca), 1)
variance.percent <- eigenvalue$variance.percent
head(eigenvalue)

# Getting relative contribution of each indicator in the PCs
load <- with(pca, unclass(rotation))
aload <- abs(load)
relative_load <- sweep(aload, 2, colSums(aload), "/")
relative_load <- as.data.frame(relative_load[,1:4])
relative_load$var_name <- rownames(relative_load)

# Getting interpretations for each variable
mapping <- read.csv("../data/pca_indicators/variables_interpretation.csv",header = T)
relative_load_explained <- merge(relative_load, mapping,by="var_name", all.x=T)
write.csv(relative_load_explained, file="../data/pca_indicators/relative_load_explained.csv")

## Preparations for plotting ###################################################
dot_alpha = 0.3
# Get coordinates
ind.coord <- as.data.frame(get_pca_ind(pca)$coord)
rownames(ind.coord) <- all_std_wide$iso3c
individual_coord = ind.coord[,1:4]
write.csv(individual_coord,file="../data/pca_indicators/indiv_coordinates_final.csv")

# Get country-continent mapping
library(countrycode)
continent_mapping <- data.frame(iso=rownames(ind.coord), continent=NA)
continent_mapping$continent <- countrycode(sourcevar=continent_mapping$iso,
                                           origin="iso3c",
                                           destination = "continent")

## Plotting ####################################################################

library(grid)

masterplot <- function(pc_x, pc_y){
  pushViewport(plotViewport(c(3,3,1,1),yscale=c(round(min(ind.coord[,pc_y]))-3,
                                                round(max(ind.coord[,pc_y]))+3),
                            xscale=c(round(min(ind.coord[,pc_x]))-3,
                                     round(max(ind.coord[,pc_x]))+3)))#)NEW
  #grid.newpage()
  grid.rect()
  grid.yaxis(gp=gpar(fontsize=8))
  grid.xaxis(gp=gpar(fontsize=8))
  grid.points(ind.coord[continent_mapping$continent=="Europe",pc_x],
              ind.coord[continent_mapping$continent=="Europe",pc_y],
              pch=16,gp=gpar(cex = 0.6,lwd=1, col="blue",alpha=dot_alpha),default.units = 'native')
  grid.points(ind.coord[continent_mapping$continent=="Asia",pc_x],
              ind.coord[continent_mapping$continent=="Asia",pc_y],
              pch=16,gp=gpar(cex = 0.6,lwd=1, col="green",alpha=dot_alpha),default.units = 'native')
  grid.points(ind.coord[continent_mapping$continent=="Oceania",pc_x],
              ind.coord[continent_mapping$continent=="Oceania",pc_y],
              pch=16,gp=gpar(cex = 0.6,lwd=1, col="#00f0d4",alpha=dot_alpha),default.units = 'native')
  
  grid.points(ind.coord[(continent_mapping$continent=="Africa")&(!continent_mapping$iso %in% coi),pc_x],
              ind.coord[(continent_mapping$continent=="Africa")&(!continent_mapping$iso %in% coi),pc_y],
              pch=16,gp=gpar(cex = 0.6,lwd=1, col="red",alpha=dot_alpha),default.units = 'native')
  grid.points(ind.coord[(continent_mapping$continent=="Africa")&(continent_mapping$iso %in% c("BDI","COD")),pc_x],
              ind.coord[(continent_mapping$continent=="Africa")&(continent_mapping$iso %in% c("BDI","COD")),pc_y],
              pch=16,gp=gpar(cex = 0.6,lwd=1, col="red"),default.units = 'native')
  grid.points(ind.coord[continent_mapping$iso=="ZWE",pc_x],
              ind.coord[continent_mapping$iso=="ZWE",pc_y],
              pch=16,gp=gpar(cex = 0.6,lwd=1, col="red"),default.units = 'native')
  grid.points(ind.coord[continent_mapping$iso=="ZWE",pc_x],
              ind.coord[continent_mapping$iso=="ZWE",pc_y],
              pch=1,gp=gpar(cex = 0.6,lwd=1, col="black"),default.units = 'native')
  grid.points(ind.coord[continent_mapping$continent=="Americas",pc_x],
              ind.coord[continent_mapping$continent=="Americas",pc_y],
              pch=16,gp=gpar(cex = 0.6,lwd=1, col="#ff992b",alpha=dot_alpha),default.units = 'native')
  
  
  grid.text('ZWE',
            y = unit((ind.coord[continent_mapping$iso=='ZWE',pc_y]-
                        (round(min(ind.coord[,pc_y]))-3))/
                       ((round(max(ind.coord[,pc_y]))+3)-(round(min(ind.coord[,pc_y]))-3))-0.04,'npc'), 
            x = unit((ind.coord[continent_mapping$iso=='ZWE',pc_x]-
                        (round(min(ind.coord[,pc_x]))-3))/
                       ((round(max(ind.coord[,pc_x]))+3)-(round(min(ind.coord[,pc_x]))-3)),'npc'),
            gp=gpar(fontsize=6,fontface='bold'))
  
  if ((pc_x==4)&(pc_y==3)){
    grid.text('BDI',
              y = unit((ind.coord[continent_mapping$iso=='BDI',pc_y]-
                          (round(min(ind.coord[,pc_y]))-3))/
                         ((round(max(ind.coord[,pc_y]))+3)-(round(min(ind.coord[,pc_y]))-3))+0.04,'npc'), 
              x = unit((ind.coord[continent_mapping$iso=='BDI',pc_x]-
                          (round(min(ind.coord[,pc_x]))-3))/
                         ((round(max(ind.coord[,pc_x]))+3)-(round(min(ind.coord[,pc_x]))-3)),'npc'),
              gp=gpar(fontsize=6,fontface='bold'))
  }
  else{
    grid.text('BDI',
              y = unit((ind.coord[continent_mapping$iso=='BDI',pc_y]-
                          (round(min(ind.coord[,pc_y]))-3))/
                         ((round(max(ind.coord[,pc_y]))+3)-(round(min(ind.coord[,pc_y]))-3))-0.04,'npc'), 
              x = unit((ind.coord[continent_mapping$iso=='BDI',pc_x]-
                          (round(min(ind.coord[,pc_x]))-3))/
                         ((round(max(ind.coord[,pc_x]))+3)-(round(min(ind.coord[,pc_x]))-3)),'npc'),
              gp=gpar(fontsize=6,fontface='bold'))
  }
  
  if ((pc_x==2)&(pc_y==1)){
    grid.text('COD',
              y = unit((ind.coord[continent_mapping$iso=='COD',pc_y]-
                          (round(min(ind.coord[,pc_y]))-3))/
                         ((round(max(ind.coord[,pc_y]))+3)-(round(min(ind.coord[,pc_y]))-3))+0.04,'npc'), 
              x = unit((ind.coord[continent_mapping$iso=='COD',pc_x]-
                          (round(min(ind.coord[,pc_x]))-3))/
                         ((round(max(ind.coord[,pc_x]))+3)-(round(min(ind.coord[,pc_x]))-3)),'npc'),
              gp=gpar(fontsize=6,fontface='bold'))
  }
  
  else if ((pc_x==3)&(pc_y==1)){
    grid.text('COD',
              y = unit((ind.coord[continent_mapping$iso=='COD',pc_y]-
                          (round(min(ind.coord[,pc_y]))-3))/
                         ((round(max(ind.coord[,pc_y]))+3)-(round(min(ind.coord[,pc_y]))-3)),'npc'), 
              x = unit((ind.coord[continent_mapping$iso=='COD',pc_x]-
                          (round(min(ind.coord[,pc_x]))-3))/
                         ((round(max(ind.coord[,pc_x]))+3)-(round(min(ind.coord[,pc_x]))-3))+0.07,'npc'),
              gp=gpar(fontsize=6,fontface='bold'))
  }
  
  else{
    grid.text('COD',
              y = unit((ind.coord[continent_mapping$iso=='COD',pc_y]-
                          (round(min(ind.coord[,pc_y]))-3))/
                         ((round(max(ind.coord[,pc_y]))+3)-(round(min(ind.coord[,pc_y]))-3))-0.04,'npc'), 
              x = unit((ind.coord[continent_mapping$iso=='COD',pc_x]-
                          (round(min(ind.coord[,pc_x]))-3))/
                         ((round(max(ind.coord[,pc_x]))+3)-(round(min(ind.coord[,pc_x]))-3)),'npc'),
              gp=gpar(fontsize=6,fontface='bold'))
  }
  
  grid.text(paste0('PC',as.character(pc_x)),y = unit(0,'npc')+unit(-2.7,'lines'),gp = gpar(fontsize=9,fontface='bold',col = 'black'))
  grid.text(paste0('PC',as.character(pc_y)),x = unit(0,'npc')+unit(-2.7,'lines'),rot=90,gp = gpar(fontsize=9,fontface='bold',col = 'black'))
  popViewport()
}

paneller=function(row,column)
{
  x_pc = column + 1
  y_pc = row
  pushViewport(viewport(layout.pos.col=column,layout.pos.row=row))
  masterplot(x_pc,y_pc)
  popViewport()
}

pc_legend=function(row,column)
{
  pushViewport(viewport(layout.pos.col=column,layout.pos.row=row))
  
  grid.text('Variance explained by each PC',
            y = unit(0.8,'npc'))
  
  grid.text('PC1\n47.9%',y = unit(0.70,'npc'), x = unit(0.35,'npc'))
  grid.text('PC2\n12.4%',y = unit(0.70,'npc'), x = unit(0.65,'npc'))
  grid.text('PC3\n10.4%',y = unit(0.55,'npc'), x = unit(0.35,'npc'))
  grid.text('PC4\n5.4%',y = unit(0.55,'npc'), x = unit(0.65,'npc'))
  
  
  grid.polygon(x = c(0.11,0.19,0.19,0.11),
               y = c(0.3,0.3,0.33,0.33),
               gp=gpar(col='red',alpha=dot_alpha,fill='red'))
  
  grid.polygon(x = c(0.40,0.48,0.48,0.40),
               y = c(0.3,0.3,0.33,0.33),
               gp=gpar(col='#ff992b',alpha=dot_alpha,fill='#ff992b'))
  
  grid.polygon(x = c(0.75,0.83,0.83,0.75),
               y = c(0.3,0.3,0.33,0.33),
               gp=gpar(col='green',alpha=dot_alpha,fill='green'))
  
  grid.polygon(x = c(0.25,0.33,0.33,0.25),
               y = c(0.20,0.20,0.23,0.23),
               gp=gpar(col='blue',alpha=dot_alpha,fill='blue'))
  
  grid.polygon(x = c(0.55,0.63,0.63,0.55),
               y = c(0.20,0.20,0.23,0.23),
               gp=gpar(col='#00f0d4',alpha=dot_alpha,fill='#00f0d4'))
  
  grid.text('Africa',y = unit(0.32,'npc'), x = unit(0.29,'npc'))
  grid.text('Americas',y = unit(0.32,'npc'), x = unit(0.61,'npc'))
  grid.text('Asia',y = unit(0.32,'npc'), x = unit(0.91,'npc'))
  grid.text('Europe',y = unit(0.22,'npc'), x = unit(0.44,'npc'))
  grid.text('Oceania',y = unit(0.22,'npc'), x = unit(0.74,'npc'))
  
  popViewport()
}

# Countries of interest
coi = c("COD","BDI","ZWE")

# Final plot
svg('../figs/Sfigs/raw/pca.svg',height=21*0.394,width=21*0.394,pointsize=10)
pushViewport(plotViewport(c(1,1,1,1)))
pushViewport(viewport(layout=grid.layout(nrow=3,ncol=3)))
paneller(1,1)
paneller(1,2)
paneller(1,3)
paneller(2,2)
paneller(2,3)
paneller(3,3)
pc_legend(3,1)
popViewport()
popViewport()
dev.off()