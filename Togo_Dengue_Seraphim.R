# install.packages("devtools"); library(devtools)
# install_github("sdellicour/seraphim/unix_OS") # for Unix systems
## install_github("sdellicour/seraphim/windows") # for Windows systems
# install.packages("diagram")

library(diagram)
require(seraphim)

loadstuff<- TRUE
if(loadstuff){

      # 1. Extracting the spatio-temporal information contained in posterior trees

      treefile<- "data/beast_files//Big_Togo_A.2_Phylogeo_1B.e01_postBurnIn.trees"

      localTreesDirectory = "Tree_extractions_Big_Togo_A.2"
      allTrees = scan(file=treefile, what="", sep="\n", quiet=T)
      burnIn = 0 
      randomSampling = FALSE
      nberOfTreesToSample = 100
       mostRecentSamplingDatum = 2024.738 
      
      coordinateAttributeName = "location"

      treeExtractions(localTreesDirectory, allTrees, burnIn, randomSampling, nberOfTreesToSample, mostRecentSamplingDatum, coordinateAttributeName)


      # 2. Extracting the spatio-temporal information embedded in the MCC tree

      treefile<- "data/beast_files/Big_Togo_A.2_Phylogeo_1B.e01_MCC.tree"

      source("mccExtractions.r")
      mcc_tre = readAnnotatedNexus(treefile)
      mcc_tab = mccTreeExtraction(mcc_tre, mostRecentSamplingDatum)
      write.csv(mcc_tab, "data/Togo_cluster_MCC.csv", row.names=F, quote=F)

      # 3. Estimating the HPD region for each time slice
      mcc_tab=read.csv("data/Togo_cluster_MCC.csv")
      nberOfExtractionFiles = nberOfTreesToSample
      prob = 0.95; precision = 0.025
      startDatum = min(mcc_tab[,"startYear"])

      polygons = suppressWarnings(spreadGraphic2(localTreesDirectory, nberOfExtractionFiles, prob, startDatum, precision))


      # 4.1 spatial boundaries

      reloadborders<-TRUE
      if(reloadborders){
          template_raster = raster("data/togo_studyarea_population_density.asc")
          borders = geodata::gadm("GADM", country="TGO", level=1)
          
      }

}

# 4.2 Defining the different colour scales

minYear = min(mcc_tab[,"startYear"]); maxYear = max(mcc_tab[,"endYear"])
endYears_indices = (((mcc_tab[,"endYear"]-minYear)/(maxYear-minYear))*100)+1

##colors have to go to max time slice - AT LEAST
n_number_colours_needed<- max(round(endYears_indices))
n_repeats_discrete<- 10
c2<- rev((brewer.pal(9,"RdBu")))
c3<-brewer.pal(9,"YlOrRd")

colours<- rev(rep(c(c2), each=n_repeats_discrete))

colour_scale<- colorRampPalette(colours)(n_number_colours_needed)

endYears_colours = colour_scale[round(endYears_indices)]
polygons_colours = rep(NA, length(polygons))
for (i in 1:length(polygons))
{
  date = as.numeric(names(polygons[[i]]))
  polygon_index = round((((date-minYear)/(maxYear-minYear))*100)+1)
  polygons_colours[i] = paste0(colour_scale[polygon_index],"15")
}

# 5. Generating the dispersal history plot

pdf('figures/Fig2D_Togo_Cluster_map_Big_A.2_clade_phylogeo.pdf',width=6, height=6.3,bg="white")

ptsize<- 0.7
pitjit<- 0.2
par(mar=c(0,0,0,0), oma=c(1.2,3.5,1,0), mgp=c(0,0.4,0), lwd=0.2, bty="o")
plot(template_raster,col = brewer.pal(9,"YlOrRd"), colNA="white",box=F, axes=F, legend=F)

plot(borders, add=T, lwd=0.7, border="black")

for (i in length(polygons):1)
{
  plot(polygons[[i]], axes=F, col=polygons_colours[i], add=T, border=NA)
}
for (i in 1:dim(mcc_tab)[1])
{
  curvedarrow(cbind(mcc_tab[i,"startLon"],mcc_tab[i,"startLat"]), cbind(mcc_tab[i,"endLon"],mcc_tab[i,"endLat"]), arr.length=0,
              arr.width=0, lwd=3.5*1.1, lty=1, lcol="grey22", arr.col=NA, arr.pos=FALSE, curve=0.4, dr=NA, endhead=F)
  curvedarrow(cbind(mcc_tab[i,"startLon"],mcc_tab[i,"startLat"]), cbind(mcc_tab[i,"endLon"],mcc_tab[i,"endLat"]), arr.length=0,
              arr.width=0, lwd=3*1.1, lty=1, lcol=endYears_colours[i], arr.col=NA, arr.pos=FALSE, curve=0.4, dr=NA, endhead=F)

}
for (i in dim(mcc_tab)[1]:1)
{
  xs<- mcc_tab[i,"startLon"]
  ys<- mcc_tab[i,"startLat"]
  xe<- jitter(mcc_tab[i,"endLon"],pitjit)
  ye<- jitter(mcc_tab[i,"endLat"],pitjit)
  if (i == 1)
  {
    points(xs, ys, pch=16, col=colour_scale[1], cex=ptsize)
    points(xs, ys, pch=1, col="gray10", cex=ptsize)
  }
  points(xe, ye, pch=16, col=endYears_colours[i], cex=ptsize)
  points(xe, ye, pch=1, col="gray10", cex=ptsize)
}

xrange<- c(xmin(template_raster), xmax(template_raster))
yrange<- c(ymin(template_raster), ymax(template_raster))
rast = raster(matrix(nrow=1, ncol=2)); rast[1] = min(mcc_tab[,"startYear"]); rast[2] = max(mcc_tab[,"endYear"])
rast1 = raster(matrix(nrow=1, ncol=2)); rast1[1] = 84; rast1[2] = 223

plot(rast, legend.only=T, add=T, col=colour_scale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.50,0.80,0.84,0.855),
     legend.args=list(text="Date of Viral Movement", cex=0.7, line=0.3, col="gray30"), horizontal=T,
     axis.args=list(cex.axis=0.6, lwd=0, lwd.tick=0.2, tck=-0.5, col.axis="gray30", line=0, mgp=c(0,-0.02,0)))
plot(rast1, legend.only=T, add=T, col=c3, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.60,0.80,0.74,0.755),
     legend.args=list(text="Population Density", cex=0.7, line=0.3, col="gray30"), horizontal=T,
     axis.args=list(cex.axis=0.6, lwd=0, lwd.tick=0.2, tck=-0.5, col.axis="gray30", line=0, mgp=c(0,-0.02,0)))

a<-dev.off()

