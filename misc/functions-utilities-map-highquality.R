


coolmaptile<-function(mapbase, fillcolor=(brewer.pal(name='Greys',7)[-1]), fillname="Habitat suitability", ...){
# Convert a map
# https://stackoverflow.com/questions/33227182/how-to-set-use-ggplot2-to-map-a-raster
# baseraste<-makebaseraster()
# Euroclim<-make_Euroclim()
# example<-Euroclim[["gs_sum"]]
example<-mapbase
test_spdf <- as(example, "SpatialPixelsDataFrame")
test_df <- as.data.frame(test_spdf)
colnames(test_df) <- c("value", "x", "y")

maptileplot<-ggplot() +
  geom_tile(data=test_df, aes(x=x, y=y, fill=value), ...)+
  # scale_color_gradientn("Polygenic score" , colours=brewer.pal(name='BrBG',5))+
  scale_fill_gradientn(fillname,  colours=fillcolor)+
  coord_equal() +
  theme_map()+  theme(legend.position = "bottom")
return(maptileplot)
}

# If you want to add more points you can do something like this
  # geom_point(aes(y=fn(degreeadapt[,'latitude']), x=fn(degreeadapt[,'longitude'])),color="black",size=1.2)+
  # geom_point(aes(y=fn(degreeadapt[,'latitude']), x=fn(degreeadapt[,'longitude'])),color="white",size=1.1)+
  # geom_point(aes(y=fn(degreeadapt[,'latitude']), x=fn(degreeadapt[,'longitude']),color=degreeadapt[,"polyscore"]),size=1)+
  # scale_color_gradientn("Polygenic score" , colours=rev(mypalettes('moispectral') ))+