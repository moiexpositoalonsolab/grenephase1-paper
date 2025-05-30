
# Broad European range
xlim = c(-10.5, +90)
ylim = c(25, +65)
Range=extent(c(xlim,ylim))
# Small range
xlim=c(-10.5,+ 53)
ylim=c(32,65)
Range=extent(c(xlim,ylim))


################################################################################
# Load climates

# now<-rbioclim::getData(name="worldclim2",var='bio', res=2.5,path = "data-external/")
# The present data is already downloaded from worldclim 2
now<-stack(list.files("data-external/wc2.1_2.5/", pattern = "\\.tif$", full.names = TRUE))
names(now)<-paste0("bio",1:19)
now<- now %>% crop(.,Range)


# fut<-getData(name="CMIP5",var='bio', res=5,model='MP', year=50, rcp=85)
# Same, future data is already downloaded from worldclim 2
fut<-stack("data-external/cmip6/2.5m/wc2.1_2.5m_bioc_MPI-ESM1-2-HR_ssp245_2021-2040.tif")
names(fut)<-paste0("bio",1:19)
fut<- fut %>% crop(.,Range)


latlon<-now[[1]]
(latlon)

