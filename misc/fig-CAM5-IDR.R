

iso1<-read.csv("data-external/cam5-iso1-IDRprediction.csv",header = F,sep = " ")
iso3<-read.csv("data-external/cam5-iso3-IDRprediction.csv",header = F,sep = " ")

# iso1 %>% head
# iso1$V1<-iso1$V1 %>% as.numeric()
# iso1$aa<-1:nrow(iso1) %>% as.numeric
dim(iso1)
dim(iso3)

figidr<-
ggplot()+
  geom_line(data=iso1, aes(x=V1,y=V2),color="grey", lwd=2)+
  geom_line(data=iso3, aes(x=V1,y=V2),color="black", lwd=2)+
  scale_y_continuous(limits = c(0,1))+
  labs(x="aa residue position",
       y="IDR prediction"
       )

save_plot("figs/fig-CAM5-IDR.pdf",
          figidr,
          base_height = 4, base_width=7
          )
save_plot("figs/fig-CAM5-IDR.png",
          figidr,
          base_height = 4, base_width=7
)
