require(RColorBrewer)
par(mfrow=c(1,1))
par(mar=c(3,3,1,1))
display.brewer.all(3, colorblindFriendly=TRUE)
cols<-brewer.pal(3,"Set2")
cols.5<-c("#66C2A5", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F")
cols.4<-c("#66C2A5",  "#E78AC3", "#A6D854", "#FFD92F")
cols.3<-c("#66C2A5", "#FC8D62", "#8DA0CB")
cols.2<-brewer.pal(3,"RdBu")
cols.2<-c("#67A9CF", "#EF8A62")
par(mfrow=c(1, 1))
display.brewer.all(4, colorblindFriendly=TRUE)
cols1<-brewer.pal(4,"PuOr")
Cols1.2<-cols1[1:2]
Cols1.3<-cols1[c(1:3)]
dark_org<-"#E66101" 
light_org<-"#FDB863"
dark_pur<-"#5E3C99"
light_pur<-"#B2ABD2"
labels.17<-c("Control", "EWS")
labels.18.19<-c("Control", "AMWT", "SMWT")


# Figure 1. Species richness graphs
# Grouped with interaction
pdf(file="22_04_14_Species_richness.pdf", width = 9, height = 6)
par(mar=c(5,5.5,3,1))
par(mfrow=c(1,1))
par(xpd=F)
spp.richness.fig<-barplot(tapply(percent.cover.no.nas$richness, list(environmental.data$TREAT, environmental.data$YEAR), mean),  
                     beside=T, space = c(0, 0.2),
                     xlab="Year", ylab="Mean plant richness per quadrat (?SE)",
                     ylim=c(0, 7), cex.lab=1.2, cex.axis = 1,
                     col=ColTREAT)
abline(h=0) 
mean.rich<-(tapply(percent.cover.no.nas$richness, list(environmental.data$TREAT, environmental.data$YEAR), mean))
se.rich<-(tapply(percent.cover.no.nas$richness, list(environmental.data$TREAT, environmental.data$YEAR), sd)
             /sqrt(tapply(percent.cover.no.nas$richness, list(environmental.data$TREAT, environmental.data$YEAR), length)))
arrows(spp.richness.fig, mean.rich+se.rich, spp.richness.fig, mean.rich-se.rich, code=3, angle=90, length=0.05)

legend(cex=1, "topright",  horiz=F, title="Alleyway treatment", bty="n",
       legend=c("Control", "Standard management", "Active management"), 
       pch=c(15, 15), pt.cex=2, col=ColTREAT)
dev.off()

pdf(file="Figure_1.pdf", width = 9, height = 12)
par(mar=c(5,5.5,3,1))
par(mfrow=c(2,1))
par(xpd=F)
# TREATMENT:
spp.richness.treat<-barplot(tapply(percent.cover.no.nas$richness, environmental.data$TREAT, mean),  
                          xlab="Alleyway treatment", ylab="Mean species richness per quadrat (?SE)",
                          ylim=c(0, 5), cex.lab=1.2, cex.axis=1,   
                          names.arg=c("Control", "Standard \nmanagement", "Active \nmanagement"))
abline(h=0) 
mean.rich<-(tapply(percent.cover.no.nas$richness, environmental.data$TREAT, mean))
se.rich<-(tapply(percent.cover.no.nas$richness, environmental.data$TREAT, sd)
          /sqrt(tapply(percent.cover.no.nas$richness, environmental.data$TREAT, length)))
arrows(spp.richness.treat, mean.rich+se.rich, spp.richness.treat, mean.rich-se.rich, 
       code=3, angle=90, length=0.05)

spp.richness.year<-barplot(tapply(percent.cover.no.nas$richness, environmental.data$YEAR, mean),
                           xlab="Study year", ylab="Mean species richness per quadrat (?SE)",
                           ylim=c(0, 5), cex.lab=1.2, cex.axis = 1, 
                           names.arg=c("Year one", "Year two", "Year three"))
abline(h=0) 
mean.rich.year<-(tapply(percent.cover.no.nas$richness, environmental.data$YEAR, mean))
se.rich.year<-(tapply(percent.cover.no.nas$richness, environmental.data$YEAR, sd)
          /sqrt(tapply(percent.cover.no.nas$richness, environmental.data$YEAR, length)))
arrows(spp.richness.year, mean.rich.year+se.rich.year, spp.richness.year, mean.rich.year-se.rich.year, 
       code=3, angle=90, length=0.05)
dev.off()

# Figure 2. Species composition ordination
par(mfrow=c(1,1))
ColTREAT<-Cols[as.numeric(environmental.grouped$Treatment)]
ColTREAT[environmental.grouped$Treatment == "control"] = "#A6CEE3"
ColTREAT[environmental.grouped$Treatment == "smwt"] = "#1F78B4"
ColTREAT[environmental.grouped$Treatment == "amwt"] = "#B2DF8A"

# identify x 10 indicator species
x=lvsplot(Tot_Percent.ord.1, biplot=T, ind.spp=10, return.vals=TRUE, jitter=T)
# SF_A_millefolium, UF_M_periflora_PC, UF_V_sativa_PC,  UF_P_lagopus_PC, UF_P_aviculare_PC, UF_M_chamomilla_PC 
# UF_F_arvensis_PC, UF_G_robertianum_PC, UF_LactucaSpp_Pc, UF_P_oleracea_PC
x$scaled.lvs
scaled.lv.coefs_df<-as.data.frame(x$scaled.lv.coefs)
z=row.names(Tot_Percent.ord.1$geweke.diag$geweke.diag$lv.coefs)
row.names(scaled.lv.coefs_df)<-row.names(Tot_Percent.ord.1$geweke.diag$geweke.diag$lv.coefs)
# create labels
scaled.lv.coefs_INDISPP<-scaled.lv.coefs_df[c("SF_A_millefolium_PC", "UF_M_periflora_PC", "UF_V_sativa_PC", 
                                              "UF_P_lagopus_PC", "UF_P_aviculare_PC", "UF_M_chamomilla_PC", 
                                              "UF_F_arvensis_PC", "UF_G_robertianum_PC", "UF_LactucaSpp_PC", 
                                              "UF_P_oleracea_PC"),]
row.names(scaled.lv.coefs_INDISPP)<-c("A_millefolium", "M_periflora", "V_sativa", 
                                      "P_lagopus", "P_aviculare", "M_chamomilla", 
                                      "F_arvensis", "G_robertianum", "Lactuca_Sp", "P_oleracea")

pdf(file="Figure_2.pdf", width = 11, height = 11)
pchr = NULL
pchr[environmental.grouped$Treatment == "amwt"] = 18
pchr[environmental.grouped$Treatment == "smwt"] = 16
pchr[environmental.grouped$Treatment == "control"] = 17
par(mfrow=c(1, 1))
par(mar=c(11, 4.5, 4, 2))
par(xpd=F)
plot(x$scaled.lvs, col=ColTREAT, pch=pchr, cex=1.2, xlim=c(-3,4), 
     cex.main=1.5, cex.lab=1.4, cex.axis=1.2, 
     xlab="Latent variable 1", ylab="Latent variable 2", main="Biplot of latent variable posterior medians")
text(scaled.lv.coefs_INDISPP, labels=row.names(scaled.lv.coefs_INDISPP), cex=0.8, col="grey30")
par(xpd=TRUE)
legend(-1.2, -3.2, title = "Alleyway treatment", legend=c("Control", "Standard \nmanagement", "Active \nmanagement"), 
       col=c("#A6CEE3", "#1F78B4", "#B2DF8A"), pch=c(17, 16, 18), cex=1.2,  horiz=TRUE)
par(xpd=FALSE)
abline(v=0, lty=2, col="grey50")
abline(h=0, lty=2, col="grey50")
dev.off()

# Figure 3 - sown species performance 
library(lattice)
latt.cols<-brewer.pal(3, "PiYG")
ss = max( abs(coef(sown.species.max)) )
colort = colorRampPalette(latt.cols)
plot.sown_effect = levelplot(t(as.matrix(coef(sown.species.max))), ylab="",
                             xlab="", col.regions=colort(1000), at=seq(-ss, ss, length=1000),
                             scales = list( x= list(rot = 45)))
pdf(file="Figure_3.pdf", width = 9, height = 4)
print(plot.sown_effect)
dev.off()

# Figure 4. Sward height and COV
Height17box <- cov.df.17 %>%
  ggplot(aes(x=Treatment, y=mean.height)) + 
  geom_boxplot(fill=dark_org, how.legend=FALSE, notch=TRUE) +
  ylab("Mean vegetation hight per plot (mm)\n") +
  ylim(0, 300) +
  # geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=12), 
        axis.title=element_text(size=14,face="bold"), 
        plot.margin = unit(c(0,0,1,1), "cm"), 
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
CoV17box <- cov.df.17 %>%
  ggplot(aes(x=Treatment, y=CoV)) + 
  geom_boxplot(fill=light_org, how.legend=FALSE, notch=TRUE) +
  ylab("Coefficient of variation\n 
       of vegetation height per plot (%)\n") +
  xlab("\nAlleyway Treatment") +
  scale_x_discrete(labels=labels.17) + 
  ylim(0, 100) +
  # geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=12), 
        axis.title=element_text(size=14,face="bold"), 
        plot.margin = unit(c(0,0,0.3,1), "cm"))
Height18.19box <- cov.df.18.19 %>%
  ggplot(aes(x=Treatment, y=mean.height)) + 
  geom_boxplot(fill=dark_pur, how.legend=FALSE, notch=TRUE) +
  ylim(0, 300) +
  # geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=12), 
        axis.title=element_text(size=14,face="bold"), 
        plot.margin = unit(c(0,0,1,0.5), "cm"), 
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
        axis.title.y=element_blank(),
        axis.text.y=element_blank())
CoV18.19box <- cov.df.18.19 %>%
  ggplot(aes(x=Treatment, y=CoV)) + 
  geom_boxplot(fill=light_pur, how.legend=FALSE, notch=TRUE) +
  xlab("\nAlleyway Treatment") +
  scale_x_discrete(labels=labels.18.19) + 
  ylim(0, 100) +
  # geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=12), 
        axis.title=element_text(size=14,face="bold"), 
        plot.margin = unit(c(0,0,0.3,0.5), "cm"), 
        axis.title.y=element_blank(),
        axis.text.y=element_blank())
library(patchwork)
pdf(file="Figure_4.pdf", width = 9, height = 9)
Height17box + Height18.19box + CoV17box +  CoV18.19box
dev.off()


# Figure 5. Resource availability
resource.long.17<-dplyr::filter(resource.long, Year == "2017")
resource.long.17$Treatment<-gsub("amwt", "ews", resource.long.17$Treatment)
resource.long.17$Treatment<-gsub("smwt", "ews", resource.long.17$Treatment)
resource.long.17<-resource.long.17 %>% drop_na() # drop rows which total 0
resource.long.17$Treatment<-as.factor(resource.long.17$Treatment)
resource.long.18<-dplyr::filter(resource.long, Year == "2018")
resource.long.19<-dplyr::filter(resource.long, Year == "2019")

pdf(file="Figure_5.pdf", width = 10.5, height = 14.85)
# 2017
nf<-layout(matrix(c(1,2,3), 
                  ncol=1, byrow=TRUE), 
           heights=c(4,4,5.1))
par(mar=c(1, 5, 1.5, 1))
resource.long.17$Treatment<-relevel(resource.long.17$Treatment, ref = "control")
fig2017<-barplot(tapply(resource.long.17$cover, list(resource.long.17$Treatment, resource.long.17$cover_type), mean),  beside = T, 
                 xlab="", ylab="Plant resource index (%)", xaxt='n', 
                 ylim=c(0, 50), cex.lab=1.5, cex.axis = 1.5, cex.names = 1.5, cex.main = 2,
                 col=c("#A6DBA0", "#33A02C"), space=c(0,0.9))
par(xpd=F)
abline(h=0, lwd = 1.5) 
mean.2017<-(tapply(resource.long.17$cover, list(resource.long.17$Treatment, resource.long.17$cover_type), mean))
se.2017<-(tapply(resource.long.17$cover, list(resource.long.17$Treatment, resource.long.17$cover_type), sd)
          /sqrt(tapply(resource.long.17$cover, list(resource.long.17$Treatment, resource.long.17$cover_type), length)))
arrows(fig2017, mean.2017+se.2017, fig2017, mean.2017-se.2017, code=3, angle=90, length=0.05, lwd = 1.5)
par(xpd=T)

#2018
resource.long.18$Treatment<-relevel(resource.long.18$Treatment, ref = "control")
fig2018<-barplot(tapply(resource.long.18$cover, list(resource.long.18$Treatment, resource.long.18$cover_type), mean),  beside = T, 
                 xlab="", ylab="Plant resource index (%)", xaxt='n', 
                 ylim=c(0, 50), cex.lab=1.5, cex.axis = 1.5, cex.names = 1.5, cex.main = 2,
                 col=Cols_3,  space=c(0,0.3))
par(xpd=F)
abline(h=0, lwd = 1.5) 
mean.2018<-(tapply(resource.long.18$cover, list(resource.long.18$Treatment, resource.long.18$cover_type), mean))
se.2018<-(tapply(resource.long.18$cover, list(resource.long.18$Treatment, resource.long.18$cover_type), sd)
          /sqrt(tapply(resource.long.18$cover, list(resource.long.18$Treatment, resource.long.18$cover_type), length)))
arrows(fig2018, mean.2018+se.2018, fig2018, mean.2018-se.2018, code=3, angle=90, length=0.05, lwd = 1.5)
par(xpd=T)

# 2019
resource.long.18$TREAT<-relevel(resource.long.18$Treatment, ref = "control")
par(mar=c(10, 5, 1.5, 1))
fig2019<-barplot(tapply(resource.long.19$cover, list(resource.long.18$Treatment, resource.long.18$cover_type), mean),  beside = T, 
                 xlab="", ylab="Plant resource index (%)",   
                 ylim=c(0, 50), cex.lab=1.5, cex.axis = 1.5, cex.names = 1.5, cex.main = 2,
                 col=Cols_3, space=c(0,0.3))
par(xpd=F)
abline(h=0, lwd = 1.5) 
mean.2019<-(tapply(resource.long.19$cover, list(resource.long.18$Treatment, resource.long.18$cover_type), mean))
SE.2019<-(tapply(resource.long.19$cover, list(resource.long.18$Treatment, resource.long.18$cover_type), sd)
          /sqrt(tapply(resource.long.19$cover, list(resource.long.18$Treatment, resource.long.18$cover_type), length)))
arrows(fig2019, mean.2019+SE.2019, fig2019, mean.2019-SE.2019, code=3, angle=90, length=0.05, lwd=1.5)
par(xpd=T)
legend("bottom", inset=c(-0,-0.2), cex=1.5, horiz=T,  
       bty = "n", text.width=c(4.7, 4.7, 4.7, 4.7),
       legend=c("Control", "Active \nmanagement", "Standard \nmanagement", "Establishing \nwildflower strips"), 
       pch=c(15, 15), pt.cex=3, col=c("#A6DBA0", "#C2A5CF",  "#7B3294", "#33A02C"))
dev.off()
