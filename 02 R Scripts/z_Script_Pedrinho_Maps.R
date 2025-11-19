install.packages("readr")


library(raster)
library(rgdal)
library(maps)
library(readr)


library(readxl)
moderadores <- read_excel("~/2. Mestrado/Mapa/moderadores.xlsx")
View(moderadores)


#Carregando os rasters
bio1<-raster(choose.files())
bio2<-raster(choose.files())
bio3<-raster(choose.files())
bio4<-raster(choose.files())
bio5<-raster(choose.files())
bio6<-raster(choose.files())
bio7<-raster(choose.files())
bio8<-raster(choose.files())
bio9<-raster(choose.files())
bio10<-raster(choose.files())
bio11<-raster(choose.files())
bio12<-raster(choose.files())
bio13<-raster(choose.files())
bio14<-raster(choose.files())
bio15<-raster(choose.files())
bio16<-raster(choose.files())
bio17<-raster(choose.files())
bio18<-raster(choose.files())
bio19<-raster(choose.files())

wc_elev<-raster(choose.files())

annual_pet<-raster(choose.files())
aridity_index<-raster(choose.files())
climatic_moisture<-raster(choose.files())
continentality<-raster(choose.files())
emberguerq<-raster(choose.files())
degdays0<-raster(choose.files())
degdays5<-raster(choose.files())
cold<-raster(choose.files())
warm<-raster(choose.files())
temp10<-raster(choose.files())
pet_cq<-raster(choose.files())
pet_dq<-raster(choose.files())
pet_season<-raster(choose.files())
pet_warm<-raster(choose.files())
pet_wet<-raster(choose.files())
therm<-raster(choose.files())

tri<-raster(choose.files())
topo<-raster(choose.files())

#Definindo as localidades
localidades<-data.frame(lon=Moderadores[,4], lat=Moderadores[,3])
coordinates(localidades)<-c("lon","lat")
View(localidades)

###################### PLOT PARA TESTAR ##########################
plot(aridity_index, main="Thornthwaite Aridity Index", xlab="Longitude",ylab="Latitude",cex.axis=1.3, cex.lab=1.4, cex.main=1.5,col=rev(heat.colors(20)))
points(localidades, pch=16)
##################################################################

#Extraindo variáveis de cada raster X para as localidades Y, e passando valores para o dataframe 'Moderadores'
Moderadores$BIO1<-extract(x=bio1, y=localidades)
Moderadores$BIO2<-extract(x=bio2, y=localidades)
Moderadores$BIO3<-extract(x=bio3, y=localidades)
Moderadores$BIO4<-extract(x=bio4, y=localidades)
Moderadores$BIO5<-extract(x=bio5, y=localidades)
Moderadores$BIO6<-extract(x=bio6, y=localidades)
Moderadores$BIO7<-extract(x=bio7, y=localidades)
Moderadores$BIO8<-extract(x=bio8, y=localidades)
Moderadores$BIO9<-extract(x=bio9, y=localidades)
Moderadores$BIO10<-extract(x=bio10, y=localidades)
Moderadores$BIO11<-extract(x=bio11, y=localidades)
Moderadores$BIO12<-extract(x=bio12, y=localidades)
Moderadores$BIO13<-extract(x=bio13, y=localidades)
Moderadores$BIO14<-extract(x=bio14, y=localidades)
Moderadores$BIO15<-extract(x=bio15, y=localidades)
Moderadores$BIO16<-extract(x=bio16, y=localidades)
Moderadores$BIO17<-extract(x=bio17, y=localidades)
Moderadores$BIO18<-extract(x=bio18, y=localidades)
Moderadores$BIO19<-extract(x=bio19, y=localidades)

Moderadores$WC_Elev<-extract(x=wc_elev, y=localidades)

Moderadores$PET<-extract(x=annual_pet, y=localidades)
Moderadores$TAI<-extract(x=aridity_index, y=localidades)
Moderadores$CliMois<-extract(x=climatic_moisture, y=localidades)
Moderadores$Cont<-extract(x=continentality, y=localidades)
Moderadores$embQ<-extract(x=emberguerq, y=localidades)
Moderadores$GDD0<-extract(x=degdays0, y=localidades)
Moderadores$GDD5<-extract(x=degdays5, y=localidades)
Moderadores$TempCold<-extract(x=cold, y=localidades)
Moderadores$TempWarm<-extract(x=warm, y=localidades)
Moderadores$Temp10<-extract(x=temp10, y=localidades)
Moderadores$PETCold<-extract(x=pet_cq, y=localidades)
Moderadores$PETDry<-extract(x=pet_dq, y=localidades)
Moderadores$PETseason<-extract(x=pet_season, y=localidades)
Moderadores$PETWarm<-extract(x=pet_warm, y=localidades)
Moderadores$PETWet<-extract(x=pet_wet, y=localidades)
Moderadores$thermInd<-extract(x=therm, y=localidades)

Moderadores$TRI<-extract(x=tri, y=localidades)
Moderadores$TopoWet<-extract(x=topo, y=localidades)

View(Moderadores)

###salvar

write.table(Moderadores, file='Moderadores_WC_ENV.csv', sep=';', dec='.', row.names=FALSE)

Moderadores <- read_delim("C:/Users/ppsfb/OneDrive/Área de Trabalho/Doutorado - LEV/Meta-análise/Análises/Moderadores_WC_ENV.txt", 
                          "\t", escape_double = FALSE, trim_ws = TRUE)
View(Moderadores)

#######matrix de correlação

library(corrplot)

corr_matrix_WC <- Moderadores[,16:34]
corr_matrix_ENV <- Moderadores[,36:51]

matrix_WC <- na.omit(corr_matrix_WC)
matrix_ENV <- na.omit(corr_matrix_ENV)

sum(is.na(matrix_WC))
sum(is.na(matrix_ENV))

corrplot(cor(matrix_WC), method = "circle")
corrplot(cor(matrix_ENV), method = "circle")

corr_matrix_Elev <- Moderadores[,c(5, 6, 35, 52, 53)]
matrix_Elev <- na.omit(corr_matrix_Elev)
sum(is.na(matrix_Elev))
elev <- cor(matrix_Elev) # Corr matrix
round(elev, 2)
corrplot(cor(matrix_Elev), method = "number")

?na.omit

corr_matrix_selec <- Moderadores[,c(16,22, 27, 30, 35, 36, 37, 48, 51)]

matrix_selec <- na.omit(corr_matrix_selec)

sum(is.na(matrix_selec))

corrplot(cor(matrix_selec), method = "number")

###

corr_matrix_total <- Moderadores[,16:53]

matrix_total <- na.omit(corr_matrix_total)

sum(is.na(matrix_total))

corrplot(cor(matrix_total), method = "circle")

###########PCAs

library(vegan)
library(ggfortify)

?biplot

round(apply(matrix_selec,2,var),2) 
pca_selec <- prcomp(matrix_selec, scale=TRUE)
screeplot(pca_selec)
biplot(pca_selec)

round(apply(matrix_ENV,2,var),2) 
pca_ENV <- prcomp(matrix_ENV, scale=TRUE)
screeplot(pca_ENV)
biplot(pca_ENV)

corr_select <- Moderadores[,c(16,22, 27, 30, 35, 36, 37, 48)]
matrix_select <- na.omit(corr_select)
corrplot(cor(matrix_select), method = "number")
pca_select <- prcomp(matrix_select, scale=TRUE)
screeplot(pca_select)
biplot(pca_select)

autoplot(pca_select, loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 4)

autoplot(pca_ENV, loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 4)
