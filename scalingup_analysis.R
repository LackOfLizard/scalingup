
#Libraries

library(xlsx)
library(stats)
library(rJava)
library(ape)
library(phytools)
library(png)
library(TreeTools)
library(geiger)
library(rlang)
library(pillar)
library(FactoMineR)
library(factoextra)
library(scatterplot3d)
library(magick)
library(devtools)
library(credentials)
library(rgl)
library(Claddis)



#Read in dataset from excel
scale <- as.matrix(read.xlsx('Dataset.xlsx',
                             sheetIndex = 1,
                             rowIndex = 1:153,
                             colIndex = 2:181,
                             header = F))
#Return taxa names
scale_taxa <- read.xlsx('Dataset.xlsx',
                        sheetIndex = 1,
                        rowIndex = 1:153,
                        colIndex = 1,
                        as.data.frame = F)
#Develop vector for denoting proper character numbers
ch <- c(1:180)
#Apply names to each section of the matrix
rownames(scale) <- scale_taxa
colnames(scale) <- ch

#Verify that taxa are identically named
name.check(pyron, scale)
#Underscores instead of spaces in our tree, data needs to be changed
rownames(scale) <- gsub(" ","_",rownames(scale))
#Now to verify
name.check(pyron, scale)

#Set up the data to use in MCA
scale.df <- as.data.frame(scale)
#Naming variables will help with interpretation later on
var_names <- (c("Rostral", "Ex.Rostral", "Spl.Rostral", "Narial.Op", "Op.Rostral", "Op.Superlab",
                "Nasal.Ros", "Supralab.Ros", "Internasal", "Ceph.Size", "Ceph.Keel", "Ceph.Elong",
                "En.Frontal", "FParietal.Con", "CirOrbital", "En.Supraoc", "Supraoc.Row", "Supracil",
                "Supraoc.Con", "Supracil.Ov", "Supralab.Orb", "En.Suboc", "Suboc.Con", "En.Interpar",
                "En.Temp", "Heat.Pit", "Mouth", "Tongue.Gr", "Eye", "Eye.Scales",
                "Eyelids", "Low.Eye.Sc", "Spectacle", "Gran.Eyelid", "Ear", "Tympanum",
                "Ear.Cover", "Dewlap", "Mental", "El.Postment", "One.Postment", "En.Med.Infralab",
                "Infralab.Sep", "Ex.Labials", "Mental.Gr", "Gular.Fold", "Annul.Body", "Varanus.Sc",
                "Homo.Dorsal", "Gen.Dorsal", "Body.Keel", "Dorsal.Ridge", "En.Middorsal", "Nuchal.Cr",
                "Gen.Lateral", "Lateral.Fold", "Gastrosteges", "Ventral.Size", "Gen.Ventral", "Belly.Keel",
                "Ventral.Edge", "Fem.Pores", "Preanal.Pores", "Tail.Tip", "Tail.Spike", "Tail.Cross",
                "Tail.Size", "Tail.Keel", "Annul.Tail", "En.Preanal", "Anal.Plate", "En.Postanal",
                "Tubercles", "One.Subcaudal", "Ex.Hindlimb", "Limb.Cross", "Toe.Pads", "Subdig.Num",
                "Subdig,Keel", "Mittens", "Claw.Sheath", "Frontal.Shape", "Frontal.Tex", "Frontal.Ov",
                "Chalcidine", "Inasal.Shape", "Inasal.Tex", "Inasal.Keel", "Inasal.K.T", "Loreal.Shape",
                "Loreal.Tex", "Loreal.Keel", "Loreal.K.T", "Loreal.Just", "Loreal.Ov", "En.Loreal",
                "Supraoc.Shape", "Supraoc.Tex", "Supraoc.Keel", "Supraoc.K.T", "Supraoc.Just", "Supraoc.Ov",
                "Postmen.Shape", "Postmen.Tex", "Postmen.Keel", "Postmen.K.T", "Postmen.Ov", "Homo.DNeck",
                "DNeck.Shape", "DNeck.Tex", "DNeck.Keel", "DNeck.K.T", "DNeck.Just", "DNeck.Ov",
                "LNeck.Shape", "LNeck.Tex", "LNeck.Keel", "LNeck.K.T", "LNeck.Just", "LNeck.Ov",
                "VNeck.Shape", "VNeck.Tex", "VNeck.Keel", "VNeck.K.T", "VNeck.Just", "VNeck.Ov",
                "Dorsal.Shape", "Dorsal.Tex", "Dorsal.Keel", "Dorsal.K.T", "Dorsal.Just", "Dorsal.Ov",
                "Lateral.Shape", "Lateral.Tex", "Lateral.Keel", "Lateral.K.T", "Lateral.Just", "Lateral.Ov",
                "Ventral.Shape", "Ventral.Tex", "Ventral.Keel", "Ventral.K.T", "Ventral.Just", "Ventral.Ov",
                "Thigh.Shape", "Thigh.Tex", "Thigh.Keel", "Thigh.K.T", "Thigh.Just", "Thigh.Ov",
                "Foot.Shape", "Foot.Tex", "Foot.Keel", "Foot.K.T", "Foot.Just", "Foot.Ov",
                "Digit.Just", "Digit.Ov", "DCaudal.Shape", "DCaudal.Tex", "DCaudal.Keel", "DCaudal.K.T",
                "DCaudal.Just", "DCaudal.Ov", "VCaudal.Shape", "VCaudal.Tex", "VCaudal.Keel", "VCaudal.K.T",
                "VCaudal.Just", "VCaudal.Ov", "Ap.Pits", "Osteoderms", "Spines", "Horns",
                "Subtympanic", "Supraoc.Disc", "Sens.Head", "Sens.Body", "Scale.Form", "Mid.Ven.Crest"))
colnames(scale.df) <- var_names

head(scale.df)

scale.body <- subset(scale.df,
                     select = -c(145:158))
scale.body <- subset(scale.body,
                     select = -c(75:81))

#Get 3 Dimensional MCA
MCA.3D <- MCA(scale.body, ncp = 3, na.method = "NA")
coords <- MCA.3D$ind$coord
coords <- as.data.frame(coords)

coords_var <- MCA.3D$var$coord
coords_var <- as.data.frame(coords_var)

var_contrib <- MCA.3D$var$contrib
var_contrib <- as.data.frame(var_contrib)
var_contrib$sum <- rowSums(var_contrib)
var_cos2 <- MCA.3D$var$cos2
var_cos2 <- as.data.frame(var_cos2)
var_cos2$sum <- rowSums(var_cos2)


#MCA to Morphospace

coord_colors <- c("#2f3d34", "#918716", rep("#ddcc77", 20), rep("#502c52", 3), rep("#d41cce", 13),
                  rep("#332288", 8), rep("#8fc2db", 8), rep("#2b95ff", 4), rep("#cc7f43", 15),
                  rep("#e66c76", 7), rep("#8a000b", 28), rep("#4ccf84", 7), rep("#117733", 38))
coords$col <- coord_colors

sphen <- coords[1,]
hfb <- subset(coords, col %in% c('#918716', '#8fc2db', '#4ccf84', '#117733'))
gecko <- coords[3:22,]
skink <- coords[23:38,]
lat <- subset(coords, col %in% c('#332288', '#2b95ff'))
ang <- coords[59:73,]
iguana <- coords[74:108,]


#Plot Morphospace


rgl.bg(color = "white")

plot3d(cube3d(scaleMatrix(1.2,1.2,1.2)), alpha = 0, front = 'lines',
       xlab = "Axis 1 (10.1%)", ylab = "Axis 2 (5.3%)", zlab = "Axis 3 (4.3%)")

shapelist3d(dodecahedron3d(), sphen$`Dim 1`, sphen$`Dim 2`, sphen$`Dim 3`,
            size = 0.05, alpha = 0.5, lit = T, shininess = 100,
            color = sphen$col)
shapelist3d(tetrahedron3d(), hfb$`Dim 1`, hfb$`Dim 2`, hfb$`Dim 3`,
            size = 0.045, alpha = 0.5, lit = T, shininess = 100,
            color = hfb$col)
shapelist3d(cube3d(), gecko$`Dim 1`, gecko$`Dim 2`, gecko$`Dim 3`,
            size = 0.04, alpha = 0.5, lit = T, shininess = 100,
            color = gecko$col)
shapelist3d(octahedron3d(), skink$`Dim 1`, skink$`Dim 2`, skink$`Dim 3`,
            size = 0.07, alpha = 0.5, lit = T, shininess = 100,
            color = skink$col)
shapelist3d(icosahedron3d(), lat$`Dim 1`, lat$`Dim 2`, lat$`Dim 3`,
            size = 0.05, alpha = 0.5, lit = T, shininess = 100,
            color = lat$col)
shapelist3d(dodecahedron3d(), ang$`Dim 1`, ang$`Dim 2`, ang$`Dim 3`,
            size = 0.05, alpha = 0.5, lit = T, shininess = 100,
            color = ang$col)
shapelist3d(cuboctahedron3d(), iguana$`Dim 1`, iguana$`Dim 2`, iguana$`Dim 3`,
            size = 0.04, alpha = 0.5, lit = T, shininess = 100,
            color = iguana$col)
rgl.lines(c(0, 0), lim(1), c(0, 0), color = "red")
rgl.lines(c(0, 0), c(0, 0), lim(1), color = "green")
rgl.lines(lim(1), c(0, 0), c(0, 0), color = "blue")

#Variable vector endpoints

shapelist3d(cube3d(), coords_var$`Dim 1`, coords_var$`Dim 2`, coords_var$`Dim 3`,
            size = 0.01, alpha = 0.3, lit = T, shininess = 70,
            color = "#000000")


#Use identify3d to isolate exact individuals and variables

identify3d(coords, labels = rownames(coords))

identify3d(coords_var, labels = rownames(coords_var))

identify3d(coords_var, labels = coords_var$`Dim 1`, adj = c(0, 2))

identify3d(coords_var, labels = coords_var$`Dim 2`, adj = c(0, 3.3))

identify3d(coords_var, labels = coords_var$`Dim 3`, adj = c(0, 4.5))


#Chronophylomorphospace, based on Zheng and Wiens phylogeny

read_nexus_matrix("Zheng_Wiens_Tree_Nodes.nex", equalize_weights = T)

install_github("manabusakamoto/evoldiver")
library(evoldiver)

chronoPTS2D(coords$`Dim 1`, coords$`Dim 2`, tree = zaw, method = "REML",
            radius = 1.5, xlab = "Axis 1", ylab = "Axis 2", zlab = "Time (MYA)",
            col = coord_colors, clade = coord_colors, node.label = F,
            tip.label = T, shadow = T, box = F)






























