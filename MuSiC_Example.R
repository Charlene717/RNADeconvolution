## MuSiC: Multi-subject Single-cell Deconvolution
## Ref: https://xuranw.github.io/MuSiC/articles/MuSiC.html

##### Presetting ######
rm(list = ls()) # Clean variable
memory.limit(150000)


# install devtools if necessary
install.packages('devtools')

# install the MuSiC package
remove.packages('rlang')
install.packages("rlang")
devtools::install_github('xuranw/MuSiC')

# load
library(MuSiC)



##### Estimation of cell type proportions #####
## Bulk data of human pancreas
# GSE50244.bulk.eset = readRDS('https://xuranw.github.io/MuSiC/data/GSE50244bulkeset.rds')
GSE50244.bulk.eset = readRDS(paste0("./Demo_Data/GSE50244bulkeset.rds"))
GSE50244.bulk.eset


##### Single cell data of human pancreas #####
# Download EMTAB single cell dataset from Github
# EMTAB.eset = readRDS('https://xuranw.github.io/MuSiC/data/EMTABesethealthy.rds')
EMTAB.eset = readRDS('./Demo_Data/EMTABesethealthy.rds')
EMTAB.eset

# Download Xin et al. single cell dataset from Github
# XinT2D.eset = readRDS('https://xuranw.github.io/MuSiC/data/XinT2Deset.rds')
XinT2D.eset = readRDS('./Demo_Data/XinT2Deset.rds')
XinT2D.eset



##### Bulk Tissue Cell Type Estimation #####
# Estimate cell type proportions
Est.prop.GSE50244 = music_prop(bulk.eset = GSE50244.bulk.eset, sc.eset = EMTAB.eset, clusters = 'cellType',
                               samples = 'sampleID', select.ct = c('alpha', 'beta', 'delta', 'gamma',
                                                                   'acinar', 'ductal'), verbose = F)
names(Est.prop.GSE50244)

# Jitter plot of estimated cell type proportions
library(reshape2)
jitter.fig = Jitter_Est(list(data.matrix(Est.prop.GSE50244$Est.prop.weighted),
                             data.matrix(Est.prop.GSE50244$Est.prop.allgene)),
                        method.name = c('MuSiC', 'NNLS'), title = 'Jitter plot of Est Proportions')

# A more sophisticated jitter plot is provided as below. We seperated the T2D subjects and normal
#subjects by their HbA1c levels.

m.prop.GSE50244 = rbind(melt(Est.prop.GSE50244$Est.prop.weighted),
                        melt(Est.prop.GSE50244$Est.prop.allgene))

colnames(m.prop.GSE50244) = c('Sub', 'CellType', 'Prop')
m.prop.GSE50244$CellType = factor(m.prop.GSE50244$CellType, levels = c('alpha', 'beta', 'delta', 'gamma', 'acinar', 'ductal'))
m.prop.GSE50244$Method = factor(rep(c('MuSiC', 'NNLS'), each = 89*6), levels = c('MuSiC', 'NNLS'))
m.prop.GSE50244$HbA1c = rep(GSE50244.bulk.eset$hba1c, 2*6)
m.prop.GSE50244 = m.prop.GSE50244[!is.na(m.prop.GSE50244$HbA1c), ]
m.prop.GSE50244$Disease = factor(c('Normal', 'T2D')[(m.prop.GSE50244$HbA1c > 6.5)+1], levels = c('Normal', 'T2D'))
m.prop.GSE50244$D = (m.prop.GSE50244$Disease == 'T2D')/5
m.prop.GSE50244 = rbind(subset(m.prop.GSE50244, Disease == 'Normal'), subset(m.prop.GSE50244, Disease != 'Normal'))

jitter.new = ggplot(m.prop.GSE50244, aes(Method, Prop)) +
  geom_point(aes(fill = Method, color = Disease, stroke = D, shape = Disease),
             size = 2, alpha = 0.7, position = position_jitter(width=0.25, height=0)) +
  facet_wrap(~ CellType, scales = 'free') + scale_colour_manual( values = c('white', "gray20")) +
  scale_shape_manual(values = c(21, 24))+ theme_minimal()

library(cowplot)
plot_grid(jitter.fig, jitter.new, labels = 'auto')


# Create dataframe for beta cell proportions and HbA1c levels
m.prop.ana = data.frame(pData(GSE50244.bulk.eset)[rep(1:89, 2), c('age', 'bmi', 'hba1c', 'gender')],
                        ct.prop = c(Est.prop.GSE50244$Est.prop.weighted[, 2],
                                    Est.prop.GSE50244$Est.prop.allgene[, 2]),
                        Method = factor(rep(c('MuSiC', 'NNLS'), each = 89),
                                        levels = c('MuSiC', 'NNLS')))
colnames(m.prop.ana)[1:4] = c('Age', 'BMI', 'HbA1c', 'Gender')
m.prop.ana = subset(m.prop.ana, !is.na(HbA1c))
m.prop.ana$Disease = factor( c('Normal', 'T2D')[(m.prop.ana$HbA1c > 6.5) + 1], c('Normal', 'T2D') )
m.prop.ana$D = (m.prop.ana$Disease == 'T2D')/5

ggplot(m.prop.ana, aes(HbA1c, ct.prop)) +
  geom_smooth(method = 'lm',  se = FALSE, col = 'black', lwd = 0.25) +
  geom_point(aes(fill = Method, color = Disease, stroke = D, shape = Disease), size = 2, alpha = 0.7) +  facet_wrap(~ Method) +
  ggtitle('HbA1c vs Beta Cell Type Proportion') + theme_minimal() +
  scale_colour_manual( values = c('white', "gray20")) +
  scale_shape_manual(values = c(21, 24))


lm.beta.MuSiC = lm(ct.prop ~ HbA1c + Age + BMI + Gender, data = subset(m.prop.ana, Method == 'MuSiC'))
lm.beta.NNLS = lm(ct.prop ~ HbA1c + Age + BMI + Gender, data = subset(m.prop.ana, Method == 'NNLS'))
summary(lm.beta.MuSiC)
summary(lm.beta.NNLS)

##### Estimation of cell type proportions with pre-grouping of cell types #####
#### Data Preparation ####
# Download Mouse bulk dataset from Github
# Mouse.bulk.eset = readRDS('https://xuranw.github.io/MuSiC/data/Mousebulkeset.rds')
Mouse.bulk.eset = readRDS('./Demo_Data/Mousebulkeset.rds')
Mouse.bulk.eset

#### Single cell data ####
# Download EMTAB single cell dataset from Github
# Mousesub.eset = readRDS('https://xuranw.github.io/MuSiC/data/Mousesubeset.rds')
Mousesub.eset = readRDS('./Demo_Data/Mousesubeset.rds')
Mousesub.eset

levels(Mousesub.eset$cellType)

## Estimation of cell type proportions
## Clustering single cell data
# Produce the first step information
Mousesub.basis = music_basis(Mousesub.eset, clusters = 'cellType', samples = 'sampleID',
                             select.ct = c('Endo', 'Podo', 'PT', 'LOH', 'DCT', 'CD-PC', 'CD-IC', 'Fib',
                                           'Macro', 'Neutro','B lymph', 'T lymph', 'NK'))

# Plot the dendrogram of design matrix and cross-subject mean of realtive abundance
par(mfrow = c(1, 2))
d <- dist(t(log(Mousesub.basis$Disgn.mtx + 1e-6)), method = "euclidean")
# Hierarchical clustering using Complete Linkage
hc1 <- hclust(d, method = "complete" )
# Plot the obtained dendrogram
plot(hc1, cex = 0.6, hang = -1, main = 'Cluster log(Design Matrix)')
d <- dist(t(log(Mousesub.basis$M.theta + 1e-8)), method = "euclidean")
# Hierarchical clustering using Complete Linkage
# hc2 <- hclust(d, method = "complete" )
hc2 <- hclust(d, method = "complete")
# Plot the obtained dendrogram
plot(hc2, cex = 0.6, hang = -1, main = 'Cluster log(Mean of RA)')

## Bulk Tissue Cell Type Estimation with Pre-grouping of Cell Types
clusters.type = list(C1 = 'Neutro', C2 = 'Podo', C3 = c('Endo', 'CD-PC', 'LOH', 'CD-IC', 'DCT', 'PT'), C4 = c('Macro', 'Fib', 'B lymph', 'NK', 'T lymph'))

cl.type = as.character(Mousesub.eset$cellType)

for(cl in 1:length(clusters.type)){
  cl.type[cl.type %in% clusters.type[[cl]]] = names(clusters.type)[cl]
}
pData(Mousesub.eset)$clusterType = factor(cl.type, levels = c(names(clusters.type), 'CD-Trans', 'Novel1', 'Novel2'))
# 13 selected cell types
s.mouse = unlist(clusters.type)
s.mouse

# load('https://xuranw.github.io/MuSiC/data/IEmarkers.RData')
load('./Demo_Data/IEmarkers.RData')
# This RData file provides two vectors of gene names Epith.marker and Immune.marker

# We now construct the list of group marker
IEmarkers = list(NULL, NULL, Epith.marker, Immune.marker)
names(IEmarkers) = c('C1', 'C2', 'C3', 'C4')
# The name of group markers should be the same as the cluster names

Est.mouse.bulk = music_prop.cluster(bulk.eset = Mouse.bulk.eset, sc.eset = Mousesub.eset, group.markers = IEmarkers, clusters = 'cellType', groups = 'clusterType', samples = 'sampleID', clusters.type = clusters.type)


##### Benchmark Evaluation #####
# Construct artificial bulk dataset. Use all 4 cell types: alpha, beta, gamma, delta
XinT2D.construct.full = bulk_construct(XinT2D.eset, clusters = 'cellType', samples = 'SubjectName')

names(XinT2D.construct.full)
# [1] "Bulk.counts" "num.real"

XinT2D.construct.full$Bulk.counts

head(XinT2D.construct.full$num.real)

# calculate cell type proportions
XinT2D.construct.full$prop.real = relative.ab(XinT2D.construct.full$num.real, by.col = FALSE)
head(XinT2D.construct.full$prop.real)

# Estimate cell type proportions of artificial bulk data
Est.prop.Xin = music_prop(bulk.eset = XinT2D.construct.full$Bulk.counts, sc.eset = EMTAB.eset,
                          clusters = 'cellType', samples = 'sampleID',
                          select.ct = c('alpha', 'beta', 'delta', 'gamma'))


# Estimation evaluation

Eval_multi(prop.real = data.matrix(XinT2D.construct.full$prop.real),
           prop.est = list(data.matrix(Est.prop.Xin$Est.prop.weighted),
                           data.matrix(Est.prop.Xin$Est.prop.allgene)),
           method.name = c('MuSiC', 'NNLS'))

#           RMSD     mAD      R
#MuSiC   0.09881 0.06357 0.9378
#NNLS    0.17161 0.11749 0.8159

library(cowplot)
prop.comp.fig = Prop_comp_multi(prop.real = data.matrix(XinT2D.construct.full$prop.real),
                                prop.est = list(data.matrix(Est.prop.Xin$Est.prop.weighted),
                                                data.matrix(Est.prop.Xin$Est.prop.allgene)),
                                method.name = c('MuSiC', 'NNLS'),
                                title = 'Heatmap of Real and Est. Prop' )

abs.diff.fig = Abs_diff_multi(prop.real = data.matrix(XinT2D.construct.full$prop.real),
                              prop.est = list(data.matrix(Est.prop.Xin$Est.prop.weighted),
                                              data.matrix(Est.prop.Xin$Est.prop.allgene)),
                              method.name = c('MuSiC', 'NNLS'),
                              title = 'Abs.Diff between Real and Est. Prop' )

plot_grid(prop.comp.fig, abs.diff.fig, labels = "auto", rel_widths = c(4,3))



