library(ggplot2)
library(igraph)
library(showtext)
library(fmsb)
library(stringr)

# Analysis for changing radius experiments

dirs <- list.dirs(path = "../data/radius_experiments", recursive=FALSE)
#Make sure that files from same replicate get merged appropriately
data <- NULL #This section adapted from stackoverflow.com answer
for (d in dirs) {
  for (d2 in list.dirs(path = d, recursive = FALSE)){
    phenotype_file_path <- paste(d2, "phenotype_count.dat", sep = "/")
    task_file_path <- paste(d2, "tasks.dat", sep = "/")
    if (file.exists(phenotype_file_path) & file.exists(task_file_path)){
      phenotypes <- read.csv(phenotype_file_path, header = F, sep = " ", na.strings = "", colClasses = "character", skip = 8)
      tasks <- read.csv(task_file_path, header = F, sep = " ", na.strings = "", colClasses = "character", skip = 15)
      tasks <- cbind(tasks[,1], tasks[,length(tasks)-1])
      dat <- merge(phenotypes, tasks, by = 1)
      dat$condition <- tail(unlist(strsplit(d, split = "/", fixed = T)), n=1)
      dat$seed <- tail(unlist(strsplit(tail(unlist(strsplit(d2, split = "/", fixed = T)), n=1), split="_", fixed=T))[3], n=1)
      data <- rbind(data, dat)
    }
  }
}

#Remove Null columns
data$V7 <- NULL

#Assign column names
names(data)[1] <- "Updates"
names(data)[2] <- "PhenotypeRichness"
names(data)[3] <- "ShannonDiversityPhenotype"
names(data)[4] <- "PhenotypeRichnessbyTaskCount"
names(data)[5] <- "ShannonDiversityPhenotypeTaskCount"
names(data)[6] <- "AverageTaskDiversity"
names(data)[7] <- "EQU"
#Include factors as factors
data$condition <- factor(data$condition)
#Convert numeric type to numeric
data$seed <- as.numeric(data$seed)
data$Updates <- as.numeric(data$Updates)
data$PhenotypeRichness <- as.numeric(data$PhenotypeRichness)
data$ShannonDiversityPhenotype <- as.numeric(data$ShannonDiversityPhenotype)
data$PhenotypeRichnessbyTaskCount <- as.numeric(data$PhenotypeRichnessbyTaskCount)
data$ShannonDiversityPhenotypeTaskCount <- as.numeric(data$ShannonDiversityPhenotypeTaskCount)
data$AverageTaskDiversity <- as.numeric(data$AverageTaskDiversity)
data$EQU <- as.numeric(as.character(data$EQU))

# Just grab final time points
data$condition <- factor(data$condition)
end_points <- subset(data, (data$Updates==100000))

# Extract radius from directory name
end_points$radius <- as.numeric(unlist(regmatches(end_points$condition, gregexpr('radius\\K\\d{1,2}', end_points$condition , perl=T))))

# Find data from time point immediately before EQU evolved
end_points$preEQUdiv <- vector(mode="numeric", length=length(end_points$seed))
for (s in unique(end_points$seed)){
  #print(seed)
  temp <- subset(data, (data$seed == s))
  #print(length(unique(temp$seed)))
  #print(temp$EQU[10000])
  for (i in order(temp$Updates)){
    if (temp$EQU[i] > 0){
      end_points$preEQUdiv[end_points$seed==s] <- temp$ShannonDiversityPhenotype[i-1]
      break()
    }
  } 
  if (end_points$preEQUdiv[end_points$seed==s] == 0){
    end_points$preEQUdiv[end_points$seed==s] <- temp$ShannonDiversityPhenotype[length(temp$ShannonDiversityPhenotype)]
  }
}

#Do stats

# Kruskal-wallis test to see if there is a difference
kruskal.test(preEQUdiv ~ condition, data=end_points)

#There is! Let's figure out which distributions are actually different with post-hoc tests
mw <- pairwise.wilcox.test(end_points$preEQUdiv, end_points$condition, "bonferroni")

# And now for a bunch of code that figures out to assign significance group labels for the figure
# Inspired by http://stackoverflow.com/questions/23681709/algorithm-for-automating-pairwise-significance-grouping-labels-in-r

n <- 10
g <- as.matrix(mw$p.value > 0.05)
g <- cbind(rbind(NA, g), NA)
g <- replace(g, is.na(g), FALSE)
g <- g + t(g)
diag(g) <- 1
rownames(g) <- 1:n
colnames(g) <- 1:n

# Load data
same <- which(g==1)
topology <- data.frame(N1=((same-1) %% n) + 1, N2=((same-1) %/% n) + 1)
topology <- topology[order(topology[[1]]),] # Get rid of loops and ensure right naming of vertices
g3 <- simplify(graph.data.frame(topology,directed = FALSE))

# Plot graph
#plot(g3)

# Calcuate the maximal cliques
res <- maximal.cliques(g3)

# Reorder given the smallest level
res <- res[order(sapply(res,function(x)paste0(sort(x),collapse=".")))]
ml<-max(sapply(res, length))

# Automatic ordering was harder than it was worth
hacky_order <- c(6, 4, 5, 2, 3, 7, 1)
res <- res[hacky_order]

lab.txt <- vector(mode="list", n)
lab <- letters[seq(res)]
for(i in seq(res)){
  for(j in res[[i]]){
    lab.txt[[j]] <- paste0(lab.txt[[j]], lab[i])
  }
}

# Data from Python simulation of entropies at each radius
# Means
ents = c(0.0, 0.056097261524557046, 0.22681635460177557, 0.505722253371525, 0.9662911524844436, 1.4452577806530846, 2.0693472932292405, 2.59354081668483, 3.158401004365323, 3.7555042724913763, 4.294790391785869, 4.760958299929938, 5.071048954538068, 5.217626813632836, 5.240720517139437, 5.336763638593836, 5.414300841234414, 5.5421317855756005, 5.54141660639101, 5.554431456990773, 5.4637396175981054, 5.460460935036288, 5.420339785053999, 5.322959184385793, 5.166554159189014, 4.922165247714721, 4.760749342247721, 4.5766803128188505, 4.451288556935073, 4.330016721344092, 4.1932342186960145, 4.04218441870403, 3.8760956534284907, 3.649000209964099, 3.429633706680144, 3.2231350132995247, 2.9392202535987515, 2.594859428797114, 2.385890622418992, 2.137912737260153, 1.890928437586362, 1.8334521716442722, 1.6407544712423345, 1.4287784878462098, 1.2731400897130387, 1.0975233088427736, 0.9616336575777238, 0.8382032927199684, 0.559187414764423, 0.3954865800942658, 0.2663181679421692, 0.19940899150338454, 0.1238052830819111, 0.09684100285096864, 0.06339390688178524, 0.05770874541487648, 0.04790915434473234, 0.021830060948708308, 0.015194599383194462, 0.003441297360887945)
# Std deviations
ent_std = c(0.0, 9.7144514654701197e-17, 2.2204460492503131e-16, 5.5511151231257827e-16, 1.2212453270876722e-15, 1.5543122344752192e-15, 3.9968028886505635e-15, 0.0027062393223679692, 0.019186727394950216, 0.039737631244691554, 0.06841635610695368, 0.096814874135081122, 0.13045583109310305, 0.16194032472843445, 0.15733907970442315, 0.16651935186426331, 0.17509468520568214, 0.18629845771344236, 0.21337111863461841, 0.18112269657180649, 0.17862593930489357, 0.16335683947304339, 0.1446867432675995, 0.12331308894186238, 0.1089748551246799, 0.090257327627271233, 0.10568863783947982, 0.1049876937542918, 0.13929314268197648, 0.18112357511176427, 0.21186649267984242, 0.24277664295307522, 0.26372307063247674, 0.24742622176258275, 0.31048141976635479, 0.36555315023907803, 0.38974451903739704, 0.36468125832545911, 0.37661034904855673, 0.36646999866726682, 0.40169749491538692, 0.38731781990933406, 0.33813253803406218, 0.35150486603969006, 0.32971851121808188, 0.30375764291170593, 0.31849864673657008, 0.25235626534673339, 0.20475865068549029, 0.183627352370793, 0.14754770259700836, 0.13879415077665963, 0.1104034581315181, 0.10928450858615596, 0.068952240189564903, 0.061319591211477065, 0.047204192647343497, 0.027240944733794788, 0.015106906681080416, 0.004372760501080074)

dists = c(0:59)
entline = data.frame(ents, ent_std, dists)

radii <- c(6, 8, 10, 15, 20, 24, 30, 36, 48, 59)
labheights <- c(1:10)
for (i in 1:10){
  labheights[i] <- max(subset(end_points$preEQUdiv, end_points$radius==radii[i])) + .25
}

# Rejigger label heights a little so they are all clearly legible
labheights[8] <- labheights[8] + .2
labheights[3] <- labheights[3] + .2

# Make actual plot
font.add("Arial", regular = "arial.ttf", bold = "arialbd.ttf", italic = "ariali.ttf", bolditalic = "arialbi.ttf")

# Save plot with correct fonts
setEPS()
cairo_ps("../figs/DiversityVsRadius.eps", width = 7.5, height= 4.8472)
showtext.begin()
ggplot() + geom_line(data=entline, color="blue", aes(x=dists, y=ents)) + 
  geom_ribbon(data = entline, fill="blue", aes(x=dists, ymin=ents - ent_std, ymax=ents+ent_std), alpha=.3) + 
  geom_boxplot(data=end_points, aes(x=radius, y=preEQUdiv, fill="red", group=condition)) + 
  theme_classic(base_size = 12, base_family = "Arial") + 
  theme(axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), legend.position="none") + 
  scale_x_continuous("Patch Radius", breaks = radii) + scale_y_continuous("Phenotypic Shannon Diversity") + 
  geom_text(aes(label=unlist(lab.txt), x = radii, y = labheights))
dev.off()

# Now let's make a bar chart showing how often EQU evolved
end_points$EQU_evolved <- as.numeric(end_points$EQU>0)

equ_props <- c(1:length(radii))
err_bars <- c(1:length(radii))
for (i in 1:length(radii)){
  equ_props[i] <- sum(subset(end_points$EQU_evolved, end_points$radius == radii[i]))/30
  err_bars[i] <- sqrt(sum(subset(end_points$EQU_evolved, end_points$radius == radii[i])))/30
}

# Do a fisher's exact test
fisher.test(matrix(c(equ_props*30, 30-equ_props*30), ncol = 2), simulate.p.value=TRUE,B=1e7)

# And now pairwise, to see which ones are different from each other
mw <- pairwise.fisher.test(equ_props*30, c(30,30,30,30,30,30,30,30,30,30), "bonferroni")

# Figuring out groups for figure
n <- 10
g <- as.matrix(mw$p.value > 0.05)
g <- cbind(rbind(NA, g), NA)
g <- replace(g, is.na(g), FALSE)
g <- g + t(g)
diag(g) <- 1
rownames(g) <- 1:n
colnames(g) <- 1:n

# Load data
same <- which(g==1)
topology <- data.frame(N1=((same-1) %% n) + 1, N2=((same-1) %/% n) + 1)
topology <- topology[order(topology[[1]]),] # Get rid of loops and ensure right naming of vertices
g3 <- simplify(graph.data.frame(topology,directed = FALSE))

# Calcuate the maximal cliques
res <- maximal.cliques(g3)
res <- res[order(sapply(res,function(x)paste0(sort(x),collapse=".")))]

ml<-max(sapply(res, length))

lab.txt <- vector(mode="list", n)
lab <- letters[seq(res)]
for(i in seq(res)){
  for(j in res[[i]]){
    lab.txt[[j]] <- paste0(lab.txt[[j]], lab[i])
  }
}

# Make figure
setEPS()
cairo_ps("../figs/RadiusEvolutionaryPotential.eps", width = 7.14, height= 4.69)
showtext.begin()
ggplot() + geom_bar(data=end_points, stat="identity", aes(x=radius, y=EQU_evolved/30, group=condition)) + 
  theme_classic(base_size = 12, base_family = "Arial") + 
  theme(axis.line.x = element_line(color = "black"),axis.line.y = element_line(color = "black"), legend.position="none") + 
  scale_x_continuous("Patch Radius", breaks = radii) + scale_y_continuous("Proportion of runs in which EQU Evolved") + 
  geom_errorbar(aes(x=radii, ymin=equ_props-err_bars, ymax=equ_props+err_bars)) + 
  geom_text(aes(label=unlist(lab.txt), x = radii, y=equ_props+err_bars+.05))
dev.off()

# Now the patch distance experiments
# Get data
dirs <- list.dirs(path = "../data/distance_experiments/", recursive=FALSE)
#we're only interested in actual test conditions, which involved changing the distance between patches
#Make sure that files from same replicate get merged appropriately
data <- NULL #This section adapted from stackoverflow.com answer
for (d in dirs) {
  for (d2 in list.dirs(path = d, recursive = FALSE)){
    phenotype_file_path <- paste(d2, "phenotype_count.dat", sep = "/")
    task_file_path <- paste(d2, "tasks.dat", sep = "/")
    if (file.exists(phenotype_file_path) & file.exists(task_file_path)){
      phenotypes <- read.csv(phenotype_file_path, header = F, sep = " ", na.strings = "", colClasses = "character", skip = 8)
      tasks <- read.csv(task_file_path, header = F, sep = " ", na.strings = "", colClasses = "character", skip = 15)
      tasks <- cbind(tasks[,1], tasks[,length(tasks)-1])
      dat <- merge(phenotypes, tasks, by = 1)
      dat$condition <- tail(unlist(strsplit(d, split = "/", fixed = T)), n=1)
      dat$seed <- tail(unlist(strsplit(tail(unlist(strsplit(d2, split = "/", fixed = T)), n=1), split="_", fixed=T))[3], n=1)
      data <- rbind(data, dat)
    }
  }
}

#Remove Null columns
data$V7 <- NULL
#Assign column names
names(data)[1] <- "Updates"
names(data)[2] <- "PhenotypeRichness"
names(data)[3] <- "ShannonDiversityPhenotype"
names(data)[4] <- "PhenotypeRichnessbyTaskCount"
names(data)[5] <- "ShannonDiversityPhenotypeTaskCount"
names(data)[6] <- "AverageTaskDiversity"
names(data)[7] <- "EQU"
#Include factors as factors
data$condition <- factor(data$condition)
#Convert numeric type to numeric
data$seed <- as.numeric(data$seed)
data$Updates <- as.numeric(data$Updates)
data$PhenotypeRichness <- as.numeric(data$PhenotypeRichness)
data$ShannonDiversityPhenotype <- as.numeric(data$ShannonDiversityPhenotype)
data$PhenotypeRichnessbyTaskCount <- as.numeric(data$PhenotypeRichnessbyTaskCount)
data$ShannonDiversityPhenotypeTaskCount <- as.numeric(data$ShannonDiversityPhenotypeTaskCount)
data$AverageTaskDiversity <- as.numeric(data$AverageTaskDiversity)
data$EQU <- as.numeric(as.character(data$EQU))

# Get data from end of experiments
end_points <- subset(data, (data$Updates==100000))
# Extract distance from condition
end_points$distance <- as.numeric(str_sub(common_end_points$condition, -2))

# Get data from before EQU evolved
end_points$preEQUdiv <- vector(mode="numeric", length=length(end_points$seed))
for (s in unique(end_points$seed)){
  #print(seed)
  temp <- subset(data, (data$seed == s))
  #print(length(unique(temp$seed)))
  #print(temp$EQU[10000])
  for (i in order(temp$Updates)){
    if (temp$EQU[i] > 0){
      end_points$preEQUdiv[end_points$seed==s] <- temp$ShannonDiversityPhenotype[i-1]
      break()
    }
  } 
  if (end_points$preEQUdiv[end_points$seed==s] == 0){
    end_points$preEQUdiv[end_points$seed==s] <- temp$ShannonDiversityPhenotype[length(temp$ShannonDiversityPhenotype)]
  }
}

#Do stats
kruskal.test(end_points$preEQUdiv ~ end_points$condition)
mw <- pairwise.wilcox.test(end_points$preEQUdiv, end_points$condition, "bonferroni")

distances = c(0,3,7,10,14,17,21,25,29)

# Figure out significance groups for figure
n <- length(distances)
g <- as.matrix(mw$p.value > 0.05)
g <- cbind(rbind(NA, g), NA)
g <- replace(g, is.na(g), FALSE)
g <- g + t(g)
diag(g) <- 1
rownames(g) <- 1:n
colnames(g) <- 1:n

# Load data
same <- which(g==1)
topology <- data.frame(N1=((same-1) %% n) + 1, N2=((same-1) %/% n) + 1)
topology <- topology[order(topology[[1]]),] # Get rid of loops and ensure right naming of vertices
g3 <- simplify(graph.data.frame(topology,directed = FALSE))

# Calcuate the maximal cliques
res <- maximal.cliques(g3)

# Reorder given the smallest level
res <- sapply(res, sort)
res <- res[order(sapply(res,function(x)paste0(sort(x),collapse=".")))]

hacky_order <- c(3, 4, 5, 2, 1)

res <- res[hacky_order]

lab.txt <- vector(mode="list", n)
lab <- letters[seq(res)]
for(i in seq(res)){
  for(j in res[[i]]){
    lab.txt[[j]] <- paste0(lab.txt[[j]], lab[i])
  }
}

labheights <- c(1:n)
for (i in 1:n){
  labheights[i] <- boxplot.stats(subset(end_points$preEQUdiv, end_points$distance==distances[i]))$stats[5] + .25
}

labheights[7] <- labheights[7] + .2
labheights[2] <- labheights[2] - .1

# Data from python simulations of entropy in environments
ents = c(0.24229218908241523, 0.432655228904054, 0.614370657736406, 0.7789529246751371, 0.9487135192210562, 1.1150109377233368, 1.289260747642993, 1.4607590939477038, 1.6324428087412246, 1.804068302337938, 1.9783643019156631, 2.1304144613528844, 2.274051144585721, 2.428125761545831, 2.5867561384160105, 2.7143408965717493, 2.8344124888634186, 2.9085001591516253, 2.970983203453481, 3.018887072187108, 3.059565686441386, 3.0651594314555695, 3.04960726821926, 3.0168640901376955, 2.969186485555546, 2.899798466822764, 2.877381768261058, 2.8652497334134766, 2.8626831892555, 2.8626831892555, 2.8600529773156, 2.8123980846515235, 2.7438324850098934, 2.6525640047635375, 2.5523241490211808)
ent_std = c(3.6082248300317588e-16, 0.0079180084774929573, 0.007850735805323342, 0.0078421504701223516, 0.0081020113386100819, 0.0094247073665004594, 0.012247854388459404, 0.019314852245278194, 0.029998409135058055, 0.037062945060039112, 0.038238834788564907, 0.037369625331275025, 0.043567382381024151, 0.042383752194369677, 0.052249868358512322, 0.047475149763152549, 0.045834898478013336, 0.037999925330217482, 0.03494166259806971, 0.030604247035814422, 0.026736597989313785, 0.020491462082110341, 0.018741525702981716, 0.012857967001057933, 0.008878140630954481, 0.0037941184790412173, 0.00089241855234213369, 0.00040780695630593725, 7.5495165674510645e-15, 7.5495165674510645e-15, 1.4335639355657646e-07, 6.7549318831457126e-05, 0.00026806375687374204, 0.00094379723705217749, 0.0016212741371206081)

dists = c(0:34)
entline = data.frame(ents, ent_std, dists)

# Make figure

font.add("Arial", regular = "arial.ttf", bold = "arialbd.ttf", italic = "ariali.ttf", bolditalic = "arialbi.ttf")

setEPS()
cairo_ps("../figs/DiversityVsDistance.eps", width = 7.5, height= 4.8472)
showtext.begin()
ggplot() + geom_line(data=entline, color="blue", aes(x=dists, y=ents)) + 
  geom_ribbon(data = entline, fill="blue", aes(x=dists, ymin=ents + ent_std*-1, ymax=ents+ent_std), alpha=.3) + 
  geom_boxplot(data=end_points, aes(x=distance, y=preEQUdiv, fill="red", group=condition)) + 
  theme_classic(base_size = 12, base_family = "Arial") + 
  theme(axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), legend.position="none") + 
  scale_x_continuous("Inter-patch Distance", breaks = distances) + 
  scale_y_continuous("Phenotypic Shannon Diversity") + 
  geom_text(aes(label=unlist(lab.txt), x = distances, y = labheights))
dev.off()

end_points$EQU_evolved <- as.numeric(end_points$EQU>0)

equ_props <- c(1:length(distances))
err_bars <- c(1:length(distances))
for (i in 1:length(distances)){
  equ_props[i] <- sum(subset(end_points$EQU_evolved, end_points$distance == distances[i]))/30
  err_bars[i] <- sqrt(sum(subset(end_points$EQU_evolved, end_points$distance == distances[i])))/30
}

fisher.test(matrix(c(equ_props*30, 30-equ_props*30), ncol = 2), simulate.p.value=TRUE,B=1e7)
mw <- pairwise.fisher.test(equ_props*30, c(30,30,30,30,30,30,30,30,30), "bonferroni")

# Make significance groups
n <- length(distances)
g <- as.matrix(mw$p.value > 0.05)
g <- cbind(rbind(NA, g), NA)
g <- replace(g, is.na(g), FALSE)
g <- g + t(g)
diag(g) <- 1
rownames(g) <- 1:n
colnames(g) <- 1:n

# Load data
same <- which(g==1)
topology <- data.frame(N1=((same-1) %% n) + 1, N2=((same-1) %/% n) + 1)
topology <- topology[order(topology[[1]]),] # Get rid of loops and ensure right naming of vertices
g3 <- simplify(graph.data.frame(topology,directed = FALSE))

# Calcuate the maximal cliques
res <- maximal.cliques(g3)

# Reorder given the smallest level
res <- res[order(sapply(res,function(x)paste0(sort(x),collapse=".")))]

ml<-max(sapply(res, length))

lab.txt <- vector(mode="list", n)
lab <- letters[seq(res)]
for(i in seq(res)){
  for(j in res[[i]]){
    lab.txt[[j]] <- paste0(lab.txt[[j]], lab[i])
  }
}

# Make figure

setEPS()
cairo_ps("../figs/DistanceEvolutionaryPotential.eps", width = 7.14, height= 4.69)
showtext.begin()
ggplot() + geom_bar(data=end_points, stat="identity", aes(x=distance, y=EQU_evolved/30, group=condition)) + 
  theme_classic(base_size = 12, base_family = "Arial") + 
  theme(axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), legend.position="none") + 
  scale_x_continuous("Inter-patch Distance", breaks = distances) + 
  scale_y_continuous("Proportion of runs in which EQU Evolved")  + 
  geom_errorbar(aes(x=distances, ymin=equ_props-err_bars, ymax=equ_props+err_bars))  + 
  geom_text(aes(label=unlist(lab.txt), x = distances, y=equ_props+err_bars+.05))
dev.off()

dirs <- list.dirs(path = "../data/niche_construction_experiments", recursive=FALSE)

# Pull out just the runs that had a period of ecology mode at the end - those are the ones we're plotting
dirs <- dirs[grepl("ecology", dirs)]

#Make sure that files from same replicate get merged appropriately
data <- data.frame() #This section adapted from stackoverflow.com answer
for (d in dirs) {
  for (d2 in list.dirs(path = d, recursive = FALSE)){
    #phenotype_file_path <- #paste(d2, "phenotype_count.dat", sep = "/")
    task_file_path <- paste(d2, "grid_task.200000.dat", sep = "/")
    if (file.exists(task_file_path)){
      #phenotypes <- read.csv(phenotype_file_path, header = F, sep = " ", na.strings = "", colClasses = "character", skip = 8)
      tasks <- read.csv(task_file_path, header = F, sep = " ", na.strings = "", colClasses = "character", skip = 15)
      EQU_Evolved <- length(grep(3, tasks)) > 0
      #print(phenotypes)
      dat <- data.frame(EQU_Evolved)
      #print(dat)
      dat$condition <- tail(unlist(strsplit(d, split = "/", fixed = T)), n=1)
     data <- rbind(data, dat)
    }
  }
}

# Extract information
data$condition <- factor(data$condition)
data$distance <- as.numeric(unlist(regmatches(data$condition, gregexpr('notXorn\\K\\d{1,2}', data$condition , perl=T))))

# Figure out proportions
distances <- unique(data$distance)
sizes <- c(164, 143, 122, 103, 84, 67, 50, 37, 24, 13, 2, 1, 0,0,0,0)
data$size <- sizes
equ_props <- c(1:length(distances))
err_bars <- c(1:length(distances))
for (i in 1:length(distances)){
  equ_props[i] <- sum(subset(data$EQU_Evolved, data$distance == distances[i]))/30
  err_bars[i] <- sqrt(sum(subset(data$EQU_Evolved, data$distance == distances[i])))/30
}

labels <- paste0(as.character(distances), "\n(", as.character(sizes),")")

# Make figure

setEPS()
cairo_ps("../figs/MinimalNiche.eps", width = 7.14, height= 4.69)
showtext.begin()

ggplot() + geom_bar(data=data, stat="identity", aes(x=distance, y=EQU_Evolved/30, group=distance)) + 
  theme_classic(base_size = 12, base_family = "Arial") + 
  theme(axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), legend.position="none") + 
  scale_x_continuous("Inter-patch Distance", breaks = distances, labels=labels) + 
  scale_y_continuous("Proportion of runs in which optimal phenotype evolved") + 
  geom_errorbar(aes(x=distances, ymin=equ_props-err_bars, ymax=equ_props+err_bars))

dev.off()
