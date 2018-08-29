#!/usr/bin/env Rscript
setwd("/home/phumbert/raw_data/RNAseq/results_DEG/")

library(limma)
library(ggplot2)
library(dplyr)
library(magrittr)
library(venneuler)
#library(grid)
#library(gridExtra)

load("matrix_DEG.RData")

# HISTOGRAMS for p.value distribution
pdf("Results_visualising.pdf")
ggplot(data=listDEGenes[[1]], aes(P.Value)) + geom_histogram(bins=30,col="green", fill="blue") + labs(title="P.Value distribution Elba1 vs YW") + labs(x="P.Value", y="Counts")
ggplot(data=listDEGenes[[2]], aes(P.Value)) + geom_histogram(bins=30,col="green", fill="blue") + labs(title="P.Value distribution Elba2 vs YW") + labs(x="P.Value", y="Counts")
ggplot(data=listDEGenes[[3]], aes(P.Value)) + geom_histogram(bins=30,col="green", fill="blue") + labs(title="P.Value distribution Elba3 vs YW") + labs(x="P.Value", y="Counts")
ggplot(data=listDEGenes[[4]], aes(P.Value)) + geom_histogram(bins=30,col="green", fill="blue") + labs(title="P.Value distribution Insv vs YW") + labs(x="P.Value", y="Counts")
### using the histograms created in 4_count_matrix.r

# SUMMARIZE #DF genes depending on different cutoffs:

### ELBA1 VS YW
sum_elba1VSyw <- c("ELBA1 VS YW")

###### adjPvalues <0.05
g05=0
for (g in listDEGenes[[1]]$adj.P.Val) {
	if (g < 0.05) g05=g05+1	
}
sum_elba1VSyw <- append(sum_elba1VSyw, c("Number of DF genes with FDR < 0.05 :", g05), after=length(sum_elba1VSyw))

###### adjPvalues <0.1
g1=0
for (g in listDEGenes[[1]]$adj.P.Val) {
	if (g < 0.1) g1=g1+1	
}
sum_elba1VSyw <- append(sum_elba1VSyw, c("Number of DF genes with FDR < 0.1 :", g1), after=length(sum_elba1VSyw))

###### Fold Change >1.3-fold
g1_3=0
for (g in listDEGenes[[1]]$logFC) {
	gbis=exp(g)
	if (gbis > 1.3) g1_3=g1_3+1	
}
sum_elba1VSyw <- append(sum_elba1VSyw, c("Number of DF genes with Fold Change > 1.3-fold :", g1_3), after=length(sum_elba1VSyw))

###### Fold Change >1.5-fold
g1_5=0
for (g in listDEGenes[[1]]$logFC) {
	gbis=exp(g)
	if (gbis > 1.5) g1_5=g1_5+1	
}
sum_elba1VSyw <- append(sum_elba1VSyw, c("Number of DF genes with Fold Change > 1.5-fold :", g1_5), after=length(sum_elba1VSyw))

###### Fold Change >2-fold
g2=0
for (g in listDEGenes[[1]]$logFC) {
	gbis=exp(g)
	if (gbis > 2) g2=g2+1	
}
sum_elba1VSyw <- append(sum_elba1VSyw, c("Number of DF genes with Fold Change > 2-fold :", g2), after=length(sum_elba1VSyw))

###### adjPvalue <0.05 and Fold Change > 1.3-fold
g05_1_3=0
length=0
for (p in listDEGenes[[1]]$adj.P.Val) {
	length=length+1
	if (p < 0.05) {
	f=listDEGenes[[1]]$logFC[length]
	fbis=exp(f)
	if (f > 1.3) g05_1_3=g05_1_3+1
	}			
}
sum_elba1VSyw <- append(sum_elba1VSyw, c("Number of DF genes with FDR < 0.05 and Fold Change > 1.3-fold :", g05_1_3), after=length(sum_elba1VSyw))

###### adjPvalue <0.05 and Fold Change > 1.5-fold
g05_1_5=0
length=0
for (p in listDEGenes[[1]]$adj.P.Val) {
	length=length+1
	if (p < 0.05) {
	f=listDEGenes[[1]]$logFC[length]
	fbis=exp(f)
	if (f > 1.5) g05_1_5=g05_1_5+1
	}			
}
sum_elba1VSyw <- append(sum_elba1VSyw, c("Number of DF genes with FDR < 0.05 and Fold Change > 1.5-fold :", g05_1_5), after=length(sum_elba1VSyw))

###### adjPvalue <0.05 and Fold Change > 2-fold
g05_2=0
length=0
for (p in listDEGenes[[1]]$adj.P.Val) {
	length=length+1
	if (p < 0.05) {
	f=listDEGenes[[1]]$logFC[length]
	fbis=exp(f)
	if (f > 2) g05_2=g05_2+1
	}			
}
sum_elba1VSyw <- append(sum_elba1VSyw, c("Number of DF genes with FDR < 0.05 and Fold Change > 2-fold :", g05_2), after=length(sum_elba1VSyw))

###### adjPvalue <0.1 and Fold Change > 1.3-fold
g1_1_3=0
length=0
for (p in listDEGenes[[1]]$adj.P.Val) {
	length=length+1
	if (p < 0.1) {
	f=listDEGenes[[1]]$logFC[length]
	fbis=exp(f)
	if (f > 1.3) g1_1_3=g1_1_3+1
	}			
}
sum_elba1VSyw <- append(sum_elba1VSyw, c("Number of DF genes with FDR < 0.1 and Fold Change > 1.3-fold :", g1_1_3), after=length(sum_elba1VSyw))

###### adjPvalue <0.1 and Fold Change > 1.5-fold
g1_1_5=0
length=0
for (p in listDEGenes[[1]]$adj.P.Val) {
	length=length+1
	if (p < 0.1) {
	f=listDEGenes[[1]]$logFC[length]
	fbis=exp(f)
	if (f > 1.5) g1_1_5=g1_1_5+1
	}			
}
sum_elba1VSyw <- append(sum_elba1VSyw, c("Number of DF genes with FDR < 0.1 and Fold Change > 1.5-fold :", g1_1_5), after=length(sum_elba1VSyw))

###### adjPvalue <0.1 and Fold Change > 2-fold
g1_2=0
length=0
for (p in listDEGenes[[1]]$adj.P.Val) {
	length=length+1
	if (p < 0.1) {
	f=listDEGenes[[1]]$logFC[length]
	fbis=exp(f)
	if (f > 2) g1_2=g1_2+1
	}			
}
sum_elba1VSyw <- append(sum_elba1VSyw, c("Number of DF genes with FDR < 0.1 and Fold Change > 2-fold :", g1_2), after=length(sum_elba1VSyw))

###### Separate up- and down-regulated

listDEGenes[[1]] <- listDEGenes[[1]]%>%mutate(threshold = ifelse(logFC >1 & adj.P.Val < 0.05, "up-regulated", ifelse(logFC<=(-1) & adj.P.Val <= 0.05, "down-regulated", "not significant")))

up=0
down=0
for (r in listDEGenes[[1]]$threshold) {
	if (r=="up-regulated") up=up+1
	if (r=="down-regulated") down=down+1
}
sum_elba1VSyw <- append(sum_elba1VSyw, c("Number of up-regulated genes:", up, "Number of down-regulated genes:", down), after=length(sum_elba1VSyw))

###### MA Plot and VOLCANO Plot 

ggplot(data=listDEGenes[[1]], aes(x=AveExpr, y=logFC)) + geom_point(aes(colour= threshold), size=1) + 
 scale_colour_manual(values=c("up-regulated"="red","down-regulated" ="green","not significant"="black")) + xlab("Average expression") + 
 ylab("Log2(Fold Change)") + labs(title="MA plot Elba1 VS Yw")

ggplot(data=listDEGenes[[1]], aes(x=logFC, y=(-log10(P.Value))))+geom_point(aes(colour= threshold), size=1) +
 scale_colour_manual(values=c("up-regulated"="red","down-regulated" ="green","not significant"="black"))+xlab("log(Fold Change)")+ ylab("-log10(p.value)") + 
 labs(title="Volcano plot: Elba1 vs YW")

### ELBA2 VS YW
sum_elba2VSyw <- c("ELBA2 VS YW")
###### adjPvalues <0.05
g05=0
for (g in listDEGenes[[2]]$adj.P.Val) {
	if (g < 0.05) g05=g05+1	
}
sum_elba2VSyw <- append(sum_elba2VSyw, c("Number of DF genes with FDR < 0.05 :", g05), after=length(sum_elba2VSyw))

###### adjPvalues <0.1
g1=0
for (g in listDEGenes[[2]]$adj.P.Val) {
	if (g < 0.1) g1=g1+1	
}
sum_elba2VSyw <- append(sum_elba2VSyw, c("Number of DF genes with FDR < 0.1 :", g1), after=length(sum_elba2VSyw))

###### Fold Change >1.3-fold
g1_3=0
for (g in listDEGenes[[2]]$logFC) {
	gbis=exp(g)
	if (gbis > 1.3) g1_3=g1_3+1	
}
sum_elba2VSyw <- append(sum_elba2VSyw, c("Number of DF genes with Fold Change > 1.3-fold :", g1_3), after=length(sum_elba2VSyw))

###### Fold Change >1.5-fold
g1_5=0
for (g in listDEGenes[[2]]$logFC) {
	gbis=exp(g)
	if (gbis > 1.5) g1_5=g1_5+1	
}
sum_elba2VSyw <- append(sum_elba2VSyw, c("Number of DF genes with Fold Change > 1.5-fold :", g1_5), after=length(sum_elba2VSyw))

###### Fold Change >2-fold
g2=0
for (g in listDEGenes[[2]]$logFC) {
	gbis=exp(g)
	if (gbis > 2) g2=g2+1	
}
sum_elba2VSyw <- append(sum_elba2VSyw, c("Number of DF genes with Fold Change > 2-fold :", g2), after=length(sum_elba2VSyw))

###### adjPvalue <0.05 and Fold Change > 1.3-fold
g05_1_3=0
length=0
for (p in listDEGenes[[2]]$adj.P.Val) {
	length=length+1
	if (p < 0.05) {
	f=listDEGenes[[2]]$logFC[length]
	fbis=exp(f)
	if (f > 1.3) g05_1_3=g05_1_3+1
	}			
}
sum_elba2VSyw <- append(sum_elba2VSyw, c("Number of DF genes with FDR < 0.05 and Fold Change > 1.3-fold :", g05_1_3), after=length(sum_elba2VSyw))

###### adjPvalue <0.05 and Fold Change > 1.5-fold
g05_1_5=0
length=0
for (p in listDEGenes[[2]]$adj.P.Val) {
	length=length+1
	if (p < 0.05) {
	f=listDEGenes[[2]]$logFC[length]
	fbis=exp(f)
	if (f > 1.5) g05_1_5=g05_1_5+1
	}			
}
sum_elba2VSyw <- append(sum_elba2VSyw, c("Number of DF genes with FDR < 0.05 and Fold Change > 1.5-fold :", g05_1_5), after=length(sum_elba2VSyw))

###### adjPvalue <0.05 and Fold Change > 2-fold
g05_2=0
length=0
for (p in listDEGenes[[2]]$adj.P.Val) {
	length=length+1
	if (p < 0.05) {
	f=listDEGenes[[2]]$logFC[length]
	fbis=exp(f)
	if (f > 2) g05_2=g05_2+1
	}			
}
sum_elba2VSyw <- append(sum_elba2VSyw, c("Number of DF genes with FDR < 0.05 and Fold Change > 2-fold :", g05_2), after=length(sum_elba2VSyw))

###### adjPvalue <0.1 and Fold Change > 1.3-fold
g1_1_3=0
length=0
for (p in listDEGenes[[2]]$adj.P.Val) {
	length=length+1
	if (p < 0.1) {
	f=listDEGenes[[2]]$logFC[length]
	fbis=exp(f)
	if (f > 1.3) g1_1_3=g1_1_3+1
	}			
}
sum_elba2VSyw <- append(sum_elba2VSyw, c("Number of DF genes with FDR < 0.1 and Fold Change > 1.3-fold :", g1_1_3), after=length(sum_elba2VSyw))

###### adjPvalue <0.1 and Fold Change > 1.5-fold
g1_1_5=0
length=0
for (p in listDEGenes[[2]]$adj.P.Val) {
	length=length+1
	if (p < 0.1) {
	f=listDEGenes[[2]]$logFC[length]
	fbis=exp(f)
	if (f > 1.5) g1_1_5=g1_1_5+1
	}			
}
sum_elba2VSyw <- append(sum_elba2VSyw, c("Number of DF genes with FDR < 0.1 and Fold Change > 1.5-fold :", g1_1_5), after=length(sum_elba2VSyw))

###### adjPvalue <0.1 and Fold Change > 2-fold
g1_2=0
length=0
for (p in listDEGenes[[2]]$adj.P.Val) {
	length=length+1
	if (p < 0.1) {
	f=listDEGenes[[2]]$logFC[length]
	fbis=exp(f)
	if (f > 2) g1_2=g1_2+1
	}			
}
sum_elba2VSyw <- append(sum_elba2VSyw, c("Number of DF genes with FDR < 0.1 and Fold Change > 2-fold :", g1_2), after=length(sum_elba2VSyw))

###### Separate up- and down-regulated

listDEGenes[[2]] <- listDEGenes[[2]]%>%mutate(threshold = ifelse(logFC >1 & adj.P.Val < 0.05, "up-regulated", ifelse(logFC<=(-1) & adj.P.Val <= 0.05, "down-regulated", "not significant")))

up=0
down=0
for (r in listDEGenes[[2]]$threshold) {
	if (r=="up-regulated") up=up+1
	if (r=="down-regulated") down=down+1
}
sum_elba2VSyw <- append(sum_elba2VSyw, c("Number of up-regulated genes:", up, "Number of down-regulated genes:", down), after=length(sum_elba2VSyw))

###### MA Plot and VOLCANO Plot 

ggplot(data=listDEGenes[[2]], aes(x=AveExpr, y=logFC)) + geom_point(aes(colour= threshold), size=1) + 
 scale_colour_manual(values=c("up-regulated"="red","down-regulated" ="green","not significant"="black")) + xlab("Average expression") + 
 ylab("Log2(Fold Change)") + labs(title="MA plot Elba2 VS Yw")

ggplot(data=listDEGenes[[2]], aes(x=logFC, y=(-log10(P.Value))))+geom_point(aes(colour= threshold), size=1) +
 scale_colour_manual(values=c("up-regulated"="red","down-regulated" ="green","not significant"="black"))+xlab("log(Fold Change)")+ ylab("-log10(p.value)") + 
 labs(title="Volcano plot: Elba2 vs YW")

### ELBA3 VS YW
sum_elba3VSyw <- c("ELBA3 VS YW")
###### adjPvalues <0.05
g05=0
for (g in listDEGenes[[3]]$adj.P.Val) {
	if (g < 0.05) g05=g05+1	
}
sum_elba3VSyw <- append(sum_elba3VSyw, c("Number of DF genes with FDR < 0.05 :", g05), after=length(sum_elba3VSyw))

###### adjPvalues <0.1
g1=0
for (g in listDEGenes[[3]]$adj.P.Val) {
	if (g < 0.1) g1=g1+1	
}
sum_elba3VSyw <- append(sum_elba3VSyw, c("Number of DF genes with FDR < 0.1 :", g1), after=length(sum_elba3VSyw))

###### Fold Change >1.3-fold
g1_3=0
for (g in listDEGenes[[3]]$logFC) {
	gbis=exp(g)
	if (gbis > 1.3) g1_3=g1_3+1	
}
sum_elba3VSyw <- append(sum_elba3VSyw, c("Number of DF genes with Fold Change > 1.3-fold :", g1_3), after=length(sum_elba3VSyw))

###### Fold Change >1.5-fold
g1_5=0
for (g in listDEGenes[[3]]$logFC) {
	gbis=exp(g)
	if (gbis > 1.5) g1_5=g1_5+1	
}
sum_elba1VSyw <- append(sum_elba1VSyw, c("Number of DF genes with Fold Change > 1.5-fold :", g1_5), after=length(sum_elba1VSyw))

###### Fold Change >2-fold
g2=0
for (g in listDEGenes[[3]]$logFC) {
	gbis=exp(g)
	if (gbis > 2) g2=g2+1	
}
sum_elba3VSyw <- append(sum_elba3VSyw, c("Number of DF genes with Fold Change > 2-fold :", g2), after=length(sum_elba3VSyw))

###### adjPvalue <0.05 and Fold Change > 1.3-fold
g05_1_3=0
length=0
for (p in listDEGenes[[3]]$adj.P.Val) {
	length=length+1
	if (p < 0.05) {
	f=listDEGenes[[3]]$logFC[length]
	fbis=exp(f)
	if (f > 1.3) g05_1_3=g05_1_3+1
	}			
}
sum_elba3VSyw <- append(sum_elba3VSyw, c("Number of DF genes with FDR < 0.05 and Fold Change > 1.3-fold :", g05_1_3), after=length(sum_elba3VSyw))

###### adjPvalue <0.05 and Fold Change > 1.5-fold
g05_1_5=0
length=0
for (p in listDEGenes[[3]]$adj.P.Val) {
	length=length+1
	if (p < 0.05) {
	f=listDEGenes[[3]]$logFC[length]
	fbis=exp(f)
	if (f > 1.5) g05_1_5=g05_1_5+1
	}			
}
sum_elba3VSyw <- append(sum_elba3VSyw, c("Number of DF genes with FDR < 0.05 and Fold Change > 1.5-fold :", g05_1_5), after=length(sum_elba3VSyw))

###### adjPvalue <0.05 and Fold Change > 2-fold
g05_2=0
length=0
for (p in listDEGenes[[3]]$adj.P.Val) {
	length=length+1
	if (p < 0.05) {
	f=listDEGenes[[3]]$logFC[length]
	fbis=exp(f)
	if (f > 2) g05_2=g05_2+1
	}			
}
sum_elba3VSyw <- append(sum_elba3VSyw, c("Number of DF genes with FDR < 0.05 and Fold Change > 2-fold :", g05_2), after=length(sum_elba3VSyw))

###### adjPvalue <0.1 and Fold Change > 1.3-fold
g1_1_3=0
length=0
for (p in listDEGenes[[3]]$adj.P.Val) {
	length=length+1
	if (p < 0.1) {
	f=listDEGenes[[3]]$logFC[length]
	fbis=exp(f)
	if (f > 1.3) g1_1_3=g1_1_3+1
	}			
}
sum_elba3VSyw <- append(sum_elba3VSyw, c("Number of DF genes with FDR < 0.1 and Fold Change > 1.3-fold :", g1_1_3), after=length(sum_elba3VSyw))

###### adjPvalue <0.1 and Fold Change > 1.5-fold
g1_1_5=0
length=0
for (p in listDEGenes[[3]]$adj.P.Val) {
	length=length+1
	if (p < 0.1) {
	f=listDEGenes[[3]]$logFC[length]
	fbis=exp(f)
	if (f > 1.5) g1_1_5=g1_1_5+1
	}			
}
sum_elba3VSyw <- append(sum_elba3VSyw, c("Number of DF genes with FDR < 0.1 and Fold Change > 1.5-fold :", g1_1_5), after=length(sum_elba3VSyw))

###### adjPvalue <0.1 and Fold Change > 2-fold
g1_2=0
length=0
for (p in listDEGenes[[3]]$adj.P.Val) {
	length=length+1
	if (p < 0.1) {
	f=listDEGenes[[3]]$logFC[length]
	fbis=exp(f)
	if (f > 2) g1_2=g1_2+1
	}			
}
sum_elba3VSyw <- append(sum_elba3VSyw, c("Number of DF genes with FDR < 0.1 and Fold Change > 2-fold :", g1_2), after=length(sum_elba3VSyw))

###### Separate up- and down-regulated

listDEGenes[[3]] <- listDEGenes[[3]]%>%mutate(threshold = ifelse(logFC >1 & adj.P.Val < 0.05, "up-regulated", ifelse(logFC<=(-1) & adj.P.Val <= 0.05, "down-regulated", "not significant")))

up=0
down=0
 for (r in listDEGenes[[3]]$threshold) {
 	if (r=="up-regulated") up=up+1
 	if (r=="down-regulated") down=down+1
 } 
sum_elba3VSyw <- append(sum_elba3VSyw, c("Number of up-regulated genes:", up, "Number of down-regulated genes:", down), after=length(sum_elba3VSyw))

###### MA Plot and VOLCANO Plot 

ggplot(data=listDEGenes[[3]], aes(x=AveExpr, y=logFC)) + geom_point(aes(colour= threshold), size=1) + 
 scale_colour_manual(values=c("up-regulated"="red","down-regulated" ="green","not significant"="black")) + xlab("Average expression") + 
 ylab("Log2(Fold Change)") + labs(title="MA plot Elba3 VS Yw")

ggplot(data=listDEGenes[[3]], aes(x=logFC, y=(-log10(P.Value))))+geom_point(aes(colour= threshold), size=1) +
 scale_colour_manual(values=c("up-regulated"="red","down-regulated" ="green","not significant"="black"))+xlab("log(Fold Change)")+ ylab("-log10(p.value)") + 
 labs(title="Volcano plot: Elba3 vs YW")

### INSV VS YW
sum_insvVSyw <- c("INSV VS YW")
###### adjPvalues <0.05
g05=0
for (g in listDEGenes[[4]]$adj.P.Val) {
	if (g < 0.05) g05=g05+1	
}
sum_insvVSyw <- append(sum_insvVSyw, c("Number of DF genes with FDR < 0.05 :", g05), after=length(sum_insvVSyw))

###### adjPvalues <0.1
g1=0
for (g in listDEGenes[[4]]$adj.P.Val) {
	if (g < 0.1) g1=g1+1	
}
sum_insvVSyw <- append(sum_insvVSyw, c("Number of DF genes with FDR < 0.1 :", g1), after=length(sum_insvVSyw))

###### Fold Change >1.3-fold
g1_3=0
for (g in listDEGenes[[4]]$logFC) {
	gbis=exp(g)
	if (gbis > 1.3) g1_3=g1_3+1	
}
sum_insvVSyw <- append(sum_insvVSyw, c("Number of DF genes with Fold Change > 1.3-fold :", g1_3), after=length(sum_insvVSyw))

###### Fold Change >1.5-fold
g1_5=0
for (g in listDEGenes[[4]]$logFC) {
	gbis=exp(g)
	if (gbis > 1.5) g1_5=g1_5+1	
}
sum_insvVSyw <- append(sum_insvVSyw, c("Number of DF genes with Fold Change > 1.5-fold :", g1_5), after=length(sum_insvVSyw))

###### Fold Change >2-fold
g2=0
for (g in listDEGenes[[4]]$logFC) {
	gbis=exp(g)
	if (gbis > 2) g2=g2+1	
}
sum_insvVSyw <- append(sum_insvVSyw, c("Number of DF genes with Fold Change > 2-fold :", g2), after=length(sum_insvVSyw))

###### adjPvalue <0.05 and Fold Change > 1.3-fold
g05_1_3=0
length=0
for (p in listDEGenes[[4]]$adj.P.Val) {
	length=length+1
	if (p < 0.05) {
	f=listDEGenes[[4]]$logFC[length]
	fbis=exp(f)
	if (f > 1.3) g05_1_3=g05_1_3+1
	}			
}
sum_insvVSyw <- append(sum_insvVSyw, c("Number of DF genes with FDR < 0.05 and Fold Change > 1.3-fold :", g05_1_3), after=length(sum_insvVSyw))

###### adjPvalue <0.05 and Fold Change > 1.5-fold
g05_1_5=0
length=0
for (p in listDEGenes[[4]]$adj.P.Val) {
	length=length+1
	if (p < 0.05) {
	f=listDEGenes[[4]]$logFC[length]
	fbis=exp(f)
	if (f > 1.5) g05_1_5=g05_1_5+1
	}			
}
sum_insvVSyw <- append(sum_insvVSyw, c("Number of DF genes with FDR < 0.05 and Fold Change > 1.5-fold :", g05_1_5), after=length(sum_insvVSyw))

###### adjPvalue <0.05 and Fold Change > 2-fold
g05_2=0
length=0
for (p in listDEGenes[[4]]$adj.P.Val) {
	length=length+1
	if (p < 0.05) {
	f=listDEGenes[[4]]$logFC[length]
	fbis=exp(f)
	if (f > 2) g05_2=g05_2+1
	}			
}
sum_insvVSyw <- append(sum_insvVSyw, c("Number of DF genes with FDR < 0.05 and Fold Change > 2-fold :", g05_2), after=length(sum_insvVSyw))

###### adjPvalue <0.1 and Fold Change > 1.3-fold
g1_1_3=0
length=0
for (p in listDEGenes[[4]]$adj.P.Val) {
	length=length+1
	if (p < 0.1) {
	f=listDEGenes[[4]]$logFC[length]
	fbis=exp(f)
	if (f > 1.3) g1_1_3=g1_1_3+1
	}			
}
sum_insvVSyw <- append(sum_insvVSyw, c("Number of DF genes with FDR < 0.1 and Fold Change > 1.3-fold :", g1_1_3), after=length(sum_insvVSyw))

###### adjPvalue <0.1 and Fold Change > 1.5-fold
g1_1_5=0
length=0
for (p in listDEGenes[[4]]$adj.P.Val) {
	length=length+1
	if (p < 0.1) {
	f=listDEGenes[[4]]$logFC[length]
	fbis=exp(f)
	if (f > 1.5) g1_1_5=g1_1_5+1
	}			
}
sum_insvVSyw <- append(sum_insvVSyw, c("Number of DF genes with FDR < 0.1 and Fold Change > 1.5-fold :", g1_1_5), after=length(sum_insvVSyw))

###### adjPvalue <0.1 and Fold Change > 2-fold
g1_2=0
length=0
for (p in listDEGenes[[4]]$adj.P.Val) {
	length=length+1
	if (p < 0.1) {
	f=listDEGenes[[4]]$logFC[length]
	fbis=exp(f)
	if (f > 2) g1_2=g1_2+1
	}			
}
sum_insvVSyw <- append(sum_insvVSyw, c("Number of DF genes with FDR < 0.1 and Fold Change > 2-fold :", g1_2), after=length(sum_insvVSyw))

###### Separate up- and down-regulated

listDEGenes[[4]] <- listDEGenes[[4]]%>%mutate(threshold = ifelse(logFC >1 & adj.P.Val < 0.05, "up-regulated", ifelse(logFC<=(-1) & adj.P.Val <= 0.05, "down-regulated", "not significant")))

 up=0
 down=0
 for (r in listDEGenes[[4]]$threshold) {
 	if (r=="up-regulated") up=up+1
 	if (r=="down-regulated") down=down+1
 } 
sum_insvVSyw <- append(sum_insvVSyw, c("Number of up-regulated genes:", up, "Number of down-regulated genes:", down), after=length(sum_insvVSyw))

###### MA Plot and VOLCANO Plot 

ggplot(data=listDEGenes[[4]], aes(x=AveExpr, y=logFC)) + geom_point(aes(colour= threshold), size=1) + 
 scale_colour_manual(values=c("up-regulated"="red","down-regulated" ="green","not significant"="black")) + xlab("Average expression") + 
 ylab("Log2(Fold Change)") + labs(title="MA plot Insv VS Yw")

ggplot(data=listDEGenes[[4]], aes(x=logFC, y=(-log10(P.Value))))+geom_point(aes(colour= threshold), size=1) +
 scale_colour_manual(values=c("up-regulated"="red","down-regulated" ="green","not significant"="black"))+xlab("log(Fold Change)")+ ylab("-log10(p.value)") + 
 labs(title="Volcano plot: Insv vs YW")

### Writing all these informations in a .txt file
summa <-file("summarize_genes.txt")
writeLines(c(sum_elba1VSyw, sum_elba2VSyw, sum_elba3VSyw, sum_insvVSyw), summa)
close(summa)


# VENN DIAGRAM: overlap DF genes among different factors

## Venn Diagram for genes up regulated
elba1_up <- listDEGenes[[1]][which(listDEGenes[[1]]$threshold=="up-regulated"),] # use subset function does the same thing (select some rows in dataframe)
elba2_up <- listDEGenes[[2]][which(listDEGenes[[2]]$threshold=="up-regulated"),]
elba3_up <- listDEGenes[[3]][which(listDEGenes[[3]]$threshold=="up-regulated"),]
insv_up <- listDEGenes[[4]][which(listDEGenes[[4]]$threshold=="up-regulated"),]
# creating dataframe with just the genes up-regulated in elba1,2,3 or insv

elba1_up_genes <- elba1_up$GeneID # create a vector with all genes up-regulated in elba1
elba1_up_set <- rep("elba1", each=length(elba1_up_genes)) # create a vector with "elba1" repeated the length of the previous vector as value 
elba2_up_genes <- elba2_up$GeneID
elba2_up_set <- rep("elba2", each=length(elba2_up_genes))
elba3_up_genes <- elba3_up$GeneID
elba3_up_set <- rep("elba3", each=length(elba3_up_genes))
insv_up_genes <- insv_up$GeneID
insv_up_set <- rep("insv", each=length(insv_up_genes))

set_up <- c(elba1_up_set, elba2_up_set, elba3_up_set, insv_up_set) # create a vector with all the values of the variables in argument
elements_up <- c(elba1_up_genes, elba2_up_genes, elba3_up_genes, insv_up_genes) # create a vector with all the values (geneID) in argument
 
up <- data.frame(elements_up,set_up) # create a data frame with 2 columns: 1 set (elba1, 2, 3 or insv), 2 elements (geneID up-regulated in set associated)

v_up <- venneuler(up) # Venn Diagram
v_up$labels <- c(paste("Elba1\n",length(elba1_up_genes)),paste("Elba2\n",length(elba2_up_genes)),paste("Elba3\n",length(elba3_up_genes)),paste("Insv\n",length(insv_up_genes))) #labels in circles= set name and number of genes in this circle
plot(v_up, main="Venn Diagramm for genes up regulated") # plot the Venn Diagramm

## Venn Diagram for genes down regulated
elba1_down <- listDEGenes[[1]][which(listDEGenes[[1]]$threshold=="down-regulated"),] # use subset function does the same thing (select some rows in dataframe)
elba2_down <- listDEGenes[[2]][which(listDEGenes[[2]]$threshold=="down-regulated"),]
elba3_down <- listDEGenes[[3]][which(listDEGenes[[3]]$threshold=="down-regulated"),]
insv_down <- listDEGenes[[4]][which(listDEGenes[[4]]$threshold=="down-regulated"),]
# creating dataframe with just the genes down-regulated in elba1,2,3 or insv

elba1_down_genes <- elba1_down$GeneID # create a vector with all genes down-regulated in elba1
elba1_down_set <- rep("elba1", each=length(elba1_down_genes)) # create a vector with "elba1" repeated the length of the previous vector as value
elba2_down_genes <- elba2_down$GeneID
elba2_down_set <- rep("elba2", each=length(elba2_down_genes))
elba3_down_genes <- elba3_down$GeneID
elba3_down_set <- rep("elba3", each=length(elba3_down_genes))
insv_down_genes <- insv_down$GeneID
insv_down_set <- rep("insv", each=length(insv_down_genes))

set_down <- c(elba1_down_set, elba2_down_set, elba3_down_set, insv_down_set) # create a vector with all the values of the variables in argument
elements_down <- c(elba1_down_genes, elba2_down_genes, elba3_down_genes, insv_down_genes) # create a vector with all the values (geneID) in argument

down <- data.frame(elements_down,set_down) # create a data frame with 2 columns: 1 set (elba1, 2, 3 or insv), 2 elements (geneID down-regulated in set associated)

v_down <- venneuler(down) # Venn Diagram
v_down$labels <- c(paste("Elba1\n",length(elba1_down_genes)),paste("Elba2\n",length(elba2_down_genes)),paste("Elba3\n",length(elba3_down_genes)),paste("Insv\n",length(insv_down_genes))) #label in circles= set name and number of genes in this circle
plot(v_down, main="Venn Diagramm for genes down regulated") # plot the Venn Diagramm


# PCA PLOT for all samples

df <- v$E
df <- as.data.frame(v$E)
df_pca <- prcomp(df)
df_out <- as.data.frame(df_pca$rotation)

df_out <- as.data.frame(df_out)
df_out$group <- c("elba1","elba1","elba1","elba2","elba2","elba2","elba3","elba3","elba3","insv","insv","insv","yw","yw","yw" )

ggplot(df_out,aes(x=PC1,y=PC2))+ geom_point(aes(colour= group), size=1) + 
 scale_colour_manual(values=c("elba1"="red","elba2" ="green","elba3"="blue","insv"="orange","yw"="purple")) + 
 labs(title="PCA plot for all samples")

dev.off()