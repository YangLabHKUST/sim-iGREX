out1 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/sparsityB0.2_sparsityG0.2_SNRy0.3_n800_4000.txt",header = T)
out2 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/sparsityB0.2_sparsityG0.5_SNRy0.3_n800_4000.txt",header = T)
out3 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/sparsityB0.2_sparsityG0.8_SNRy0.3_n800_4000.txt",header = T)
out4 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/sparsityB0.2_sparsityG0.2_SNRy0.3_n1000_4000.txt",header = T)
out5 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/sparsityB0.2_sparsityG0.5_SNRy0.3_n1000_4000.txt",header = T)
out6 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/sparsityB0.2_sparsityG0.8_SNRy0.3_n1000_4000.txt",header = T)
out7 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/sparsityB0.2_sparsityG0.2_SNRy0.3_n2000_4000.txt",header = T)
out8 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/sparsityB0.2_sparsityG0.5_SNRy0.3_n2000_4000.txt",header = T)
out9 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/sparsityB0.2_sparsityG0.8_SNRy0.3_n2000_4000.txt",header = T)

sparsityG <- c(replicate(4,c(0.2,0.5,0.8)))
n1 <- c(sapply(c(800,1000,2000),rep,3))
n_rep <- nrow(out1)

out_all <- data.frame()
for(i in 1:9){
  out <- get(paste("out",i,sep=""))
  out <- melt(out[,1:6],id.vars =NULL)
  out <- cbind(out,method=sapply(substring(out[,1],6),switch,REML="REML",MoM="MoM",ss="SummaryStat"),
               PVE=ifelse(substr(out[,1],4,4)=="g","gene","SNP"),sparsityG=sparsityG[i],n1=n1[i])
  out <- out[,-1]
  out_all <- rbind(out_all,out)
}

names(out_all) <- c("PVE","method","level","sparsity","n1")
out_all$method <- factor(out_all$method,levels = c("REML","MoM","SummaryStat"))
levels(out_all$level) <- c("GREX","Alternative")
pdf("sparsityG.pdf",width=13,height=10)
P  <- ggplot(data=out_all,aes(x=level,y=PVE,fill=method)) +
  geom_boxplot() +
  geom_hline(yintercept = 0.2,linetype="dashed",color="blue") +
  geom_hline(yintercept = 0.3,linetype="dashed",color="red") +
  facet_grid(sparsity~n1, labeller = label_bquote(cols="n"[1]*"="*.(n1),rows=pi[G]*"="*.(sparsity))) +
  scale_fill_brewer(palette="RdBu") +
  theme(text = element_text(size=20))
P
dev.off()



out1 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/sparsityB0.2_sparsityG0.2_SNRy0.3_n1000_4000.txt",header = T)
out2 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/sparsityB0.2_sparsityG0.5_SNRy0.3_n1000_4000.txt",header = T)
out3 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/sparsityB0.2_sparsityG0.8_SNRy0.3_n1000_4000.txt",header = T)
out4 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/sparsityB0.5_sparsityG0.2_SNRy0.3_n1000_4000.txt",header = T)
out5 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/sparsityB0.5_sparsityG0.5_SNRy0.3_n1000_4000.txt",header = T)
out6 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/sparsityB0.5_sparsityG0.8_SNRy0.3_n1000_4000.txt",header = T)
out7 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/sparsityB0.8_sparsityG0.2_SNRy0.3_n1000_4000.txt",header = T)
out8 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/sparsityB0.8_sparsityG0.5_SNRy0.3_n1000_4000.txt",header = T)
out9 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/sparsityB0.8_sparsityG0.8_SNRy0.3_n1000_4000.txt",header = T)

sparsityG <- c(replicate(4,c(0.2,0.5,0.8)))
sparsityB <- c(sapply(c(0.2,0.5,0.8),rep,3))
n_rep <- nrow(out1)

out_all <- data.frame()
for(i in 1:9){
  out <- get(paste("out",i,sep=""))
  out <- melt(out[,1:6],id.vars =NULL)
  out <- cbind(out,method=sapply(substring(out[,1],6),switch,REML="REML",MoM="MoM",ss="SummaryStat"),
               PVE=ifelse(substr(out[,1],4,4)=="g","gene","SNP"),sparsityG=sparsityG[i],sparsityB=sparsityB[i])
  out <- out[,-1]
  out_all <- rbind(out_all,out)
}

names(out_all) <- c("PVE","method","level","sparsityG","sparsityB")
out_all$method <- factor(out_all$method,levels = c("REML","MoM","SummaryStat"))
levels(out_all$level) <- c("GREX","Alternative")
pdf("sparsityG_B.pdf",width=13,height=10)
P  <- ggplot(data=out_all,aes(x=level,y=PVE,fill=method)) +
  geom_boxplot() +
  geom_hline(yintercept = 0.2,linetype="dashed",color="blue") +
  geom_hline(yintercept = 0.3,linetype="dashed",color="red") +
  facet_grid(sparsityG~sparsityB, labeller = label_bquote(cols=pi[B]*"="*.(sparsityB),rows=pi[G]*"="*.(sparsityG))) +
  scale_fill_brewer(palette="RdBu") +
  theme(text = element_text(size=20))
P
dev.off()



library(reshape2)
library(ggplot2)
out1 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/sparsityB0.2_SNRy0.1_n800_4000.txt",header = T)
out2 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/sparsityB0.5_SNRy0.1_n800_4000.txt",header = T)
out3 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/sparsityB0.8_SNRy0.1_n800_4000.txt",header = T)
out4 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/sparsityB0.2_SNRy0.2_n800_4000.txt",header = T)
out5 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/sparsityB0.5_SNRy0.2_n800_4000.txt",header = T)
out6 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/sparsityB0.8_SNRy0.2_n800_4000.txt",header = T)
out7 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/sparsityB0.2_SNRy0.3_n800_4000.txt",header = T)
out8 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/sparsityB0.5_SNRy0.3_n800_4000.txt",header = T)
out9 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/sparsityB0.8_SNRy0.3_n800_4000.txt",header = T)
out10 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/sparsityB0.2_SNRy0.5_n800_4000.txt",header = T)
out11 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/sparsityB0.5_SNRy0.5_n800_4000.txt",header = T)
out12 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/sparsityB0.8_SNRy0.5_n800_4000.txt",header = T)

sparsity <- c(replicate(4,c(0.2,0.5,0.8)))
SNRy <- c(sapply(c(0.1,0.2,0.3,0.5),rep,3))
n_rep <- nrow(out1)

out_all <- data.frame()
for(i in 1:12){
  out <- get(paste("out",i,sep=""))
  out <- melt(out[,1:6],id.vars =NULL)
  out <- cbind(out,method=sapply(substring(out[,1],6),switch,REML="REML",MoM="MoM",ss="SummaryStat"),
               PVE=ifelse(substr(out[,1],4,4)=="g","gene","SNP"),sparsity=sparsity[i],SNRy=SNRy[i])
  out <- out[,-1]
  out_all <- rbind(out_all,out)
}

names(out_all) <- c("PVE","method","level","sparsity","SNRy")
out_all$method <- factor(out_all$method,levels = c("REML","MoM","SummaryStat"))
levels(out_all$level) <- c("GREX","Alternative")
pdf("/Users/cmx/Desktop/presentation/iGREX/image/simu2.pdf",width=13,height=10)
P  <- ggplot(data=out_all,aes(x=level,y=PVE,fill=method)) +
  geom_boxplot() +
  geom_hline(yintercept = 0.2,linetype="dashed",color="blue") +
  geom_hline(yintercept = 0.3,linetype="dashed",color="red") +
  facet_grid(sparsity~SNRy, labeller = label_bquote(cols="SNR"[y]*"="*.(SNRy),rows=pi[1]*"="*.(sparsity))) +
  scale_fill_brewer(palette="RdBu") +
  theme(text = element_text(size=20))
P
dev.off()



out1 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/5prop23_SNRy0.3_n700_4000.txt",header = T)
out2 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/5prop23_SNRy0.3_n800_4000.txt",header = T)
out3 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/5prop23_SNRy0.3_n1000_4000.txt",header = T)
out4 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/5prop23_SNRy0.3_n2000_4000.txt",header = T)
out5 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/5prop23_SNRy0.2_n700_4000.txt",header = T)
out6 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/5prop23_SNRy0.2_n800_4000.txt",header = T)
out7 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/5prop23_SNRy0.2_n1000_4000.txt",header = T)
out8 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/5prop23_SNRy0.2_n2000_4000.txt",header = T)
out9 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/5prop23_SNRy0.1_n700_4000.txt",header = T)
out10 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/5prop23_SNRy0.1_n800_4000.txt",header = T)
out11 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/5prop23_SNRy0.1_n1000_4000.txt",header = T)
out12 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/5prop23_SNRy0.1_n2000_4000.txt",header = T)

n1 <- c(replicate(3,c(700,800,1000,2000)))
SNRy <- c(sapply(c(0.3,0.2,0.1),rep,4))
n_rep <- nrow(out1)

out_all <- data.frame()
for(i in 1:12){
  out <- get(paste("out",i,sep=""))
  out <- melt(out[,1:10],id.vars =NULL)
  out <- cbind(out,method=sapply(substring(out[,1],6),switch,REML="REML",REML0="REML0",MoM="MoM",MoM0="MoM0",ss="SummaryStat"),
                level=ifelse(substr(out[,1],4,4)=="g","gene","SNP"),n1=n1[i],SNRy=SNRy[i])
  out <- out[,-1]
  out_all <- rbind(out_all,out)
}

names(out_all) <- c("PVE","method","level","n1","SNRy")
out_all$method <- factor(out_all$method,levels = c("REML0","MoM0","REML","MoM","SummaryStat"))
levels(out_all$level) <- c("GREX","Alternative")


pdf("/Users/cmx/Desktop/presentation/iGREX/image/simu1.pdf",width=13,height=10)
P  <- ggplot(data=out_all,aes(x=level,y=PVE,fill=method)) +
  geom_boxplot() +
  geom_hline(yintercept = 0.2,linetype="dashed",color="blue") +
  geom_hline(yintercept = 0.3,linetype="dashed",color="red") +
  facet_grid(SNRy~n1, labeller = label_bquote(cols="n"[1]*"="*.(n1),rows="SNR"[y]*"="*.(SNRy))) +
  scale_fill_brewer(palette="RdBu") +
  theme(text = element_text(size=20))
P
dev.off()


out1 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/5prop14_SNRy0.3_n800_4000.txt",header = T)
out2 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/5prop14_SNRy0.3_n1000_4000.txt",header = T)
out3 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/5prop14_SNRy0.3_n2000_4000.txt",header = T)
out4 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/5prop23_SNRy0.3_n800_4000.txt",header = T)
out5 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/5prop23_SNRy0.3_n1000_4000.txt",header = T)
out6 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/5prop23_SNRy0.3_n2000_4000.txt",header = T)
out7 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/5prop32_SNRy0.3_n800_4000.txt",header = T)
out8 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/5prop32_SNRy0.3_n1000_4000.txt",header = T)
out9 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/5prop32_SNRy0.3_n2000_4000.txt",header = T)
out10 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/5prop41_SNRy0.3_n800_4000.txt",header = T)
out11 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/5prop41_SNRy0.3_n1000_4000.txt",header = T)
out12 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/5prop41_SNRy0.3_n2000_4000.txt",header = T)

n1 <- c(replicate(4,c(800,1000,2000)))
PVEg <- c(sapply(c(0.1,0.2,0.3,0.4),rep,3))
n_rep <- nrow(out1)

out_all <- data.frame()
for(i in 1:12){
  out <- get(paste("out",i,sep=""))
  out <- melt(out[,1:10],id.vars =NULL)
  out <- cbind(out,method=sapply(substring(out[,1],6),switch,REML="REML",REML0="REML0",MoM="MoM",MoM0="MoM0",ss="SummaryStat"),
               level=ifelse(substr(out[,1],4,4)=="g","gene","SNP"),n1=n1[i],PVEg=PVEg[i])
  out <- out[,-1]
  out_all <- rbind(out_all,out)
}

names(out_all) <- c("PVE","method","level","n1","PVEg")
out_all$method <- factor(out_all$method,levels = c("REML0","MoM0","REML","MoM","SummaryStat"))
levels(out_all$level) <- c("GREX","Alternative")


pdf("/Users/cmx/Desktop/presentation/iGREX/image/simu3.pdf",width=13,height=10)
P  <- ggplot(data=out_all,aes(x=level,y=PVE,fill=method)) +
  geom_boxplot() +
  geom_hline(aes(yintercept = PVEg),linetype="dashed",color="blue") +
  geom_hline(aes(yintercept = 0.5-PVEg),linetype="dashed",color="red") +
  facet_grid(PVEg~n1, labeller = label_bquote(cols="n"[1]*"="*.(n1),rows="PVE"[G]*"="*.(PVEg))) +
  scale_fill_brewer(palette="RdBu") +
  theme(text = element_text(size=20))
P
dev.off()





################################################ genotype 012 matrix ################################################


out1 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/geno_5prop23_SNRy0.3_n800_4000.txt",header = T)
out2 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/geno_5prop23_SNRy0.3_n1000_4000.txt",header = T)
out3 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/geno_5prop23_SNRy0.3_n2000_4000.txt",header = T)
out4 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/geno_5prop23_SNRy0.3_n5000_4000.txt",header = T)
out5 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/geno_5prop23_SNRy0.2_n800_4000.txt",header = T)
out6 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/geno_5prop23_SNRy0.2_n1000_4000.txt",header = T)
out7 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/geno_5prop23_SNRy0.2_n2000_4000.txt",header = T)
out8 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/geno_5prop23_SNRy0.2_n5000_4000.txt",header = T)
out9 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/geno_5prop23_SNRy0.1_n800_4000.txt",header = T)
out10 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/geno_5prop23_SNRy0.1_n1000_4000.txt",header = T)
out11 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/geno_5prop23_SNRy0.1_n2000_4000.txt",header = T)
out12 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/geno_5prop23_SNRy0.1_n5000_4000.txt",header = T)

n1 <- c(replicate(3,c(800,1000,2000,5000)))
SNRy <- c(sapply(c(0.3,0.2,0.1),rep,4))
n_rep <- nrow(out1)

out_all <- data.frame()
for(i in 1:12){
  out <- get(paste("out",i,sep=""))
  out <- melt(out[,1:10],id.vars =NULL)
  out <- cbind(out,method=sapply(substring(out[,1],6),switch,REML="REML",REML0="REML0",MoM="MoM",MoM0="MoM0",ss="SummaryStat"),
               level=ifelse(substr(out[,1],4,4)=="g","gene","SNP"),n1=n1[i],SNRy=SNRy[i])
  out <- out[,-1]
  out_all <- rbind(out_all,out)
}

names(out_all) <- c("PVE","method","level","n1","SNRy")
out_all$method <- factor(out_all$method,levels = c("REML0","MoM0","REML","MoM","SummaryStat"))
levels(out_all$level) <- c("GREX","Alternative")

pdf("geno_SNRy_n1.pdf",width=13,height=10)
P  <- ggplot(data=out_all,aes(x=level,y=PVE,fill=method)) +
  geom_boxplot() +
  geom_hline(yintercept = 0.2,linetype="dashed",color="blue") +
  geom_hline(yintercept = 0.3,linetype="dashed",color="red") +
  facet_grid(SNRy~n1, labeller = label_bquote(cols="n"[1]*"="*.(n1),rows="SNR"[y]*"="*.(SNRy))) +
  scale_fill_brewer(palette="RdBu") +
  theme(text = element_text(size=20))
P
dev.off()


out1 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/geno_rho0.1_5prop23_SNRy0.3_n1000_4000.txt",header = T)
out2 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/geno_rho0.3_5prop23_SNRy0.3_n1000_4000.txt",header = T)
out3 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/geno_rho0.5_5prop23_SNRy0.3_n1000_4000.txt",header = T)
out4 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/geno_rho0.8_5prop23_SNRy0.3_n1000_4000.txt",header = T)

out1 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/geno_rho0.1_5prop23_SNRy0.3_n5000_40002.txt",header = T)
out2 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/geno_rho0.3_5prop23_SNRy0.3_n5000_40002.txt",header = T)
out3 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/geno_rho0.5_5prop23_SNRy0.3_n5000_40002.txt",header = T)
out4 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/geno_rho0.8_5prop23_SNRy0.3_n5000_40002.txt",header = T)

# rho <- c(replicate(3,c(0.1,0.3,0.5,0.8)))
rho <- c(0.1,0.3,0.5,0.8)
n_rep <- nrow(out1)

out_all <- data.frame()
for(i in 1:4){
  out <- get(paste("out",i,sep=""))
  out <- melt(out[,1:10],id.vars =NULL)
  out <- cbind(out,method=sapply(substring(out[,1],6),switch,REML="REML",REML0="REML0",MoM="MoM",MoM0="MoM0",ss="SummaryStat"),
               level=ifelse(substr(out[,1],4,4)=="g","gene","SNP"),rho=rho[i])
  out <- out[,-1]
  out_all <- rbind(out_all,out)
}

names(out_all) <- c("PVE","method","level","rho")
out_all$method <- factor(out_all$method,levels = c("REML0","MoM0","REML","MoM","SummaryStat"))
levels(out_all$level) <- c("GREX","Alternative")

pdf("geno_SNRy.pdf",width=13,height=10)
P  <- ggplot(data=out_all,aes(x=level,y=PVE,fill=method)) +
  geom_boxplot() +
  geom_hline(yintercept = 0.2,linetype="dashed",color="blue") +
  geom_hline(yintercept = 0.3,linetype="dashed",color="red") +
  facet_grid(.~rho, labeller = label_bquote(cols="rho"*"="*.(rho))) +
  scale_fill_brewer(palette="RdBu") +
  coord_cartesian(ylim = c(0.1,0.45)) +
  theme(text = element_text(size=20))
P
dev.off()




out1 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/MoM_subsample100_prop23_SNRy0.3_n800_4000.txt",header = T)
out2 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/MoM_subsample300_prop23_SNRy0.3_n800_4000.txt",header = T)
out3 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/MoM_subsample500_prop23_SNRy0.3_n800_4000.txt",header = T)
out4 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/MoM_subsample1000_prop23_SNRy0.3_n800_4000.txt",header = T)
out5 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/MoM_subsample100_prop23_SNRy0.3_n1000_4000.txt",header = T)
out6 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/MoM_subsample300_prop23_SNRy0.3_n1000_4000.txt",header = T)
out7 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/MoM_subsample500_prop23_SNRy0.3_n1000_4000.txt",header = T)
out8 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/MoM_subsample1000_prop23_SNRy0.3_n1000_4000.txt",header = T)
out9 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/MoM_subsample100_prop23_SNRy0.3_n2000_4000.txt",header = T)
out10 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/MoM_subsample300_prop23_SNRy0.3_n2000_4000.txt",header = T)
out11 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/MoM_subsample500_prop23_SNRy0.3_n2000_4000.txt",header = T)
out12 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/MoM_subsample1000_prop23_SNRy0.3_n2000_4000.txt",header = T)

subsample <- c(replicate(3,c(100,300,500,1000)))
n1 <- c(sapply(c(800,1000,2000),rep,4))
n_rep <- nrow(out1)

out_all <- data.frame()
for(i in 1:12){
  out <- get(paste("out",i,sep=""))
  out <- melt(out[,1:2],id.vars =NULL)
  out <- cbind(out,method=sapply(substring(out[,1],6),switch,REML="REML",REML0="REML0",MoM="MoM",MoM0="MoM0",ss="SummaryStat"),
               level=ifelse(substr(out[,1],4,4)=="g","gene","SNP"),subsample=subsample[i],n1=n1[i])
  out <- out[,-1]
  out_all <- rbind(out_all,out)
}

names(out_all) <- c("PVE","method","level","subsample","n1")
out_all$method <- factor(out_all$method,levels = c("REML0","MoM0","REML","MoM","SummaryStat"))
levels(out_all$level) <- c("GREX","Alternative")

pdf("geno_subsample.pdf",width=13,height=10)
P  <- ggplot(data=out_all,aes(x=level,y=PVE,fill=level)) +
  geom_boxplot() +
  geom_hline(yintercept = 0.2,linetype="dashed",color="blue") +
  geom_hline(yintercept = 0.3,linetype="dashed",color="red") +
  facet_grid(n1~subsample, labeller = label_bquote(cols="m"*"="*.(subsample),rows="n"[1]*"="*.(n1))) +
  scale_fill_brewer(palette="RdBu") +
  theme(text = element_text(size=20))
P
dev.off()


out1 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/MoM_subsample100_prop23_SNRy0.3_n5000_8000.txt",header = T)
out2 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/MoM_subsample300_prop23_SNRy0.3_n5000_8000.txt",header = T)
out3 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/MoM_subsample500_prop23_SNRy0.3_n5000_8000.txt",header = T)
out4 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/MoM_subsample1000_prop23_SNRy0.3_n5000_8000.txt",header = T)

subsample <- c(replicate(3,c(100,300,500,1000)))
n_rep <- nrow(out1)

out_all <- data.frame()
for(i in 1:4){
  out <- get(paste("out",i,sep=""))
  out <- melt(out[,1:2],id.vars =NULL)
  out <- cbind(out,method=sapply(substring(out[,1],6),switch,REML="REML",REML0="REML0",MoM="MoM",MoM0="MoM0",ss="SummaryStat"),
               level=ifelse(substr(out[,1],4,4)=="g","gene","SNP"),subsample=subsample[i])
  out <- out[,-1]
  out_all <- rbind(out_all,out)
}

names(out_all) <- c("PVE","method","level","subsample")
out_all$method <- factor(out_all$method,levels = c("REML0","MoM0","REML","MoM","SummaryStat"))
levels(out_all$level) <- c("GREX","Alternative")

pdf("geno_subsample.pdf",width=13,height=10)
P  <- ggplot(data=out_all,aes(x=level,y=PVE,fill=level)) +
  geom_boxplot() +
  geom_hline(yintercept = 0.2,linetype="dashed",color="blue") +
  geom_hline(yintercept = 0.3,linetype="dashed",color="red") +
  facet_grid(.~subsample, labeller = label_bquote(cols="m"*"="*.(subsample))) +
  # scale_fill_brewer(palette="RdBu") +
  theme(text = element_text(size=20))
P
dev.off()




out1 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/geno5_rho0.1sparsityB0.2_sparsityG0.2_SNRy0.3_n1000_4000.txt",header = T)
out2 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/geno5_rho0.3sparsityB0.2_sparsityG0.2_SNRy0.3_n1000_4000.txt",header = T)
out3 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/geno5_rho0.5sparsityB0.2_sparsityG0.2_SNRy0.3_n1000_4000.txt",header = T)
out4 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/geno5_rho0.8sparsityB0.2_sparsityG0.2_SNRy0.3_n1000_4000.txt",header = T)

rho <- c(replicate(3,c(0.1,0.3,0.5,0.8)))
n_rep <- nrow(out1)

out_all <- data.frame()
for(i in 1:4){
  out <- get(paste("out",i,sep=""))
  out <- melt(out[,1:10],id.vars =NULL)
  out <- cbind(out,method=sapply(substring(out[,1],6),switch,REML="REML",REML0="REML_0",MoM="MoM",MoM0="MoM_0",ss="SummaryStat"),
               level=ifelse(substr(out[,1],4,4)=="g","gene","SNP"),rho=rho[i])
  out <- out[,-1]
  out_all <- rbind(out_all,out)
}

names(out_all) <- c("PVE","method","level","rho")
out_all$method <- factor(out_all$method,levels = c("REML_0","MoM_0","REML","MoM","SummaryStat"))
levels(out_all$level) <- c("GREX","Alternative")

pdf("autocorr_sparsity.pdf",width=13,height=10)
P  <- ggplot(data=out_all,aes(x=level,y=PVE,fill=method)) +
  geom_boxplot() +
  geom_hline(yintercept = 0.2,linetype="dashed",color="blue") +
  geom_hline(yintercept = 0.3,linetype="dashed",color="red") +
  facet_grid(.~rho, labeller = label_bquote(cols=rho*"="*.(rho))) +
  scale_fill_brewer(palette="RdBu",labels=c(bquote(REML[0]),bquote(MoM[0]),"REML","MoM","SummaryStat")) +
  theme(text = element_text(size=20))
P
dev.off()



library(ggplot2)
library(reshape2)
out1 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/iGREXvsRHOGE_SNRz0.2_SNRy0.1_n800_4000.txt",header = T)
out2 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/iGREXvsRHOGE_SNRz0.2_SNRy0.3_n800_4000.txt",header = T)
out3 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/iGREXvsRHOGE_SNRz0.2_SNRy0.5_n800_4000.txt",header = T)
out4 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/iGREXvsRHOGE_SNRz0.2_SNRy0.7_n800_4000.txt",header = T)
out5 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/iGREXvsRHOGE_SNRz0.2_SNRy0.9_n800_4000.txt",header = T)
out6 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/iGREXvsRHOGE_SNRz0.2_SNRy0.1_n1000_4000.txt",header = T)
out7 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/iGREXvsRHOGE_SNRz0.2_SNRy0.3_n1000_4000.txt",header = T)
out8 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/iGREXvsRHOGE_SNRz0.2_SNRy0.5_n1000_4000.txt",header = T)
out9 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/iGREXvsRHOGE_SNRz0.2_SNRy0.7_n1000_4000.txt",header = T)
out10 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/iGREXvsRHOGE_SNRz0.2_SNRy0.9_n1000_4000.txt",header = T)
out11 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/iGREXvsRHOGE_SNRz0.2_SNRy0.1_n2000_4000.txt",header = T)
out12 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/iGREXvsRHOGE_SNRz0.2_SNRy0.3_n2000_4000.txt",header = T)
out13 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/iGREXvsRHOGE_SNRz0.2_SNRy0.5_n2000_4000.txt",header = T)
out14 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/iGREXvsRHOGE_SNRz0.2_SNRy0.7_n2000_4000.txt",header = T)
out15 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/iGREXvsRHOGE_SNRz0.2_SNRy0.9_n2000_4000.txt",header = T)
out16 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/iGREXvsRHOGE_SNRz0.5_SNRy0.1_n800_4000.txt",header = T)
out17 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/iGREXvsRHOGE_SNRz0.5_SNRy0.3_n800_4000.txt",header = T)
out18 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/iGREXvsRHOGE_SNRz0.5_SNRy0.5_n800_4000.txt",header = T)
out19 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/iGREXvsRHOGE_SNRz0.5_SNRy0.7_n800_4000.txt",header = T)
out20 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/iGREXvsRHOGE_SNRz0.5_SNRy0.9_n800_4000.txt",header = T)
out21 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/iGREXvsRHOGE_SNRz0.5_SNRy0.1_n1000_4000.txt",header = T)
out22 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/iGREXvsRHOGE_SNRz0.5_SNRy0.3_n1000_4000.txt",header = T)
out23 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/iGREXvsRHOGE_SNRz0.5_SNRy0.5_n1000_4000.txt",header = T)
out24 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/iGREXvsRHOGE_SNRz0.5_SNRy0.7_n1000_4000.txt",header = T)
out25 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/iGREXvsRHOGE_SNRz0.5_SNRy0.9_n1000_4000.txt",header = T)
out26 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/iGREXvsRHOGE_SNRz0.5_SNRy0.1_n2000_4000.txt",header = T)
out27 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/iGREXvsRHOGE_SNRz0.5_SNRy0.3_n2000_4000.txt",header = T)
out28 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/iGREXvsRHOGE_SNRz0.5_SNRy0.5_n2000_4000.txt",header = T)
out29 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/iGREXvsRHOGE_SNRz0.5_SNRy0.7_n2000_4000.txt",header = T)
out30 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/iGREXvsRHOGE_SNRz0.5_SNRy0.9_n2000_4000.txt",header = T)

SNRy <- c(replicate(6,c(0.1,0.3,0.5,0.7,0.9)))
SNRz <- c(sapply(c(0.2,0.5),rep,15))
n1 <- c(replicate(2,c(sapply(c(800,1000,2000),rep,5))))
n_rep <- nrow(out1)

out_all <- data.frame()
for(i in 1:30){
  out <- get(paste("out",i,sep=""))
  out <- melt(out[,c(1,2,3,5)],id.vars =NULL)
  out <- cbind(out,SNRz=SNRz[i],SNRy=SNRy[i],n1=n1[i])
  out_all <- rbind(out_all,out)
}
names(out_all) <- c("method","PVE","SNRz","SNRy","n1")
out_all$SNRy <- as.factor(out_all$SNRy)

pdf("RHOGE_n1_SNRz.pdf",width=13,height=10)
P  <- ggplot(data=out_all,aes(x=SNRy,y=PVE,fill=method)) +
  geom_boxplot() +
  geom_hline(aes(yintercept = SNRz),linetype="dashed",color="blue") +
  facet_grid(SNRz~n1,scales="free", labeller = label_bquote(cols="n"[1]*"="*.(n1),rows=SNR[z]*"="*.(SNRz))) +
  scale_fill_brewer(palette="RdBu")
P
dev.off()

