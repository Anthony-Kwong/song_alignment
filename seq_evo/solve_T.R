#OPT approach to estimating note transitions on the population level

X = readRDS(file = "./results/seq_evo/sf_mat.rds")
Y = readRDS(file = "./results/seq_evo/sons_mat.rds")

library("lpSolve")
readr::write_csv(X, file = "./results/seq_evo/sf_mat.csv")
write.table(X, file="./results/seq_evo/sf_mat.txt", row.names=F, col.names=F)
write.table(Y, file="./results/seq_evo/sons_mat.txt", row.names=F, col.names=F)

source("./seq_evo/EstimateTransmissionMat.R")

res = estimateTransmissionMatrix(father.mat = X, son.mat = Y, cost.mat=NULL, tol=1.0e-10)
opt_T = matrix(res$par, nrow = 16)
rownames(opt_T) = LETTERS[1:16]
colnames(opt_T) = LETTERS[1:16]

T_star = as.data.frame(opt_T)
ggplot(data = T_star, aes(x = Var2, y = Var1, fill = as.factor(matrix))) +
  geom_tile() +
  scale_fill_manual(values = rev(heat.colors(256))) +
  labs(x = "Column", y = "Row", title = "Heatmap of Matrix")

