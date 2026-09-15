## In-sample diagnostics for the fitted heterogeneous NB-INGARCH model.
## Pearson residuals use the one-step-ahead NB predictive distribution;
## PIT values are randomized within the discrete probability mass at y.
source("aux_data_manage4.R")
.libPaths(c("Rlib", .libPaths()))
library(RTMB)
source("funs_INGARCH.R")

W_train <- !is.na(Y_train)
obj <- make_state_nb_objective()
fit <- optim(obj$par, obj$fn, obj$gr, method="L-BFGS-B",
             lower=c(rep(-Inf, p), 0, -Inf),
             upper=c(rep( Inf, p), .999, Inf),
             control=list(maxit=5000, factr=1e7))
u <- unpack_state_nb(fit$par)
lambda <- exp(train_predictor(u$beta))
f <- gamma_filter(Y_train, W_train, lambda, u$delta, u$shape)
A <- f$A[, 1:ncol(Y_train)]
B <- f$B[, 1:ncol(Y_train)]
mu <- A*lambda/B
y <- Y_train
ok <- W_train & is.finite(mu) & mu > 0 & A > 0

pearson <- (y - mu) / sqrt(mu + mu^2/A)
set.seed(20260912)
pit <- rep(NA_real_, length(y))
pit[ok] <- pnbinom(y[ok] - 1, size=A[ok], mu=mu[ok]) +
  runif(sum(ok)) * dnbinom(y[ok], size=A[ok], mu=mu[ok])

diag <- data.frame(policy=rep(seq_len(nrow(y)), ncol(y)),
                   year=rep(2006:2010, each=nrow(y)),
                   y=as.vector(y), mean=as.vector(mu),
                   shape=as.vector(A), pearson=as.vector(pearson),
                   pit=pit)
diag <- diag[ok, ]
write.csv(diag, "heterogeneous_LGPIF_diagnostics.csv", row.names=FALSE)

png("heterogeneous_LGPIF_diagnostics.png", width=2400, height=850, res=160)
par(mfrow=c(1,3), mar=c(4.2,4.2,3,1))
hist(diag$pit, breaks=seq(0, 1, length.out=11), col="#70AD47", border="white",
     main="Randomized PIT", xlab="PIT", xlim=c(0,1))
abline(h=nrow(diag)/10, lty=2, lwd=1.2)

pearson_breaks <- seq(-1.2, 12, by=0.1)
pearson_plot <- pmin(diag$pearson, 12)
pearson_hist <- hist(pearson_plot, breaks=pearson_breaks, plot=FALSE)
hist_height <- max(pearson_hist$counts)
plot(pearson_hist, col="#4472C4", border="white",
     xlim=range(pearson_breaks), ylim=c(-0.09*hist_height, 1.08*hist_height),
     main="Pearson residuals", xlab="Standardized Pearson residual",
     ylab="Frequency")
set.seed(20260913)
strip_y <- runif(nrow(diag), -0.075*hist_height, -0.015*hist_height)
strip_extreme <- abs(diag$pearson) > 5
points(pearson_plot[!strip_extreme], strip_y[!strip_extreme], pch=16, cex=0.22,
       col=adjustcolor("#4472C4", alpha.f=0.28))
points(pearson_plot[strip_extreme], strip_y[strip_extreme], pch=16, cex=0.38,
       col="#C00000")
abline(h=0, col="grey50", lwd=0.8)
outlier_n <- sum(diag$pearson > 12)
abline(v=0, lty=2, lwd=1.2)
signed_log_resid <- sign(diag$pearson) * log10(1 + abs(diag$pearson))
is_extreme <- abs(diag$pearson) > 5
plot(diag$mean, signed_log_resid, log="x", pch=16, cex=0.45,
     col=ifelse(is_extreme, "#C00000", adjustcolor("#4472C4", alpha.f=0.35)),
     xlab="Fitted mean (log scale)", ylab="Signed log Pearson residual",
     main="Residuals vs fitted mean")
abline(h=0, lty=2, lwd=1.2)
abline(h=c(-log10(6), log10(6)), lty=3, col="grey60")
legend("topright", legend=c("|Pearson residual| > 5", "other"),
       col=c("#C00000", "#4472C4"), pch=16, pt.cex=0.65,
       bty="n", cex=0.75)

dev.off()

sink("heterogeneous_LGPIF_diagnostics_summary.txt")
cat("convergence:", fit$convergence, "\n")
cat("n observations:", nrow(diag), "\n")
cat("observations above the histogram display limit:", outlier_n, "\n")
cat("parameter estimates:\n")
print(c(u$beta, delta=u$delta, gamma_1_0=u$shape))
cat("Pearson mean/sd:", mean(diag$pearson), sd(diag$pearson), "\n")
cat("PIT mean/sd:", mean(diag$pit), sd(diag$pit), "\n")
sink()
