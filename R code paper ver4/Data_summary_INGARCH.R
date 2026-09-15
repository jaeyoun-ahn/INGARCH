source("aux_data_manage2.R")


data.train <- data.train[, c("PolicyNum", "Year", "TypeCity", "TypeCounty", "TypeMisc", "TypeSchool", "TypeTown", "TypeVillage", 
                             "CoverageIM", "lnDeductIM", "NoClaimCreditIM", "FreqIM", "col.Cov", "col.Cov.Idx")]
data.valid <- data.valid[, c("PolicyNum", "Year", "TypeCity", "TypeCounty", "TypeMisc", "TypeSchool", "TypeTown", "TypeVillage", 
                             "CoverageIM", "lnDeductIM", "NoClaimCreditIM", "FreqIM", "col.Cov", "col.Cov.Idx")]

data.train[, c("PolicyNum", "Year", "CoverageIM")]
#do not use NoClaimCredit - this is just an indicator for if 2 years consecutively 0 claims
#allow data to naturally indicate how previous 0's interact with next years

library(data.table)

data.train <- as.data.table(data.train)
data.valid <- as.data.table(data.valid)
data.train[, CoverageIM1 := CoverageIM[1], by = PolicyNum]


get.summary <- function(data){
  data[, .(years = length(unique(Year)), policies = length(unique(PolicyNum)), meanCoverage = mean(CoverageIM), obs = .N, 
                 AvgFreq = mean(FreqIM), NoClaims = sum(FreqIM))]
}
data.train.summary <- get.summary(data.train)
data.valid.summary <- get.summary(data.valid)
makefactor <- function(data){
  data[, Type := factor(max.col(cbind(TypeCity, TypeCounty, TypeMisc, TypeSchool, TypeTown, TypeVillage)))]
  levels(data$Type) <- c("City", "County", "Misc", "School", "Town", "Village")
  data
}
data.train <- makefactor(data.train)
data.valid <- makefactor(data.valid)

data.train.Type <- data.train[, .(policies = length(unique(PolicyNum)), obs = .N, 
                              AvgFreq = round(mean(FreqIM), 3), NoClaims = sum(FreqIM)), 
                              by = Type]
data.train.Type[, ':='(policies_pct = round(100*policies/length(unique(data.train$PolicyNum)), 1), 
                    obs_pct = round(100*obs/nrow(data.train), 1))]

data.train.Idx <- data.train[, .(obs = .N, 
                                 AvgFreq = round(mean(FreqIM), 3), NoClaims = sum(FreqIM)), 
                             by = col.Cov.Idx]

setkeyv(data.train.Idx, "col.Cov.Idx")
data.train.Idx[, col.Cov.Idx := as.numeric(col.Cov.Idx)]
names(data.train.Idx)[1] <- "Coverage tercile"
data.train.Type
data.train.Idx

coverage1 <- unique(data.train[, c("PolicyNum", "CoverageIM1")], by = "PolicyNum")
data.valid <- coverage1[data.valid, on = "PolicyNum"]

library(ggplot2)

coverage_plot_data <- data.frame(
  value=c(data.train$CoverageIM, log(data.train$CoverageIM)),
  scale=factor(rep(c("Raw coverage", "Log coverage"), each=nrow(data.train)),
               levels=c("Raw coverage", "Log coverage"))
)

coverage_plot <- ggplot(coverage_plot_data, aes(value, fill=scale)) +
  geom_histogram(aes(y=after_stat(density)), bins=40,
                 color="white", linewidth=0.2) +
  facet_wrap(~scale, scales="free", nrow=1) +
  scale_fill_manual(values=c("Raw coverage"="#4472C4",
                             "Log coverage"="#70AD47")) +
  labs(title="Distribution of insurance coverage",
#       subtitle="The log transformation substantially reduces right skewness",
       x=NULL, y="Density") +
  theme_minimal(base_size=12) +
  theme(legend.position="none", panel.grid.minor=element_blank(),
        plot.title=element_text(face="bold"),
        strip.text=element_text(face="bold"))

ggsave("CoverageIM_transformation.png", coverage_plot,
       width=9, height=4.5, dpi=300)


# Coverage changes for continuously observed policies ------------------

coverage <- data.train[, .(CoverageIM=CoverageIM[1]), by=.(PolicyNum, Year)]
setorder(coverage, PolicyNum, Year)
complete_ids <- coverage[, .(complete=all(2006:2010 %in% Year)),
                         by=PolicyNum][complete==TRUE, PolicyNum]
complete_coverage <- coverage[PolicyNum %in% complete_ids]
complete_coverage[, indexed_log_coverage :=
  log(CoverageIM/CoverageIM[Year==2006]), by=PolicyNum]

sample_coverage_paths <- function(k, seed=101){
  set.seed(seed)
  complete_coverage[PolicyNum %in% sample(complete_ids, min(k, length(complete_ids)))]
}

coverage_change_figure <- function(k=10, seed=101){
  sampled <- sample_coverage_paths(k, seed)
  ribbon <- complete_coverage[Year>2006,
    .(lower=quantile(indexed_log_coverage, .25),
      median=median(indexed_log_coverage),
      upper=quantile(indexed_log_coverage, .75)), by=Year]

  p1 <- ggplot(sampled, aes(Year, indexed_log_coverage, group=PolicyNum)) +
    geom_line(color="#4472C4", alpha=.35, linewidth=.3) +
    geom_ribbon(data=ribbon, aes(Year, ymin=lower, ymax=upper, group=NULL),
                inherit.aes=FALSE, fill="#4472C4", alpha=.25) +
    geom_line(data=ribbon, aes(Year, median, group=NULL),
              inherit.aes=FALSE, color="#1F4E79", linewidth=1) +
    geom_hline(yintercept=0, color="grey45", linetype="dashed") +
    scale_x_continuous(breaks=2006:2010) +
    labs(title="Indexed coverage trajectories",
         subtitle=paste(k, "randomly selected complete policies"),
         x=NULL, y="Log coverage relative to 2006") +
    theme_minimal(base_size=11)

  changes <- copy(coverage)
  changes[, `:=`(previous=shift(CoverageIM), previous_year=shift(Year)),
          by=PolicyNum]
  changes <- changes[Year==previous_year+1]
  changes[, `:=`(
    annual_change=log(CoverageIM/previous),
    transition=factor(paste0(previous_year, "\u2013", Year),
                      levels=paste0(2006:2009, "\u2013", 2007:2010)))]

  p2 <- ggplot(changes, aes(transition, annual_change)) +
    geom_violin(fill="#70AD47", color=NA, alpha=.45, trim=FALSE) +
    geom_boxplot(width=.16, outlier.shape=NA, fill="white", linewidth=.35) +
    geom_hline(yintercept=0, color="grey45", linetype="dashed") +
    coord_cartesian(ylim=quantile(changes$annual_change, c(.01, .99))) +
    labs(title="Annual changes in log coverage",
         subtitle="All endpoint-complete policies; central 98% displayed",
         x=NULL, y="Annual log-coverage change") +
    theme_minimal(base_size=11)

  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout=grid::grid.layout(1, 2)))
  print(p1, vp=grid::viewport(layout.pos.row=1, layout.pos.col=1))
  print(p2, vp=grid::viewport(layout.pos.row=1, layout.pos.col=2))
  invisible(list(trajectories=p1, annual_changes=p2))
}

png("CoverageIM_changes.png", width=3000, height=1400, res=300)
coverage_change_plots <- coverage_change_figure(k=20, seed = 500)
dev.off()






glm.cov <- glm(FreqIM~log(CoverageIM), 
               family = poisson(link="log"), 
               data = data.train)

glm.cov2 <- glm(FreqIM~factor(col.Cov.Idx), 
                family = poisson(link="log"), 
                data = data.train)

glm.cov3 <- glm(FreqIM~log(CoverageIM1), 
                family = poisson(link="log"), 
                data = data.train)
preds_cov <- predict(glm.cov, newdata = data.valid, type = "response")
preds_cov3 <- predict(glm.cov3, newdata = data.valid, type = "response")
sum(dpois(data.valid$FreqIM, preds_cov, log=TRUE))
sum(dpois(data.valid$FreqIM, preds_cov3, log=TRUE))

glm.basic <- glm(FreqIM~TypeMisc + TypeCounty + TypeSchool + TypeTown + TypeVillage + 
                   log(CoverageIM),
                 family = quasipoisson(link = "log"), 
                 data = data.train)
glm.basic2 <- glm(FreqIM~TypeMisc + TypeCounty + TypeSchool + TypeTown + TypeVillage + 
                   log(CoverageIM1),
                 family = quasipoisson(link = "log"), 
                 data = data.train)

install.packages("lme4")
library(lme4)
re.model <- lmer(log(CoverageIM)~1|PolicyNum, data = data.train)
groupmeans <- data.train[, mean(log(CoverageIM)), by = PolicyNum]
groupvars <- data.train[, var(log(CoverageIM)), by = PolicyNum]
