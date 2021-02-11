# likelihood spot check
wd <- "~/2021_SeedDispersal/"
setwd(wd)
dir("res_ouwie/")
Rsaves <- paste0(wd, "res_ouwie/Rosaceae/", dir("res_ouwie/Rosaceae/"))
Rsaves <- Rsaves[grep("se", Rsaves)]
Rsaves <- Rsaves[grep("temp", Rsaves)]
load("~/2021_SeedDispersal/simmaps/Rosaceae-dat.temp.se-21_02_09-simmap-1.Rsave")


load(Rsaves[4])
lapply(obj, function(x) x$solution[3,])


require(OUwie)
Tmax <- max(branching.times(simmaps[[1]]))
test.dummy <- obj[[1]]
phy <- simmaps[[1]]
data <- data.frame(sp = as.character(rownames(test.dummy$data)), reg = test.dummy$data[,1], d = test.dummy$data[,2], ms = test.dummy$data[,3])
alpha <- test.dummy$solution[1,]
sigma <- test.dummy$solution[2,]
theta <- test.dummy$solution[3,]

OUwie.fixed(phy = phy, data = data, model = test.dummy$model, simmap.tree = TRUE, scaleHeight = TRUE, alpha = alpha, sigma.sq = sigma, theta = theta, mserr = "known", algorithm = "three.point")

OUwie.fixed(phy = phy, data = data, model = test.dummy$model, simmap.tree = TRUE, scaleHeight = TRUE, alpha = alpha, sigma.sq = sigma, theta = theta, mserr = "known", algorithm = "invert")



test <- OUwie(phy = phy, data = data, model = "OUM", simmap.tree = TRUE, scaleHeight = TRUE, mserr = "known", algorithm = "three.point", ub = 10)

alpha <- test$solution[1,]
sigma <- test$solution[2,]
theta <- test$solution[3,]

OUwie.fixed(phy = phy, data = data, model = test$model, simmap.tree = TRUE, scaleHeight = TRUE, alpha = alpha, sigma.sq = sigma, theta = theta, mserr = "known", algorithm = "three.point")

OUwie.fixed(phy = phy, data = data, model = test$model, simmap.tree = TRUE, scaleHeight = TRUE, alpha = alpha, sigma.sq = sigma, theta = theta, mserr = "known", algorithm = "invert")
