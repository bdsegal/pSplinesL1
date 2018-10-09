# psplinesl1
This package fits additive mixed models with P-splines and an l1 penalty using alternating direction method of multipliers and cross validation (Segal et al., to appear)

## Installation
```{r}
library(devtools)
install_github("bdsegal/psplinesl1")
```

## Examples
```{r}
library(psplinesl1)
data(simData)

# setup p-spline matrices
X <- list(ps(x = "x", data = simData, 
             norder = 2, k = 1, width = 0.05,
             center = TRUE))

# setup random effect matrices
rand <- re(x = "x", id = "id", data = simData,
           randomCurves = FALSE)

# run cross validation and view paths
cvOut <- cv(y = "y", id = "id", X = X, rand = rand,
                  K = 5,
                  pathLength = 20,
                  data = simData)
plot(cvOut)

# fit model with all data
a1 <- admm(y = "y", id = "id", X = X, rand = rand,
            lambda = cvOut$smoothOpt[2:(length(X)+1)],
            lmeUpdate = TRUE,
            rho = min(5, max(cvOut$smoothOpt)),
            data = simData)

# get and plot fitted model with confidence bands
CI <- ci(model = a1, alpha = 0.05)
plot(CI)

# extract values from ci object for custom plotting
CIpoly <- data.frame(x = c(CI[[1]]$x, rev(CI[[1]]$x)), 
                     y = c(CI[[1]]$lower, rev(CI[[1]]$upper)))

ggplot(aes(x = x, y = y), data = trueMean)+
  geom_polygon(data = CIpoly, fill = "grey")+
  geom_line(color = "black") +
  geom_line(aes(y = smooth), data = CI[[1]], color = "red")
```

## References
Segal, B. D., Elliott, M. R., Braun, T., Jiang, H. (to appear). P-splines with an l1 penalty for repeated measures. Electronic Journal of Statistics.