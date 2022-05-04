library(bigstatsr)
library(bigsnpr)
library(dplyr)
library(ggplot2)
library(tibble)
library(magrittr)
library(profvis)



file.sources = list.files(c("R/"), 
                          pattern="*.R$", full.names=TRUE, 
                          ignore.case=TRUE)

sapply(file.sources,source,.GlobalEnv)

n <- 1000
causal <- numeric(n)
causal[1:n/10] <- 1

aids <- simulate_disease(n, 0.05, 0.8,causal = causal)
p1 <- simulate_people(n,aids,"p1", T)
p2 <- simulate_people(n,aids,"p2", T)
child_sibs <- simulate_child_with_sibs(p1, p2, aids, "child_sibs1")

g_liabs <- LTFH(child_sibs)


x <- GWAS(child_sibs, child_sibs$FAM$Status)
Manhattan_plot(x)

y <- GWAX(child_sibs) 
Manhattan_plot(y)

z <- GWAS(child_sibs, g_liabs)
Manhattan_plot(z)