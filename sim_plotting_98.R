#####
# These functions will create visuals based on
# the stat 98 sim study.
# Goldberg, Mehta, Wassermann
#####

library(ggplot2)

# load simulated data
sim1 <- read.csv('../sim_1.csv')
sim2 <- read.csv('../sim_2.csv')
sim3 <- read.csv('../sim_3.csv')
sim4 <- read.csv('../sim_4.csv')

# create vector for looping
sim.vect <- list(sim1, sim2, sim3, sim4)

# create vector for names
sim.name <- c('Normal 1', 'Normal 2', 'Normal 3', 'Uniform')


####################
# BP Test plotting #
####################
extract.prop.bp <- function(sim.tbl){
  # extract metrics on type II error of BP hyp test
  return(dim(sim.tbl[sim.tbl$bp.test == 'accept',])[1] / dim(sim.tbl)[1])
}

# populate data on Type II errors
bp.type2.vect <- sapply(sim.vect, extract.prop.bp)

bp.t2 <- data.frame(
  Distribution = factor(sim.name, levels=sim.name),
  Type_II_Error = bp.type2.vect
)

# plot bargraph of TII errors
ggplot(data=bp.t2, aes(x=Distribution, y=Type_II_Error, fill=Distribution, label=Type_II_Error)) +
  geom_bar(colour="black", stat="identity") + 
  guides(fill=FALSE) +
  geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  ggtitle('Type II Errors For Breusch-Pagan Test')

ggsave(file="../TII.png")

###############
# B0 Coverage #
###############

extract.prop.b0 <- function(sim.tbl){
  # extract metrics on type II error of BP hyp test
  return(dim(sim.tbl[sim.tbl$b0.cover.bool== FALSE,])[1] / dim(sim.tbl)[1])
}

b0.cov.vect <- sapply(sim.vect, extract.prop.b0)
cov.probs.b0 <- data.frame(
  Distribution = factor(sim.name, levels=sim.name),
  b0_Coverage = b0.cov.vect
)
ggplot(data=cov.probs.b0, aes(x=Distribution, y=b0_Coverage, fill=Distribution, label=b0_Coverage)) +
  geom_bar(colour="black", stat="identity") + 
  guides(fill=FALSE) +
  geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  ggtitle('Beta 0 Confidence Interval Coverage Proportions')

ggsave(file="../B0coverage.png")

###############
# B1 Coverage #
###############

extract.prop.b1 <- function(sim.tbl){
  # extract metrics on type II error of BP hyp test
  return(dim(sim.tbl[sim.tbl$b1.cover.bool== FALSE,])[1] / dim(sim.tbl)[1])
}

b1.cov.vect <- sapply(sim.vect, extract.prop.b1)
cov.probs.b1 <- data.frame(
  Distribution = factor(sim.name, levels=sim.name),
  b1_Coverage = b1.cov.vect
)
ggplot(data=cov.probs.b1, aes(x=Distribution, y=b1_Coverage, fill=Distribution, label=b1_Coverage)) +
  geom_bar(colour="black", stat="identity") + 
  guides(fill=FALSE) +
  geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  ggtitle('Beta 1 Confidence Interval Coverage Proportions')

ggsave(file="../B1coverage.png")

