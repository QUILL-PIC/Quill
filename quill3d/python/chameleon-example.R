library(ggplot2)
library(grid)
source("chameleon.R")

foo(2)
configure("../results/log", FALSE)
get("ne")

ne <- read2d("../results/rho20", "xy")
ne$z <- (-ne$z)
p <- ggplot(aes(x = x, y = y, z = z), data = ne) + geom_raster(aes(fill = z))
ggsave("example-R.png")
