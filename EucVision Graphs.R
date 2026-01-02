# Libraries
library(ggplot2)
library(dplyr)
library(tinyplot)
library(viridis)



# create a dataset
data <- data.frame(
  name=c( rep("A",500), rep("B",500), rep("B",500), rep("C",20), rep('D', 100)  ),
  value=c( rnorm(500, 10, 5), rnorm(500, 13, 1), rnorm(500, 18, 1), rnorm(20, 25, 4), rnorm(100, 12, 1) )
)

# sample size
sample_size = data %>% group_by(name) %>% summarize(num=n())

# Plot
data %>%
  left_join(sample_size) %>%
  mutate(myaxis = paste0(name, "\n", "n=", num)) %>%
  ggplot( aes(x=myaxis, y=value, fill=name)) +
  geom_violin(width=1.4) +
  geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  scale_fill_viridis(discrete = TRUE) +
  ggtitle("A Violin wrapping a boxplot") +
  xlab("")






# "violin" type convenience string
tinyplot(count ~ spray, data = InsectSprays, type = "violin")

# aside: to match the defaults of `ggplot2::geom_violin()`, use `trim = TRUE`
# and `joint.bw = FALSE`
tinyplot(count ~ spray, data = InsectSprays, type = "violin",
         trim = TRUE, joint.bw = FALSE)

# use flip = TRUE to reorient the axes
tinyplot(count ~ spray, data = InsectSprays, type = "violin", flip = TRUE)

# for flipped plots with long group labels, it's better to use a theme for
# dynamic plot resizing
tinytheme("clean")
tinyplot(weight ~ feed, data = chickwts, type = "violin", flip = TRUE)

tinytheme("clean")
# you can group by the x var to add colour (here with the original orientation)
tinyplot(weight ~ feed | feed, data = chickwts, type = "violin", legend = FALSE)

# dodged grouped violin plot example (different dataset)
tinyplot(len ~ dose | supp, data = ToothGrowth, type = "violin", fill = 0.2)

# note: above we relied on `...` argument passing alongside the "violin"
# type convenience string. But this won't work for `width`, since it will
# clash with the top-level `tinyplot(..., width = <width>)` arg. To ensure
# correct arg passing, it's safer to use the formal `type_violin()` option.
tinyplot(len ~ dose | supp, data = ToothGrowth, fill = 0.2,
         type = type_violin(width = 0.8))

# reset theme
tinytheme()