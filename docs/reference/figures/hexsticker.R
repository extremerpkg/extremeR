library(hexSticker)
library(ggplot2)

plt <- ggplot(data.frame(x = c(-4, 4)), aes(x = x)) +
  stat_function(fun = dnorm, color = "white") +
  labs(x = "", y = "") +
  xlim(0,4) + theme(axis.ticks = element_blank(),
                    axis.text = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    axis.line = element_line(colour = "black"),
                    panel.background = element_rect(fill = "transparent"),
                    plot.background = element_rect(fill = "transparent", color = NA)) +
  geom_area(stat = "function",
            fun = dnorm, fill = "#00998a",
            xlim = c(2, 4))

sticker(plt, package="extremeR", p_size=20, s_x=1.05, s_y=.8, s_width=1, s_height=1,
        h_fill="#373737", h_color="grey", filename = "man/figures/hex.png")
