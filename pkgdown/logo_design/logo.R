library(ggplot2)
library(hexSticker)
#remotes::install_github('coolbutuseless/cssparser')
#remotes::install_github('coolbutuseless/svgparser')
tiger_df <- svgparser::read_svg("pkgdown/noun-mosquito-2711297.svg", obj_type = 'data.frame')
temp <- ggplot(tiger_df) +
  geom_path(aes(x, y,colour = col, group = interaction(elem_idx, path_idx)),size=1) +
  scale_color_manual(values = c("white"))+
  scale_y_reverse()+ theme_void()+ggeasy::easy_remove_legend()

#s<-sticker(temp, package="plasmoRUtils", p_size=14, s_x=1.05, s_y=.8, s_width=1.5, s_height=1.0,
#           h_fill="#F96167", h_color="#F9E795", filename="inst/figures/lattice.png",h_size=4)

ts <- sticker(temp,
              package="plasmoRUtils",
              p_size=26,
              s_x=1.05,
              s_y=.8,
              s_width=1.5,
              s_height=1.0,
              h_fill="#66A5AD",
              h_color="#C4DFE6",
              filename="inst/figures/lattice.png",
              h_size=5,
              url="Rohit-Satyam/plasmoRUtils",
              u_size = 7,
              u_color = "white",white_around_sticker = FALSE,
              asp = 1,
              dpi = 300)

ts

ggsave('inst/figures/lattice.png', ts, bg='transparent')

#ggsave(file="test.svg", plot=ts, bg='transparent')
