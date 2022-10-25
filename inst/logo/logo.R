sysfonts::font_add_google("Zilla Slab", "pf", regular.wt = 500)

hexSticker::sticker(
  subplot = ~ plot.new(), s_x = 1, s_y = 1, s_width = 0.1, s_height = 0.1,
  package = "TempStable", p_x = 1, p_y = 1, p_size = 20, h_size = 1.2, p_family = "pf",
  p_color = "#4F76B6", h_fill = "#f3edd2", h_color = "#18357A",
  dpi = 320, filename = "man/figures/logo.png"
)

magick::image_read("man/figures/logo.png")

rstudioapi::restartSession()
