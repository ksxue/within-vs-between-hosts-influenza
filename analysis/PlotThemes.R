# Set common palette for all plots.
library(RColorBrewer)

DARK2PALETTE <- brewer.pal(n=8, name="Dark2")

SET2PALETTE <- c("gray30", # dark grey
             "darkorange1", # orange
             "mediumseagreen", # light bluish-green
             "plum3", # light purple
             "palevioletred1", # light pink
             "olivedrab3", # light yellow-green
             "goldenrod1", # golden yellow
             "dodgerblue3", # light blue
             "wheat", # beige
             "lightgray", # light grey
             "black" #black
)

DESERTPALETTE <- c("gray30", # dark grey
                   "coral3", # reddish-orange
                   "thistle3", #light purple
                   "darkseagreen", #light green
                   "lightgoldenrod1", #light yellow
                   "rosybrown2", #pink
                   "sandybrown", #light orange
                   "lightgray" #light gray
)

PALETTE <- DESERTPALETTE

SHAPES <- c(16,15,17)

# Set common formatting for all plots.
THEME_ALL <- theme(
  text=element_text(size=7),
  axis.title=element_text(size=7),
  axis.text=element_text(size=5.5),
  strip.text.x=element_text(margin = margin(3,0,3,0), size=7),
  strip.text.y=element_text(margin = margin(0,3,0,3), size=7),
  strip.background=element_blank())

THEME_PRES <- theme(
  text=element_text(size=14),
  axis.title=element_text(size=14, face="bold"),
  axis.text=element_text(size=8),
  strip.text.x=element_text(size=10, margin = margin(3,0,3,0)),
  strip.text.y=element_text(size=10, margin = margin(0,3,0,3)),
  strip.background=element_blank())
