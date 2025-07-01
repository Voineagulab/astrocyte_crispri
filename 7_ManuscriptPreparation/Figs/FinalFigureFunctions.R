## Palettes
  # pal_iv_discrete <- c("#fe9696", # pink
  #                      "#7a1639", # magenta
  #                      "#aee680", # light-green
  #                      "#619264", # mid-green
  #                      "#354a69", # dark-grey-blue 
  #                      "#f4ead4", # cream
  #                      "#000000") # black
  # 
  # pal_iv_discrete2 <- c("#f3ccde", # pink
  #                       "#dd5d2b", # dark-orange
  #                       "#6f083d", # magenta
  #                       "#fc9d77", # light-orange
  #                       "#3a5f9a", # mid-blue
  #                       "#2c2c54", # dark-blue-grey
  #                       "#ddf4f8") # light-blue
  # 
  # pal_iv_discrete2_outline <- c("#da97b6", # pink
  #                               "#b64519", # dark-orange
  #                               "#39041f", # magenta
  #                               "#e4784e", # light-orange
  #                               "#2d4976", # mid-blue
  #                               "#1c1c35", # dark-blue-grey
  #                               "#a7d2da") # light-blue
  # 
  # # subsets to use for separating n unrelated groups
  # pal_single <- "#354a69" # dark-grey-blue
  # 
  # 
  # # continuous 
  # pal_grn2orng <- c("#b51a00", "#ee5900", "#ff9d68", "#feceb8", "#ffe2d6", "#d9daab", "#84c053", "#4c852e", "#003f01")
  # 
  # # meaningful group separation
  # pal_hits <- c("#3A5F9A", "#F29D6C")
  # pal_hits_outline <- c("#2d4976", "#e47636")
  # # pal_iv_candidates
  
  # function to visualise
  showColours <- function(colours) { 
    if (any(names(pals) %in% colours)) colours <- pals[[colours]]
    
    ggplot(data.frame(id = seq_along(colours), colour = colours)) + 
      geom_tile(aes(1, id, fill = colours)) + 
      geom_text(aes(1, id, label = colours)) + 
      scale_fill_identity()
  }
  
## Latest palettes
  pals <- list()
  
  # hits
  pals$Hits <- c("#3A5F9A", "#F29D6C")
  pals$Hits_Darker <- c("#2d4976", "#e47636")
  
  # continuous
  pals$grn2orng <- c("#b51a00", "#ee5900", "#ff9d68", "#feceb8", "#ffe2d6", "#d9daab", "#84c053", "#4c852e", "#003f01")
  
  # sets of neutral colour
  pals$One <- "#354a69" # dark-grey-blue
  pals$Two <- c("#6F083D", "#8DCACE")
  
  # main discrete palette
  pals$Primary <- c("#fc9d77", "#3d8f80", "#2c2c54", "#f3ccde", "#8dcace", "#3a5f9a", "#6f083d", "#dd5d2b")
  pals$Primary_Darker <- c("#E4784E", "#317266", "#1C1C35", "#DA97B6", "#43979D", "#2D4976", "#39041F", "#B64519")
  
  # a gradient for ggplot2
  gradientInc_fill <- scale_fill_gradient2(low = "white", mid = "grey97", high = "#b51a00")
  # gradientInc_fill_log <- scale_fill_gradient2(low = "white", mid = "grey97", high = "#b51a00", trans = "log2")
  gradientInc_col <- scale_colour_gradient2(low = "white", mid = "grey97", high = "#b51a00")
  # gradientInc_col_log <- scale_colour_gradient2(low = "white", mid = "grey97", high = "#b51a00", trans = "log2")
  
  
## Themes
  invis <- element_blank()
  ts <- function(x) element_text(size = x)
  theme_resizeText <- theme_bw() + theme(axis.text = ts(6), legend.text = ts(6),
                                         axis.title = ts(7))
  standard_theme <- theme_bw() + theme(panel.grid = element_blank(), text = ts(6), strip.text = ts(7),
                                       axis.title = ts(7), axis.line = element_line(), legend.position = "bottom")
  
## Sizing
  # a4
  h_a4 <- 29.7 / 2.54
  w_a4 <- 21.0 / 2.54
  
  h_margin <- (29.7 - (2.54)) / 2.54
  w_margin <- (21 - (2.54)) / 2.54
  
## Text
  text90 <- element_text(angle = 90, hjust = 0.5, vjust = 1)
  