theme(
  # Legend title and text labels
  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  # Title font color size and face
  legend.title = element_text(color, size, face),
  # Title alignment. Number from 0 (left) to 1 (right)
  legend.title.align = NULL,             
  # Text label font color size and face
  legend.text = element_text(color, size, face), 
  # Text label alignment. Number from 0 (left) to 1 (right)
  legend.text.align = NULL,
  
  # Legend position, margin and background
  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  # Legend position: right, left, bottom, top, none
  legend.position = "right", 
  # Margin around each legend
  legend.margin = margin(0.2, 0.2, 0.2, 0.2, "cm"),
  # Legend background
  legend.background = element_rect(fill, color, size, linetype),
  
  # Legend direction and justification
  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  # Layout of items in legends ("horizontal" or "vertical")
  legend.direction = NULL, 
  # Positioning legend inside or outside plot 
  # ("center" or two-element numeric vector) 
  legend.justification = "center", 
  
  # Background underneath legend keys
  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  legend.key = element_rect(fill, color),  # Key background
  legend.key.size = unit(1.2, "lines"),    # key size (unit)
  legend.key.height = NULL,                # key height (unit)
  legend.key.width = NULL,                 # key width (unit)
  
  # Spacing between legends. 
  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  legend.spacing = unit(0.4, "cm"), 
  legend.spacing.x = NULL,                 # Horizontal spacing
  legend.spacing.y = NULL,                 # Vertical spacing
  
  # Legend box
  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  # Arrangement of multiple legends ("horizontal" or "vertical")
  legend.box = NULL, 
  # Margins around the full legend area
  legend.box.margin = margin(0, 0, 0, 0, "cm"), 
  # Background of legend area: element_rect()
  legend.box.background = element_blank(), 
  # The spacing between the plotting area and the legend box
  legend.box.spacing = unit(0.4, "cm")
)