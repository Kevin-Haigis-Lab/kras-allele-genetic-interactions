+++
# An estimate the timeline for this project from it's current place to first submission.

widget = "blank"  # See https://sourcethemes.com/academic/docs/page-builder/
headless = true  # This file represents a page section.
active = true  # Activate this widget? true/false
weight = 80  # Order that this section will appear.

title = "Estimated Timeline"
subtitle = ""

[design]
  # Choose how many columns the section has. Valid values: 1 or 2.
  columns = "1"

[design.background]
  # Apply a background color, gradient, or image.
  #   Uncomment (by removing `#`) an option to apply it.
  #   Choose a light or dark text color by setting `text_color_light`.
  #   Any HTML color name or Hex value is valid.

  # Background color.
  # color = "#e8f4ff"
  
  # Background gradient.
  # gradient_start = "#d9d9ff"
  # gradient_end = "#d9edff"
  
  # Background image.
  # image = "headers/bubbles-wide.jpg"  # Name of image in `static/img/`.
  # image_darken = 0.6  # Darken the image? Range 0-1 where 0 is transparent and 1 is opaque.

  # Text color (true=light or false=dark).
  text_color_light = false

[design.spacing]
  # Customize the section spacing. Order is top, right, bottom, left.
  padding = ["20px", "0", "20px", "0"]

[advanced]
 # Custom CSS. 
 css_style = ""
 
 # CSS class.
 css_class = ""
+++

## To-Do

1. modeling of allele-specific synthetic lethality (2 weeks)
	1. linear modeling
	2. latent factor analysis
	3. machine learning techniques
2. anlyze the RC-test for comutation and mutual exclusivity results (1 week)
	1. compare against Fisher's extact test results
	2. check overlap with known oncogenes (COSMIC, etc.)
	3. functional enrichment
	4. overlap with synthetic lethal results
3. more detailed analysis of genes-of-interest from the above two steps (2 weeks)
4. write (2 weeks)
	1. figure out storyline
	2. re-make figures

Raw total: 7 weeks  
Accounting for other (TAing, lab retreat, Thanksgiving, misc. nonsense): + 2 weeks  
**Estimated time needed: 9 weeks (~ 2.5 months)**  
**Target date for submission: Friday, Dec. 6, 2019**
