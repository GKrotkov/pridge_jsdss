# Visualizing bootstrap CIs
rm(list = ls())

library(tidyverse)
library(gt)
library(forestplot)

blair_red <- "#a7000a"
load("data/bootstrap_results.rda")

###############
#### Table ####
###############

bounds_gt <- gt(bounds_tbl) |>
    tab_header(title = "pRidge Mean % Improvement against Baseline",
               subtitle = "Bootstrap 95% CIs (2.5th, 50th, 97.5th)") |>
    tab_style(
        style = cell_borders(
            sides = c("left", "right", "bottom", "top"),
            color = "black",
            weight = px(2)
        ),
        locations = cells_body(columns = 1:4)
    ) |>
    tab_style(
        style = cell_borders(
            sides = c("t"),
            color = blair_red,
            weight = px(4)
        ),
        locations = list(cells_title())
    ) |>
    tab_style(
        style = cell_borders(
            sides = c("t", "b", "l", "r"),
            color = "black",
            weight = px(3)
        ),
        locations = list(cells_column_labels())
    )

gtsave(bounds_gt, "output/bootstrap_cis.png")

#####################
#### Forest Plot ####
#####################

bounds_epa$baseline <- "EPA"
bounds_opr$baseline <- "OPR"

viz <- rbind(bounds_opr, bounds_epa) |>
    rename(labeltext = year, mean = median) |>
    mutate(lower = round(lower, digits = 1),
           mean = round(mean, digits = 1),
           upper = round(upper, digits = 1))

viz |>
    group_by(baseline) |>
    forestplot(ci.vertices = TRUE,
               ci.vertices.height = 0.1,
               boxsize = 0.2,
               # cutting out the high outlier 2018 performance for visibility
               clip = c(-5, 50),
               xlab = "% improvement over baseline",
               txt_gp = fpTxtGp(
                   xlab = gpar(fontsize = 18),
                   ticks = gpar(fontsize = 18)
               )) |>
    fp_add_lines("steelblue") |>
    fp_add_header("Subgroup") |>
    fp_set_style(box = c(blair_red, "forestgreen"),
                 default = gpar(verticies = TRUE),
                 line = "black") |>
    fp_set_zebra_style("#EFEFEF") |>
    fp_set_favors(low = "Favors baseline", high = "Favors pRidge",
                  txt_gp = gpar(cex = 0.6))
