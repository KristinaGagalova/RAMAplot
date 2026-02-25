#!/usr/bin/env R

## Author: Kristina K. Gagalova

library(dplyr)
library(ggplot2)

# Read one RAMAplot reference file (phi psi density) from .bz2
read_rama_ref_long <- function(path) {
  stopifnot(file.exists(path))
  con <- bzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  
  dat <- read.table(con, header = TRUE, comment.char = "#", sep = "",
                    stringsAsFactors = FALSE)
  names(dat) <- tolower(names(dat))
  
  if (!all(c("phi", "psi", "density") %in% names(dat))) {
    stop("Expected columns phi, psi, density; got: ", paste(names(dat), collapse = ", "))
  }
  
  dat %>%
    transmute(
      phi = as.numeric(phi),
      psi = as.numeric(psi),
      density = as.numeric(density)
    )
}

# Apply the SAME transform RAMAplot uses:
# pivot(phi x psi) -> transpose -> vertical flip ([::-1])
to_rama_imshow_df <- function(dat_long) {
  phi_vals <- sort(unique(dat_long$phi))
  psi_vals <- sort(unique(dat_long$psi))
  
  # Build matrix M[phi, psi]
  M <- matrix(NA_real_, nrow = length(phi_vals), ncol = length(psi_vals),
              dimnames = list(phi = phi_vals, psi = psi_vals))
  i_phi <- match(dat_long$phi, phi_vals)
  i_psi <- match(dat_long$psi, psi_vals)
  M[cbind(i_phi, i_psi)] <- dat_long$density
  
  # RAMAplot: rama_ref.transpose()
  H <- t(M)  # now rows=psi, cols=phi
  
  # RAMAplot: imshow(ref_obj.histo2d[::-1]) (flip vertically)
  H <- H[nrow(H):1, , drop = FALSE]
  
  # Convert back to a data.frame for ggplot (x=phi, y=psi)
  psi_plot <- as.numeric(rownames(H))
  phi_plot <- as.numeric(colnames(H))
  
  expand.grid(psi = psi_plot, phi = phi_plot) %>%
    mutate(density = as.vector(H)) %>%
    transmute(phi = phi, psi = psi, density = density)
}

# RAMAplot bounds + colors
rama_styles <- list(
  Gen  = list(bounds = c(0, 0.0005, 0.02, 1), cols = c("#FFFFFF", "#B3E8FF", "#26C3FF")),
  Gly  = list(bounds = c(0, 0.0020, 0.02, 1), cols = c("#FFFFFF", "#FFE8C5", "#FED479")),
  PreP = list(bounds = c(0, 0.0020, 0.02, 1), cols = c("#FFFFFF", "#FFD7DA", "#F3ABB0")),
  Pro  = list(bounds = c(0, 0.0020, 0.02, 1), cols = c("#FFFFFF", "#D0FFC5", "#7EE589"))
)

# Turn density into 3-level factor using the same bounds
bin_density <- function(density, bounds) {
  cut(
    density,
    breaks = bounds,
    include.lowest = TRUE,
    right = TRUE,
    labels = c("low", "mid", "high")
  )
}

# Load all four from user paths and rename panel labels
load_rama_refs_from_paths <- function(gen_path, gly_path, prepro_path, pro_path) {
  paths <- list(Gen = gen_path, Gly = gly_path, PreP = prepro_path, Pro = pro_path)
  
  bind_rows(lapply(names(paths), function(k) {
    dat <- read_rama_ref_long(paths[[k]])
    imgdf <- to_rama_imshow_df(dat)
    style <- rama_styles[[k]]
    
    imgdf %>%
      mutate(
        panel = recode(k,
                       "Gen"  = "General",
                       "Gly"  = "Glycine",
                       "PreP" = "Pre-proline",
                       "Pro"  = "Proline"),
        level = bin_density(density, style$bounds)
      )
  }))
}

# Plot (background only) in a 2x2 layout
plot_rama_backgrounds <- function(bg_df) {
  
  # Control facet order
  bg_df$panel <- factor(
    bg_df$panel,
    levels = c("General", "Glycin", "Pre-prolin", "Prolin")
  )
  
  # Key for per-panel fills
  bg_df <- bg_df %>%
    mutate(key = paste(panel, level, sep = ":"))
  
  # Fill map updated to renamed panels
  fill_map <- c(
    "General:low"      = rama_styles$Gen$cols[1],
    "General:mid"      = rama_styles$Gen$cols[2],
    "General:high"     = rama_styles$Gen$cols[3],
    
    "Glycine:low"       = rama_styles$Gly$cols[1],
    "Glycine:mid"       = rama_styles$Gly$cols[2],
    "Glycine:high"      = rama_styles$Gly$cols[3],
    
    "Pre-proline:low"   = rama_styles$PreP$cols[1],
    "Pre-proline:mid"   = rama_styles$PreP$cols[2],
    "Pre-proline:high"  = rama_styles$PreP$cols[3],
    
    "Proline:low"       = rama_styles$Pro$cols[1],
    "Proline:mid"       = rama_styles$Pro$cols[2],
    "Proline:high"      = rama_styles$Pro$cols[3]
  )
  
  ggplot(bg_df, aes(phi, psi)) +
    geom_raster(aes(fill = key), interpolate = FALSE) +
    scale_fill_manual(values = fill_map, guide = "none") +
    facet_wrap(~panel, nrow = 2) +
    coord_fixed(xlim = c(-180, 180), ylim = c(-180, 180), expand = FALSE) +
    
    # slightly thicker axis lines
    geom_hline(yintercept = 0, linewidth = 0.6) +
    geom_vline(xintercept = 0, linewidth = 0.6) +
    
    scale_x_continuous(breaks = seq(-180, 180, 60)) +
    scale_y_continuous(breaks = seq(-180, 180, 60)) +
    
    labs(x = expression(phi), y = expression(psi)) +
    
    theme_classic(base_size = 16) +   # increase global size
    
    theme(
      # Bigger axis tick labels
      axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 14),
      
      # Bigger axis titles
      axis.title.x = element_text(size = 25, face = "bold"),
      axis.title.y = element_text(size = 25, face = "bold"),
      
      # Bigger facet labels
      strip.text = element_text(size = 18, face = "bold"),
      
      # More space between panels
      panel.spacing = unit(1.8, "lines"),
      
      # ← THIS draws the box around each facet
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    )
}

# ---- Run ----
bg <- load_rama_refs_from_paths(
  gen_path    = "pyrama_data/pref_general.data.bz2",
  gly_path    = "pyrama_data/pref_glycine.data.bz2",
  prepro_path = "pyrama_data/pref_preproline.data.bz2",
  pro_path    = "pyrama_data/pref_proline.data.bz2"
)

plot_rama_backgrounds(bg)
