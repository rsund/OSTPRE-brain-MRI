# Statistical analyses for the OSTPRE-brain-MRI paper
# Reijo Sund 3/2025

library(dplyr)
library(ggplot2)

bids <- "P:/ostpre_bids" # Path to data

# Load data
qs::qload(file = file.path(bids, "clinical/mem_dates.qs")) # Diagnosis dates for Dementia, MCI and SMC
qs::qload(file = file.path(bids, "nifti/tags.qs")) # DICOM tags from the NIfTI converted scans

# Number of individuals after NIfTI conversion
tags |>
  distinct(lomno1) |>
  count()

# Number of separate examination visits after NIfTI conversion
tags |>
  distinct(lomno1, ses) |>
  count()

# Number of scans after NIfTI conversion
tags |> count()


# Extraction interesting parameters from the DICOM tags
reltags <- tags |>
  mutate(
    StudyID = as.character(StudyInstanceUID),
    SeriesID = as.character(SeriesInstanceUID),
    MF = as.character(Manufacturer),
    SN = as.character(DeviceSerialNumber),
    MagStre = as.character(MagneticFieldStrength)
  ) |>
  select(lomno1, ses, StudyID, SeriesID, MagneticFieldStrength, Manufacturer, ManufacturersModelName, DeviceSerialNumber, StationName, MF, SN, MagStre) |>
  distinct()

# Load CAT12 summary measures
imgdata_in <- readxl::read_excel(file.path(bids, "bids/derivatives/summary_measures/ostpre_CAT12_summary_measures.xlsx")) |>
  mutate(
    lomno1 = as.numeric(gsub("sub-", "", subj)),
    ses = as.numeric(gsub("ses-", "", session))
  )

# Results of the outlier analysis comparing parameters from CAT12 and SynthSeg
outliers <- readxl::read_excel(file.path(bids, "bids/derivatives/summary_measures/ostpre_outliers.xlsx"))

# How many outliers?
outliers |> count()

# How many individuals with outliers?
outliers |>
  distinct(subj) |>
  count()


# CAT12 summary measures (outliers excluded)
imgdata <- imgdata_in |>
  anti_join(outliers, by = c("subj", "session"))

# Number of scans for which there are CAT12 summary measures
imgdata |> count()

# Number of individuals with scans for which there are CAT12 summary measures
imgdata |>
  distinct(lomno1) |>
  count()


# Combine data for visits with detected T1-weighted scans
sessions <- readr::read_delim(file.path(bids, "bids/sessions.csv"), delim = ";", col_names = FALSE, show_col_types = FALSE) |>
  mutate(
    lomno1 = as.numeric(substr(X6, 4, 8)),
    ses = as.numeric(gsub("_", "", substr(X6, 14, 15))),
    sesdate = as.Date(X1)
  ) |>
  arrange(lomno1, ses) |>
  select(lomno1, ses, sesdate, StudyID = X2, SeriesID = X3, SeriesName = X4) |>
  left_join(imgdata, by = c("lomno1", "ses")) |>
  left_join(mem_dates |> rename(lomno1 = LOMNO1), by = c("lomno1")) |>
  left_join(reltags, by = c("lomno1", "SeriesID")) |>
  mutate(
    memtype = case_when(
      dem - 366 < sesdate ~ "Dementia",
      atc - 366 < sesdate ~ "ATC-Dementia",
      mci - 366 < sesdate ~ "MCI", # Mild Congnitive Impairment
      mem - 366 < sesdate ~ "SMC", # Subjective Memory Complaint
      dem - 366 < mci ~ "Pre MCI-Dementia",
      dem - 366 < mem ~ "Pre SMC-Dementia",
      is.na(dem) & !is.na(mci) ~ "Pre MCI",
      is.na(dem) & !is.na(mem) ~ "Pre SMC",
      !is.na(dem) ~ "Pre Dementia",
      is.na(dem) & is.na(mem) & is.na(atc) ~ "No memory concerns",
      is.na(dem) & is.na(mem) & !is.na(atc) ~ "Pre ATC-Dementia"
    ),
    Age = lubridate::interval(spvm, sesdate) / lubridate::dyears(1),
    mt = factor(case_when(
      memtype %in% c("Dementia", "ATC-Dementia") ~ "Dementia",
      grepl("^Pre|^No memory concerns", memtype) ~ "No memory concerns",
      TRUE ~ memtype
    ), levels = c("Dementia", "MCI", "SMC", "No memory concerns")),
    FS = case_when(
      MagStre == 3 ~ "3Tesla",
      TRUE ~ "1.5Tesla"
    )
  )

sessions |> count() # Number of examination visits with detected T1-weighted scans
sessions |>
  distinct(lomno1) |>
  count() # Number of individuals who had visit with detected T1-weighted scans


# Load relevant data for ADNI cohort
adni <- readxl::read_excel(file.path(bids, "bids/derivatives/summary_measures/ADNI_CAT12_ostprematched_summary_measures_v2.xlsx")) |>
  filter(DX_bl != "NA") |>
  mutate(
    mtyp = factor(case_when(
      DX_bl %in% c("AD") ~ "Dementia",
      DX_bl %in% c("EMCI", "LMCI") ~ "MCI",
      DX_bl %in% c("SMC") ~ "SMC",
      TRUE ~ "No memory concerns"
    ), levels = c("Dementia", "MCI", "SMC", "No memory concerns")),
    mt = factor(case_when(
      DX_bl %in% c("AD") ~ "Dementia - ADNI",
      DX_bl %in% c("EMCI", "LMCI") ~ "MCI - ADNI",
      DX_bl %in% c("SMC") ~ "SMC - ADNI",
      TRUE ~ "No memory concerns - ADNI"
    ), levels = c("Dementia - ADNI", "MCI - ADNI", "SMC - ADNI", "No memory concerns - ADNI")),
    MF = case_when(
      MANUFACTURER == 1 ~ "Siemens",
      MANUFACTURER == 2 ~ "Philips",
      MANUFACTURER == 3 ~ "GE",
    ),
    FS = case_when(
      grepl("^3", FLDSTRENG) ~ "3Tesla",
      TRUE ~ "1.5Tesla"
    )
  ) |>
  rename(adni_subj = subj) |>
  filter(!is.na(Age)) # Remove the participant with missing age

adni |>
  filter(quality_percent <= 68) |>
  count() # How many scans below the QC threshold in the ADNI data


# Combine ADNI and OSTPRE data and include only scans with quality score over 68% threshold
yhd <- adni |>
  bind_rows(sessions) |>
  filter(quality_percent > 68)


# Create Figure 2
yls <- yhd |>
  tidyr::pivot_longer(
    cols = quality_percent:ventricle
  ) |>
  select(value, mt, name) |>
  filter(name %in% c("hippocampus", "ventricle", "GM_volume"))

graph_titles <- c(
  GM_volume = "Gray Matter",
  hippocampus = "Hippocampus",
  ventricle = "Ventricle"
)

grob_adni <- grid::grobTree(grid::textGrob("ADNI", x = 0.17, y = 0.98, hjust = 0, gp = grid::gpar(col = "#525252", fontsize = 10, fontface = "italic")))
grob_ostpre <- grid::grobTree(grid::textGrob("OSTPRE", x = 0.6, y = 0.98, hjust = 0, gp = grid::gpar(col = "#525252", fontsize = 10, fontface = "italic")))

xpla <- c(0, 30, 55, 80)
polys <- tribble(
  ~group, ~x, ~y,
  "Dementia", xpla[1] + 0, 10,
  "Dementia", xpla[1] + 10, 0,
  "Dementia", xpla[1] + 10, 10,
  "Dementia - ADNI", xpla[1] + 0, 0,
  "Dementia - ADNI", xpla[1] + 10, 0,
  "Dementia - ADNI", xpla[1] + 0, 10,
  "MCI", xpla[2] + 0, 10,
  "MCI", xpla[2] + 10, 0,
  "MCI", xpla[2] + 10, 10,
  "MCI - ADNI", xpla[2] + 0, 0,
  "MCI - ADNI", xpla[2] + 10, 0,
  "MCI - ADNI", xpla[2] + 0, 10,
  "SMC", xpla[3] + 0, 10,
  "SMC", xpla[3] + 10, 0,
  "SMC", xpla[3] + 10, 10,
  "SMC - ADNI", xpla[3] + 0, 0,
  "SMC - ADNI", xpla[3] + 10, 0,
  "SMC - ADNI", xpla[3] + 0, 10,
  "No memory concerns", xpla[4] + 0, 10,
  "No memory concerns", xpla[4] + 10, 0,
  "No memory concerns", xpla[4] + 10, 10,
  "No memory concerns - ADNI", xpla[4] + 0, 0,
  "No memory concerns - ADNI", xpla[4] + 10, 0,
  "No memory concerns - ADNI", xpla[4] + 0, 10,
)

texts <- tribble(
  ~group, ~x, ~y,
  "Dementia", xpla[1] + 12, 5,
  "MCI", xpla[2] + 12, 5,
  "SMC", xpla[3] + 12, 5,
  "No memory complaints", xpla[4] + 12, 5,
)

viol_colors <- c(
  "Dementia" = "#CD5C5C",
  "MCI" = "#FFA07A",
  "SMC" = "#87CEFA",
  "No memory concerns" = "#98FB98",
  "Dementia - ADNI" = "#8B0000",
  "MCI - ADNI" = "#FF8C00",
  "SMC - ADNI" = "#4682B4",
  "No memory concerns - ADNI" = "#228B22"
)

leg <- ggplot(polys, aes(x = x, y = y)) +
  geom_polygon(aes(fill = group, group = group), show.legend = FALSE) +
  geom_text(data = texts, aes(label = group, x = x, y = y, hjust = "left")) +
  scale_fill_manual(
    name = "group",
    values = viol_colors
  ) +
  theme_void() +
  expand_limits(x = c(0, 140)) # , y=c(0, 200))

viol <- yls |>
  ggplot(
    aes(x = "", y = value)
  ) +
  geom_violin(aes(fill = mt), position = position_dodge(.9), show.legend = FALSE) +
  facet_wrap(~name, scale = "free", labeller = labeller(name = graph_titles)) +
  labs(x = "", y = "Volume, ml") +
  scale_fill_manual(
    name = "",
    values = viol_colors,
    breaks = c(
      "Dementia", "MCI", "SMC", "No memory concerns"
    )
  ) +
  theme(legend.position = "top") +
  theme(
    panel.grid.major.x = element_line(color = "blue", linewidth = 0.5, linetype = 2),
  ) +
  annotation_custom(grob_adni) +
  annotation_custom(grob_ostpre)

cowplot::ggdraw() +
  cowplot::draw_plot(leg, x = 0.15, y = 0.90, width = 0.85, height = 0.070) +
  cowplot::draw_plot(viol, x = 0, y = 0, height = 0.9)

ggsave("../images/Figure 2.pdf", width = 7, height = 5)


# Select and harmonize variables to be used in the statistical analyses
yhd_ana <- yhd |>
  mutate(
    subj = factor(case_when(
      is.na(subj) ~ sprintf("ADNI-%05d", adni_subj),
      TRUE ~ sub("sub", "OSTPRE", subj)
    )),
    gr = factor(case_when(
      grepl("^ADNI", subj) ~ "ADNI",
      TRUE ~ "OSTPRE"
    )),
    mtyp = case_when(
      gr == "OSTPRE" ~ mt,
      TRUE ~ mtyp
    ),
    mpf = factor(case_when(
      grepl("^Dementia", mtyp) ~ "3-Dementia",
      grepl("^MCI", mtyp) ~ "2-MCI",
      grepl("^SMC", mtyp) ~ "1-SMC",
      grepl("^No memory concerns", mtyp) ~ "0-No memory concerns"
    ), labels = c("No memory complaints", "SMC", "MCI", "Dementia")),
    cage = Age - 75,
    MF = factor(MF),
    FS = factor(FS)
  ) |>
  select(gr, subj, session, Age, cage, mpf, MF, FS, TIV:ventricle)


# Data for Supplementary Figure 2
memdist_ana <- yhd_ana |>
  mutate(age = floor(Age)) |>
  count(gr, age, mpf) |>
  tidyr::pivot_wider(
    id_cols = c(gr, age),
    names_from = "mpf",
    values_from = "n",
    values_fill = 0
  )
# openxlsx::write.xlsx(memdist_ana,"../memdist_ana.xlsx")


# Data for Table 1

# N and age
yhd_ana |>
  group_by(gr) |>
  summarise(
    n = n(),
    nsub = n_distinct(subj),
    age = mean(Age, na.rm = TRUE),
    agesd = sd(Age, na.rm = TRUE),
  )

# Field Strength
yhd_ana |>
  count(gr, FS) |>
  group_by(gr) |>
  mutate(np = n / sum(n) * 100)

# Manufacturer
yhd_ana |>
  count(gr, MF) |>
  group_by(gr) |>
  mutate(np = n / sum(n) * 100)

# Cognitive Status
yhd_ana |>
  count(gr, mpf) |>
  group_by(gr) |>
  mutate(np = n / sum(n) * 100)


# Estimate the models

library(ggeffects)
library(marginaleffects)
library(clubSandwich)

melist <- c("hippocampus", "ventricle", "GM_volume")

# Combined cohort
mdpl <- list()
for (i in 1:length(melist)) {
  me <- melist[i]
  md <- yhd_ana |>
    select(gr, subj, Age, cage, mpf, MF, FS, TIV, meas = {{ me }}) |>
    filter(!is.na(meas) & !is.na(Age)) |>
    mutate(
      smeas = scale(meas)[, 1],
      nmeas = meas - mean(meas),
      sTIV = scale(TIV)[, 1]
    )

  m1 <- lm(smeas ~ splines::bs(Age) + splines::bs(sTIV) + gr * mpf + MF + FS, data = md)

  mdp1 <- ggeffects::predict_response(m1, terms = c("Age [65:90]"), vcov_fun = "vcovCR", vcov_type = "CR0", vcov_args = list(cluster = md$subj))
  mdp2 <- ggeffects::predict_response(m1, terms = c("sTIV [-3,-2.75,-2.5,-2.25,-2,-1.75,-1.5,-1.25,-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5,2.75,3]"), vcov_fun = "vcovCR", vcov_type = "CR0", vcov_args = list(cluster = md$subj))
  mdp <- ggeffects::predict_response(m1, terms = c("mpf", "gr"), vcov_fun = "vcovCR", vcov_type = "CR0", vcov_args = list(cluster = md$subj))

  mdpl[[me]] <- list(m1 = m1, mdp1 = mdp1, mdp2 = mdp2, mdp = mdp)
}

# Cohorts separately
mdpl_sep <- list()
for (i in 1:length(melist)) {
  me <- melist[i]
  md <- yhd_ana |>
    select(gr, subj, Age, cage, mpf, MF, FS, TIV, meas = {{ me }}) |>
    filter(!is.na(meas) & !is.na(Age)) |>
    mutate(
      smeas = scale(meas)[, 1],
      nmeas = meas - mean(meas),
      sTIV = scale(TIV)[, 1]
    )

  md1 <- md |> filter(gr == "ADNI")
  md2 <- md |> filter(gr == "OSTPRE")

  m1 <- lm(smeas ~ splines::bs(Age) + splines::bs(sTIV) + mpf + MF + FS, data = md1)
  m2 <- lm(smeas ~ splines::bs(Age) + splines::bs(sTIV) + mpf + MF + FS, data = md2)

  mdp1a <- ggeffects::predict_response(m1, terms = c("Age [65:89]"), vcov_fun = "vcovCR", vcov_type = "CR0", vcov_args = list(cluster = md1$subj))
  mdp1b <- ggeffects::predict_response(m2, terms = c("Age [65:89]"), vcov_fun = "vcovCR", vcov_type = "CR0", vcov_args = list(cluster = md2$subj))
  mdp2a <- ggeffects::predict_response(m1, terms = c("sTIV [-2.5,-2.25,-2,-1.75,-1.5,-1.25,-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5,2.75,3]"), vcov_fun = "vcovCR", vcov_type = "CR0", vcov_args = list(cluster = md1$subj))
  mdp2b <- ggeffects::predict_response(m2, terms = c("sTIV [-2.5,-2.25,-2,-1.75,-1.5,-1.25,-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5,2.75,3]"), vcov_fun = "vcovCR", vcov_type = "CR0", vcov_args = list(cluster = md2$subj))

  mdpl_sep[[me]] <- list(m1 = m1, m2 = m2, mdp1a = mdp1a, mdp1b = mdp1b, mdp2a = mdp2a, mdp2b = mdp2b)
}


# Create Figure 3

mdp_1 <- as.data.frame(mdpl[["GM_volume"]]$mdp) |> mutate(para = "Gray Matter Volume")
mdp_2 <- as.data.frame(mdpl[["hippocampus"]]$mdp) |> mutate(para = "Hippocampus Volume")
mdp_3 <- as.data.frame(mdpl[["ventricle"]]$mdp) |> mutate(para = "Ventricle Volume")

mdp_panel <- mdp_1 |>
  bind_rows(mdp_2) |>
  bind_rows(mdp_3) |>
  mutate(
    status = factor(case_when(
      x == "No memory complaints" ~ "NMC",
      TRUE ~ x
    ), levels = c("NMC", "SMC", "MCI", "Dementia"))
  )

named_colors <- c(
  "NMC ADNI" = "#228B22",
  "NMC OSTPRE" = "#98FB98",
  "SMC ADNI" = "#4682B4",
  "SMC OSTPRE" = "#87CEFA",
  "MCI ADNI" = "#FF8C00",
  "MCI OSTPRE" = "#FFA07A",
  "Dementia ADNI" = "#8B0000",
  "Dementia OSTPRE" = "#CD5C5C"
)

named_shapes <- c(
  "NMC ADNI" = 21,
  "NMC OSTPRE" = 23,
  "SMC ADNI" = 21,
  "SMC OSTPRE" = 23,
  "MCI ADNI" = 21,
  "MCI OSTPRE" = 23,
  "Dementia ADNI" = 21,
  "Dementia OSTPRE" = 23
)

mdp_panel |>
  mutate(
    sg = paste(status, group)
  ) |>
  ggplot(
    aes(
      y = predicted,
      x = status,
      ymin = conf.low,
      ymax = conf.high,
      color = sg,
      shape = group,
      fill = sg,
      group = sg
    )
  ) +
  facet_wrap(
    facets = ~para
  ) +
  geom_linerange(
    linewidth = 1,
    position = position_dodge(width = 0.3)
  ) +
  geom_point(
    size = 3,
    shape = rep(named_shapes, times = 3),
    color = "white",
    stroke = 0.5,
    position = position_dodge(width = 0.3)
  ) +
  scale_color_manual(
    values = named_colors,
    limits = c("NMC ADNI", "NMC OSTPRE", "SMC ADNI", "SMC OSTPRE", "MCI ADNI", "MCI OSTPRE", "Dementia ADNI", "Dementia OSTPRE")
  ) +
  scale_fill_manual(
    values = named_colors,
    limits = c("NMC ADNI", "NMC OSTPRE", "SMC ADNI", "SMC OSTPRE", "MCI ADNI", "MCI OSTPRE", "Dementia ADNI", "Dementia OSTPRE")
  ) +
  scale_x_discrete(
    name = ""
  ) +
  scale_y_continuous(
    name = "Standardized measurement",
    limits = c(-1.02, 0.7)
  ) +
  theme(
    legend.position = "top",
    legend.title = element_blank()
  ) +
  guides(
    color = guide_legend(
      override.aes = list(
        shape = named_shapes,
        fill = named_colors
      )
    )
  )

ggsave("../images/Figure 3.pdf", height = 7, width = 7)


# Table 2

# Grey matter volume
m1p <- ggeffects::predict_response(mdpl[["GM_volume"]]$m1, terms = c("mpf", "gr"), vcov_fun = "vcovCR", vcov_type = "CR0", vcov_args = list(cluster = yhd_ana$subj))
ggeffects::test_predictions(m1p)
ggeffects::test_predictions(
  m1p,
  test = c(
    "(b1 - b2) = (b5 - b6)", # None vs SMC difference between ADNI and OSTPRE
    "(b1 - b3) = (b5 - b7)", # None vs MCI difference between ADNI and OSTPRE
    "(b1 - b4) = (b5 - b8)", # None vs Dementia difference between ADNI and OSTPRE
    "(b3 - b4) = (b7 - b8)" # MCI vs Dementia difference between ADNI and OSTPRE
  )
)

# Hippocampus volume
m2p <- ggeffects::predict_response(mdpl[["hippocampus"]]$m1, terms = c("mpf", "gr"), vcov_fun = "vcovCR", vcov_type = "CR0", vcov_args = list(cluster = yhd_ana$subj))
ggeffects::test_predictions(m2p)
ggeffects::test_predictions(
  m2p,
  test = c(
    "(b1 - b2) = (b5 - b6)", # None vs SMC difference between ADNI and OSTPRE
    "(b1 - b3) = (b5 - b7)", # None vs MCI difference between ADNI and OSTPRE
    "(b1 - b4) = (b5 - b8)", # None vs Dementia difference between ADNI and OSTPRE
    "(b3 - b4) = (b7 - b8)" # MCI vs Dementia difference between ADNI and OSTPRE
  )
)

# Ventricle volume
m3p <- ggeffects::predict_response(mdpl[["ventricle"]]$m1, terms = c("mpf", "gr"), vcov_fun = "vcovCR", vcov_type = "CR0", vcov_args = list(cluster = yhd_ana$subj))
ggeffects::test_predictions(m3p)
ggeffects::test_predictions(
  m3p,
  test = c(
    "(b1 - b2) = (b5 - b6)", # None vs SMC difference between ADNI and OSTPRE
    "(b1 - b3) = (b5 - b7)", # None vs MCI difference between ADNI and OSTPRE
    "(b1 - b4) = (b5 - b8)", # None vs Dementia difference between ADNI and OSTPRE
    "(b3 - b4) = (b7 - b8)" # MCI vs Dementia difference between ADNI and OSTPRE
  )
)


# Create Figure 4

mod_age <- NULL |>
  bind_rows(as_tibble(c(mdpl[["GM_volume"]]$mdp1, gr = "Gray Matter Volume"))) |>
  bind_rows(as_tibble(c(mdpl[["hippocampus"]]$mdp1, gr = "Hippocampus Volume"))) |>
  bind_rows(as_tibble(c(mdpl[["ventricle"]]$mdp1, gr = "Ventricle Volume"))) |>
  mutate(para = "Age")

mod_tiv <- NULL |>
  bind_rows(as_tibble(c(mdpl[["GM_volume"]]$mdp2, gr = "Gray Matter Volume"))) |>
  bind_rows(as_tibble(c(mdpl[["hippocampus"]]$mdp2, gr = "Hippocampus Volume"))) |>
  bind_rows(as_tibble(c(mdpl[["ventricle"]]$mdp2, gr = "Ventricle Volume"))) |>
  mutate(para = "Standardised TIV")

mod_yhd <- mod_age |>
  bind_rows(mod_tiv)

mod_yhd |>
  ggplot(aes(x = x, y = predicted, colour = gr)) +
  facet_wrap(~para, scale = "free_x") +
  geom_line() +
  geom_ribbon(data = mod_yhd, aes(ymin = conf.low, ymax = conf.high), alpha = 0.1) +
  theme(legend.position = "top", legend.title = element_blank()) +
  labs(x = "", y = "Standardized measurement")

ggsave("../images/Figure 4.pdf", width = 5, height = 5)


# Create Supplementary Figure 4
mod_age_sep <- NULL |>
  bind_rows(as_tibble(c(mdpl_sep[["GM_volume"]]$mdp1a, gr = "Gray Matter Volume", grp = "ADNI"))) |>
  bind_rows(as_tibble(c(mdpl_sep[["GM_volume"]]$mdp1b, gr = "Gray Matter Volume", grp = "OSTPRE"))) |>
  bind_rows(as_tibble(c(mdpl_sep[["hippocampus"]]$mdp1a, gr = "Hippocampus Volume", grp = "ADNI"))) |>
  bind_rows(as_tibble(c(mdpl_sep[["hippocampus"]]$mdp1b, gr = "Hippocampus Volume", grp = "OSTPRE"))) |>
  bind_rows(as_tibble(c(mdpl_sep[["ventricle"]]$mdp1a, gr = "Ventricle Volume", grp = "ADNI"))) |>
  bind_rows(as_tibble(c(mdpl_sep[["ventricle"]]$mdp1b, gr = "Ventricle Volume", grp = "OSTPRE"))) |>
  mutate(para = "Age")

mod_tiv_sep <- NULL |>
  bind_rows(as_tibble(c(mdpl_sep[["GM_volume"]]$mdp2a, gr = "Gray Matter Volume", grp = "ADNI"))) |>
  bind_rows(as_tibble(c(mdpl_sep[["GM_volume"]]$mdp2b, gr = "Gray Matter Volume", grp = "OSTPRE"))) |>
  bind_rows(as_tibble(c(mdpl_sep[["hippocampus"]]$mdp2a, gr = "Hippocampus Volume", grp = "ADNI"))) |>
  bind_rows(as_tibble(c(mdpl_sep[["hippocampus"]]$mdp2b, gr = "Hippocampus Volume", grp = "OSTPRE"))) |>
  bind_rows(as_tibble(c(mdpl_sep[["ventricle"]]$mdp2a, gr = "Ventricle Volume", grp = "ADNI"))) |>
  bind_rows(as_tibble(c(mdpl_sep[["ventricle"]]$mdp2b, gr = "Ventricle Volume", grp = "OSTPRE"))) |>
  mutate(para = "Standardised TIV")

mod_yhd_sep <- mod_age_sep |>
  bind_rows(mod_tiv_sep)

mod_yhd_sep |>
  ggplot(aes(x = x, y = predicted, color = gr)) +
  facet_wrap(~ para + grp, scale = "free_x") +
  geom_line() +
  geom_ribbon(data = mod_yhd_sep, aes(ymin = conf.low, ymax = conf.high), alpha = 0.1) +
  theme(legend.position = "top", legend.title = element_blank()) +
  labs(x = "", y = "Standardized measurement")

ggsave("../images/Supplementary Figure 4.svg", width = 10, height = 15)
