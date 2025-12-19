

library(tidyverse)
library(stringr)

parse_cjs_txt <- function(path) {
  stopifnot(file.exists(path))
  txt <- readLines(path)
  fn  <- basename(path)
  
  # --- 0) Parse metadata from filename: YYYY_species_hatchery_releasegroup.ext ---
  m <- str_match(fn, "^(\\d{4})_([A-Za-z]+)_([A-Za-z0-9]+)_([A-Za-z0-9]+)\\.[A-Za-z]+$")
  if (any(is.na(m))) {
    stop("Filename does not match 'YYYY_species_hatchery_releasegroup.ext' pattern: ", fn)
  }
  year          <- as.integer(m[2])
  species       <- m[3]
  hatchery      <- m[4]
  release_group <- m[5]
  
  # --- helpers ---
  extract_table <- function(lines, marker) {
    i <- str_which(lines, fixed(marker))
    if (length(i) == 0) stop(sprintf("Section not found: %s", marker))
    tibble(header = lines[i + 1], data = lines[i + 2])
  }
  
  parse_est_se <- function(header_line, data_line) {
    h <- str_split(header_line, "\t")[[1]]
    d <- str_split(data_line,   "\t")[[1]]
    
    cols <- h[-1]                       # drop "Population"
    vals <- d[-1]                       # drop "webrequest"
    
    mat <- str_match(vals, "^\\s*([0-9.]+)\\s*\\(([^)]+)\\)\\s*$")
    tibble(
      name     = cols,
      estimate = as.numeric(mat[, 2]),
      se       = as.numeric(mat[, 3])
    )
  }
  
  # --- 1) Survival Estimates ---
  surv_raw <- extract_table(txt, "Survival Estimates")
  survival_tbl <- parse_est_se(surv_raw$header, surv_raw$data) |>
    mutate(
      metric = if_else(str_detect(name, "^Reach \\d+$"), "reach", "overall"),
      reach  = if_else(metric == "reach", as.integer(str_match(name, "Reach (\\d+)")[,2]), NA_integer_)
    )
  
  # legend labels for reaches
  legend_start_reach <- str_which(txt, fixed("Survival Reaches:"))
  reach_labels <- if (length(legend_start_reach)) {
    lines <- txt[(legend_start_reach + 1):(legend_start_reach + 7)]
    tibble(
      reach       = as.integer(str_match(lines, "^\\s*(\\d+):")[,2]),
      reach_label = str_trim(str_match(lines, "^\\s*\\d+:\\s*(.+)$")[,2])
    )
  } else {
    tibble(reach = integer(), reach_label = character())
  }
  
  survival_tbl <- survival_tbl |>
    left_join(reach_labels, by = "reach") |>
    transmute(
      year, species, hatchery, release_group,
      metric, reach, reach_label, estimate, se
    )
  
  # --- 2) Capture Probability Estimates ---
  cap_raw <- extract_table(txt, "Capture Probability Estimates")
  capture_tbl <- parse_est_se(cap_raw$header, cap_raw$data) |>
    mutate(
      metric = if_else(str_detect(name, "^Site \\d+$"), "site", "final"),
      site   = if_else(metric == "site", as.integer(str_match(name, "Site (\\d+)")[,2]), NA_integer_)
    )
  
  # legend labels for sites
  legend_start_site <- str_which(txt, fixed("Capture Sites:"))
  site_labels <- if (length(legend_start_site)) {
    lines <- txt[(legend_start_site + 1):(legend_start_site + 7)]
    tibble(
      site       = as.integer(str_match(lines, "^\\s*(\\d+):")[,2]),
      site_label = str_trim(str_match(lines, "^\\s*\\d+:\\s*(.+)$")[,2])
    )
  } else {
    tibble(site = integer(), site_label = character())
  }
  
  capture_tbl <- capture_tbl |>
    left_join(site_labels, by = "site") |>
    transmute(
      year, species, hatchery, release_group,
      metric, site, site_label, estimate, se
    )
  
  list(
    survival = survival_tbl,
    capture  = capture_tbl
  )
}


out <- parse_cjs_txt("outputs/dart_outputs/2025_chinook_CLWH_CLEARC.txt")

dart.files <- c("outputs/dart_outputs/2025_chinook_CLWH_CLEARC.txt",
                "outputs/dart_outputs/2025_sthd_DNFH_CLEARC.txt",
                "outputs/dart_outputs/2025_sthd_DNFH_CW.txt")



survival_all <- map_dfr(dart.files, ~ parse_cjs_txt(.x)$survival, .id = "file_index") |>
  mutate(source_file = dart.files[file_index]) |>
  select(-file_index)



capture_all  <- map_dfr(dart.files, ~ parse_cjs_txt(.x)$capture, .id = "file_index") |>
  mutate(source_file = dart.files[file_index]) |> 
  select(-file_index)
  

# key to translate this to RMark and Nimble summary
# format

survival_key <- tibble(site_label=c("Release-LGR","LGR-LGS","LGS-LMN","LMN-IHR",
                                   "IHR-MCN","MCN-JDA","JDA-BON","overall_phi"),
                      descriptor=c("release_lgr_phi","lgr_lgs_phi","lgs_lomo_phi","lomo_iha_phi",
                                   "iha_mcn_phi","mnc_jda_phi","jda_bon_phi","overall_phi"))


survival_estimates <- survival_all |> 
  select(-c(metric,reach)) |>
  rename(site_label=reach_label) |> 
  mutate(site_label=replace_na(site_label,"overall_phi")) |> 
  left_join(survival_key,by="site_label") |> 
  mutate(category="phi",
         species=ifelse(species=="chinook",
                        "SCS","SST"),
         timing_group=NA,
         lcl=as.numeric(NA),
         ucl=as.numeric(NA),
         estimate_source="DART") |> 
  select(hatchery,species,release_year=year,timing_group,
         release_group,descriptor,category,estimate,se,
         lcl,ucl,estimate_source)
  

# key to translate this to RMark and Nimble summary
# format

capture_key <- tibble(site_label=c("LGR","LGS","LMN","IHR",
                                   "MCN","JDA","BON","final_joint"),
                      descriptor=c("lgr_p","lgs_p","lomo_p","iha_p",
                                   "mcn_p","jda_p","bon_p","final_joint"))

capture_estimates  <- capture_all |> 
  select(-site) |> 
  mutate(site_label=replace_na(site_label,"final_joint")) |> 
  left_join(capture_key,by="site_label") |> 
  mutate(category=ifelse(descriptor=="final_joint",
                         "final","p"),
         species=ifelse(species=="chinook",
                        "SCS","SST"),
         timing_group=NA,
         lcl=as.numeric(NA),
         ucl=as.numeric(NA),
         estimate_source="DART") |> 
  select(hatchery,species,release_year=year,timing_group,
         release_group,descriptor,category,estimate,se,
         lcl,ucl,estimate_source)

estimates.bind <- survival_estimates |> 
  bind_rows(capture_estimates)
