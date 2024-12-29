library(readxl)
library(tidyverse)

sourceFilePath <- "data/Mandai_data_20231204.xlsx"

years <- c(2011, 2012, 2013, 2014, 2015)

# Generate trees object ----

trees <- read_excel(sourceFilePath,
                   sheet = "tree") %>%
  mutate_at(.vars = c("DBH (2011)",
                      "DBH (2012)",
                      "DBH (2013)",
                      "DBH (2014)",
                      "DBH (2015)"),
            .funs = list(~ case_match(.,
                                      "can't find" ~ "cnf",
                                      "cnf" ~ "cnf",
                                      "didn't measure" ~ "missing",
                                      "dead" ~ "dead",
                                      "<1" ~ "0.9",
                                      .default = as.character(round(as.numeric(.), 1))
                                      )
                         )
            ) %>%
  rename(plotID = `Plot No.`,
         treeID = Tag,
         species = Species)

trees <- lapply(years, function(x) {
  trees %>%
    select(plotID, treeID, species,
           paste0("DBH (", x, ")")) %>%
    rename(DBH = paste0("DBH (", x, ")")) %>%
    filter(!is.na(DBH))
})
names(trees) <- years

# Generate plot object ----

plots_temp <- excel_sheets(sourceFilePath)[c(1:5)] %>%
  set_names() %>%
  map(read_excel,
      path = sourceFilePath,
      col_names = FALSE) %>%
  lapply(function(x) {
    names(x) <- apply(x[1:2,], 2, function(x) {
      paste(x, collapse = "-") %>%
        str_replace_all("NA-", "") %>%
        str_replace_all("-NA", "")
    })
    
    x[-c(1:2),]
  })
names(plots_temp) <- years

plots <- lapply(plots_temp, function(x) {
    x %>%
      select(`Plot No.`, Date, `Map forest type`, `Map damage type`, N, P, K) %>%
      rename(plotID = `Plot No.`,
             date = Date,
             type_forest = `Map forest type`,
             type_damage = `Map damage type`) %>%
      mutate(N = as.numeric(N)/1000,
             P = as.numeric(P),
             K = as.numeric(K)) %>%
      arrange(plotID)
  })
plots$`2011`$date <- as.Date(as.numeric(plots$`2011`$date), origin = "1899-12-30")
plots$`2012`$date <- as.Date(as.numeric(plots$`2012`$date), origin = "1899-12-30")
plots$`2013`$date <- as.Date(plots$`2013`$date, format = "%d/%m/%Y")
plots$`2014`$date <- as.Date(plots$`2014`$date, format = "%d/%m/%Y")
plots$`2015`$date <- as.Date(plots$`2015`$date, format = "%d/%m/%Y")


# Generae quadrants object ----

quadrants <- vector("list", length(years))
names(quadrants) <- years

quadrants$`2011` <- plots_temp$`2011`[,c(1,8:11)] %>%
  rename(plotID = `Plot No.`,
         NE = `Leaf Litter-NE`) %>%
  pivot_longer(cols = -plotID,
               names_to = "quadrant",
               values_to = "leafLitterDepth") %>%
  full_join(
    plots_temp$`2011`[,c(1,12:15)] %>%
      rename(plotID = `Plot No.`,
             NE = `Canopy-NE`) %>%
      pivot_longer(cols = -plotID,
                   names_to = "quadrant",
                   values_to = "canopyCover") %>%
      mutate(canopyCover = as.numeric(canopyCover)/96*100),
    by = c("plotID", "quadrant")
  ) %>%
  mutate(quadrantID = paste0(plotID, quadrant)) %>%
  select(quadrantID, plotID, leafLitterDepth, canopyCover)

quadrants$`2012` <- plots_temp$`2012`[,c(1,8:11)] %>%
  rename(plotID = `Plot No.`,
         NE = `Leaf Litter-NE`) %>%
  pivot_longer(cols = -plotID,
               names_to = "quadrant",
               values_to = "leafLitterDepth") %>%
  full_join(
    plots_temp$`2012`[,c(1,12:27)] %>%
      rename(plotID = `Plot No.`,
             NE1 = `Canopy (dark)-NE1`) %>%
      mutate_at(-1, as.numeric) %>%
      mutate(NE = rowMeans(select(., starts_with("NE")), na.rm = TRUE),
             NW = rowMeans(select(., starts_with("NW")), na.rm = TRUE),
             SE = rowMeans(select(., starts_with("SE")), na.rm = TRUE),
             SW = rowMeans(select(., starts_with("SW")), na.rm = TRUE)) %>%
      select(plotID, NE, NW, SE, SW) %>%
      pivot_longer(cols = -plotID,
                   names_to = "quadrant",
                   values_to = "canopyCover"),
    by = c("plotID", "quadrant")
  ) %>%
  mutate(quadrantID = paste0(plotID, quadrant)) %>%
  select(quadrantID, plotID, leafLitterDepth, canopyCover)

quadrants$`2013` <- plots_temp$`2013`[,c(1,8:23)] %>%
  rename(plotID = `Plot No.`,
         NE = `Leaf Litter-NE1`) %>%
  mutate_at(-1, as.numeric) %>%
  mutate(NE = rowMeans(select(., starts_with("NE")), na.rm = TRUE),
         NW = rowMeans(select(., starts_with("NW")), na.rm = TRUE),
         SE = rowMeans(select(., starts_with("SE")), na.rm = TRUE),
         SW = rowMeans(select(., starts_with("SW")), na.rm = TRUE)) %>%
  select(plotID, NE, NW, SE, SW) %>%
  pivot_longer(cols = -plotID,
               names_to = "quadrant",
               values_to = "leafLitterDepth") %>%
  full_join(
    plots_temp$`2013`[,c(1,24:39)] %>%
      rename(plotID = `Plot No.`,
             NE1 = `Canopy (dark)-NE1`) %>%
      mutate_at(-1, as.numeric) %>%
      mutate(NE = rowMeans(select(., starts_with("NE")), na.rm = TRUE),
             NW = rowMeans(select(., starts_with("NW")), na.rm = TRUE),
             SE = rowMeans(select(., starts_with("SE")), na.rm = TRUE),
             SW = rowMeans(select(., starts_with("SW")), na.rm = TRUE)) %>%
      select(plotID, NE, NW, SE, SW) %>%
      pivot_longer(cols = -plotID,
                   names_to = "quadrant",
                   values_to = "canopyCover"),
    by = c("plotID", "quadrant")
  ) %>%
  mutate(quadrantID = paste0(plotID, quadrant)) %>%
  select(quadrantID, plotID, leafLitterDepth, canopyCover)

quadrants$`2014` <- plots_temp$`2014`[,c(1,8:23)] %>%
  rename(plotID = `Plot No.`,
         NE = `Leaf Litter-NE1`) %>%
  mutate_at(-1, as.numeric) %>%
  mutate(NE = rowMeans(select(., starts_with("NE")), na.rm = TRUE),
         NW = rowMeans(select(., starts_with("NW")), na.rm = TRUE),
         SE = rowMeans(select(., starts_with("SE")), na.rm = TRUE),
         SW = rowMeans(select(., starts_with("SW")), na.rm = TRUE)) %>%
  select(plotID, NE, NW, SE, SW) %>%
  pivot_longer(cols = -plotID,
               names_to = "quadrant",
               values_to = "leafLitterDepth") %>%
  full_join(
    plots_temp$`2014`[,c(1,24:39)] %>%
      rename(plotID = `Plot No.`,
             NE1 = `Canopy (dark)-NE1`) %>%
      mutate_at(-1, as.numeric) %>%
      mutate(NE = rowMeans(select(., starts_with("NE")), na.rm = TRUE),
             NW = rowMeans(select(., starts_with("NW")), na.rm = TRUE),
             SE = rowMeans(select(., starts_with("SE")), na.rm = TRUE),
             SW = rowMeans(select(., starts_with("SW")), na.rm = TRUE)) %>%
      select(plotID, NE, NW, SE, SW) %>%
      pivot_longer(cols = -plotID,
                   names_to = "quadrant",
                   values_to = "canopyCover"),
    by = c("plotID", "quadrant")
  ) %>%
  mutate(quadrantID = paste0(plotID, quadrant)) %>%
  select(quadrantID, plotID, leafLitterDepth, canopyCover)

quadrants$`2015` <- plots_temp$`2015`[,c(1,8:23)] %>%
  rename(plotID = `Plot No.`,
         NE = `Leaf Litter-NE1`) %>%
  mutate_at(-1, as.numeric) %>%
  mutate(NE = rowMeans(select(., starts_with("NE")), na.rm = TRUE),
         NW = rowMeans(select(., starts_with("NW")), na.rm = TRUE),
         SE = rowMeans(select(., starts_with("SE")), na.rm = TRUE),
         SW = rowMeans(select(., starts_with("SW")), na.rm = TRUE)) %>%
  select(plotID, NE, NW, SE, SW) %>%
  pivot_longer(cols = -plotID,
               names_to = "quadrant",
               values_to = "leafLitterDepth") %>%
  full_join(
    plots_temp$`2015`[,c(1,24:39)] %>%
      rename(plotID = `Plot No.`,
             NE1 = `Canopy (dark)-NE1`) %>%
      mutate_at(-1, ~ str_replace(., pattern = "L", replacement = "")) %>%
      mutate_at(-1, as.numeric) %>%
      mutate(NE = rowMeans(select(., starts_with("NE")), na.rm = TRUE),
             NW = rowMeans(select(., starts_with("NW")), na.rm = TRUE),
             SE = rowMeans(select(., starts_with("SE")), na.rm = TRUE),
             SW = rowMeans(select(., starts_with("SW")), na.rm = TRUE)) %>%
      select(plotID, NE, NW, SE, SW) %>%
      pivot_longer(cols = -plotID,
                   names_to = "quadrant",
                   values_to = "canopyCover") %>%
      mutate(canopyCover = 100 - canopyCover/96*100),
    by = c("plotID", "quadrant")
  ) %>%
  mutate(quadrantID = paste0(plotID, quadrant)) %>%
  select(quadrantID, plotID, leafLitterDepth, canopyCover)

# Take the chance to get rid of ? in undergrowth sheet

# Make data file ----

save(trees, plots, quadrants,
     file = "data/Mandai_data.RData")

# Mandai <- new.env()
# load(file = "data/Mandai_data.RData", envir = Mandai)
# Mandai <- as.list(Mandai)