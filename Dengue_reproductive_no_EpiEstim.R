## Load packages
library(readxl)
library(readr)
library(dplyr)
library(lubridate)
library(tidyr)
library(EpiEstim)
library(ggplot2)


## Read in the Excel file
Togo_cases_all_tests <- read_excel("data/Togo_cases_all_tests.xlsx")
View(Togo_cases_all_tests)


## Filter for confirmed cases and prepare time series

# Ensure date column is correctly formatted
Togo_cases_all_tests$`DATE DE PRELEVEMENT` <- as.Date(Togo_cases_all_tests$`DATE DE PRELEVEMENT`)

summary(Togo_cases_all_tests$`DATE DE PRELEVEMENT`)

# Filter for confirmed cases only
confirmed_cases <- Togo_cases_all_tests %>% 
  filter(Classification == "Confirmé")

View(confirmed_cases)

summary(confirmed_cases$`DATE DE PRELEVEMENT`)


# Aggregate cases by collection date
daily_cases <- confirmed_cases %>%
  group_by(`DATE DE PRELEVEMENT`) %>%
  count(`DATE DE PRELEVEMENT`, name = "cases") %>%
  rename(dates = `DATE DE PRELEVEMENT`)

# Plot time series
ggplot(daily_cases, aes(x = dates, y = cases)) +
  geom_bar(stat = "identity") +
  labs(x = "Date", y = "Number of Cases (daily)") +
  theme_minimal()

# Save as .png
ggsave("output/Daily_cases.png")

# Rename 'cases' column to 'I' as required by EpiEstim
incidence_data <- daily_cases %>%
  rename(I = cases)

# Create a full set of dates from min to max (required by EpiEstim)
daily_cases_alldates <- data.frame(
  dates = seq(min(daily_cases$dates), max(daily_cases$dates), by = "day")
)

# Join with your existing data, fill missing days with 0
incidence_data <- daily_cases_alldates %>%
  left_join(daily_cases, by = "dates") %>%
  mutate(I = ifelse(is.na(cases), 0, cases)) %>%
  select(dates, I)


## Estimate R(t) using EpiEstim
res <- estimate_R(
  incid = incidence_data,
  method = "parametric_si",
  config = make_config(list(
    mean_si = 17,   # From Codeco et al. 2018 (https://10.1016/j.epidem.2018.05.011)
    std_si = 8
  ))
)

# Plot the estimated R over time
plot(res)


## Estimate R(t) using EpiEstim with restricted date range 
## (because the early cases are distorting estimates)
# Date range
start_date <- as.Date("2024-05-07")
end_date <- as.Date("2024-11-25")

# Filter the aggregated time series
incidence_data_less_outliers <- incidence_data %>%
  filter(dates >= start_date & dates <= end_date)

res_less_outliers <- estimate_R(
  incid = incidence_data_less_outliers,
  method = "parametric_si",
  config = make_config(list(
    mean_si = 17,
    std_si = 8
  ))
)

plot(res_less_outliers, legend = F)

## Export the Rt values
write.csv(res_less_outliers[["R"]], file = "output/Rt_estimates_Togo_17daySI.csv", row.names = FALSE)


################################################################################


### Regional analysis ###

# Aggregate cases by collection date and region
daily_cases_region <- confirmed_cases %>%
  group_by(`DATE DE PRELEVEMENT`, REGION) %>%
  count(`DATE DE PRELEVEMENT`, name = "cases") %>%
  rename(dates = `DATE DE PRELEVEMENT`)

# Plot: time series by region (overlaid)
ggplot(daily_cases_region, aes(x = dates, y = cases, fill = REGION)) +
  geom_bar(stat = "identity") +
  labs(title = "Daily Cases by Region",
       x = "Date", y = "Number of Cases") +
  theme_minimal()

# Plot: time series by region (individual)
ggplot(daily_cases_region, aes(x = dates, y = cases, fill = REGION)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ REGION, scales = "free_y") +
  labs(title = "Daily Cases by Region",
       x = "Date", y = "Number of Cases") +
  theme_minimal()

# Save as .png
ggsave("output/Daily_cases_by_region.png")


## GRAND LOME ##

# Rename 'cases' column to 'I' as required by EpiEstim
daily_cases_GL <- daily_cases_region %>%
  filter(REGION == 'GRAND LOME')

# Create a full set of dates from min to max (required by EpiEstim)
daily_cases_alldates <- data.frame(
  dates = seq(min(daily_cases_GL$dates), max(daily_cases_GL$dates), by = "day")
)

# Join with your existing data, fill missing days with 0
incidence_data_GL <- daily_cases_alldates %>%
  left_join(daily_cases_GL, by = "dates") %>%
  mutate(I = ifelse(is.na(cases), 0, cases)) %>%
  select(dates, I)


## Estimate R(t) using EpiEstim
res_GL <- estimate_R(
  incid = incidence_data_GL,
  method = "parametric_si",
  config = make_config(list(
    mean_si = 17,   # From Codeco et al. 2018 (https://10.1016/j.epidem.2018.05.011)
    std_si = 8
  ))
)

# Plot the estimated R over time
plot(res_GL, legend = F)

## Export the Rt values
write.csv(res_GL[["R"]], file = "output/Rt_estimates_GL_17daySI.csv", row.names = FALSE)


## SAVANES ##

# Rename 'cases' column to 'I' as required by EpiEstim
daily_cases_Savanes <- daily_cases_region %>%
  filter(REGION == 'SAVANES') %>%
  filter(dates >= as.Date("2024-06-01"))

# Create a full set of dates from min to max (required by EpiEstim)
daily_cases_alldates <- data.frame(
  dates = seq(min(daily_cases_Savanes$dates), max(daily_cases_Savanes$dates), by = "day")
)

# Join with your existing data, fill missing days with 0
incidence_data_Savanes <- daily_cases_alldates %>%
  left_join(daily_cases_Savanes, by = "dates") %>%
  mutate(I = ifelse(is.na(cases), 0, cases)) %>%
  select(dates, I)


## Estimate R(t) using EpiEstim
res_Savanes <- estimate_R(
  incid = incidence_data_Savanes,
  method = "parametric_si",
  config = make_config(list(
    mean_si = 17,   # From Codeco et al. 2018 (https://10.1016/j.epidem.2018.05.011)
    std_si = 8
  ))
)

# Plot the estimated R over time
plot(res_Savanes, legend = F)

## Export the Rt values
write.csv(res_Savanes[["R"]], file = "output/Rt_estimates_Savanes_17daySI.csv", row.names = FALSE)



## KARA ##

# Rename 'cases' column to 'I' as required by EpiEstim
daily_cases_Kara <- daily_cases_region %>%
  filter(REGION == 'KARA')

# Create a full set of dates from min to max (required by EpiEstim)
daily_cases_alldates <- data.frame(
  dates = seq(min(daily_cases_Kara$dates), max(daily_cases_Kara$dates), by = "day")
)

# Join with your existing data, fill missing days with 0
incidence_data_Kara <- daily_cases_alldates %>%
  left_join(daily_cases_Kara, by = "dates") %>%
  mutate(I = ifelse(is.na(cases), 0, cases)) %>%
  select(dates, I)


## Estimate R(t) using EpiEstim
res_Kara <- estimate_R(
  incid = incidence_data_Kara,
  method = "parametric_si",
  config = make_config(list(
    mean_si = 17,   # From Codeco et al. 2018 (https://10.1016/j.epidem.2018.05.011)
    std_si = 8
  ))
)

# Plot the estimated R over time
plot(res_Kara, legend = F)

## Export the Rt values
write.csv(res_Kara[["R"]], file = "output/Rt_estimates_Kara_17daySI.csv", row.names = FALSE)



################################################################################

### Epidemiological analysis (Figures for manuscript) ###

library(dplyr)
library(lubridate)
library(readr)
library(ggplot2)
library(grid)   # for unit()


### National Rt estimates overlaid on weekly cases ###

# Rt window used when you created the CSV
start_date <- as.Date("2024-05-07")
end_date   <- as.Date("2024-11-25")

# Weekly cases over the full data
weekly_cases_full <- daily_cases %>%
  mutate(week = floor_date(dates, unit = "week", week_start = 7)) %>%
  group_by(week) %>%
  summarise(cases = sum(cases), .groups = "drop")

# Rt (map t_end -> calendar date = window end)
rt <- read_csv("output/Rt_estimates_Togo_17daySI.csv", show_col_types = FALSE) %>%
  mutate(
    R_median = dplyr::coalesce(`Median(R)`, `Mean(R)`),
    R_low    = `Quantile.0.025(R)`,
    R_high   = `Quantile.0.975(R)`
  )
dates_seq <- seq.Date(start_date, end_date, by = "day")
rt <- rt %>%
  mutate(date = dates_seq[t_end]) %>%
  select(date, R_median, R_low, R_high)

write_csv(rt, "output/Rt_estimates_Togo_clean.csv")

# Dual-axis scaling
scale_factor_full <- max(weekly_cases_full$cases, na.rm = TRUE) /
  max(rt$R_high, na.rm = TRUE)
teal <- "#008B8B"

# Create plot
p_full <- ggplot() +
  geom_col(data = weekly_cases_full,
           aes(x = week, y = cases),
           width = 7, fill = "grey70", color = "grey30", linewidth = 0.2) +
  geom_ribbon(data = rt,
              aes(x = date, ymin = R_low * scale_factor_full, ymax = R_high * scale_factor_full),
              fill = teal, alpha = 0.18) +
  geom_line(data = rt,
            aes(x = date, y = R_median * scale_factor_full),
            color = teal, linewidth = 1) +
  geom_hline(yintercept = 1 * scale_factor_full,
             linetype = "dashed", linewidth = 0.6, color = teal) +
  scale_y_continuous(
    name = "Weekly cases",
    sec.axis = sec_axis(~ . / scale_factor_full, name = "Rt")
  ) +
  labs(x = "Date") +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_line(color = "grey82", linewidth = 0.5),
    panel.grid.major.y = element_line(color = "grey82", linewidth = 0.5),
    panel.grid.minor   = element_line(color = "grey90", linewidth = 0.3),
    panel.border = element_rect(color = "grey35", fill = NA, linewidth = 0.7),
    axis.ticks = element_line(color = "grey35"),
    axis.ticks.length = unit(3, "pt"),
    plot.title = element_blank()
  )

plot(p_full)

# Save as PNG and PDF
ggsave("output/Rt_cases_full.png", plot = p_full, dpi = 300)
ggsave("output/Rt_cases_full.pdf", plot = p_full)


## Restricted range ##

# Weekly cases (restricted)
weekly_cases_res <- daily_cases %>%
  filter(dates >= start_date, dates <= end_date) %>%
  mutate(week = floor_date(dates, unit = "week", week_start = 7)) %>%
  group_by(week) %>%
  summarise(cases = sum(cases), .groups = "drop")

scale_factor_res <- max(weekly_cases_res$cases, na.rm = TRUE) /
  max(rt$R_high, na.rm = TRUE)

p_res <- ggplot() +
  geom_col(data = weekly_cases_res,
           aes(x = week, y = cases),
           width = 7, fill = "grey70", color = "grey30", linewidth = 0.2) +
  geom_ribbon(data = rt,
              aes(x = date, ymin = R_low * scale_factor_res, ymax = R_high * scale_factor_res),
              fill = teal, alpha = 0.18) +
  geom_line(data = rt,
            aes(x = date, y = R_median * scale_factor_res),
            color = teal, linewidth = 1) +
  geom_hline(yintercept = 1 * scale_factor_res,
             linetype = "dashed", linewidth = 0.6, color = teal) +
  scale_y_continuous(
    name = "Weekly cases",
    sec.axis = sec_axis(~ . / scale_factor_res, name = "Rt")
  ) +
  labs(x = "Date") +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_line(color = "grey82", linewidth = 0.5),
    panel.grid.major.y = element_line(color = "grey82", linewidth = 0.5),
    panel.grid.minor   = element_line(color = "grey90", linewidth = 0.3),
    panel.border = element_rect(color = "grey35", fill = NA, linewidth = 0.7),
    axis.ticks = element_line(color = "grey35"),
    axis.ticks.length = unit(3, "pt"),
    plot.title = element_blank()
  )

plot(p_res)

# Save as PNG and PDF
ggsave("output/Rt_cases_restricted.png", plot = p_res, dpi = 300)
ggsave("output/Rt_cases_restricted.pdf", plot = p_res)







