library(data.table)
library(ggplot2)
library(zoo)
library(mgcv)
library(lubridate)
library(stringr)
library(cowplot)
library(qs)
library(ogwrangler)

# Load Google Mobility data and interface with ogwrangler
gm = qread("./data/google_mobility_uk.qs");
CreateCache();
gm[, name := ifelse(sub_region_2 != "", paste0(sub_region_2, ", ", sub_region_1), sub_region_1)];
gm = gm[name != ""];
gm_match = data.table(code = ogcode("*", "gmcty"));
gm_match[, name := ogwhat(code, "name")];
gm_match[, rgn := ogwhat(code, "rgn")];
gm_match[, country := str_sub(rgn, 1, 1)];
gm_match[, rgn := NULL];
gm_match[, pop2019 := ogwhat(code, "pop2019")];
gm_match[, nhs_cd := ogwhat(code, "nhser")];
gm_match[, lad_cd := ogwhat(code, "lad")];
gm_match[, lad_nm := ogwhat(lad_cd)];
gm_match[nhs_cd %like% "^E", nhs_nm := ogwhat(nhs_cd)];
gm_match[nhs_cd %like% "^N", nhs_nm := "Northern Ireland"];
gm_match[nhs_cd %like% "^S", nhs_nm := "Scotland"];
gm_match[nhs_cd %like% "^W", nhs_nm := "Wales"];

gm = merge(gm, gm_match, by = "name");
gm = gm[, .SD, .SDcols = c(1, 16:22, 9:15)];


# Format google mobility data
gm_melt = function(gm) {
    g = melt(gm, id.vars = 1:9)
    g[, variable := str_remove_all(variable, "_percent_change_from_baseline")]
    return (g)
}
ggg = gm_melt(gm)
qsave(ggg, "./data/gm_for_analysis-2020-01-18.qs")
fwrite(ggg, "./data/gm_for_analysis-2020-01-18.csv")
