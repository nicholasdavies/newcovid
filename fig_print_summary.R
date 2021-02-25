library(data.table)
library(ggplot2)
library(cowplot)
library(colourpal)
library(zoo)
library(lubridate)
library(maps)
library(mapdata)
library(maptools)
library(rgdal)
library(ggmap)
library(ggplot2)
library(rgeos)
library(broom)
library(plyr)
library(pals)

# Load map of England

shapefile = readOGR(dsn = "~/Documents/uk_covid_data_sensitive/Local_Authority_Districts_(December_2019)_Boundaries_UK_BGC-shp/", 
    layer = "Local_Authority_Districts_(December_2019)_Boundaries_UK_BGC")

mapdata = tidy(shapefile, region = "lad19cd")
setDT(mapdata)
mapdata = mapdata[id %like% "^E"]

# https://github.com/wch/ggplot2/wiki/New-theme-system
new_theme_empty <- theme_bw()
new_theme_empty$line <- element_blank()
new_theme_empty$rect <- element_blank()
new_theme_empty$strip.text <- element_blank()
new_theme_empty$axis.text <- element_blank()
new_theme_empty$axis.title <- element_blank()
new_theme_empty$plot.title <- element_text(hjust = 0.5, size = 10, face = "bold")
new_theme_empty$plot.margin <- structure(c(0, 0, -1, -1), unit = "lines", valid.unit = 3L, class = "unit")

# load SGTF data by month

ll = fread("~/Documents/uk_covid_data_sensitive/phe/20210208/Anonymised Combined Line List 20210208.csv")
sgtf = fread("~/Documents/uk_covid_data_sensitive/phe/20210208/SGTF_linelist_20210208.csv")
ll[, specimen_date := dmy(specimen_date)]
sgtf[, specimen_date := ymd(specimen_date)]
d = merge(ll, sgtf, by = c("FINALID", "specimen_date"), all = TRUE)
d2 = d[pillar == "Pillar 2" & specimen_date >= "2020-10-01" & specimen_date <= "2021-01-31", 
    .(sgtf = sum(sgtf_under30CT == 1, na.rm = T), other = sum(sgtf_under30CT == 0, na.rm = T)), by = .(mo = month(specimen_date), LTLA_code)]

plot_map = function(mapdata, d2)
{
    mapdata = merge(mapdata, d2, by.x = "id", by.y = "LTLA_code", all = TRUE)

    ggplot(mapdata) + 
        geom_polygon(aes(x = long, y = lat, group = group, fill = sgtf / (sgtf + other)), colour = NA) +
        coord_fixed(1) +
        scale_fill_gradientn(colours = rev(ocean.matter(20)), limits = c(0, 1)) +
        new_theme_empty
}

p_oct
p_oct = plot_map(mapdata, d2[mo == 10]) + labs(title = "October 2020", fill = NULL) + theme(legend.position = c(0.23, 0.6))
p_nov = plot_map(mapdata, d2[mo == 11]) + labs(title = "November 2020", fill = NULL) + theme(legend.position = "none")
p_dec = plot_map(mapdata, d2[mo == 12]) + labs(title = "December 2020", fill = NULL) + theme(legend.position = "none")
p_jan = plot_map(mapdata, d2[mo == 1]) + labs(title = "January 2021", fill = NULL) + theme(legend.position = "none")

p_spread = plot_grid(p_oct, p_nov, p_dec, p_jan, nrow = 2)


# Results data
r_results = fread(
"ModelType	ModelNo	Geography	RCentral	RLo	RHi
GLMM	1a	UK	77	73	81
GLMM	1b	UK	67	65	69
GLMM	2a	UK	90	85	95
GLMM	2h	England	83	81	84
Rt regression	4a	England	43	38	48
Rt regression	4b	England	57	52	62
Transmission model	5a	England	82	43	130
GLMM	3a	Denmark	55	45	66
GLMM	3b	Switzerland	74	66	82
GLMM	3c	USA	59	56	83")

r_results[, Geography := factor(Geography, unique(Geography))]
r_results[, X := c(1, 2, 3, 6, 7, 8, 9, 12, 15, 18)]

pr = ggplot(r_results) + 
    geom_pointrange(aes(x = X, y = RCentral, ymin = RLo, ymax = RHi, colour = ModelType, group = Geography), 
        position = position_dodge(width = 0.5), fatten = 1) +
    geom_hline(aes(yintercept = 0), linetype = "33") +
    scale_x_continuous(breaks = c(2, 7.5, 12, 15, 18), labels = c("UK", "EN", "DK", "CH", "US")) +
    labs(x = "Country", y = "Estimated change\nin reproduction number (%)", colour = "Model type") +
    theme(legend.position = c(0.1, 0.25)) +
    ylim(-100, NA) +
    scale_colour_manual(values = c("#f4b141", "#97c1c9", "#b6a195"))

deaths = fread(
"Vaccination	NPI	Central	Lo	Hi
No vaccination	Moderate	216000	205000	227000
No vaccination	High (schools open)	146000	138000	152000
No vaccination	High (schools closed)	147000	139000	155000
No vaccination	Very high	149000	140000	157000
200,000 vaccinations per week	Moderate	202000	192000	213000
200,000 vaccinations per week	High (schools open)	137000	130000	143000
200,000 vaccinations per week	High (schools closed)	129000	123000	135000
200,000 vaccinations per week	Very high	119000	112000	125000
2 million vaccinations per week	Moderate	140000	133000	146000
2 million vaccinations per week	High (schools open)	98900	94600	103000
2 million vaccinations per week	High (schools closed)	81000	77600	84200
2 million vaccinations per week	Very high	58200	56100	60300")

deaths[, NPI := str_replace_all(NPI, " \\(schools ", "\n(schools\n")]
deaths[, Vaccination := factor(Vaccination, levels = unique(Vaccination))]
deaths[, NPI := factor(NPI, levels = unique(NPI))]

pd = ggplot(deaths) +
    geom_pointrange(aes(x = NPI, y = Central/1000, ymin = Lo/1000, ymax = Hi/1000, colour = Vaccination), 
        fatten = 1, position = position_dodge(width = 0.8)) +
    ylim(0, NA) +
    labs(x = "Stringency of non-pharmaceutical interventions", y = "Estimated additional\nCOVID-19 deaths (thousands)") +
    theme(legend.position = c(0.1, 0.25)) +
    scale_colour_manual(values = c("#881133", "#44aa55", "#55aaff"))

pg = plot_grid(p_spread, plot_grid(pr, pd, nrow = 2, labels = c("B", "C"), label_size = 10, align = "v", axis = "bottom"), labels = c("A", ""), label_size = 10)
ggsave("./output/fig0.pdf", pg, width = 20, height = 14, units = "cm", useDingbats = FALSE)

