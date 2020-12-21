suppressPackageStartupMessages({
  require(data.table)
  require(ggplot2)
})

.debug <- "~/Downloads"
.args <- if (interactive()) sprintf(c(
  "%s/nextstrain_groups_ngs-sa_COVID19-ZA-2020.12.17_metadata.tsv",
  "%s/newlin.csv",
  "comparison.png"
), .debug) else commandArgs(trailingOnly = TRUE)

rsa.dt <- fread(.args[1])[Country == "South Africa"][order(`Collection Data`)]
uk.dt <- fread(.args[2])[order(sample_date)]

plot.dt <- rbind(
  rsa.dt[,
    .(.N, iso3c = "ZAF"),
    by = .(date=`Collection Data`, newvariant = `Clade`=="20C")
  ],
  uk.dt[,
    .(.N, iso3c="GBR"),
    by=.(date=sample_date, newvariant = var2)
  ]
)

plot.dt[, total := sum(N), by=.(date, iso3c) ]

plot.dt[,
  c("lo95","hi95") := 
    as.data.table(t(mapply(
      function(x, n, p=x/n) binom.test(x, n, p, conf.level = .95)$conf.int,
      x = N, n = total
    )))
][,
  c("lo50","hi50") := 
    as.data.table(t(mapply(
      function(x, n, p=x/n) binom.test(x, n, p, conf.level = .50)$conf.int,
      x = N, n = total
    )))
]

p <- ggplot(plot.dt[newvariant == TRUE]) + aes(date) +
  facet_grid(. ~ iso3c) +
  geom_ribbon(aes(ymin = lo95, ymax = hi95), alpha = 0.1) +
  geom_ribbon(aes(ymin = lo50, ymax = hi50), alpha = 0.2) +
  geom_line(aes(y=N/total)) +
  scale_x_date(
    name = NULL,
    date_breaks = "months", date_minor_breaks = "weeks",
    date_labels = "%b"
  ) +
  scale_y_continuous("Novel Variant Fraction") +
  coord_cartesian(ylim = c(0, 1), xlim = c(as.Date("2020-10-01"), NA)) +
  theme_minimal()

ggsave(tail(.args, 1), p, width = 6, height = 3, units = "in", dpi = 300)