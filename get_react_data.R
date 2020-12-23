

library(data.table)

react_Rt_1 <- 
  data.table(time_from = as.Date(c("2020-05-01", "2020-06-19", "2020-07-24",
                               "2020-08-20", "2020-09-18", "2020-10-16", "2020-10-26")),
         creation_date = as.Date(c("2020-06-01", "2020-07-07", "2020-08-11",
                                   "2020-09-08", "2020-10-05", "2020-10-25", "2020-11-02")),
         `0.5` = c(0.6, 0.5, 1.3, 1.7, 1.16, 1.56, 0.85),
         `0.05` = c(0.4, 0.4, 0.9, 1.4, 1.05, 1.27, 0.73),
         `0.95` = c(0.7, 0.8, 1.8, 2.0, 1.27, 1.89, 0.99),
         estimate = "per_round")
react_Rt_2 <-
  data.table(time_from = react_Rt_1$time_from[1:(nrow(react_Rt_1) - 2)],
         creation_date = react_Rt_1$creation_date[c(2:5,7)],
         `0.5` = c(0.9, 0.8, 1.3, 1.39, 1.19),
         `0.05` = c(0.8, 0.8, 1.2, 1.34, 1.17),
         `0.95` = c(0.9, 0.9, 1.4, 1.43, 1.21),
         estimate = "two_rounds")
react_Rt <- rbind(react_Rt_1, react_Rt_2)

