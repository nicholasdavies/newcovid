# Specification of tiers and lockdown effect

tier2 = list(res = 0.55, wor = -0.36, gro = -1.41, ret = -1.99, tra = -1.32)
tier3 = list(res = 1.44, wor = -2.72, gro = -1.98, ret = -8.83, tra = -4.62)

lockdown0 = list(res = 0, wor = 0, gro = 0, ret = 0, tra = 0)

# SE of tiers:
#                     indic tier      change         sd
#  1:  Grocery and pharmacy  T 1   0.0000000 0.00000000
#  2:  Grocery and pharmacy  T 2  -1.4081018 0.17815559
#  3:  Grocery and pharmacy  T 3  -1.9777151 0.28459703
#  4:                 Parks  T 1   0.0000000 0.00000000
#  5:                 Parks  T 2  -0.3708141 0.77230785
#  6:                 Parks  T 3 -13.3090377 2.97390548
#  7:           Residential  T 1   0.0000000 0.00000000
#  8:           Residential  T 2   0.5468330 0.05782779
#  9:           Residential  T 3   1.4384973 0.11679388
# 10: Retail and recreation  T 1   0.0000000 0.00000000
# 11: Retail and recreation  T 2  -1.9850383 0.22942847
# 12: Retail and recreation  T 3  -8.8287906 0.55393636
# 13:      Transit stations  T 1   0.0000000 0.00000000
# 14:      Transit stations  T 2  -1.3173915 0.28510438
# 15:      Transit stations  T 3  -4.6199864 0.68732933
# 16:            Workplaces  T 1   0.0000000 0.00000000
# 17:            Workplaces  T 2  -0.3615680 0.19479747
# 18:            Workplaces  T 3  -2.7212377 0.31536569


# IMPACT OF LOCKDOWN FOR WALES:
#                 variable       mean       se
# 1:  grocery_and_pharmacy -16.131957 1.330679
# 2:                 parks -39.064812 6.721604
# 3:           residential   6.616149 0.452793
# 4: retail_and_recreation -37.521544 1.460017
# 5:      transit_stations -18.956221 1.154021
# 6:            workplaces -18.517148 1.565395
lockdownW = list(res = 6.62, wor = -18.52, gro = -16.13, ret = -37.52, tra = -18.96)
seW = c(1.57, 1.33, 1.46, 1.15, 0.19, 0.18, 0.23, 0.29, 0.32, 0.28, 0.55, 0.69)
# for all se vectors, the definition is c(ld_wplc, ld_groc, ld_rtrc, ld_trns, t2_wplc, t2_groc, t2_rtrc, t2_trns, t3_wplc, t3_groc, t3_rtrc, t3_trns)

# IMPACT OF LOCKDOWN FOR NORTHERN IRELAND:
#                 variable        mean        se
# 1:  grocery_and_pharmacy  -1.3130269 0.9487168
# 2:                 parks   0.7009911 4.6486876
# 3:           residential   4.6109637 0.3439340
# 4: retail_and_recreation -14.9799870 1.4119360
# 5:      transit_stations -10.5919078 1.1598166
# 6:            workplaces -14.8809006 0.6219085
lockdownN = list(res = 4.61, wor = -14.88, gro = -1.31, ret = -14.98, tra = -10.59)
seN = c(0.62, 0.95, 1.41, 1.16, 0.19, 0.18, 0.23, 0.29, 0.32, 0.28, 0.55, 0.69)
# for all se vectors, the definition is c(ld_wplc, ld_groc, ld_rtrc, ld_trns, t2_wplc, t2_groc, t2_rtrc, t2_trns, t3_wplc, t3_groc, t3_rtrc, t3_trns)

# IMPACT OF LOCKDOWN FOR SCOTLAND:
#                 variable       mean        se
# 1:  grocery_and_pharmacy  0.7392551 0.6692967
# 2:                 parks  6.3261325 7.3463237
# 3:           residential  1.9880707 0.4309069
# 4: retail_and_recreation -6.4052156 0.4742608
# 5:      transit_stations -5.7244764 1.1456064
# 6:            workplaces -8.4104283 1.1061221
lockdownS = list(res = 1.99, wor = -8.41, gro = 0.74, ret = -6.41, tra = -5.72)
seS = c(1.11, 0.67, 0.47, 1.15, 0.19, 0.18, 0.23, 0.29, 0.32, 0.28, 0.55, 0.69)
# for all se vectors, the definition is c(ld_wplc, ld_groc, ld_rtrc, ld_trns, t2_wplc, t2_groc, t2_rtrc, t2_trns, t3_wplc, t3_groc, t3_rtrc, t3_trns)

# IMPACT OF LOCKDOWN FOR ENGLAND:
#                 variable       mean        se
# 1:  grocery_and_pharmacy  -4.782386 0.5578814
# 2:                 parks -19.335939 2.8611438
# 3:           residential   5.055462 0.2409448
# 4: retail_and_recreation -27.749227 0.8168795
# 5:      transit_stations -14.496543 1.1498847
# 6:            workplaces -10.295991 1.6405690
lockdownE = list(res = 5.05, wor = -10.30, gro = -4.78, ret = -27.75, tra = -14.50)
seE = c(1.64, 0.56, 0.82, 1.15, 0.19, 0.18, 0.23, 0.29, 0.32, 0.28, 0.55, 0.69)
# for all se vectors, the definition is c(ld_wplc, ld_groc, ld_rtrc, ld_trns, t2_wplc, t2_groc, t2_rtrc, t2_trns, t3_wplc, t3_groc, t3_rtrc, t3_trns)

# NEW:
# IMPACT OF LOCKDOWN FOR ENGLAND, WITH REGIONAL VARIATION
#                 variable       mean         sd
# 1:  Grocery and pharmacy  -4.042386  2.4651887
# 2:                 Parks -21.319373 13.2545134
# 3:           Residential   4.851414  0.9116758
# 4: Retail and recreation -26.643174  5.0953206
# 5:      Transit stations -14.259708  4.8450944
# 6:            Workplaces  -9.740738  1.8109323

lockdownER = list(res = 4.85, wor = -9.74, gro = -4.04, ret = -26.64, tra = -14.26)
seER = c(1.81, 2.47, 5.10, 4.85, 0.19, 0.18, 0.23, 0.29, 0.32, 0.28, 0.55, 0.69)
# for all se vectors, the definition is c(ld_wplc, ld_groc, ld_rtrc, ld_trns, t2_wplc, t2_groc, t2_rtrc, t2_trns, t3_wplc, t3_groc, t3_rtrc, t3_trns)
