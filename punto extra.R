# punto extra

metadata_db <- data.frame(
  individual = c(rep(1,2), rep(8,2), rep(14, 2), rep(12, 2), rep(17, 2), rep(18,2)),
  sample = c("S001", "S003", "S006", "S008", "S016", "S018", "S011", "S013", "S021", "S023", "S026", "S028"),
  condition = c(rep("SLE_R", 6), rep("SLE_NR", 6)),
  reponder = c(rep("R", 6), rep("NR", 6)),
  timepoint = c(rep(c(0, 7), 6)),
  timelabel = c(rep(c("d0", "d7"), 6)),
  identifiers = c("p1.S001.day0", "p1.S003.day7", "p8.S006.day0", "p8.S008.day7", "p14.S016.day0", "p14.S018.day7", "p12.S011.day0", "p12.S013.day7", "p17.S021.day0", "p17.S023.day7", "p18.S026.day0", "p18.S028.day7"),
  BioSample = c("SAMN38793003", "SAMN38793001", "SAMN38792999", "SAMN38792997", "SAMN38792990", "SAMN38792988", "SAMN38792994", "SAMN38792993", "SAMN38792985", "SAMN38792983", "SAMN38792980", "SAMN38792978")
)