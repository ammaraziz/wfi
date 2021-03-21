files = list.files(path = ".", pattern = "*.txt")
lowCovHA = list()
lowCovPos = list()
for (f in files) {
  print(f)
  depth = read.table(f, sep = "\t", header = T)
  
  pos_min_cov = min(which(depth$Coverage.Depth > 20))
  if (pos_min_cov > 20) {
    lowCovHA = append(lowCovHA, f)
    lowCovPos = append(lowCovPos, pos_min_cov)}

}

tmp = data.frame(Sample = unlist(lowCovHA), PosDepthAt20 = unlist(lowCovPos))

write.csv(tmp, file = "iSeq22_HA_issues.txt", row.names = F, quote = F)
