
library(randomPlantedForest)

randomPlantedForest:::start_profiler("/tmp/profile.out")
rpf(Species ~ ., iris, ntrees = 50)
randomPlantedForest:::stop_profiler()

system("~/go/bin/pprof -svg -nodefraction 0.01 pkg/src/pkg.so /tmp/profile.
out")

rsvg::rsvg_png("profile001.svg", "profile001.png", width = 2000)
