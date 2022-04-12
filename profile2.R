
library(randomPlantedForest)

randomPlantedForest:::start_profiler("/tmp/profile.out")
rpf(Species ~ ., iris, ntrees = 50)
randomPlantedForest:::stop_profiler()

# Run from terminal to actually write to svg file
# TODO: Fix to run this from R?
system("~/go/bin/pprof -svg -nodefraction 0.01 src/randomPlantedForest.so /tmp/profile.out")

rsvg::rsvg_png("profile001.svg", "profile001.png", width = 2000)
