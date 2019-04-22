rehaping<-function(loc){

###format as restab
ResTabC <- format_as_restab(out)
# Make an empty res_tab2 for comparison to res_tab2_1
library(dplyr)
res_tab2_2 <- make_empty_res_tab2()
# convert to a matrix with integer values where levels would be
res_tab2_2 %<>% mutate_if(is.factor, as.integer) %>% as.matrix

# Import the C++ Reshaper
library(inline)
# We use readr to ensure that the UTF8 encoding of the .cpp file is preserved,
# newlines are interpreted properly, etc. Other approaches, such as using readLines
# don't automatically respect newlines and tabs properly.
cpp_reshaper <- cxxfunction(
  signature(ResTab='numeric', ResTabus='numeric', ResTabfb='numeric', res_tab2 = 'numeric'),
  plugin='Rcpp',
  body=readr::read_file(
    system.file('inline_cpp/format_restab2.cpp', package='tabus')))

res_tab2_2 <- cpp_reshaper(ResTabC$ResTab, ResTabC$ResTabus, ResTabC$ResTabfb, res_tab2_2)
}
