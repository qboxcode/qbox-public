#
# mkvext.R: generate an arbitrary potential in a vext.xml file
#
# todo: insert function describing potential and XML header and trailer
library("caTools")
x <- rnorm(10)
s <- base64encode(x)
