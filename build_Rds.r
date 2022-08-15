wd <- getwd()
setwd("man")
Rds <- system2("ls", stdout = TRUE)
for (Rd in Rds) {
    message(Rd)
    parsed <- tools::parse_Rd(Rd)
    name <- substr(Rd, 1, nchar(Rd) - 3)
    tools::Rd2HTML(parsed, out = paste0("../documentation/", name, ".html"))
}
setwd(wd)