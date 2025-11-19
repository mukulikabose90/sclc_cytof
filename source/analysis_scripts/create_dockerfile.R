library(dockerfiler)

the_lockfile <- "renv.lock"

my_dock <- dock_from_renv(
  lockfile = the_lockfile,
  distro = "focal",
  FROM = "rocker/verse"
)

my_dock
