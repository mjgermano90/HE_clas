# Forzar a usar R 4.2.3 solo si existe en esta PC exacta
ruta_r <- "C:/Program Files/R/R-4.2.3/bin/x64/R.exe"

if (file.exists(ruta_r) && (R.version$major != "4" || package_version(R.version$minor) >= "4.0")) {
  message("🔄 Redirigiendo a R 4.2.3...")
  system(ruta_r, wait = FALSE)
  q("no")
}

source("renv/activate.R")
