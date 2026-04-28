# HE_clas: Análisis de Expresión Diferencial y Supervivencia

Este repositorio contiene los scripts y el entorno virtual para el análisis genómico y clínico del proyecto **HE_clas**. Para garantizar la reproducibilidad y evitar conflictos de compatibilidad de R, el proyecto está empaquetado con `renv`.

## ⚙️ Especificaciones del Entorno
- **Versión de R:** 4.2.3 (2023-03-15)
- **Sistema Operativo:** Windows 10
- **Gestión de paquetes:** `renv`
- **Configuración de Repositorios:** Se utilizó una "cápsula del tiempo" de CRAN y Bioconductor para fijar la compatibilidad:
  ```r
  options(repos = c(CRAN = "[https://packagemanager.posit.co/cran/2023-06-01](https://packagemanager.posit.co/cran/2023-06-01)"))

Bioconductor 3.16: Se ancló estrictamente esta versión mediante BiocManager porque es la única totalmente compatible con R 4.2.3. Si se intenta usar una versión superior, se generarán conflictos letales con el motor interno de Matrix.