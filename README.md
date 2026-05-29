# HE_clas: Análisis de Expresión Diferencial y Supervivencia

Este repositorio contiene los scripts y el entorno virtual para el análisis genómico y clínico del proyecto **HE_clas**. 

Para garantizar la reproducibilidad científica y evitar conflictos de compatibilidad en R (especialmente con paquetes de Bioconductor), el proyecto utiliza un entorno virtual administrado con `renv` y anclado a una fecha específica de repositorios (cápsula del tiempo).

---

## ⚙️ Especificaciones del Entorno
* **Versión de R requerida:** `4.2.3` (totalmente compatible con Bioconductor 3.16).
* **Gestión de paquetes:** `renv` (la configuración fija repositorios a la fecha `2023-06-01` en Posit Package Manager).
* **Sistema de archivos:** Las rutas se gestionan de forma relativa al directorio del proyecto (`getwd()`).

---

## 🚀 Instrucciones de Instalación y Configuración

Sigue estos pasos para clonar e inicializar el proyecto en cualquier computadora nueva sin conflictos:

1. **Clonar el repositorio:**
   ```bash
   git clone <URL_DEL_REPOSITORIO>
   cd HE_clas_v2
   ```

2. **Instalar R 4.2.3:**
   Asegúrate de tener instalada la versión exacta de R 4.2.3. Puedes descargarla del [sitio oficial del CRAN](https://cran.r-project.org/bin/windows/base/old/4.2.3/).

3. **Abrir el proyecto:**
   Abre el archivo [HE_clas.Rproj](./HE_clas.Rproj) en RStudio. Esto automáticamente activará el entorno `renv` a través del archivo `.Rprofile`.

4. **Restaurar el entorno virtual:**
   En la consola de R de RStudio, ejecuta el siguiente comando para descargar e instalar todas las dependencias exactas grabadas en `renv.lock`:
   ```r
   renv::restore()
   ```

---

## 📂 Dependencias de Datos (Acción Requerida)

Por motivos de tamaño y almacenamiento, las carpetas de datos (`Datos/`) y resultados (`Resultados/`) se encuentran ignoradas en Git (excepto por marcadores de estructura `.gitkeep`). 

Para que los scripts funcionen, debes descargar y colocar manualmente los siguientes archivos dentro de la carpeta `Datos/` en la raíz del proyecto:

### Datos Principales (en `Datos/`):
* **RNA-Seq de tumores (RDA):** `TCGA_BRCA_RNAseq_tumor_female_unique_immune.rda` *(Matriz de expresión génica de tumores)*
* **RNA-Seq tejido normal (RDA):** `TCGA_BRCA_RNAseq_normal_female.rda`
* **Datos unificados (RDA):** `TCGA_tumor_normal_luminal_HE_RNAseq_without_treatment.rda`
* **Mutaciones somáticas (TXT):** `maf_BRCA.txt` *(Archivo MAF de mutaciones de TCGA BRCA)*
* **Datos H&E (CSV):** `brca_H&E.csv`
* **Estado clínico (CSV):** `clinical_status.csv`
* **Tipos celulares HE (XLSX):** `immune_cell_HE.xlsx`
* **Base de datos Immport (TXT):** `ImmuneGeneList.txt`
* **Checkpoints inmunológicos (XLSX):** `Immune_checkpoint_genes_stimulatory_inhibitory_The_Immune_Landscape_Cancer_listo.xlsx`

### Pathways y GMTs (en `Datos/Pathways/`):
* Archivos de firmas genéticas en formato `.gmt` descargados de MSigDB:
  * `c2.cp.kegg_legacy.v2023.2.Hs.symbols.gmt`
  * `c5.go.bp.v2023.2.Hs.symbols.gmt`
  * `c5.go.cc.v2023.2.Hs.symbols.gmt`
  * `c5.go.mf.v2023.2.Hs.symbols.gmt`
  * `c6.all.v2023.2.Hs.symbols.gmt`
  * `h.all.v2024.1.Hs.symbols.gmt`
  * `immunologic_signature_gene_sets.gmt`
* Archivos adicionales de categorías:
  * `hallmark_category.csv`
  * `filtered_kegg_pathways_by_category.csv`

---

## 📊 Estructura de Salida
Los scripts guardan automáticamente todos los archivos PDF, SVG, Excel y CSV generados en la carpeta `Resultados/HE clas Results/`. El script inicializa dinámicamente este directorio si no existe previamente en la nueva computadora.