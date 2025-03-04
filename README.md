# Análisis de datos de RNA-Seq 👾

Instructora: Dra. Evelia Coss, Posdoc de la Dra. Alejandra Medina, LIIGH-UNAM

Clases para los alumnos de Ciencias Genomicas de 4to semestre de la ENES, Juriquilla (25 y 27 Febrero, 4 y 6 Marzo 2025). Formando parte de la clase de Bioinformática y Estadística 2. 

## Descripción

El módulo consta de sesiones teóricas y prácticas impartidas de forma presencial, que cubrirán aspectos básicos del tópico como:

- Calidad y limpieza de archivos fastq
- Alineamiento y ensamblaje con el genoma de referencia usando STAR
- Generación del archivo de cuentas crudas
- Importar datos en R
- Normalización y corrección por batch
- Expresión diferencial con DESEq2 y edgeR
- Análisis funcional de los genes detectados
- Visualización grafica de los resultados

Se ofrecerán presentaciones detalladas sobre el uso de programas clave, todos de código abierto, utilizando datos extraídos de bases de datos. Además, se introducirá el uso de scripts sencillos en Bash y R, con el objetivo de aprender los conceptos básicos de estos lenguajes para el análisis de datos transcriptómicos.

## Contenido 📌

- Dia 1. Aspectos generales de RNA-Seq / Control de calidad de los datos
- Dia 2. Diversos pipeline para Alineamiento, ensamblaje y conteo de reads
- Dia 3. Importar datos en R, Normalización y Corrección por batch / DEG con DESeq2
- Dia 4. GSEA - Análisis funcional

### Dia 1. Aspectos generales de RNA-Seq / Control de calidad de los datos

- Fecha: martes 25 de febrero 2025
- Presentación:
    - [Aspectos generales de RNA-Seq](https://eveliacoss.github.io/RNAseq_classFEB2025/Presentaciones/Dia1_AspectosGenerales.html#1)
    - [Control de calidad de los datos](https://eveliacoss.github.io/RNAseq_classFEB2025/Presentaciones/Dia1_AspectosGenerales.html#43)
    - [Mis primeros pasos en Bash](https://eveliacoss.github.io/RNAseq_classFEB2025/Presentaciones/Dia1_AspectosGenerales.html#68)
- Datos: `/mnt/atgc-d1/bioinfoII/rnaseq/BioProject_2025/`
- Tarea para el **jueves 27 de feb**: Elegir en equipos los transcriptomas que emplearán en su proyecto - [Información diapositiva 60 y Google Classroom ENES](https://eveliacoss.github.io/RNAseq_classFEB2025/Presentaciones/Dia1_AspectosGenerales.html#64)

Ejemplo de entregable: [Reporte](https://eveliacoss.github.io/RNAseq_classFEB2025/Entregable_ejemplo1.html)

> Bioproject que usaron anteriormente (NO SE PUEDEN USAR): PRJNA821620 (At), PRJNA256121 (At), PRJNA858106 (Hs), PRJNA826506 (Hs), PRJNA649971 (Mm), PRJNA743296 (Hs), PRJNA739157 (Hs), PRJNA983389 (Hs), PRJNA759864 (Hs).
>
> Bioproject dados en cursos (NO SE PUEDEN USAR): PRJNA649971 (Mm) - [RNASeq_Workshop_Nov2023](https://github.com/EveliaCoss/RNASeq_Workshop_Nov2023) y [RNAseq_classFEB2024](https://github.com/EveliaCoss/RNAseq_classFEB2024), PRJNA821620 (At) - [RNAseq_classFEB2023](https://github.com/EveliaCoss/RNAseq_classFEB2023)

- Lecturas y cursos recomendados:
    - [Mi Clase de 2023](https://github.com/EveliaCoss/RNASeq_Workshop_Nov2023)
    - [Introduction to RNAseq Methods - Presentacion](https://bioinformatics-core-shared-training.github.io/Bulk_RNAseq_Course_Nov22/Bulk_RNAseq_Course_Base/Markdowns/01_Introduction_to_RNAseq_Methods.pdf)
    - [Intro-to-rnaseq-hpc-O2](https://github.com/hbctraining/Intro-to-rnaseq-hpc-O2/tree/master/lessons)
    - [RNA-seq technology and analysis overview](https://github.com/mdozmorov/presentations/tree/master/RNA-seq)


### Dia 2. Diversos pipeline para Alineamiento, ensamblaje y conteo de reads

- Fecha: jueves 27 de febrero 2025
- Presentación:
    - [Diversos pipeline para Alineamiento, ensamblaje y conteo de reads](https://eveliacoss.github.io/RNAseq_classFEB2025/Presentaciones/Dia2_QCAlineamiento.html#1)
- Lecturas y cursos recomendados:
    - Alineamiento *de novo* - [Trinity](https://github.com/trinityrnaseq/trinityrnaseq)
    - Alineamiento con el genoma de referencia - [STAR](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)
    - *Pseudoalineamiento* con [Kallisto](https://pachterlab.github.io/kallisto/manual)
    - [*Pseudoalineamiento* con Kallisto - practica](https://github.com/EveliaCoss/RNAseq_classFEB2023/tree/main/RNA_seq#practica2)

- Tarea para el **martes 4 de marzo**: 
    + Descripción de los datos (entregable 1, fecha límite: 25 de febrero).
    + Explicación de las gráficas generadas por MultiQC.
    + Conclusión sobre los datos: 
      - ¿Son viables para continuar con el análisis? 
      - ¿Qué pasos deben seguirse para mejorar la calidad de los datos?

Ejemplo de entregable: [Reporte](https://eveliacoss.github.io/RNAseq_classFEB2025/Entregable_ejemplo2.html)

**Subir en la tarea de Google Classroom**

### Dia 3. Importar datos en R, Normalización y Corrección por batch effect 🪲 / DEG con DESeq2

- Fecha: martes 4 de marzo 2025
- Presentación:
    - [Importar datos en R, Normalización y Corrección por batch / DEG con DESeq2](https://eveliacoss.github.io/RNAseq_classFEB2025/Presentaciones/Dia3_ImportarDatos.html#1)
- Scripts: https://github.com/EveliaCoss/RNAseq_classFEB2025/tree/main/Practica_Dia3/scripts/
- Lecturas y cursos recomendados:
    - [Metodos de normalizacion](https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html#2-create-deseq2-object)
    - [RNA-seq workflow: gene-level exploratory analysis and differential expression](https://www.bioconductor.org/packages/devel/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#pca-plot-using-generalized-pca)
    - [End-to-end RNA-Seq workflow](https://www.bioconductor.org/help/course-materials/2015/CSAMA2015/lab/rnaseqCSAMA.html)
    - [Transformation, Normalization, and Batch Effect Removal](https://bio-protocol.org/exchange/protocoldetail?id=4462&type=1)


### Dia 4. GSEA - Análisis funcional

- Fecha: jueves 6 de marzo 2025
- Presentación:
   - [GSEA - Análisis funcional](https://eveliacoss.github.io/RNAseq_classFEB2025/Presentaciones/Dia4_GSEA.html#1)
- Lecturas y cursos recomendados:
    - Base de datos [Gene Ontology Resource](http://geneontology.org/)
    - Base de datos [AmiGo2](https://amigo.geneontology.org/amigo/landing)
    - [Reducir terminos con REVIGO](http://revigo.irb.hr/)
    - Heatmap con [ComplexHeatmap -  Github](https://github.com/jokergoo/ComplexHeatmap)
    - Heatmap con [ComplexHeatmap -  manual](https://jokergoo.github.io/ComplexHeatmap-reference/book/)

## Requisitos

- Contar con una terminal en tu sistema operativo
- Los paquetes que emplearemos en R v4.0.2, se encuentran presentes en el cluster DNA, por lo que, no es necesario instalar nada en nuestras computadoras.
- Nodo de prueba (qlogin)

## Pipeline ⚡
### Pasos a seguir para el análisis de los datos de **RNA-Seq**

- Script [`load_data_inR.R`](https://github.com/EveliaCoss/RNAseq_classFEB2024/blob/main/Practica_Dia3/scripts/load_data_inR.R):

  **1)** Importar datos en R (archivo de cuentas) + metadatos y **2)** Crear una matriz de cuentas con todos los transcriptomas

- Script [`DEG_analysis.R`](https://github.com/EveliaCoss/RNAseq_classFEB2024/blob/main/Practica_Dia3/scripts/DEG_analysis.R):

  **3)** Crear el archivo `dds` con `DESeq2`, **4)** Correr el análisis de Expresión Diferencial de los Genes (DEG), **5)** Normalización de los datos, **6)** Detección de batch effect y **7)** Obtener los resultados de los contraste de DEG

- Script [`VisualizacionDatos.R`](https://github.com/EveliaCoss/RNAseq_classFEB2024/blob/main/Practica_Dia3/scripts/VisualizacionDatos.R):

  **8)** Visualización de los datos

- Script [`GOterms_analysis.R`](https://github.com/EveliaCoss/RNAseq_classFEB2024/blob/main/Practica_Dia4/scripts/GOterms_analysis.R):

  **9)** Análisis de Terminos funcionales (GO terms)

## Clase para mejorar tus skills

- Crear llaves y alias

Puedes crear llaves (ssh-keygen) y alias para acceder a los servidores de una manera segura y rápida: https://github.com/EveliaCoss/keygen

## Clase de retroalimentación (Fecha pendiente)


## Cursos para practicar 📕

- [VieRnes de Bioinformática 2023](https://github.com/EveliaCoss/ViernesBioinfo2023)
- [VieRnes de Bioinformática 2024](https://github.com/EveliaCoss/ViernesBioinfo2024)

## Referencias 📚
- [GOterms en S. cereviacae](https://www.yeastgenome.org/goSlimMapper)
- [Go Term finder](https://go.princeton.edu/cgi-bin/GOTermFinder?)
- [REVIGO - pagina principal](http://revigo.irb.hr/FAQ)
- [REVIGO - Reducir terminos GO, ejemplos](https://www.bioconductor.org/packages/release/bioc/vignettes/rrvgo/inst/doc/rrvgo.html)
- [ggprofiler2](https://cran.r-project.org/web/packages/gprofiler2/vignettes/gprofiler2.html)
- [Pathway enrichment analysis and visualization of omics data](https://cytoscape.org/cytoscape-tutorials/protocols/enrichmentmap-pipeline/#/)
- [Biomedical knowledge mining using GOSemSim and Clusterprofiler](https://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-kegg.html)
- [Pathview - Pagina principal](https://pathview.r-forge.r-project.org/)
- [Pathview - Manual](https://pathview.r-forge.r-project.org/pathview.pdf)
- [KEGG - Pathway ID](https://www.genome.jp/kegg/pathway.html)
- [AnnotationDbi](https://hbctraining.github.io/DGE_workshop_salmon_online/lessons/AnnotationDbi_lesson.html)
- [Gene Ontology enrichment analysis - Uso de varias bases de datos en R](https://davetang.org/muse/2010/11/10/gene-ontology-enrichment-analysis/)
