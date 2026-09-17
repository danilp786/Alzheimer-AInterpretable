# 🧬 Análisis del microbioma oral en la enfermedad de Alzheimer mediante Machine Learning

> **Trabajo Fin de Grado — Grado en Bioinformática**
> **Universidad Católica de Ávila**
> **Autor:** Daniel de Lara Pérez

---

## 📌 Descripción

Este repositorio contiene el código y la infraestructura desarrollados para el **análisis del microbioma oral como fuente potencial de biomarcadores de la enfermedad de Alzheimer (EA)** mediante técnicas de **aprendizaje automático, clustering y explicabilidad de modelos**.

El proyecto parte de datos públicos de secuenciación del gen **16S rRNA**, correspondientes a muestras de microbioma oral de pacientes con enfermedad de Alzheimer y sujetos control.

El objetivo principal es explorar si las características microbianas pueden proporcionar información útil para:

* 🧠 Clasificar muestras asociadas a enfermedad de Alzheimer.
* 🦠 Identificar características microbianas potencialmente relevantes.
* 🤖 Evaluar diferentes técnicas de aprendizaje automático.
* 🔎 Analizar la interpretabilidad de los modelos mediante técnicas de explicabilidad.
* 🧪 Establecer una **infraestructura reproducible** que pueda servir como base para futuros estudios con cohortes de mayor tamaño.

> ⚠️ **Importante:** este proyecto constituye una **prueba de concepto**. Los resultados obtenidos no permiten establecer biomarcadores clínicamente validados ni confirmar la capacidad del microbioma oral para diagnosticar la enfermedad de Alzheimer.

---

## 🔬 Flujo general del proyecto

```text
                    DATOS PÚBLICOS
                         │
                         ▼
              ┌─────────────────────┐
              │   NCBI / SRA        │
              │   PRJNA1243403      │
              └──────────┬──────────┘
                         │
                         ▼
                 FASTQ paired-end
                         │
                         ▼
              ┌─────────────────────┐
              │      QIIME 2        │
              ├─────────────────────┤
              │ Quality control      │
              │ Cutadapt             │
              │ V3-V4               │
              │ Join pairs           │
              │ Deblur               │
              │ ASVs                 │
              │ SILVA                │
              └──────────┬──────────┘
                         │
                         ▼
                   Tabla de ASVs
                         │
                         ▼
              Preprocesamiento Python
                         │
              ┌──────────┼───────────┐
              │          │           │
              ▼          ▼           ▼
          GBT / SVC   Clustering    MGM
              │          │           │
              ▼          ▼           ▼
           SHAP /     PCA /         Transformer
         BorutaSHAP   K-means
              │          │           │
              └──────────┼───────────┘
                         ▼
                 Interpretación
                  de resultados
```

---

## 📊 Datos utilizados

Los datos proceden del **Sequence Read Archive (SRA) del National Center for Biotechnology Information (NCBI)** y pertenecen al proyecto:

**PRJNA1243403**

El análisis se centró exclusivamente en las muestras correspondientes a la región **V3-V4 del gen 16S rRNA**.

### Cohorte

| Característica           |                                        Valor |
| ------------------------ | -------------------------------------------: |
| Individuos               |                                           51 |
| Pacientes con EA         |                                           26 |
| Controles                |                                           25 |
| Parejas paciente-control |                                           25 |
| Paciente sin pareja      |                                            1 |
| Muestras microbiológicas |                                          101 |
| Región analizada         |                               16S rRNA V3-V4 |
| Tipos de muestra         | Saliva no estimulada / exudado supragingival |

Una característica especialmente importante de la cohorte es que los controles están emparejados con los pacientes. Por este motivo, la estructura de las parejas se tuvo en cuenta durante la validación para reducir el riesgo de **data leakage**.

---

# 🧫 Pipeline bioinformático

El procesamiento de las secuencias se implementó mediante una pipeline reproducible desarrollada en **Bash** y basada en **QIIME 2**.

El objetivo es transformar las lecturas FASTQ originales en una tabla de **Amplicon Sequence Variants (ASVs)** y posteriormente realizar su clasificación taxonómica.

### Flujo de procesamiento

```text
FASTQ paired-end
       │
       ▼
Importación QIIME 2
       │
       ▼
Control de calidad
       │
       ▼
Eliminación de cebadores
       │
       ▼
Unión R1 + R2
       │
       ▼
Deblur
       │
       ▼
ASVs
       │
       ▼
Clasificación taxonómica
       │
       ▼
SILVA 138
```

### Cebadores utilizados

| Región  | Cebador                        |
| ------- | ------------------------------ |
| Forward | `341F — CCTACGGGNGGCWGCAG`     |
| Reverse | `805R — GACTACHVGGGTATCTAATCC` |

### Principales herramientas

* **QIIME 2** — procesamiento de microbioma.
* **Cutadapt** — eliminación de cebadores.
* **VSEARCH** — unión de lecturas paired-end.
* **Deblur** — generación de ASVs.
* **SILVA 138** — clasificación taxonómica.
* **Python** — procesamiento y análisis.
* **Bash** — automatización de la pipeline.

La pipeline utiliza además:

```bash
set -euo pipefail
```

para detener la ejecución ante errores y evitar continuar el procesamiento cuando alguna etapa falla.

---

# 🤖 Machine Learning

Una vez obtenidas las tablas de ASVs, los datos se procesaron en Python para construir modelos capaces de diferenciar entre:

```text
0 → Control
1 → Enfermedad de Alzheimer
```

Antes del entrenamiento se realizó:

1. Filtrado de ASVs por prevalencia.
2. Transformación **Centered Log-Ratio (CLR)**.
3. Integración con los metadatos.
4. Codificación *one-hot* del medio de obtención de la muestra.
5. Validación cruzada agrupada por pareja.

Se estudiaron principalmente dos enfoques de clasificación:

### 🌳 Gradient Boosted Trees

Se utilizó un modelo de árboles potenciados y se realizó optimización de hiperparámetros mediante **Optuna**.

El espacio de búsqueda incluyó parámetros como:

* `n_estimators`
* `max_depth`
* `learning_rate`
* `subsample`
* `colsample_bytree`
* `min_child_weight`
* `gamma`
* `reg_alpha`
* `reg_lambda`

La optimización utilizó como métrica principal **ROC-AUC**.

### 📐 Support Vector Classification

Se evaluó también un modelo **SVC**, estudiando dos configuraciones diferentes de filtrado de ASVs:

* Prevalencia del 5 %
* Prevalencia del 20 %

Para la optimización se exploraron parámetros como `C`, `gamma` y `class_weight`.

---

# 🔎 Explicabilidad

Una de las partes principales del proyecto consiste en estudiar **qué características utilizan los modelos para realizar sus predicciones**.

Para ello se utilizaron diferentes aproximaciones.

### SHAP

Se empleó **SHAP (SHapley Additive exPlanations)** para estudiar la contribución de las variables a las predicciones.

Para los modelos basados en árboles se utilizó:

```text
TreeExplainer
```

mientras que para SVC se empleó:

```text
KernelExplainer
```

El análisis se realizó mediante valores **out-of-fold**, respetando la estructura de las parejas.

### BorutaSHAP

En el caso del modelo GBT también se exploró **BorutaSHAP** para seleccionar variables relevantes de forma más robusta entre diferentes particiones.

Posteriormente, las ASVs identificadas se relacionaron con su clasificación taxonómica obtenida mediante SILVA.

> ⚠️ Debido a las limitaciones computacionales, algunos análisis BorutaSHAP no pudieron completarse en todos los folds. Por ello, las características identificadas deben considerarse **exploratorias** y no biomarcadores validados.

---

# 🧩 Clustering

Además del aprendizaje supervisado, se realizaron análisis de aprendizaje no supervisado para explorar la estructura del conjunto de datos.

Se utilizaron:

* **PCA**
* **Clustering jerárquico**
* **Dendrogramas**
* **Heatmaps**
* **K-means**
* **Silhouette Score**

Para K-means se evaluaron diferentes valores de `k`.

El análisis del coeficiente de silueta mostró que:

```text
k = 2
```

presentaba el mejor valor dentro del rango evaluado.

El coeficiente obtenido fue aproximadamente:

```text
Silhouette Score ≈ 0.4
```

lo que sugiere una estructura de agrupamiento observable, aunque con una separación entre grupos relativamente difusa.

---

# 🧠 Microbial General Model (MGM)

Como segunda aproximación se exploró un modelo Transformer basado en **Microbial General Model (MGM)**.

El procedimiento se estructuró en dos etapas:

```text
MGM preentrenado
       │
       ▼
MicroCorpus-260K
       │
       ▼
Domain-Adaptive Pretraining
       │
       ▼
NHANES Oral Microbiome
       │
       ▼
MGM-NHANES-DAPT
       │
       ▼
Fine-tuning supervisado
       │
       ▼
Cohorte Alzheimer
       │
       ▼
Clasificación binaria
```

### Domain-Adaptive Pretraining

Debido a que el corpus original del MGM presenta una mayor representación de microbiota intestinal, se realizó una adaptación al dominio de la microbiota oral utilizando datos de **NHANES Oral Microbiome 2009–2012**.

La adaptación permitió obtener un modelo inicial específico para el dominio oral antes del ajuste sobre la cohorte de Alzheimer.

### Fine-tuning

Posteriormente se realizó un ajuste supervisado sobre la cohorte de Alzheimer utilizando:

* `StratifiedGroupKFold`
* Agrupación por `pair_id`
* Optimización mediante Optuna
* ROC-AUC como métrica
* Early stopping

---

# 👁️ Interpretabilidad del Transformer

Se exploró inicialmente el uso de SHAP sobre el Transformer.

Sin embargo, **KernelExplainer presentó unos requerimientos computacionales demasiado elevados** para el entorno disponible, especialmente en memoria GPU.

Por este motivo, el análisis SHAP del MGM no pudo completarse para el conjunto de muestras.

Como alternativa exploratoria se analizaron los **pesos de atención del Transformer**.

Se estudiaron las atenciones de la última capa y se calcularon las diferencias de atención entre:

```text
Pacientes con EA
        vs.
Controles
```

Los tokens representan géneros bacterianos, permitiendo explorar qué taxones reciben una atención diferencial por parte del modelo.

> Este análisis debe interpretarse como una aproximación exploratoria a la interpretabilidad del modelo y no como una demostración de causalidad biológica.

---

# 🔐 Validación y prevención de data leakage

Uno de los aspectos metodológicos más importantes del proyecto es la estructura emparejada de los datos.

Cada paciente con Alzheimer está asociado a un control, por lo que dividir las muestras de forma aleatoria podría provocar que información relacionada con una pareja apareciese simultáneamente en entrenamiento y prueba.

Para evitarlo se utilizó:

```text
Nested Stratified Group Cross-Validation
```

con:

* **5 folds externos**
* **5 folds internos**
* Agrupación mediante `pair_id`
* Estratificación de las clases

De esta forma, las muestras pertenecientes a una misma pareja permanecen siempre dentro de la misma partición.

Esta estrategia busca reducir el riesgo de obtener estimaciones artificialmente optimistas debido a similitudes compartidas entre los miembros de una pareja.

---

# ⚠️ Limitaciones

Los resultados de este proyecto deben interpretarse teniendo muy presentes sus limitaciones.

### Tamaño de la cohorte

El conjunto contiene únicamente:

```text
101 muestras
```

frente a varios miles de características derivadas principalmente de las ASVs.

Esta elevada relación:

```text
nº de variables  >>>  nº de observaciones
```

incrementa considerablemente el riesgo de **overfitting** y dificulta la generalización de los modelos.

### Calidad de las secuencias

El número de ASVs observado por muestra fue elevado, lo que podría indicar que los controles de calidad empleados no fueron suficientemente restrictivos y que pueden haberse incorporado artefactos de secuenciación.

El propio trabajo plantea la necesidad de repetir y optimizar el procesamiento con diferentes parámetros y estrategias de control de calidad.

### Recursos computacionales

Los análisis se realizaron con recursos computacionales limitados.

Esto afectó especialmente a:

* Optimización de hiperparámetros.
* BorutaSHAP.
* KernelExplainer.
* Entrenamiento del Transformer.
* Análisis de atención.

En consecuencia, algunos experimentos tuvieron que ejecutarse de manera conservadora.

### Interpretabilidad

Los resultados de interpretabilidad presentaron una elevada variabilidad entre modelos.

Por ello, **no se recomienda interpretar las ASVs identificadas mediante SHAP como biomarcadores confirmados**. El propio trabajo señala la necesidad de realizar análisis adicionales con mayor potencia computacional y cohortes más amplias.

---

# 📈 Resultados y conclusiones

El proyecto permitió construir una infraestructura reproducible capaz de llevar a cabo el flujo completo:

```text
Secuencias
   ↓
QIIME 2
   ↓
ASVs
   ↓
Taxonomía SILVA
   ↓
Preprocesamiento
   ↓
Machine Learning
   ↓
Clustering
   ↓
Explicabilidad
   ↓
Transformer MGM
```

El pipeline generó correctamente la tabla de ASVs y la clasificación taxonómica mediante SILVA.

Sin embargo, las limitaciones relacionadas con el número de muestras, la elevada dimensionalidad y los recursos computacionales condicionaron el rendimiento y estabilidad de los modelos.

Por tanto, **la hipótesis del trabajo no puede confirmarse ni rechazarse con los datos disponibles**.

El principal resultado del proyecto es, por tanto, la construcción de una **infraestructura reproducible de bioinformática + machine learning + explicabilidad** que puede servir como punto de partida para futuros estudios con cohortes mayores y recursos computacionales más amplios.

---

# 🗂️ Estructura del proyecto

Una estructura recomendada para el repositorio es:

```text
.
├── README.md
│
├── data/
│   ├── raw/
│   ├── metadata/
│   └── processed/
│
├── pipeline/
│   ├── qiime2/
│   ├── cutadapt/
│   ├── deblur/
│   └── taxonomy/
│
├── notebooks/
│   ├── preprocessing/
│   ├── classification/
│   ├── clustering/
│   ├── shap/
│   └── mgm/
│
├── src/
│   ├── preprocessing/
│   ├── models/
│   ├── explainability/
│   └── utils/
│
├── results/
│   ├── figures/
│   ├── models/
│   └── tables/
│
├── requirements.txt
└── LICENSE
```

> La estructura anterior representa una organización recomendada del repositorio; debe adaptarse a la estructura real de los archivos incluidos en el proyecto.

---

# 🛠️ Tecnologías

| Tecnología                    | Uso                             |
| ----------------------------- | ------------------------------- |
| 🧬 **QIIME 2**                | Procesamiento de microbioma     |
| 🦠 **Deblur**                 | Generación de ASVs              |
| 🧪 **SILVA**                  | Clasificación taxonómica        |
| ✂️ **Cutadapt**               | Eliminación de cebadores        |
| 🐍 **Python**                 | Análisis y Machine Learning     |
| 💻 **Bash**                   | Automatización de la pipeline   |
| 🌳 **Gradient Boosted Trees** | Clasificación                   |
| 📐 **SVC**                    | Clasificación                   |
| 🔎 **SHAP**                   | Explicabilidad                  |
| 🧬 **BorutaSHAP**             | Selección de variables          |
| 📊 **PCA / K-means**          | Clustering                      |
| 🤗 **Transformer / MGM**      | Modelado de microbioma          |
| ⚙️ **Optuna**                 | Optimización de hiperparámetros |
| 📓 **Jupyter**                | Desarrollo y experimentación    |

---

# 🔁 Reproducibilidad

Uno de los objetivos fundamentales del proyecto fue desarrollar una infraestructura reproducible.

La pipeline mantiene separadas las diferentes etapas del procesamiento y conserva los resultados intermedios, permitiendo:

* Repetir el procesamiento de las secuencias.
* Modificar parámetros de control de calidad.
* Reentrenar los modelos.
* Comparar diferentes configuraciones.
* Incorporar nuevas cohortes.
* Extender los análisis de explicabilidad.

La automatización mediante Bash permite mantener constantes los parámetros utilizados y reducir errores derivados de la ejecución manual.

---

# 🚀 Trabajo futuro

Los resultados obtenidos permiten plantear varias líneas de trabajo futuro:

* 📚 Incorporar cohortes significativamente mayores.
* 🧬 Repetir el pipeline con controles de calidad más estrictos.
* 🔬 Evaluar diferentes parámetros de filtrado y generación de ASVs.
* 🌍 Validar los modelos en cohortes externas.
* ⚖️ Estudiar posibles factores de confusión.
* 🧪 Analizar por separado saliva y exudado supragingival.
* 🤖 Realizar una optimización más exhaustiva de los modelos.
* 🔎 Repetir los análisis SHAP y BorutaSHAP con mayor capacidad computacional.
* 🧠 Profundizar en la interpretabilidad del MGM.
* 🧬 Investigar la estabilidad de las características microbianas identificadas entre cohortes.
* 🔁 Evaluar la reproducibilidad de los posibles candidatos a biomarcadores.

---

# 📚 Referencias principales

El trabajo se fundamenta en literatura relacionada con:

* Enfermedad de Alzheimer y demencia.
* Biomarcadores.
* Microbioma oral.
* Secuenciación 16S rRNA.
* QIIME 2.
* ASVs y Deblur.
* Clasificación taxonómica mediante SILVA.
* Machine Learning aplicado a datos biológicos.
* SHAP y métodos de explicabilidad.
* BorutaSHAP.
* Microbial General Model (MGM).
* Domain-Adaptive Pretraining.

Las referencias bibliográficas completas se encuentran en la memoria del Trabajo Fin de Grado.

---

# ⚖️ Aviso científico

Este repositorio tiene una finalidad **académica y de investigación exploratoria**.

Los resultados obtenidos en este proyecto **no constituyen una herramienta diagnóstica**, no deben utilizarse para diagnosticar la enfermedad de Alzheimer y las características microbianas identificadas no deben considerarse biomarcadores clínicamente validados.

El proyecto debe entenderse como una **prueba de concepto metodológica** destinada a evaluar la viabilidad de integrar bioinformática, microbioma, aprendizaje automático y explicabilidad.

---

## 👨‍💻 Autor

**Daniel de Lara Pérez**

Grado en Bioinformática
Universidad Católica de Ávila

---

<p align="center">
  <b>🧬 Oral Microbiome × 🤖 Machine Learning × 🧠 Alzheimer's Disease</b>
</p>

<p align="center">
  <i>Explorando el potencial del microbioma oral mediante bioinformática e inteligencia artificial.</i>
</p>
