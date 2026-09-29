# Metodología paso a paso — simulación HELIOS++ Ruil_MLS

Proyecto: rodal **Ruil**, escaneo **MLS** (Stonex X120GO + VLP-16 simulado).  
Objetivo: reproducir la nube MLS con HELIOS++ y comparar métricas por árbol (3001–3079) con terreno y campo.

Repositorio de la tesis: [patricioerazo-create/VLS_HELIOS-](https://github.com/patricioerazo-create/VLS_HELIOS-)

---

## Convención de rutas

El código y los datos de la simulación viven en la carpeta **`HELIOS++/`** dentro del repositorio.  
Todas las rutas de este documento son **relativas a `HELIOS++/`** (ahí están `data/`, `Metricas/`, `output/`, etc.).

```bash
git clone https://github.com/patricioerazo-create/VLS_HELIOS-.git
cd VLS_HELIOS-
cd "HELIOS++"
export HELIOS_ROOT="$(pwd)"   # opcional; útil en scripts y Slurm
```

En los ejemplos se asume que ya estás dentro de **`HELIOS++/`**.

---

## 1. Requisitos

| Componente | Detalle |
|------------|---------|
| **HELIOS++** | Entorno Conda `helios` con `pyhelios` y binario `helios++` |
| **Python** | `laspy`, `numpy` (scripts en `Metricas/`) |
| **Cluster** | Slurm (opcional), partición según tu cluster, ~512 GB RAM y 16 CPUs en el job de referencia |
| **Herramientas LAZ** | `lasmerge`, `las2las` (LASzip) en el mismo entorno |

---

## 2. Datos de entrada

### 2.1 Inventario y centros

- **Centros XY por árbol:** `Metricas/arboles_plot_metrics_xy.csv`
- **Métricas de campo (referencia):** `Metricas/Metricas terreno.csv` (DAP, altura, radios N/E/S/O)

### 2.2 Terreno

- **DEM 2 m:** `data/sceneparts/tif/DEM_Ruiles_2m.tif` (+ `DEM_Ruiles_2m.prj`)
- Material suelo en escena: `data/sceneparts/basic/groundplane/groundplane.mtl`

### 2.3 Nubes MLS por árbol (pre-simulación)

- **Originales segmentadas (MLS / Stonex):** `Metricas/Nubes Terreno/` — archivos `{3001..3079}.las` o `.laz`
  - Si las nubes están en otra subcarpeta de `HELIOS++/`, indicá la ruta con `--input` en el paso 3.

### 2.4 Referencia de trayectoria (documentación)

- GeoJSON de campo: `data/trayectorias/Panoscanner - Trayectoria 1.geojson`  
  **No** lo lee HELIOS; los waypoints oficiales están en el survey XML (paso 5).

---

## 3. Alinear nubes MLS al DEM

**Script:** `Metricas/alinear_las_dem.py`

**Qué hace:** Para cada árbol, estima la base (XY de puntos cerca de `z_min`) y desplaza **Z** para que `z_min` coincida con la cota del DEM 2 m (interpolación bilinear).

**Salida:** `Metricas/Nubes Reales/{id}.las`  
**Registro:** `Metricas/las_dem_alignment.csv`

```bash
cd "HELIOS++"
python3 Metricas/alinear_las_dem.py \
  --input "Metricas/Nubes Terreno" \
  --output "Metricas/Nubes Reales" \
  --dem data/sceneparts/tif/DEM_Ruiles_2m.tif
```

Opcional: `--dry-run` para revisar deltas sin escribir LAS.

---

## 4. Voxelizar árboles para HELIOS++

**Script:** `Metricas/Voxelizacion.py`

**Qué hace:** Discretiza cada nube alineada en voxels 3D y escribe formato `VOXEL SPACE` con `PadBVTotal=1.0` (voxel opaco).

| Parámetro producción | Valor |
|----------------------|-------|
| Tamaño voxel (`--vsize`) | **0.035 m** (3.5 cm) |
| Opacidad (`--pad`) | **1.0** |
| Entrada (default) | `Metricas/Nubes Reales/` |
| Salida (default) | `data/sceneparts/vox/{id}.vox` |

```bash
cd "HELIOS++"
python3 Metricas/Voxelizacion.py --vsize 0.035 --pad 1.0
# Un árbol de prueba:
python3 Metricas/Voxelizacion.py --vsize 0.035 --only 3001
```

**Notas:**

- IDs **3000** y **3032** no forman parte del inventario estándar de simulación (3032 = segmentación; omitir en CSVs si aplica).
- La escena carga todos los vox con patrón `30[0-9][0-9].vox` en `data/sceneparts/vox/`.

---

## 5. Definir escena y survey HELIOS++

### 5.1 Escena — `data/scenes/Ruil_MLS.xml`

- **Geotiffloader:** DEM 2 m (`data/sceneparts/tif/...`, rutas relativas a `--assets`)
- **detailedvoxels:** `data/sceneparts/vox/30[0-9][0-9].vox`

### 5.2 Survey — `data/surveys/Ruil_MLS.xml`

| Elemento | Configuración |
|----------|----------------|
| Plataforma | `data/platforms.xml#x120go_handheld` (montaje **45°** sobre eje X) |
| Escáner | `data/scanners_tls.xml#vlp16` |
| Escena | `data/scenes/Ruil_MLS.xml#Ruil_MLS` |
| **Legs** | **13** segmentos de trayectoria |
| Altura sensor | **z = 1.3 m** (`onGround="true"`) |
| Velocidad | **0.5 m/s** (`movePerSec_m`) |
| Template escáner | `x120go_v3`: 18 750 Hz, rotación 360°/s, muestreo trayectoria 0.05 s |

Coordenadas de inicio de cada leg: UTM en el XML (véase `platformSettings x/y`).

### 5.3 Catálogos

- `data/platforms.xml` — plataforma Stonex handheld
- `data/scanners_tls.xml` — modelo VLP-16

---

## 6. Ejecutar la simulación

HELIOS++ resuelve rutas del XML y de `--assets` respecto al directorio **`HELIOS++/`** (raíz del proyecto de simulación).

### 6.1 Línea de comando (equivalente al job Slurm)

```bash
cd "HELIOS++"
helios++ data/surveys/Ruil_MLS.xml \
  --lasOutput \
  --assets "$(pwd)" \
  --output "$(pwd)/output" \
  --rebuildScene \
  --njobs 16
```

Si definiste `HELIOS_ROOT`:

```bash
cd "HELIOS++"
helios++ data/surveys/Ruil_MLS.xml \
  --lasOutput \
  --assets "$HELIOS_ROOT" \
  --output "$HELIOS_ROOT/output" \
  --rebuildScene \
  --njobs 16
```

### 6.2 Job Slurm — `run_Ruil_MLS.slurm`

Antes de enviar el job, editá en el script las rutas de Conda y binarios según tu máquina. La variable `HELIOS_ROOT` debe apuntar a **tu copia de `HELIOS++/`**.

Flujo del job de referencia:

1. Carga Conda y ejecuta `helios++` con `--assets` = raíz `HELIOS++/` y `--output` = `HELIOS++/output`.
2. Localiza la carpeta de corrida bajo `output/MLS/20*/`.
3. Fusiona `leg*_points.las` → `output/Ruil_MLS.las` con `lasmerge`.
4. Comprime → **`output/Ruil_MLS.laz`** (LAS 1.4) con `las2las`.
5. Elimina subcarpetas temporales `output/MLS/`.

```bash
cd "HELIOS++"
sbatch run_Ruil_MLS.slurm
```

**Salida principal:** `output/Ruil_MLS.laz` (nube rodal simulada; suele estar en `.gitignore`).

---

## 7. Post-proceso: métricas por árbol (comparación, no parte del motor HELIOS)

Tras **recortar o exportar** una nube por árbol desde `output/Ruil_MLS.laz`, guardá los archivos por ID en `Metricas/Nubes HELIOS/` (o en `Metricas/Nubes Terreno/` para MLS medido).

### 7.1 Altura

```bash
cd "HELIOS++"
python3 Metricas/altura.py --set helios
python3 Metricas/altura.py --set terreno
```

Por defecto leen las carpetas `Metricas/Nubes HELIOS` y `Metricas/Nubes Terreno` respectivamente.

### 7.2 Radios sectoriales y superficie de copa

**Script:** `Metricas/Metricas.py`

- Centro: `Metricas/arboles_plot_metrics_xy.csv`
- Radios en azimuts **40°, 130°, 220°, 310°** (sectores ±25°)
- Superficie de copa: proyección XY en **rejilla 10×10 cm**

```bash
cd "HELIOS++"
python3 Metricas/Metricas.py --set terreno
python3 Metricas/Metricas.py --set helios
python3 Metricas/Metricas.py --set terreno --input "Metricas/Nubes Terreno"
python3 Metricas/Metricas.py --set helios --input "Metricas/Nubes HELIOS"
```

**CSV generados:**

- `Metricas/terreno/metricas.csv`
- `Metricas/helios/metricas.csv`

Figuras opcionales en `Metricas/terreno/figuras/` y `Metricas/helios/figuras/`.

### 7.3 Suelo / 3DFin (si aplica)

```bash
cd "HELIOS++"
python3 Metricas/agregar_suelo.py --set helios
python3 Metricas/agregar_suelo.py --set terreno
```

---

## 8. Flujo resumido (diagrama)

Rutas relativas a **`HELIOS++/`**.

```mermaid
flowchart LR
  A[Metricas/Nubes Terreno] --> B[alinear_las_dem.py]
  B --> C[Metricas/Nubes Reales]
  C --> D[Voxelizacion.py 3.5cm]
  D --> E[data/sceneparts/vox]
  F[data/sceneparts/tif DEM] --> G[data/scenes/Ruil_MLS.xml]
  E --> G
  H[data/surveys/Ruil_MLS.xml] --> I[helios++]
  G --> I
  I --> J[output/Ruil_MLS.laz]
  J --> K[Metricas/Nubes HELIOS]
  K --> L[Metricas.py / altura.py]
  M[Metricas/Metricas terreno.csv] --> N[Comparación]
  L --> N
```

---

## 9. Checklist de reproducibilidad

- [ ] Clonar [VLS_HELIOS-](https://github.com/patricioerazo-create/VLS_HELIOS-) y `cd "HELIOS++"`
- [ ] DEM + `.prj` en `data/sceneparts/tif/`
- [ ] Nubes terreno por ID en `Metricas/Nubes Terreno/` (o `--input` equivalente)
- [ ] Nubes Reales en `Metricas/Nubes Reales/` y `Metricas/las_dem_alignment.csv`
- [ ] Archivos `.vox` en `data/sceneparts/vox/` (79 o subset)
- [ ] Survey en `data/surveys/Ruil_MLS.xml` (13 legs, z=1.3 m, 0.5 m/s)
- [ ] Corrida HELIOS++ → `output/Ruil_MLS.laz`
- [ ] Centros en `Metricas/arboles_plot_metrics_xy.csv`
- [ ] Métricas en `Metricas/terreno/` y `Metricas/helios/` vs campo

---

## 10. Artefactos que pueden estar fuera del repo

Tras una limpieza para transferencia, pueden faltar en el disco pero documentarse en anexos:

- `Metricas/Nubes Terreno`, `Metricas/Nubes Reales`, `Metricas/Nubes HELIOS`
- `output/Ruil_MLS.laz` (grande; conservar copia aparte)

Los **`.vox`**, XML de escena/survey y scripts en `Metricas/` permiten **re-simular** sin re-voxelizar; re-voxelizar exige de nuevo las LAS en `Metricas/Nubes Terreno/`.

---

*Documento alineado con `run_Ruil_MLS.slurm`, `data/scenes/Ruil_MLS.xml`, `data/surveys/Ruil_MLS.xml` y scripts en `Metricas/` (dentro de `HELIOS++/`).*
