# VLS_HELIOS-

Repositorio de la tesis **Simulación LiDAR a partir de nube de puntos de *Nothofagus alessandrii*** — Estación Experimental Frutillar, X Región de Los Lagos, Chile.

Aquí está el material para **reproducir la simulación MLS del rodal Ruil** con [HELIOS++](https://github.com/3dgeo-heidelberg/helios): preparación de vox, escena, survey, ejecución y métricas comparativas (campo / MLS / simulado).

---

## Estructura

```
VLS_HELIOS-/
└── HELIOS++/          ← todo el proyecto (datos, scripts, XML, vox)
    ├── docs/          ← metodología paso a paso
    ├── data/          ← DEM, vox, escenas, surveys
    ├── Metricas/      ← pipeline Python
    └── run_Ruil_MLS.slurm
```

---

## Clonar y empezar

```bash
git clone https://github.com/patricioerazo-create/VLS_HELIOS-.git
cd VLS_HELIOS-
cd "HELIOS++"
```

**Guía completa:** [`HELIOS++/docs/METODOLOGIA_SIMULACION_Ruil_MLS.md`](HELIOS%2B%2B/docs/METODOLOGIA_SIMULACION_Ruil_MLS.md)

---

## Requisitos (resumen)

- Entorno Conda con **HELIOS++** (`helios++`, `pyhelios`)
- Python: `laspy`, `numpy` (scripts en `Metricas/`)
- Opcional: cluster Slurm para la corrida rodal (`run_Ruil_MLS.slurm`)

---

## Licencia

HELIOS++ es software de [3dgeo-heidelberg/helios](https://github.com/3dgeo-heidelberg/helios) (LGPL). Los archivos en `HELIOS++/` incluyen las licencias del upstream (`COPYING`, `LICENSE.md`).
