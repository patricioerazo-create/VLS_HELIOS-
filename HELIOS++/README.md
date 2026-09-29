# HELIOS++ — simulación Ruil_MLS (tesis)

Proyecto de simulación LiDAR **MLS** en el rodal **Ruil** (Estación Experimental Frutillar): escena HELIOS++, survey Stonex/VLP-16, voxelización 3.5 cm y pipeline de métricas por árbol.

## Documentación

| Recurso | Ubicación |
|---------|-----------|
| **Metodología paso a paso** | [docs/METODOLOGIA_SIMULACION_Ruil_MLS.md](docs/METODOLOGIA_SIMULACION_Ruil_MLS.md) |
| Escena | `data/scenes/Ruil_MLS.xml` |
| Survey (13 legs) | `data/surveys/Ruil_MLS.xml` |
| Job Slurm | `run_Ruil_MLS.slurm` |
| Scripts métricas | `Metricas/` |

## Inicio rápido

```bash
cd "HELIOS++"
export HELIOS_ROOT="$(pwd)"
mkdir -p output
```

Instalá HELIOS++ (Conda recomendado): [3dgeo-heidelberg/helios](https://github.com/3dgeo-heidelberg/helios). Luego seguí la guía en `docs/METODOLOGIA_SIMULACION_Ruil_MLS.md`.

## Licencia

Motor HELIOS++: [3dgeo-heidelberg/helios](https://github.com/3dgeo-heidelberg/helios) (LGPL). Ver `COPYING`, `COPYING.LESSER` y `LICENSE.md` en este directorio.
