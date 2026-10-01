# Woher diese Kopie stammt

Code und Doku in diesem Ordner (`analysis/`, `scripts/`, `interface/README.md`,
`calibrations/README.md`, `report.md`) sind eine **Kopie** aus AKHExperiment
(`Akhelab/AKHExperiment/lina/`). Bearbeitet wird dort, nicht hier.
`tools/measure_tau_aux.py` importiert `analysis/lina_ofdr.py` aus dieser Kopie.

| Stand | Quelle | Was dazukam |
|---|---|---|
| 2026-09-01 | `carolina-ofdr` | erste Kopie, Mk1 LINEAR |
| 2026-10-01 | `carolina-ofdr` @ `ea2fde5` **plus nicht committete Aenderungen** im Arbeitsverzeichnis | Mk2-LOG-Umstellung; Aux-Achse in `lina_wl_cal.py`; neu `analysis/lina_detector.py`, `scripts/lina_chip_measure.py`, `scripts/lina_voa_series.py`, `scripts/lina_npz_to_json.py` |

Bewusst **nicht** uebernommen:

- `.gitignore`: hier werden die Messdaten in `data/` eingecheckt, drueben nicht.
- `interface/config/TunableFilter_config.json`: drueben nur andere Zeilenenden.

Abgleich pruefen:

    diff -r --strip-trailing-cr logs/lina/analysis <AKHExperiment>/lina/analysis
