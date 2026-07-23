# Dependencies

This document lists the dependencies of the `gws_microbial_genomics` brick: the Constellab bricks it requires, and the external bioinformatics tools used by its tasks.

## Constellab bricks

Declared in [`settings.json`](./settings.json):

| Brick | Version |
|---|---|
| `gws_core` | >= 0.19.3 |

## Task environments

`settings.json` declares no pip or git packages: each task that wraps an external bioinformatics tool creates its own isolated conda environment at runtime (via `MambaShellProxy`/`CondaShellProxy`), defined by an env file next to the task's source file.

| Task | Env file | Key packages |
|---|---|---|
| Build/Update Bakta Database, Prokaryotes Genome Annotation (`bakta/bakta_db_task.py`, `bakta/bakta.py`) | `bakta/Bakta_env.yml` | `bakta==1.11.4`, `pycirclize==1.10.0`, `matplotlib==3.10.6`, `biopython==1.85`, `mpld3==0.5.11` |
| eggNOG Mapper, eggNOG DB download (`eggnog-mapper/eggnog_emapper.py`, `eggnog-mapper/eggnog_db_download.py`) | `eggnog-mapper/Eggnog_env.yml` | `eggnog-mapper==2.1.12` |
| Bactopia_run, Bactopia_metadata (`Bactopia/bactopia_run.py`, `Bactopia/bactopiaMetadata.py`) | `Bactopia/bactopia_env.yml` | `bactopia==3.2.0` |

`Bactopia_explore` (`Bactopia/bactopia_explore.py`) only parses the output folder produced by `Bactopia_run` and does not create its own environment.
