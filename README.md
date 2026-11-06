# IRP mono-véhicule et faits stylisés de la distribution de carburant en Guinée

Code, instances et résultats de l'article :

> M. Maomou, B. M. Touré, B. Mansaré, A. Sonko. *Un modèle de programmation linéaire en nombres entiers pour le problème de tournées avec gestion des stocks : validation numérique sur des instances calibrées par faits stylisés de la distribution de carburant en Guinée.* (soumis)

Affiliations : Département de Mathématiques, Université Gamal Abdel Nasser de Conakry ; Département de Mathématiques, Université de N'Zérékoré (Guinée).

## Contenu
| Dossier / fichier | Rôle |
|---|---|
| `code/gen.py` | Générateur d'instances calibrées par faits stylisés (Basse-Guinée, Guinée forestière, familles témoins) |
| `code/irp.py` | Modèle PLNE et branch-and-cut (SCIP), inégalités de fenêtre, variante MTZ |
| `code/enum_check.py`, `code/validate.py` | Validation par énumération exhaustive (120 petites instances) |
| `code/irp_mv.py`, `code/validate_bench.py` | Validation sur les instances de référence d'Archetti et al. (2007), version à 2 camions |
| `code/run_experiments.py`, `code/solve_one.py` | Expérimentations (288 instances, inégalités de fenêtre, sensibilité) |
| `code/analyze.py` | Tableaux LaTeX (`tables/`) et figure (`figures/`) |
| `results_*.jsonl` | Résultats bruts (une ligne JSON par instance résolue) |

## Installation
Python ≥ 3.10, puis :

    pip install pyscipopt networkx numpy scipy matplotlib openpyxl

Le solveur SCIP est libre ; aucune licence commerciale n'est nécessaire.

## Reproduire les résultats
Depuis le dossier `code/` :

    python3 validate.py ../results_validation.jsonl
    git clone https://github.com/vidalt/HGS-IRP.git ../HGS-IRP
    python3 validate_bench.py ../HGS-IRP ../results_bench.jsonl
    python3 run_experiments.py main ../results_main.jsonl 2
    python3 run_experiments.py windows ../results_windows.jsonl 2
    python3 run_experiments.py sens ../results_sens.jsonl 2
    python3 analyze.py

Limite de 600 s par instance, un fil d'exécution (Intel Xeon 2,1 GHz, 8 Go de mémoire).

Résoudre une seule instance :

    python3 solve_one.py '{"fam":"BG","N":10,"H":3,"eta":"L","rep":1,"which":"sty","windows":true,"TL":60}'

## Sources des données
- Populations : Institut National de la Statistique de Guinée, RGPH-3 (2014).
- Coordonnées des localités : GeoNames.
- Instances de référence : Archetti et al. (2007), fichiers redistribués par le dépôt HGS-IRP ; valeurs de référence du DIMACS Implementation Challenge (2022).

## Licence
Code sous licence MIT (voir `LICENSE`).
