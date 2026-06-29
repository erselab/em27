# EM27/SUN ↔ GERT

GERT's first real-data example project. It compares **GERT** radiative-transfer
simulations and retrievals against real **EM27/SUN** (COCCON) measurements and
the **PROFFAST L2** product, closing the full loop:

```
real .BIN spectra ──► gert.ForwardModel (TransmissionSolver) ──► simulated spectrum
                                                                      │
        PROFFAST invparms (XCO2/XCH4/XCO) ◄── gert.GERTRetrieval ◄────┘   compare
```

This project **depends on `gert` as an installed library** — it ships only its
own data readers, the EM27 instrument definition, and the scene builder. All RT
and retrieval physics lives in `gert`.

## Data

`data/GSFC_SN245_260406/` — EM27/SUN SN245 at NASA GSFC, 2026-04-06: calibrated
spectra (`260406_spectra/cal/*.BIN`), a-priori `.map` profiles, ground pressure,
ILS (`ils_list.csv`), and the PROFFAST L2 `comb_invparms_*.csv`.

## Setup

```bash
pip install -e .            # installs em27gert + pulls in gert
# gert's Fortran extension: build once in the gert repo (cd ../../gert/fortran && make ext)
```

You supply your own ABSCO and solar tables (paths via the gert loaders).

## Status

See [PLAN.md](PLAN.md) for the milestone plan (M0 ingest → M4 closed-loop
time series → M5 RT-fidelity experiments).

- `em27gert/proffast_bin.py` — PROFFAST/OPUS `.BIN` reader (seeded from gert).
- `notebooks/` — seeded from the GERT TCCON notebooks; being adapted to EM27.