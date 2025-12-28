# QE → Allegro → LAMMPS (Li/LLZO-like plumbing demo)

Reproducible end-to-end demo:
1) QE (short sampling / AIMD)
2) Export dataset to extxyz (E/F)
3) Allegro training
4) LAMMPS MD using the trained model

## Quickstart (goal)
- make quickstart produces:
  - dataset/*.extxyz
  - allegro/runs/<timestamp>/
  - lammps/runs/<timestamp>/
  - artifacts/<timestamp>/

This public demo is intentionally small.
Scaling (atoms/ns/sweeps/on-prem/air-gapped) is covered via Founding 5 / support.

## Status (important)
- Current: **synthetic extxyz** is used to validate the end-to-end plumbing (dataset → train → compile → LAMMPS run).
- Next: replace dataset generation with **QE-derived Li–LLZO small interface** extxyz.

