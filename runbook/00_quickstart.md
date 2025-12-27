# Quickstart

## Goal
Run QE → Allegro → LAMMPS end-to-end with a small Li/LLZO case.

## One-command run
make quickstart

## Success criteria
- dataset/ contains .extxyz
- allegro/runs/<timestamp>/ contains a trained model
- lammps/runs/<timestamp>/ contains logs/trajectory

> NOTE: This quickstart currently uses **synthetic extxyz** for plumbing validation. Next step is to swap in **QE-derived Li–LLZO small interface** extxyz.
