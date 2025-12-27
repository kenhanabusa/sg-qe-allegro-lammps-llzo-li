.PHONY: quickstart qe dataset train deploy lammps

quickstart: qe dataset train deploy lammps

qe:
	@bash scripts/10_run_qe.sh

dataset:
	@bash scripts/20_make_dataset.sh

train:
	@bash scripts/30_train_allegro.sh

deploy:
	@bash scripts/40_deploy_lammps.sh

lammps:
	@bash scripts/50_run_lammps.sh
