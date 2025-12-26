import sys
from pathlib import Path
import yaml

src = Path(sys.argv[1])
dst = Path(sys.argv[2])
dataset_path = sys.argv[3]
type_names = sys.argv[4].split(",")

cfg = yaml.safe_load(src.read_text())

# よくあるキーに対して上書き（見つかったものだけ更新）
def set_if_exists(d, key, val):
    if isinstance(d, dict) and key in d:
        d[key] = val
        return True
    return False

# run出力先をこのrepo内に寄せたい
for k in ["root", "output_dir", "workdir"]:
    set_if_exists(cfg, k, "allegro/runs")

# dataset file名の候補
updated = False
for k in ["dataset_file_name", "dataset_filename", "dataset_path", "dataset"]:
    if set_if_exists(cfg, k, dataset_path):
        updated = True

# type_namesの候補（dict階層にもあるので再帰）
def replace_type_names(obj):
    if isinstance(obj, dict):
        for k,v in obj.items():
            if k == "type_names" and isinstance(v, (list, tuple)):
                obj[k] = type_names
            else:
                replace_type_names(v)
    elif isinstance(obj, list):
        for x in obj:
            replace_type_names(x)

replace_type_names(cfg)

# 軽量化：epoch系があれば減らす
for k in ["max_epochs", "epochs", "n_epochs"]:
    set_if_exists(cfg, k, 5)

# もしdatasetキーが見つからなかったら、ユーザが手で直すためにメモを残す
cfg.setdefault("_sg_note", {})
cfg["_sg_note"]["dataset_path_expected"] = dataset_path
cfg["_sg_note"]["type_names"] = type_names
cfg["_sg_note"]["edited"] = str(src)

dst.write_text(yaml.safe_dump(cfg, sort_keys=False))
print(f"Wrote {dst}")
print(f"NOTE: if training fails due to dataset key mismatch, open {dst} and set dataset file key to {dataset_path}")
