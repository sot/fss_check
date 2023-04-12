import os
from pathlib import Path
from ska_helpers.utils import LazyDict
import yaml

SKA = Path(os.environ["SKA"])
CONFIG_DIR = "."


def _load_config_file():
    """Load fss_check_config.yaml from current dir or $SKA/data/fss_check"""
    path = Path(CONFIG_DIR) / "fss_check_config.yml"
    with open(path) as fh:
        config = yaml.safe_load(fh)

    config['config_path'] = str(path.absolute())

    return config


CONFIG = LazyDict(load_func=_load_config_file)
