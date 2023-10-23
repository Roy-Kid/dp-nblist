import logging
import pytest
import tempfile
import tarfile
from pathlib import Path
from typing import Optional
from urllib.request import urlretrieve
log = logging.getLogger(__name__)

def fetch(self, url, filename, dir: Optional[Path | str] = None) -> Path:
    if dir is None:
        dir = self.data_dir

    fpath = Path(dir) / Path(filename)
    if not fpath.exists():
        log.info(f"downloading from {url}... to {fpath}")
        urlretrieve(url, fpath)
    else:
        log.info(f"{fpath} already exists.")
    return fpath

@pytest.fixture(scope="session")
def qm9():
    
    data_dir = Path(tempfile.gettempdir()) / Path("qm9")
    if not data_dir.exists():
        data_dir.mkdir()
        url = "https://ndownloader.figshare.com/files/3195389"
        tar_path = fetch(url, "gdb9.tar.gz")
        raw_path = data_dir / Path("gdb9_xyz")

        if not raw_path.exists():
            log.info("Extracting files...")
            tar = tarfile.open(tar_path)
            tar.extractall(raw_path)
            tar.close()
            log.info("Done.")

        log.info("Parse xyz files...")
        ordered_files = sorted(
            raw_path.rglob("*.xyz"),
            key=lambda x: (int(re.sub(r"\D", "", str(x))), str(x)),
        )  # sort by index in filename
        return ordered_files
    