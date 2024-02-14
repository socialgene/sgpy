import requests
from rich.progress import (
    BarColumn,
    DownloadColumn,
    Progress,
    TextColumn,
    TimeRemainingColumn,
    TransferSpeedColumn,
)

progress = Progress(
    TextColumn("[bold blue]{task.fields[filename]}", justify="right"),
    BarColumn(bar_width=None),
    "[progress.percentage]{task.percentage:>3.1f}%",
    "•",
    DownloadColumn(),
    "•",
    TransferSpeedColumn(),
    "•",
    TimeRemainingColumn(),
    transient=True,
)


def download(url, filepath):
    res = requests.get(
        "https://www.npatlas.org/static/downloads/NPAtlas_download.json", stream=True
    )
    if "Content-length" in res.headers:
        total = int(res.headers["Content-length"])
    else:
        total = 0
    tt = 0
    with progress:
        task_id = progress.add_task(
            "download", filename="NPAtlas_download.json", total=total
        )
        with open(filepath, "wb") as f:
            for i in res.iter_content(8192):
                _ = f.write(i)
                progress.update(task_id, advance=len(i))
                tt += len(i)
