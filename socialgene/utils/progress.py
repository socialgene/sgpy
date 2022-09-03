from rich.progress import Progress, BarColumn, DownloadColumn
import time


class MyProgress(Progress):
    def get_renderables(self):
        for task in self.tasks:
            if task.fields.get("progress_type") == "mygreenbar":
                self.columns = ("[green]Rich is awesome!", BarColumn())
            if task.fields.get("progress_type") == "mybluebar":
                self.columns = (
                    "[blue]Another bar with a different layout",
                    BarColumn(bar_width=None),
                    "•",
                    DownloadColumn(),
                )
            yield self.make_tasks_table([task])


with MyProgress() as progress:
    taskone = progress.add_task("taskone", progress_type="mygreenbar")
    tasktwo = progress.add_task("tasktwo", progress_type="mybluebar")


def a():
    aa = range(0, 10)
    with MyProgress() as progress:
        task = progress.add_task("taskone", progress_type="mygreenbar")
        for i in aa:
            time.sleep(0.5)
            progress.update(task, advance=1)


def b():
    aa = range(0, 10)
    with Progress(transient=True) as pg:
        task = pg.add_task("Progress...", total=len(aa))
        for i in aa:
            time.sleep(0.5)
            pg.update(task, advance=1)


with MyProgress() as progress:
    taskone = progress.add_task("taskone", progress_type="mygreenbar")
    tasktwo = progress.add_task("tasktwo", progress_type="mybluebar")
