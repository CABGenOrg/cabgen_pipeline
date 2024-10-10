import schedule
from time import sleep
from typing import List
from os import getenv, path
from dotenv import load_dotenv
from src.utils.handle_log import logging_conf
from src.utils.handle_errors import fatal_error
from concurrent.futures import ProcessPoolExecutor, wait
from src.models.CabgenPipeline import CabgenPipeline
from src.utils.handle_tasks import get_fastqc_tasks, get_complete_tasks, \
    get_genomic_tasks

load_dotenv()


def process_task(task: dict, mode: str):
    try:
        output_path = getenv("UPLOADED_SEQUENCES_PATH") or ""
        if not output_path:
            raise ValueError("There is no output path.")

        log_path = getenv("LOG_PATH") or ""
        if not log_path:
            raise ValueError("There is no log path.")

        sample = int(task.get("_id", 0))
        read1 = task.get("arquivofastqr1", "")
        read2 = task.get("arquivofastqr2", "")
        output = path.join(output_path, f"output_{sample}")

        logger = logging_conf(sample, log_path)
        if not logger:
            print("Logger could not be started")

        pipe = CabgenPipeline(sample, read1, read2,
                              output, logger)  # type: ignore

        if mode == "fastqc":
            pipe.run(only_fastqc=True)
        elif mode == "complete":
            pipe.run(complete=True)
        elif mode == "genomic":
            pipe.run(only_genomic=True)
        else:
            raise ValueError("Invalid mode provided for task processing.")

    except Exception as e:
        print(f"Failed to process task {sample}.\n\n{e}")


def process_tasks_in_parallel(tasks: List[dict], mode: str):
    try:
        with ProcessPoolExecutor(max_workers=4) as executor:
            futures = [executor.submit(process_task, task, mode)
                       for task in tasks]
            wait(futures)
    except Exception as e:
        print(f"Failed to process tasks in parallel.\n\n{e}")


def pipeline_job():
    try:
        fastqc_tasks = get_fastqc_tasks()
        complete_tasks = get_complete_tasks()
        genomic_tasks = get_genomic_tasks()

        if fastqc_tasks:
            print(f"Processing {len(fastqc_tasks)} FastQC tasks...")
            process_tasks_in_parallel(fastqc_tasks, "fastqc")

        if complete_tasks:
            print(f"Processing {len(complete_tasks)} Complete tasks...")
            process_tasks_in_parallel(complete_tasks, "complete")

        if genomic_tasks:
            print(f"Processing {len(genomic_tasks)} Genomic tasks...")
            process_tasks_in_parallel(genomic_tasks, "genomic")
    except Exception as e:
        print(f"Failed to run pipeline_job.\n\n{e}")


def main():
    try:
        timeout = 5
        schedule.every(timeout).minutes.do(pipeline_job)

        while True:
            schedule.run_pending()
            sleep(20)
    except Exception as e:
        fatal_error(f"Failed to run main function.\n\n{e}")


if __name__ == "__main__":
    main()
