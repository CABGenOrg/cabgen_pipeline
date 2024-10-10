import os
import logging
from typing import Union


def logging_conf(task_id: int, log_dir: str,
                 log_level: int = logging.INFO) -> Union[logging.Logger, None]:
    """
    Configures logging for a specific task.

    Args:
        task_id (int): ID of the task for which logging is configured.
        log_dir (str): Directory where log files should be stored.
        log_level (int): Logging level (default is INFO).
    """
    try:
        if not os.path.exists(log_dir):
            os.mkdir(log_dir)

        logger = logging.getLogger(f"tarefas_{task_id}")
        logger.setLevel(log_level)

        log_file = os.path.join(log_dir, f"tarefas_{task_id}.log")

        handler = logging.FileHandler(log_file, mode='w')

        formatter = logging.Formatter(
            '%(asctime)s - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        handler.setFormatter(formatter)

        # Add the handler to the logger (check if handler exists to avoid
        # duplication)
        if not logger.hasHandlers():
            logger.addHandler(handler)

        return logger

    except Exception as e:
        print(f"Error configuring logging for task {task_id}: {e}.")
        return None
