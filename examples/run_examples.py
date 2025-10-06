"""Test examples"""

from __future__ import annotations

import argparse
import multiprocessing
import os
import shutil
import subprocess
import sys
from pathlib import Path

GREEN = "\033[32m"
RED = "\033[31m"
RESET = "\033[0m"
BOLD = "\033[1m"
UNDERLINE = "\033[4m"


def get_example_paths(exclude_dirs: list[str]) -> list[str]:
    """Get all the example paths in the examples directory.
    Exclude the directories in exclude_dirs."""
    paths: list[str] = []
    for root, _, files in os.walk(Path(__file__).parent):
        if root not in exclude_dirs:
            paths.extend(
                [str(Path(root) / file) for file in files if file.endswith(".py")]
            )

    paths.remove(str(Path(__file__)))
    return paths


def dashed_line(string: str | None = "", char: str = "-") -> str:
    """Create a dashed line with a string in the middle."""
    terminal_width = shutil.get_terminal_size().columns
    title = f" {string} " if string else ""
    left_padding = (terminal_width - len(title)) // 2
    right_padding = terminal_width - left_padding - len(title)

    return char * left_padding + title + char * right_padding


def launch_test_examples(
    exclude_dirs: list[str],
    n_procs: int = 1,
) -> tuple[list[str], list[str]]:
    """Launch all the examples and return the paths of the examples that failed."""
    examples = get_example_paths(exclude_dirs=exclude_dirs)
    failed: list[str] = []
    if n_procs == 1:
        for i, example in enumerate(examples):
            print(dashed_line(f"[{i}/{len(examples)}] {example}"))
            cmd = [sys.executable, example]
            process = subprocess.run(cmd, check=True)
            if process.returncode != 0:
                failed.append(example)

            print(dashed_line())
            print()
    elif n_procs > 1:
        results: dict[str, multiprocessing.pool.AsyncResult] = {}
        with multiprocessing.Pool(processes=n_procs) as pool:
            for example in examples:
                cmd = [sys.executable, example]
                result = pool.apply_async(subprocess.run, (cmd,), {"check": True})
                results[example] = result

            pool.close()
            pool.join()

        for example, result in results.items():
            if result.successful():
                print(f"{GREEN}.{RESET}", end="")
            else:
                failed.append(example)
                print(f"{RED}X{RESET}", end="")
        print()
    else:
        raise ValueError("nprocs must be greater than 0")
    return examples, failed


def display_results(examples: list[str], failed: list[str]):
    """Display the results of the test examples.
    Raise an error if there are failed examples."""
    print(dashed_line(f"{UNDERLINE}{BOLD}RESULTS{RESET}"))
    print(
        f"{UNDERLINE}{BOLD}{GREEN}\
            Successful: {len(examples) - len(failed)}/ {len(examples)}\
                {RESET}"
    )
    if failed:
        error_str = f"{UNDERLINE}{BOLD}{RED}Failed: {len(failed)}{RESET}\n"
        for example in failed:
            error_str += f"{RED}X\t{example}{RESET}\n"
        print(error_str)
        raise RuntimeError(f"{len(failed)} examples failed")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run all examples")
    parser.add_argument(
        "--exclude",
        nargs="*",
        help="Exclude directories from the examples",
        default=[],
    )
    parser.add_argument(
        "--no-mmg",
        action="store_true",
        help="Exclude MMG examples",
    )
    parser.add_argument(
        "-n",
        "--nprocs",
        help="Number of processes to use, \
            'auto' to use all available cores. Default is 1.",
        default="1",
    )
    str_nprocs = parser.parse_args().nprocs
    if str_nprocs == "auto":
        nprocs = multiprocessing.cpu_count()
    elif str_nprocs.isdecimal():
        number = int(str_nprocs)
        if number > 0:
            nprocs = number
    else:
        raise ValueError("nprocs must be a positive integer")

    print(f"Running {nprocs} processes")

    PARENT_DIR = Path(__file__).parent

    EXCLUDE = parser.parse_args().exclude
    if parser.parse_args().no_mmg:
        EXCLUDE.extend(
            [
                "Mesh/mmg",
                "Mesh/mmg-voro",
                "Mesh/remesh",
                "Mesh/gyroid",
            ]
        )
    EXCLUDE.append("Fibers")
    EXCLUDE = [str(PARENT_DIR / exclude) for exclude in EXCLUDE]

    display_results(
        *launch_test_examples(
            exclude_dirs=EXCLUDE,
            n_procs=nprocs,
        )
    )
