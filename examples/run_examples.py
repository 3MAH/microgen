"""Test examples"""

import argparse
import os
import shutil
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Tuple


def get_example_paths(exclude_dirs: List[str]) -> List[str]:
    """Get all the example paths in the examples directory.
    Exclude the directories in exclude_dirs."""
    paths: List[str] = []
    for root, _, files in os.walk(Path(__file__).parent):
        if root not in exclude_dirs:
            paths.extend(
                [str(Path(root) / file) for file in files if file.endswith(".py")]
            )

    paths.remove(str(Path(__file__)))
    return paths


def dashed_line(string: Optional[str] = "", char: str = "-") -> str:
    """Create a dashed line with a string in the middle."""
    terminal_width = shutil.get_terminal_size().columns
    title = f" {string} " if string else ""
    left_padding = (terminal_width - len(title)) // 2
    right_padding = terminal_width - left_padding - len(title)

    return char * left_padding + title + char * right_padding


def launch_test_examples(
    exclude_dirs: List[str],
    nprocs: int = 1,
) -> Tuple[List[str], Dict[str, str]]:
    """Launch all the examples and return the paths of the examples that failed."""
    examples = get_example_paths(exclude_dirs=exclude_dirs)
    failed: Dict[str, str] = {}
    if nprocs == 1:
        for i, example in enumerate(examples):
            print(dashed_line(f"[{i}/{len(examples)}] {example}"))
            cmd = ["python", example]
            try:
                subprocess.check_output(cmd, stderr=subprocess.STDOUT)
            except subprocess.CalledProcessError as e:
                failed[example] = e.output.decode("utf-8")
                print(e.output.decode("utf-8"))

            print(dashed_line())
            print()
    elif nprocs > 1:
        # import multiprocessing

        raise NotImplementedError("nprocs > 1 not implemented")
    else:
        raise ValueError("nprocs must be greater than 0")
    return examples, failed


GREEN = "\033[32m"
RED = "\033[31m"
RESET = "\033[0m"
BOLD = "\033[1m"
UNDERLINE = "\033[4m"


def display_results(examples: List[str], failed: Dict[str, str]):
    """Display the results of the test examples.
    Raise an error if there are failed examples."""
    print(dashed_line(f"{UNDERLINE}{BOLD}RESULTS{RESET}"))
    print(
        f"{UNDERLINE}{BOLD}{GREEN}Successful: {len(examples) - len(failed)}/ {len(examples)}{RESET}"
    )
    if len(failed) > 0:
        error_str = f"{UNDERLINE}{BOLD}{RED}Failed: {len(failed)}{RESET}\n"
        for example, error in failed.items():
            error_str += f"{RED}X\t{example}{RESET}\n"
            error_str += error + "\n"
        raise RuntimeError(error_str)


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
    type=int,
    help="Number of processes to use",
    default=1,
)

PARENT_DIR = Path(__file__).parent

EXCLUDE = parser.parse_args().exclude
if parser.parse_args().no_mmg:
    EXCLUDE.extend(["Mesh/mmg", "Mesh/mmg-voro"])
EXCLUDE.append("Fibers")
EXCLUDE = [str(PARENT_DIR / exclude) for exclude in EXCLUDE]

display_results(
    *launch_test_examples(
        exclude_dirs=EXCLUDE,
        nprocs=parser.parse_args().nprocs,
    )
)
