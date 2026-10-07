"""
Build the two NuFast libraries used by the benchmark from nufast_bridge.cpp:
    harmonised - Hamiltonian parameter set and conventions
    native - original NuFast parameter set and conventions
"""

from pathlib import Path
import shutil
import subprocess
import sys


HERE = Path(__file__).resolve().parent
NUFAST = HERE.parent / "NuFast-LBL-main" / "Benchmarks"
SRC = NUFAST / "src" / "NuFast_LBL.cpp"
INCLUDE = NUFAST / "include"
BRIDGE = HERE / "nufast_bridge.cpp"


GPP = shutil.which("g++")
if GPP is None:
    raise RuntimeError("g++ was not found in PATH. For MSYS2 UCRT64 on Windows, add this directory to PATH: C:\\msys64\\ucrt64\\bin")


if sys.platform.startswith("win"):
    PREFIX = ""
    SUFFIX = ".dll"
    PLATFORM_FLAGS = ["-shared"]

elif sys.platform == "darwin":
    PREFIX = "lib"
    SUFFIX = ".dylib"
    PLATFORM_FLAGS = ["-dynamiclib"]

else:
    PREFIX = "lib"
    SUFFIX = ".so"
    PLATFORM_FLAGS = ["-shared", "-fPIC"]


def build(mode):
    output = HERE / f"{PREFIX}nufast_{mode}{SUFFIX}"
    command = [GPP, "-O3", "-std=gnu++17", *PLATFORM_FLAGS, "-I", str(INCLUDE)]
    if mode == "harmonised":
        command.append("-DNUFAST_HARMONISED=1")
    command += [str(BRIDGE), str(SRC), "-o", str(output)]
    print(f"Building {output.name} ...")
    subprocess.run(command, check=True)
    print(f"Created {output.name}")


for mode in ("harmonised", "native"):
    build(mode)
