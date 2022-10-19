import subprocess

func_names = [
    "LoopRowMajorOrdering",
    "LoopRowMajorBlocking",
    "LoopRowMajorPacking",
    "RecursionRowMajorOrdering",
    "RecursionRowMajorBlocking",
    "RecursionRowMajorPacking",
    "RecursionZmortonOrdering",
    "RecursionZmortonPacking"
]

m, n, k = 256, 256, 256
block_size = 64
file_name = "./result.txt"
toolset = "g++"

subprocess.call("make clean", shell=True)
subprocess.call(f"make test toolset={toolset}", shell=True)
subprocess.call(f"rm -f {file_name}", shell=True)

for func_no in range(len(func_names)):
    subprocess.call(f"./test {m} {n} {k} {func_no} {block_size} >> {file_name}", shell=True)
