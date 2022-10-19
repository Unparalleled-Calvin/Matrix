import subprocess
import pandas as pd

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

def run(m = 256, n = 256, k = 256, block_size = 64, file_name = "./result.txt", toolset = "g++"):
    subprocess.call("make clean", shell=True)
    subprocess.call(f"make test toolset={toolset}", shell=True)
    subprocess.call(f"rm -f {file_name}", shell=True)

    for func_no in [0, 1, 2, 4, 5, 7]: # remove naive func
        subprocess.call(f"./test {m} {n} {k} {func_no} {block_size} >> {file_name}", shell=True)

def parse(file_name = "result.txt"):
    funcs = []
    times = []
    with open(file_name, "r") as f1:
        lines = [line.strip() for line in f1.readlines()]
        for i in range(0, len(lines), 5):
            func_name = lines[i + 1].split(": ")[1]
            time = lines[i + 3].split(": ")[1].split()[0]
            times.append(time)
            funcs.append(func_name)
    return times, funcs

def test(lb = 8, ub = 8, file_name = "result.txt", save = True, save_type = "table",  toolset = "g++"):
    print(f"{1 << lb} <= matrix size <= {1 << ub}")
    if save:
        print(f"{save_type} type results will be saved in {file_name}")
    print()
    result = {}
    names = []
    for size_pow in range(lb, ub + 1):
        for block_pow in range(min(5, lb, max(size_pow - 2, 0)), max(size_pow - 2, 0) + 1):
            m = n = k = 1 << size_pow
            block_size = 1 << block_pow
            print(f"start running a test with m={m} n={n} k={k} block_size={block_size}.")
            run(m = m, n = n, k = k, block_size = block_size, file_name = file_name, toolset = toolset)
            print("done.\n")
            result[f"{m}-{block_size}"], names = parse(file_name = file_name)
    result["func"] = names
    result = pd.DataFrame(data = result).set_index("func")
    if save:
        if save_type == "csv":
            result.to_csv(file_name)
        elif save_type == "table":
            result.to_string(file_name)
    print("All tests are done.")
    return result

if __name__ == "__main__":
    test(8, 8)