### Usage

#### Test by Python

Dependencies:

- Linux (for compilation)
- Python 3.6+
- Pandas

Specify lower bound and upper bound of the matrix size (logarithm to the base 2, 8 and 8 by default) in `test.py`.

You can get csv type file (by default) or table like text file in `result.txt` (by deafult).

```shell
# for foreground process 
python3 test.py
# for background process
nohup python3 test.py 1>log.txt 2>&1 </dev/null &
```

### Compilation

`test.py` calls a subprocess to compile source files.

You can compile and test the code manually by editting `makefile`.

The recommend toolset is `g++` (by default) or `clang++` and the `-O2` flag helps a lot.

```shell
make clean
# for g++
make test toolset=g++
# for clang++
make test toolset=clang++
```

