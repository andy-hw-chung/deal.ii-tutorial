# Example 1 -  Building projects

1. Make sure the environment variable `DEAL_II_DIR` is set in the terminal environment (Hint: export the variable automatically each time you start the terminal by editing your `.bashrc` file. Add the line
```bash
export DEAL_II_DIR=path
```
where `path` is the path to your deal.ii installation

2. Make sure that there is a `CMakeLists.txt` present and that it is configured properly. Do
```
cmake .
```
to build the Make files

3. To build the application, do either
```
make debug
```
or
```
make release
```
Use `make debug` while writing your application to check more thouroughly for errors during compilation and runtime. Use `make release` when the code has been completely debugged to make your application run more quickly. Alternatively, do
```
make
```
to build with the previous used build type. This default build type is `debug`.

4. Do
```
make clean
```
to cleanup uneeded files between builds.

5. If switching between different deal.ii installations (e.g. using the same code on a different computer) make sure to delete the `CMakeCache.txt` file and re-run `cmake .` to ensure the Make files point to the correct `DEAL_II_DIR` path.
