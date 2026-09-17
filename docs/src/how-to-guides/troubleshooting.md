# Troubleshooting
This page may be useful if you are having problems.

### Julia version is incompatible / Errors about OpenSSL library
Julia 1.12 shipped an OpenSSL library that was incompatible with Python, which previously required pinning Obliqua to Julia 1.11. Obliqua now targets Julia 1.13 by default; whether the same OpenSSL/Python conflict resurfaces under 1.13 has not been separately verified here, since it depends on how your Python environment is set up. If you see OpenSSL-related errors when running Obliqua alongside Python on 1.13, fall back to Julia 1.11 using the `juliaup` command:
```bash
juliaup status        # Show the versions you have installed
juliaup add 1.11      # Install Julia 1.11
juliaup default 1.11  # Make 1.11 your default version of Julia
julia --version       # This should say 1.11.*
```

After following these steps, try installing Obliqua again.

### Julia errors on start, potentially referencing the CURL library
It is important that the shell environment variable `LD_LIBRARY_PATH` is
not set when running Obliqua. This will cause Julia to use the wrong libraries, which will causes problems. You can unset this variable or reset using either of the following commands
```bash
unset LD_LIBRARY_PATH
export LD_LIBRARY_PATH=""
```
If this does not help, it's possible that you are using a Julia distribution provided by your system package manager. It's important that you only use Julia distributed from the official website.
