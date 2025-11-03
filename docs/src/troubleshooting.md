# Troubleshooting
This page may be useful if you are having problems. However, I would suggest that you also double check that you followed all of the [Getting started](@ref) instructions.

## Julia version is incompatible / Errors about OpenSSL library
You must use Julia version 1.11 because there are incompatibilities between the OpenSSL library required by Julia1.12 and Python. 

Switch Julia versions using the `juliaup` command. E.g:
```bash
juliaup status        # Show the versions you have installed
juliaup add 1.11      # Install Julia 1.11
juliaup default 1.11  # Make 1.11 your default version of Julia
julia --version       # This should say 1.11.something
```

After following these steps, try installing AGNI again.

## Wget is not installed
You need to install [wget](https://www.gnu.org/software/wget/). This is a tool for tranferring files over a network. Wget is used by AGNI to obtain lookup data files from Zenodo. Most Linux distributions come with wget; otherwise see [this page](https://www.tecmint.com/install-wget-in-linux/).

To install wget on MacOS:
```bash
brew install wget
```

## Unzip is not installed
You need to install [unzip](https://www.gnu.org/software/wget/). This command is used by AGNI to extract some data files once downloaded. Most computers come with this command; otherwise see [this page](https://ioflood.com/blog/install-unzip-command-linux/).

To install unzip on Ubuntu/Debian:
```bash
sudo apt-get install unzip
```

## Julia errors on start, potentially referencing the CURL library
It is important that the shell environment variable `LD_LIBRARY_PATH` is
not set when running AGNI. This will cause Julia to use the wrong libraries,
which will causes problems. You can unset this variable or reset using either of the
following commands
```bash
unset LD_LIBRARY_PATH
export LD_LIBRARY_PATH=""
```
If this does not help, it's possible that you are using a Julia distribution provided by
your system package manager. It's important that you only use Julia distributed from the
official website.


