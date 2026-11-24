# Installing cdid

These pages describe `cdid` 1.4.0.0 for Stata.

## Requirements

`cdid` requires Stata 19.5 (StataNow) or later. It has been tested with StataNow 19.5 on macOS (Apple silicon). The estimation commands are written in Stata and Mata; no other software is needed.

Chained mode, the default, has no other dependency. Ordinary mode (`mode(ordinary)`, or `rcs`) uses the `csdid` package, version 2.0.2, by Brantly Callaway, Fernando Rios-Avila and Pedro H. C. Sant'Anna. It is distributed with `cdid` under its MIT license but is installed separately, so that an existing installation of `csdid` is never replaced without your consent.

## Install from GitHub

```stata
net install cdid, from("https://raw.githubusercontent.com/dvbn/cdid-stata/main")
mata: mata mlib index
```

To use ordinary mode, also install `csdid` 2.0.2:

```stata
net install csdid, from("https://raw.githubusercontent.com/dvbn/cdid-stata/main/csdid")
mata: mata mlib index
csdid version
```

If `csdid` is already installed, check its version with `csdid version` and its location with `which csdid, all` before installing. Ordinary mode stops with an error if the version found is not 2.0.2. To keep another version for other projects, install `cdid` and `csdid` 2.0.2 into a project folder as described below.

## Install into a project folder

The following commands, typed in a new Stata session started in the project folder, install both packages into subfolders of the project. They change the search path of the current session only; no profile file is modified. See `help sysdir` and `help net`.

```stata
mkdir "./ado_plus"
mkdir "./ado_personal"
sysdir set PLUS "./ado_plus"
sysdir set PERSONAL "./ado_personal"
net install cdid, from("https://raw.githubusercontent.com/dvbn/cdid-stata/main")
net install csdid, from("https://raw.githubusercontent.com/dvbn/cdid-stata/main/csdid")
mata: mata mlib index
```

Repeat the `sysdir set` commands at the start of each session that uses this installation. To install from a local copy of the repository, replace the URLs with the path of that copy, for example `net install cdid, from("/path/to/cdid-stata")`.

## Example files

The example do-files and the documentation are ancillary files. To copy them into the current working folder, type

```stata
net get cdid, from("https://raw.githubusercontent.com/dvbn/cdid-stata/main")
```

This copies `cdid_hands_on.do` and `cdid_custom_estimators.do`, among other files; see [Hands-on guide](hands-on.md) and [Custom estimators](custom-estimators.md). `cdid_hands_on.do` begins with a block of settings that you can edit before running it; it clears the data in memory only if you allow it to, and it keeps output files when you give it an output folder.

## Check the installation

```stata
which cdid
which cdid_pretest
which cdid_diagnostics
which cdid_plot
which cdid_export
which cdid_load
which cdid_pretest_load
which cdid_nuisance
which cdid_chain
findfile lcdid_native.mlib
```

Each command should be found in the folder where you installed `cdid`.

## Building the Mata library

The installed library `lcdid_native.mlib` is compiled from the Mata source files in the repository's `mata/` folder. To rebuild it from a copy of the repository, start Stata in the repository folder and type

```stata
mkdir "./build_lib"
do tools/build_native.do "." "." "./build_lib"
```

The new library is written to `build_lib/lcdid_native.mlib`. Building is never required to use `cdid`.
