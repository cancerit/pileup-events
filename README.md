# pileup-events

Count alleles and alignment events per position for a specified genomic region.

## Installation

You will need `cmake>=3.16` and `htslib>=1.14` installed,
and preferably `pkg-config` if you want to install
the easy way. If you want to create python and/or R
bindings you will also need `SWIG>=4.0`.
If you don't already have these dependencies,
it would be best to install via your package
manager. e.g. on Ubuntu
```bash
  apt install cmake swig pkg-config libhts-dev
```
or on mac
```bash
  brew install cmake swig pkg-config htslib
```

Once you have these dependencies, to
create the pileup-events binary and/or R and python bindings:
```bash
  # git clone <this repo>
  # cd <repo folder>
  mkdir build
  cd build
  cmake .. \
    -DMAKE_EXE=ON \ # create binary
    -DMAKE_PY_BINDS=ON \ # make python bindings
    -DMAKE_R_BINDS=ON # R
  cmake --build .
```

You can then invoke the binary like `./path/to/pileup-events --help`.
Add the binary path to PATH to be able to call `pileup-events` from anywhere.
For usage of the R and python bindings see subsequent sections.

## Usage

See `pileup-events --help`.
Then pick a location to investigate,
for example a variant from a VCF or a region from a BED file,
and run on that location to see the event counts.

e.g.:
```bash
  pileup-events \
    -i 2 \  # only count reads in proper pairs
    ~/path/to/sample.bam \
    chr1:1000  # single location

  pileup-events \
    ~/path/to/sample.bam \
    chr1:1000-2000  # 1001bp range
```

The region string is 1-indexed, end-inclusive, i.e. identical to `samtools view` -
excepting the fact that `pileup-events` allows a series of shorthands such as `<chr>:<pos>` 
for a single location. See the helptext for more details.
Assuming compilation against a recent version of htslib,
both .bam and .cram are in principle supported.
No testing on cram has been done as of yet.

## R & Python Bindings

Bindings for usage of pileup events in both R and Python can be generated during compilation. This is enabled by [SWIG](https://www.swig.org/).
As described in the installation section above, python bindings may be created with the `-DMAKE_PY_BINDS=ON` option and R with `-DMAKE_R_BINDS=ON`. Doing so will result in the generation of `build/python` and `build/r` respectively.

At present bindings are somewhat rudimentary - this may be improved upon in the future. Usage is as follows:

**R:**
```R
  build_dir <- "</absolute/path/to/build>/r"  # replace between <>
  dyn.load(paste0(build_dir, "/pileup_events.so"))
  source(paste0(build_dir, "/pileup_events.R"))
  count_events(
    "absolute/path/to/bam",
    "chrX:150"
    # no_overlaps=False,  # optional parameters with default values
    # min_mapq=25,
    # min_baseq=30,
    # include_flag=0,
    # exclude_flag=3844,
    # max_depth=1000000,
    # clip_bound=0
  )
```

**Python:**
```bash
  # setup in BASH
  pev_bind_path=</path/to/build>/python   # replace between <>
  export PYTHONPATH=${PYTHONPATH}:${pev_bind_path}  # add to python path (not permanent)
```
```python
  import pileup_events as pev
  help(pev.count_events)
  pev.count_events(
    "absolute/path/to/bam",
    "chrX:100-200"
    # no_overlaps=False,  # optional parameters with default values
    # min_mapq=25,
    # min_baseq=30,
    # include_flag=0,
    # exclude_flag=3844,
    # max_depth=1000000,
    # clip_bound=0
  )
```

In either case the return value of the count events function is a 1D vector, where each of the genomic positions counted is a block of 24 cells in the vector. It is currently left to the user to transform this strucutre into any desirable alternative.

The compiled bindings directories (`python/` and `r/`) can be renamed, and moved anywhere appropriate on the system. They are not dependent on other build artefacts. Do not modify or rename any of the files within these directories. Note that for the python bindings if you do move/rename the `python/` directory you will need to add the new location to PYTHONPATH.

## Development & Testing

To enable compliation with debug flags and compilation of the test binary, add the following options to the `cmake ..` step:

```bash
  cmake .. \
    -DENABLE_DEBUG_FLAGS=ON \
    -DMAKE_TEST=ON # \
    # other options
```

Compliation of the test binary will produce an additional artefact, `build/test-pev`. Execution of this artefact will run the test suite. The test suite is currently quite brief, and may be expanded upon in the future.

## Authors & Acknowledgements

`pileup-events` is the work of Alex Byrne (alex@blex.bio) & Luca Barbon of CASM Informatics, Wellcome Sanger Institute.

This tool uses [htslib](https://github.com/samtools/htslib) for accessing sequence data, and [cxxopts](https://github.com/jarro2783/cxxopts) for the CLI interface. Testing is made using [Catch2](https://github.com/catchorg/Catch2)

<!-- NOTE: I will add that is an upgrade over the original having proved it actually works under load! -->
The function of the tool is to provide an standalone version of the `bam2R()` functionality found in the R package [deepSNV](https://github.com/gerstung-lab/deepSNV). It is a complete rewrite of the concepts found therein.

### Citations
htslib:
> HTSlib: C library for reading/writing high-throughput sequencing data
> James K Bonfield, John Marshall, Petr Danecek, Heng Li, Valeriu Ohan, Andrew Whitwham, Thomas Keane, Robert M Davies
> GigaScience, Volume 10, Issue 2, February 2021, giab007, https://doi.org/10.1093/gigascience/giab007  

deepSNV:  
> Gerstung M, Beisel C, Rechsteiner M, Wild P, Schraml P, Moch H, Beerenwinkel N (2012). “Reliable detection of subclonal single-nucleotide variants in tumor cell populations.” Nat Commun, 3, 811.  
>  Gerstung M, Papaemmanuil E, Campbell PJ (2014). “Subclonal variant calling with multiple samples and prior knowledge.” Bioinformatics, 30, 1198-1204.

