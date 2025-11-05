# pileup-events

Count alignment events per position for a specified genomic region.

## Installation

You will need cmake and htslib installed,
and preferably pkg-config if you want to install
the easy way.
If you don't already have these dependencies,
it would be best to install via your package
manager. e.g. on Ubuntu
```bash
  apt install cmake cxxopts htslib pkg-config
```
or on mac
```bash
  brew install cmake cxxopts htslib pkg-config
```

Once you have these dependencies, to
create the pileup-events binary:
```bash
  mkdir build
  cd build
  cmake ..
  cmake --build .  # create binary
```

You can invoke the binary like `./path/to/pileup-events --help`.
Add the binary path to PATH to be able to do `pileup-events` from anywhere.

## Usage

See `./pileup-events --help`.
Then pick a location to investigate,
for example a variant from a VCF or a region from a BED file,
and run on that location to see the event counts.

e.g.:
```bash
  pileup-events ~/path/to/sample.bam chr1:1000  # single location
  pileup-events ~/path/to/sample.bam chr1:1000-2000  # 1001bp range
```

The region string is 1-indexed, end-inclusive, i.e. identical to `samtools view`.
Both .bam and .cram are in principle supported,
though no testing on cram has been done as yet.

# Acknowledgements

This tool uses [htslib](https://github.com/samtools/htslib) by the samtools team, and [cxxopts](https://github.com/jarro2783/cxxopts) by Jarryd Beck & contributors.

<!-- NOTE: I will add that is an upgrade over the original having proved it actually works under load! -->
The function of the tool is to provide an standalone version of the `bam2R()` functionality found in the R package [deepSNV](https://github.com/gerstung-lab/deepSNV). It is a complete rewrite of the concepts found therein. The citations for `deepSNV` are as follows:
> Gerstung M, Beisel C, Rechsteiner M, Wild P, Schraml P, Moch H, Beerenwinkel N (2012). “Reliable detection of subclonal single-nucleotide variants in tumor cell populations.” Nat Commun, 3, 811.
>  Gerstung M, Papaemmanuil E, Campbell PJ (2014). “Subclonal variant calling with multiple samples and prior knowledge.” Bioinformatics, 30, 1198-1204.

