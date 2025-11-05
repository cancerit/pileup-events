# pileup-events

## Installation

You will need cxxopts and htslib installed,
and preferably pkgconfig.
Any standard install method should work for
these dependecies, e.g. `brew`, `apt`.  

To create the pileup-events binary:
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

