TelSeq is a software that estimates telomere length from
whole genome sequencing data (BAMs).

The most current development version is available from our
git repository:
[git://github.com/zd1/telseq.git](git://github.com/zd1/telseq.git)

The software is implemented in C++.

Citation:

_Estimating telomere length from whole genome sequence data_

Zhihao Ding; Massimo Mangino; Abraham Aviv; Tim Spector; Richard Durbin.
Nucleic Acids Research 2014; doi: 10.1093/nar/gku181
[http://nar.oxfordjournals.org/content/42/9/e75](http://nar.oxfordjournals.org/content/42/9/e75)


## Compile TelSeq

### TelSeq dependency:
- the bamtools library (https://github.com/pezmaster31/bamtools)

- A modern version of GCC (version 4.8 or above)
This can been seen by "gcc --version".
If multiple GCCs are installed in your system, please set environmental
variables pointing to the one of version 4.8 or above. e.g. in bash,

```
export CXX=/path/to/gcc/gcc-4.8.1/bin/g++
export CC=/path/to/gcc/gcc-4.8.1/bin/gcc
```

One easy way to install a new GCC is to use homebrew,

```
# install homebrew if you don't have it
ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

# install GCC
brew install gcc
```


### Compile
Go to the src directory and run autogen.sh from the src directory to generate the configure file
`./autogen.sh`
Then run
```
./configure
make
```

The executable binary will be at src/Telseq/telseq.
If bamtools are installed not at the system location, you can
specify their location by

`./configure  --with-bamtools=/path/to/bamtools`

The /path/to/bamtools directory is the directory that contains 'lib' and 'include' sub directories.

## Run TelSeq

##### Read options and usage information
`telseq`

##### Analyse one or more BAMs by specifying BAM file path as command line arguments.
`telseq a.bam b.bam`

##### Analyse a list of BAMs whose paths are specified in a 'bamlist' file.
bamlist should contain only 1 column with each row the path of a BAM. i.e.

```
/path/to/a.bam
/path/to/b.bam
```
`telseq -f bamlist`

##### BAM file path can also be provided by piping in a 'bamlist', whose format must be same as above
`cat bamlist | telseq`


##### output
By default the result will be printed out to stdout. To change it to a file, use '-o'
option to specify a file path. i.e.

`telseq -o /path/to/output a.bam b.bam c.bam`

This can also be achived by just direct the output to a file using '>', i.e.

`telseq a.bam b.bam c.bam > /path/to/output`

The software will print out running status to stderr as well. To separate them from stdout, one
could direct log to a file, ie.

`telseq a.bam b.bam c.bam 2>outputlog`

Merge results from read groups by taking a weighted mean. However, it is benetifical run without
-m to output the result per lane, so to have an idea about inter-lane variation. The merging
can be done afterwards.
`telseq -m a.bam > output`

#### Output file format

|  Column | Definitions |
| -------------|----------------------------------------------|
| ReadGroup | read group, Defined by the RG tag in BAM header. |
| Library   | sequencing library that the read group belongs to.|
| Sample    | defined by the SM tag in BAM header. |
| Total     | total number of reads in this read group. |
| Mapped    | total number of mapped reads, SAM flag 0x4. |
| Duplicates | total number of duplicate reads, SAM flag 0x400. |
| LENGH_ESTIMATE | estimated telomere length. |
| TEL0 | read counts for reads containing no TTAGGG/CCCTAA repeats. |
| TEL1 | read counts for reads containing only 1 TTAGGG/CCCTAA repeats. |
| TELn | read counts for reads containing only n TTAGGG/CCCTAA repeats. |
| TEL16 | read counts for reads containing 16 TTAGGG/CCCTAA repeats. |
| GC0 | read counts for reads with GC between 40%-42%. |
| GC1 | read counts for reads with GC between 42%-44%. |
| GCn | read counts for reads with GC between (40%+n*2%)-(42%+(n+1)*2%). |
| GC9 | read counts for reads with GC between 58%-60%.  |

By default for each BAM a header line will be printed out. This can be suppressed by using the '-H' option. It is useful when one has multiple BAMs to scan and wish the output to be merged together. i.e.

`telseq -H a.bam b.bam c.bam > myresult`

To just print out the header, use '-h' option. i.e.

`telseq -h`

## Docker

### Install Docker

Please refer to the official website for installing Docker
[https://docs.docker.com/engine/installation/](https://docs.docker.com/engine/installation/)


### Build telseq Docker image

```
docker build -t telseq-docker github.com/zd1/telseq
```

### Run telseq

```
docker run telseq-docker
```

## Contact

zhihao.ding at gmail.com
