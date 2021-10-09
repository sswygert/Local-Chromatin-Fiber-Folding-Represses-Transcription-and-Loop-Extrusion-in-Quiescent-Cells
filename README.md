# Local Chromatin Fiber Folding Represses Transcription and Loop Extrusion in Quiescent Cells

## Micro-C data processing and analysis

The tools include a set of C++ source code, python scripts and nextflow scripts. On a high level, the included nextflow scripts are used to drive the C++ and python programs to process and analyse the micro-C data. The user is expected to put the micro-c data in a directory structure as "/path/to/parentdir/micro-c/MySampleName/*.fastq". To begin, clonethis repo to a directory, e.g., "/path/to/repo", and make a subdirectory "/path/to/repo/build" and cd into that directory and run "cmake .. && make" to compile the C++ programs. Then add "/path/to/repo/build" to the environment variable PATH, e.g., "export PATH=/path/to/repo/build:$PATH"

### Example: Process Micro-C fastq files

The scripts processSample.nf calls various nextflow workflows in utils.nf and workflows.nf to process micro-C fastq files into *.hic and *.mcool files. As mentioned above, if the user has put the micro-C data in "/path/to/parentdir/micro-c/MySampleName/*.fastq", they can run "processSample.nf --input /path/to/parentdir/micro-c/MySampleName/*.fastq" to get the output files in "/path/to/parentdir/micro-c/MySampleName/".
