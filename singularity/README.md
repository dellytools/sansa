You can build a [sansa](https://github.com/dellytools/sansa) singularity container (SIF file) using

`sudo singularity build sansa.sif sansa.def`

Once you have built the container you can run analysis using

`singularity exec sansa.sif sansa --help`
