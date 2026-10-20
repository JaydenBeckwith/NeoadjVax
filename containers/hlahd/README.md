# HLA-HD container

HLA_TYPING runs HLA-HD from an image you build yourself, because HLA-HD
is distributed under an academic licence after registration.

1. Register and download HLA-HD from
   <https://w3.genome.med.kyoto-u.ac.jp/HLA-HD/> (the tarball is named
   like `hlahd.1.7.1.tar.gz`). Put it in this directory. Do not commit it.
2. Build on a machine with Docker (not Gadi):

   ```bash
   docker build --build-arg HLAHD_TARBALL=hlahd.1.7.1.tar.gz \
       -t neoadjvax/hlahd:1.7.1 containers/hlahd
   ```

   The build fails if `install.sh` did not produce `bin/hlahd.sh`, or if
   `dictionary/`, `freq_data/` or `HLA_gene.split.txt` are missing from
   the tarball's layout.
3. Convert to Singularity and copy to Gadi:

   ```bash
   singularity build hlahd_1.7.1.sif docker-daemon://neoadjvax/hlahd:1.7.1
   scp hlahd_1.7.1.sif gadi.nci.org.au:/scratch/jo11/singularity_images/
   ```

4. Point the pipeline at it, e.g. in the params file:

   ```yaml
   hlahd_container: /scratch/jo11/singularity_images/hlahd_1.7.1.sif
   ```

The pipeline calls `hlahd.sh` with `freq_data/`, `HLA_gene.split.txt`
and `dictionary/` under `/opt/hlahd` (`--hlahd_install_dir`). To update
the IMGT dictionary, run HLA-HD's `update.dictionary.sh` during the build
(it needs network access) before converting.

Before the first real run, check the image on a Gadi login node:

```bash
singularity exec hlahd_1.7.1.sif bash -c 'hlahd.sh 2>&1 | head; ls /opt/hlahd'
```
