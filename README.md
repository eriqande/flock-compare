# flock-compare

This is an archive of scripts and data to reproduce the simulations we did for the paper
"Interpreting the FLOCK algorithm from a statistical perspective" to be published 
in Molecular Ecology Resources.


## Reproducing the simulations

### Flockture versus structure
The R script 'Flock2Structure' will rerun the analyses comparing Structure and FLOCKTURE from the 
manuscript.  In order to make this work you ought to be on a sane Unix system (you will need `rsync`,
`awk`, `make` and other such utilities.  Easiest, probably, is going to
be if you are running Mac OS Mavericks or higher with the Developer Tools installed.

The code is divided into 5 main parts: 

1. Simulate data
2. Run FLOCKTURE
3. Run STRUCTURE
4. Run CLUMPP and DISTRUCT
5. Compute zero-one loss and graph output.  

In order to do this you will need to clone a few git repositories: `flock-compare`, `flockture` and 
`slg_pipe`.  (Note that the current state of `flockture` will be archived with Dryad as well).
With the latter two there are some further setup things that must be done as described
in their READMEs.  See the README section at 
[https://github.com/eriqande/flockture](https://github.com/eriqande/flockture) and
[https://github.com/eriqande/slg_pipe](https://github.com/eriqande/slg_pipe)

When you do this, make sure that you clone all the repositories so that the directories
that get cloned all appear together at the same directory level.  When I do this it looks like
so:

```sh
# first get the three repositories
git clone https://github.com/eriqande/flock-compare.git
git clone https://github.com/eriqande/flockture.git
git clone https://github.com/eriqande/slg_pipe.git

# notice how these are all together in the same directory level:
ls

## output is: flock-compare/ flockture/     slg_pipe/


# Now make flockture:
cd flockture/src/
make
cd ../../   # return to the directory that holds the cloned repos


# Then get the necessary binaries for slg_pipe and place them where needed
curl -o slg_pipe_binaries.tar.gz  https://dl.dropboxusercontent.com/u/19274778/slg_pipe_binaries.tar.gz
gunzip slg_pipe_binaries.tar.gz 
tar -xvf slg_pipe_binaries.tar 
rsync -avh slg_pipe_binaries/* slg_pipe


```
Please note that the executables in the slg_pipe are for Mavericks OS and should be changed if you have a different OS. Follow the instructions in each of the readme files on how to:

1. compile the FLOCKTURE code, and
2. set up slg_pipe.

Once FLOCKTURE and the slg_pipe are setup, open the FLOCKTUREvSTRUCTURE.R script in your favorite editor such as Rstudio, Tinn-R, Vim (with VimR), etc. The code should be executed for each of the 5 main parts separately. There are a few options that can be changed:

1. Marker type (line10): Set to 1 for microsatellites and 2 for SNPs
2. Reps (line33): How many times do you want to run all models. Can be set from 1-9.
3. MGrate (lines37/40): This is a vector of migration rates for the coalescent simulation. These values determine the amount of population differentiation among the five populations simulated. 
4. Seedset (line32): Set from 1-3. These are seeds for simulating data. In the manuscript we simulated 16 datasets with varying migration rates using the same random seeds, and then changed the random seeds twice to simulate the remaining 32 for each marker type. Setting this to 3 will run the full set of simulations.

#### The data sets
Gzipped versions of the data sets that get simulated in the script are available in the files
inside `SNPuSAT_InFiles/SNPs_Simulations` and `SNPuSAT_InFiles/uSat_Simulations`.

### Flockture versus Flock

The comparison made between flock and flockture can be re-created using the R
script `FLOCKvFLOCKture.R`.  There is a little bit of  hand-curation that has 
to go on because FLOCK is not scriptable. The data needed for these, and the
FLOCK output are in `FLOCKvFLOCKTURE_data`.  The FLOCK output, in particular, is
a bit of a beast since it is all stored in Microsoft's xlsx format that doesn't appear
to be very efficient.


## Terms 

As a work partially of the United States Government, this package is in the
public domain within the United States. 

See TERMS.md for more information.


