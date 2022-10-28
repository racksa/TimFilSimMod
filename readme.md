This repository contains a snapshot of my code as it was on 28th October 2022.

# Compiling the code
Each machine the code is designed to run on has a custom entry in the makefile of the form `cilia_XXX` that should be used to generate the appropriate executable. In general, the code should be compiled by calling
```
make cilia_clean cilia_XXX
```

Note: To ensure you have access to mobility options based on code in other projects (e.g. FCM implemented using UAMMD), the code should have been downloaded using
```
git clone https://github.com/timwestwood/TimFilSim --recursive && git submodule update --init --recursive
```
and whenever you pull an updated version of the code you should use
```
git pull --recurse && git submodule update --recursive
```

# Running the code
**On a free-use machine**: The GPUs on which the code should run must be set using the environment variable `CUDA_VISIBLE_DEVICES` in `run_cilia_sim.sh`. The GPU identifiers do not need to be in ascending order. The code can then be run by calling `./run_cilia_sim.sh`, although it may need to be assigned the correct permissions first (as discussed inside the script itself).

**On Imperial HPC using the PBS queue system**: The execution will be allocated specific GPUs to use in accordance with the requested resources for the job; these GPUs will be defined by `CUDA_VISIBLE_DEVICES`. The request for resources and eventual execution of the code are done by modifying the bash script `run_cilia_sim.pbs` and then running the script `queue_cilia_sim.sh`. If no argument is supplied the job will be placed on the normal/free queue, or if the argument `express` is passed it will go to the express queue (n.b. the script must be modified to replace `XXXXX` with the express account number). `queue_cilia_sim.sh` will create a directory, named according to the job's ID in the PBS queue, into which it will automatically move the executable (if any other input data is required, the script should be modified to move this too), as well as all data produced by the simulation. *Note*: For the free queue, you can provide a realistic estimate for the resources required and the queue system will allocate you a node with at least enough resourses, but for the express queue you must request resources that match one of classes they provide; e.g. to run on a 2-GPU node, you must also ask for 8 CPUs and 48GB of memory, whether you need it or not. This is how the free queue used to work.
