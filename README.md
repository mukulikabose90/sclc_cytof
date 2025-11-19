# CyTOF-based profiling of small-cell lung cancer circulating tumor cells reveals treatment-associated subtypes and phenotypes in patient liquid biopsies

This project is fully reproducible using Docker and `renv`. 

---

## Prerequisites

- [Docker Desktop](https://www.docker.com/products/docker-desktop) installed on your system  
- Git (to clone this repository)

## Project Structure

sclc_cytof/    
├─ Dockerfile  
├─ Makefile
├─ renv.lock  
├─ run_all_scripts.R # Master script to reproduce all results.    
├─ source/final_scripts # scripts to be run.   
├─ data/ # Input data files.   
└─ figures/ # Generated output.   

- `run_all_scripts.R` executes all analysis steps in the correct order.  
- `renv.lock` locks the exact package versions used in the project.  
- `Dockfile` contains environment specifications
- `Makefile` allows for simple reproducibility


## Quick Start

1. **Clone the repository**

```bash     
git clone https://github.com/coleruoff/sclc_cytof.git      
cd sclc_cytof
```

2. **Build the environment**

```bash 
make build
```

3. **Reproduce results**

```bash 
make build
```