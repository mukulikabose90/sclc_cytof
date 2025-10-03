# CyTOF-based profiling of circulating tumor cells reveals treatment-associated phenotypes and subtype switching in small cell lung cancer liquid biopsies.

This project is fully reproducible using Docker and `renv`. 

---

## Prerequisites

- [Docker Desktop](https://www.docker.com/products/docker-desktop) installed on your system  
- Git (to clone this repository)

## Project Structure

sclc_cytof/    
├─ Dockerfile  
├─ renv.lock  
├─ run_all_scripts.R # Master script to reproduce all results.    
├─ source/final_scripts # scripts to be run.   
├─ data/ # Input data files.   
└─ figures/ # Generated output.   

- `run_all_scripts.R` executes all analysis steps in the correct order.  
- `renv.lock` locks the exact package versions used in the project.  

## Quick Start

1. **Clone the repository**

```bash     
git clone https://github.com/coleruoff/sclc_cytof.git      
cd sclc_cytof

2. 