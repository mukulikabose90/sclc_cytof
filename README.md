# CyTOF-based profiling of circulating tumor cells reveals treatment-associated subtypes and phenotypes in small-cell lung cancer liquid biopsies.

This project is fully reproducible using Docker and `renv`. 

---

## Prerequisites

- [Docker Desktop](https://www.docker.com/products/docker-desktop) installed on your system  
- Git (to clone this repository)

## Project Structure

```
sclc_cytof/   
├─ Dockerfile  
├─ Makefile
├─ renv.lock 
├─ renv/ 
├─ source/
│   └─ final_scripts/             # Scripts to be run   
│       └─ run_all_scripts.R      # Master script to reproduce all results 
├─ data/                          # Input data files and intermediate files          
│   ├─ cytof_raw_data/                  (available to download below)
│   ├─ run_stats/
│   ├─ cytof_panel_info.csv
│   └─ cytof_metadata.csv                                 
└─ figures/                       # Generated figures 
```  


- `Dockfile` contains environment specifications
- `Makefile` allows for simple reproducibility
- `renv.lock` and `renv/` locks the exact package versions used in the project 

## Data download

- All data needed to reproduce results available at [LINK](https://zenodo.org/records/17651946?preview=1&token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6IjBkNzE3YjlhLTc3YTQtNGI5ZS1hMjNkLWY0OWVmODNjYWQzNiIsImRhdGEiOnt9LCJyYW5kb20iOiIwMTZkZGFhNmM5NDZiMWI4MDU1ZWE1MGRhYzkyOTIxNSJ9.RUBmxpc3x5K9aS8Oc7kJc0Ls6x4WwL_azRQ4uMoUpQoz5Bw09DMFt3milQVPtckT8ylXKpuWZ6d5Aa9PaEZyqw)  

- Download all data into `data/` in project directory

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
make reproduce
```