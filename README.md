# Alliance of Genome Resources Pre-processing
Take source DQM/ferret files from FMS to generate and upload new files to FMS for the agr_loader to run on.

## Running Pre-processing
- command line run with :
  - python3 src/aggregate_preprocessor.py

- docker run with :
  - docker build -t preprocess . && docker run -e ALLIANCE_RELEASE=3.1.0 -e FMS_API_URL=https://fmsdev.alliancegenome.org -e API_KEY=<key> --entrypoint /usr/bin/python3 preprocess src/aggregate_preprocessor.py

## Editing variables
- Set alliance release version and api url by passing environment variables to docker, or edit src/default_env_vars.yml if running from command line for src/data_manager/data_file_manager.py  to get latest files for each datatype + dataSubType in the config src/config/default.yml 
  - FMS_API_URL: "https://fmsdev.alliancegenome.org"
  - ALLIANCE_RELEASE: "3.1.0"

## Caveats
- Cannot use the release snapshot, because this generates files and uploads to FMS, which need to be in the release snapshot. (Could change this later if snapshots get taken after agr_ferret runs and before agr_preprocess runs, and then another release snapshot gets taken for agr_loader to run)
- Processing of interaction-source files must be done sequentially because (publication-taxon1-taxon2) entries from earlier sources take precedence over those from later sources.
