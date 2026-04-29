#!/bin/bash

o2-analysis-propagationservice -b  --configuration json:///Users/jimun/Analysis/JetMC/new-my-config.json --aod-memory-rate-limit 1000000000 --shm-segment-size 8000000000 --min-failure-level error  | \
o2-analysis-tracks-extra-v002-converter -b  --configuration json:///Users/jimun/Analysis/JetMC/new-my-config.json --aod-memory-rate-limit 1000000000 --shm-segment-size 8000000000 --min-failure-level error  | \
o2-analysis-event-selection-service -b  --configuration json:///Users/jimun/Analysis/JetMC/new-my-config.json --aod-memory-rate-limit 1000000000 --shm-segment-size 8000000000 --min-failure-level error  | \
o2-analysis-mccollision-converter -b  --configuration json:///Users/jimun/Analysis/JetMC/new-my-config.json --aod-memory-rate-limit 1000000000 --shm-segment-size 8000000000 --min-failure-level error  | \
o2-analysis-multcenttable -b  --configuration json:///Users/jimun/Analysis/JetMC/new-my-config.json --aod-memory-rate-limit 1000000000 --shm-segment-size 8000000000 --min-failure-level error  | \
o2-analysis-je-jet-deriveddata-producer -b  --configuration json:///Users/jimun/Analysis/JetMC/new-my-config.json --aod-memory-rate-limit 1000000000 --shm-segment-size 8000000000 --min-failure-level error  | \
o2-analysis-lf-mc-centrality -b  --configuration json:///Users/jimun/Analysis/JetMC/new-my-config.json --aod-memory-rate-limit 1000000000 --shm-segment-size 8000000000 --min-failure-level error  | \
o2-analysis-trackselection -b  --configuration json:///Users/jimun/Analysis/JetMC/new-my-config.json --aod-memory-rate-limit 1000000000 --shm-segment-size 8000000000 --min-failure-level error  | \
o2-analysis-pid-tpc-service -b  --configuration json:///Users/jimun/Analysis/JetMC/new-my-config.json --aod-memory-rate-limit 1000000000 --shm-segment-size 8000000000 --min-failure-level error  | \
o2-analysis-pid-tof-full -b  --configuration json:///Users/jimun/Analysis/JetMC/new-my-config.json --aod-memory-rate-limit 1000000000 --shm-segment-size 8000000000 --min-failure-level error  | \
o2-analysis-pid-tof-base -b  --configuration json:///Users/jimun/Analysis/JetMC/new-my-config.json --aod-memory-rate-limit 1000000000 --shm-segment-size 8000000000 --min-failure-level error  | \
o2-analysis-je-jet-finder-mcd-charged -b  --configuration json:///Users/jimun/Analysis/JetMC/new-my-config.json --aod-memory-rate-limit 1000000000 --shm-segment-size 8000000000 --min-failure-level error  | \
o2-analysis-je-jet-finder-mcp-charged -b  --configuration json:///Users/jimun/Analysis/JetMC/new-my-config.json --aod-memory-rate-limit 1000000000 --shm-segment-size 8000000000 --min-failure-level error  | \
o2-analysis-lf-kstarinoo -b  --configuration json:///Users/jimun/Analysis/JetMC/new-my-config.json --aod-memory-rate-limit 1000000000 --shm-segment-size 8000000000 --min-failure-level error
