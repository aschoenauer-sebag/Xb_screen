#!/bin/bash
setenv PYTHONPATH /cbio/donnees/twalter/workspace/cecog/pysrc:/cbio/donnees/twalter/data/migration/src
cd /cbio/donnees/twalter/data/migration/src/tracking/simulator
python extract_features.py --what "full displacements" 

