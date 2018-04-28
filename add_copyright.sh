#!/bin/bash

#https://github.com/cloudposse/copyright-header


copyright-header --license GPL3 --copyright-software "Canadian Hydrological Model"  --copyright-software-description "The Canadian Hydrological Model (CHM) is a novel modular unstructured mesh based approach for hydrological modelling" --copyright-holder "Christopher Marsh" --copyright-year 2018 --add-path src/:tools/--remove-path ./src/modules/snowpack:./wiki/:./doc/ -o .
