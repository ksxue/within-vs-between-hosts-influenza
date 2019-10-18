# within-vs-between-hosts-influenza
Linking influenza virus evolution within and between human hosts

This repository contains code and small intermediate data files associated with a comparison of influenza's evolutionary rates within and between hosts. All scripts are meant to be run from the top level of the repository.

Data availability: Note that due to data-sharing restrictions, sequences downloaded from GISAID are not made public in this repository.

Overview
--------

The directory is organized as follows:

    project
    |- README			# Overall description of repository.
	|
	|- analysis			# Custom workflows to perform analyses in the manuscript.
	|					# Each analysis directory contains a "Run.sh" script that performs all analyses.
    |  |- DownloadDataCallVariants	# Download and organize raw data, map to influenza genomes, and summarize reads at each site.
	|  |- CallWithinHostVariants		# Call within-host variants in each sample.
	|  |- CalculateAcuteRates			# Calculate rates of within-host evolution in acute infections.
	|  |- CalculateChronicRates					# Calculate rates of within-host evolution in chronic infections.
	|  |- CalculateGlobalRates				# Calculate global rates of evolution.
	|  |- CompareVariantCallingThresholds				# Analyze the distribution of within-host diversity at different variant-frequency thresholds.
	|  |- InferGlobalVariants				# Calculate the proportion of NS variants in the global influenza population.
	|  |- figures-WvB					# Generate the figures for the accompanying manuscript.
	|- data				# Small data files, primarily associated with data organization.
	|  |- metadata					# Contains metadata and organizational tracking information for each study.
	|- pipelines		# Standardized computational pipelines for alignment and analysis of sequencing data.
	|- reference		# Reference genomes for H3N2 influenza.
	|- scripts			# Custom C++, R, and shell scripts for calling and annotating variants.

    
