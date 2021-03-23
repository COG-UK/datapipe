# Datapipe

Pipeline to process SARS-CoV-2 sequences and metadata, clean up irregularities, align and variant call then publish matched subsets of FASTA sequences and metadata for groups with different access to sensitive data. 

Runs weekly on global sequences downloaded from GISAID.

Runs daily on COG-UK sequences, and combines with non-UK GISAID sequences.

### Pipeline Overview

#### GISAID processing

1. Parse GISAID dump (`export.json`) and extract FASTA of sequences and associated metadata. 

   - Excludes known problematic sequences listed in `gisaid_omissions.txt`

   - Excludes sequences where `covv_host.lower() != 'human'` 
   - Excludes sequences where malformed (not `YYYY-MM-DD`) or impossible (earlier than `2019-11-30` or later than today) date in `covv_collection_date`
   - Reformat FASTA header
   - Add `epi-week` and `epi-day` columns to metadata

2. Run `pangolin` (https://github.com/cov-lineages/pangolin) on all new sequences. If new release of `pangolin` run on all sequences.
3. Calculate the `unmapped_genome_completeness` as the proportion of sequence length which is unambiguous (not `N`)
4. Deduplicate by date, keeping the earliest example
5. Align to the reference (`Wuhan/WH04/2020`) with `minimap2`
6. Variant call using `gofasta` and type specific mutations of interest listed in `AAs.csv` and `dels.csv`
7. Filter out low quality sequences with mapped completeness < 93%, and trim and pad alignment outside of reference coordinates `265:29674`
8. Calculate distance to reference and exclude sequences with distance to more than 4.0 epi-week std devs.

#### COG-UK processing

1. Parse matched FASTA and metadata TSV output by Elan/Majora

   - Reformats header and unaligns sequences which have already been aligned to the reference

   - Manual date correction for samples listed in `date_corrections.csv`
   - Excludes early sequences which have been resequenced as listed in `resequencing_omissions.txt`
   - Adds GISAID accession if recently submitted

   - Excludes sequences where malformed (not `YYYY-MM-DD`) or impossible (earlier than `2019-11-30` or later than today) date in `covv_collection_date`
   - Add `epi-week` and `epi-day`, `source_id` and `pillar_2` columns to metadata

2. Run `pangolin` (https://github.com/cov-lineages/pangolin) on all new sequences. If new release of `pangolin` run on all sequences.
3. Calculate the `unmapped_genome_completeness` as the proportion of sequence length which is unambiguous (not `N`)
4. Deduplicate COG-ID by completeness and label samples with duplicate `source_id`
5. Align to the reference (`Wuhan/WH04/2020`) with `minimap2`
6. Variant call using `gofasta` and type specific mutations of interest listed in `AAs.csv` and `dels.csv`
7. Filter out low quality sequences with mapped completeness < 93%, and trim and pad alignment outside of reference coordinates `265:29674`
8. Clean up geographical metadata (https://github.com/COG-UK/geography_cleaning)
9. Combine COG-UK sequences and metadata with non-UK GISAID sequences and metadata
10. Publish subsets of the data as described in `publish_cog_global_recipes.json`

### What is grapevine?

`grapevine` (https://github.com/COG-UK/grapevine) was the name of the original pipeline which did all of the above, made phylogenetic trees and more. As the number of sequences has grown the tree building steps take increasingly long to complete. As the majority of users only interact with the alignments and cleaned metadata, it was decided that a robust implementation of the alignment and metadata processing steps run daily would be more useful and that is what is provided here.
