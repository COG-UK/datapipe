# Summary of published datapipe outputs

### Alignments

- `cog_<date>_all.fa` : all unaligned sequences after deduplication
- `cog_<date>_all_alignment.fa` : all aligned sequences after deduplication
- `cog_<date>_all_metadata.csv` : all corresponding metadata
- `cog_<date>_alignment.fa` : filtered, trimmed alignment with sequences matching those in the corresponding metadata
- `cog_<date>_metadata.csv` : corresponding metadata for filtered, trimmed alignment

### Cog

- `cog.insertions.tsv` and `cog.deletions.tsv` containing all found insertions and deletions for the UK sequences

### Metadata

- `cog_global_<date>_geography.csv` : metadata containing the following columns `"central_sample_id","sequence_name","sample_date","epi_week","country","adm1","adm2","outer_postcode","adm2_raw","adm2_source","NUTS1","region","latitude","longitude","location"`
- `cog_global_<date>_mutations.csv` : metadata containing the following columns `"sequence_name", "sample_date", "lineage","lineages_version"` and additionally columns for specifically typed mutations of interest
- `cog_global_<date>_public.csv` : metadata containing the following columns `"sequence_name","cog_id","gisaid_id","sample_date","epi_week","country","adm1","pillar_2","is_surveillance","is_travel_history","travel_history","lineage","lineages_version"`
- `cog_global_<date>_consortium.csv` : metadata containing all columns as in the public metadata, extended with the following columns `"submission_org_code","root_sample_id","adm2","outer_postcode","adm2_raw","adm2_source","NUTS1","region","latitude","longitude","location","source_age","source_sex","sample_type_collected","sample_type_received","swab_site","ct_n_ct_value","ct_n_test_kit","ct_n_test_platform","ct_n_test_target"`
- `cog_<date>_unlinked.csv` : shuffled metadata with no ids containing the following columns `"sample_date","epi_week", "location","lineage","lineages_version","is_surveillance", "collection_pillar", "pillar_2"`

### Public

- `cog_<date>_all.fa` : all unaligned sequences after deduplication
- `cog_<date>_unmasked_alignment.fa` : all aligned sequences
- `cog_<date>_alignment.fa` : filtered, trimmed alignment with sequences matching those in the corresponding metadata
- `cog_<date>_metadata.csv` : corresponding metadata for filtered, trimmed alignment with the following columns `"sequence_name", "country","adm1","pillar_2","sample_date", "epi_week","lineage","lineages_version"`



