# Summary of published datapipe outputs

### Alignments

- `cog_<date>_all.fa` : all unaligned sequences after deduplication
- `cog_<date>_all_alignment.fa` : all aligned sequences after deduplication
- `cog_<date>_all_metadata.csv` : all corresponding metadata
- `cog_<date>_alignment.fa` : filtered, trimmed alignment with sequences matching those in the corresponding metadata
- `cog_<date>_metadata.csv` : corresponding metadata for filtered, trimmed alignment

### Cog

- `cog.insertions.tsv` and `cog.deletions.tsv` containing all found insertions and deletions for the UK sequences
- `UTLA_genome_counts_<date>.csv` containing counts of delta sequences by date and UTLA

### Metadata

- `cog_global_<date>_geography.csv` : metadata containing the following columns `"central_sample_id","sequence_name","sample_date","epi_week","country","adm1","adm2","outer_postcode","adm2_raw","adm2_source","NUTS1","region","latitude","longitude","location"`
- `cog_global_<date>_mutations.csv` : metadata containing the following columns `"sequence_name", "sample_date", "lineage","lineages_version"` and additionally columns for specifically typed mutations of interest
- `cog_global_<date>_public.csv` : metadata containing the following columns `"sequence_name","cog_id","gisaid_id","sample_date","epi_week","country","adm1","is_pillar_2","is_surveillance","is_travel_history","travel_history","lineage","lineages_version"`
- `cog_global_<date>_consortium.csv` : metadata containing all columns as in the public metadata, extended with the following columns `"received_date","collection_date","published_date","sequencing_org_code","submission_org_code","submission_user","root_sample_id","adm2","outer_postcode","adm2_raw","adm2_source","NUTS1","region","latitude","longitude","location","utla","utla_code","suggested_adm2_grouping","source_age","source_sex","sample_type_collected","sample_type_received","swab_site","ct_n_ct_value","ct_n_test_kit","ct_n_test_platform","ct_n_test_target"`
- `cog_<date>_unlinked.csv` : shuffled metadata with no ids containing the following columns `"safe_sample_date","epi_week", "location","lineage","lineages_version","is_surveillance", "collection_pillar", "is_pillar_2"`
- `cog_global_<date>_epidemiology.csv` : metadata containing the following columns `"sequence_name","cog_id","gisaid_id","sample_date","epi_week","collection_date","received_date","sequencing_submission_date","sequencing_org_code","root_sample_id","biosample_source_id","country","adm1","adm2","utla","utla_code","outer_postcode","NUTS1","latitude","longitude","location","source_age","source_sex","collection_pillar","is_pillar_2","is_surveillance","is_travel_history","travel_history","lineage","lineage_support","lineages_version","scorpio_call","scorpio_support","ambiguity_count"`

### Public

- `cog_<date>_all.fa` : all unaligned sequences after deduplication
- `cog_<date>_unmasked_alignment.fa` : all aligned sequences
- `cog_<date>_alignment.fa` : filtered, trimmed alignment with sequences matching those in the corresponding metadata
- `cog_<date>_metadata.csv` : corresponding metadata for filtered, trimmed alignment with the following columns `"sequence_name", "country","adm1","is_pillar_2","sample_date", "epi_week","lineage","lineages_version"`

### Civet3
- `cog_global_<date>_private_alignment.fa` : masked, trimmed, filtered alignment of COG and GLOBAL sequences
- `cog_global_<date>_private_metadata.csv` : corresponding metadata with the following columns `"sequence_name","gisaid_id","cog_id","source_id","sample_date","epi_week","country","adm1","adm2","suggested_adm2_grouping","outer_postcode","is_surveillance","is_travel_history","travel_history","is_pillar_2","collection_pillar","lineage","lineages_version","scorpio_call"`
- `cog_global_<date>_mutations.csv` : metadata file produced by gofasta updown list, providing information about nucleotide mutations and ambiguous regions in aligned sequences

