# DATA-SPECIFIC INFORMATION FOR: [Supplementary_Table3.tsv]

* Number of variables: 15
* Number of cases/rows: 136964
* Variable List: variable name, description, unit, and value labels/ranges
| Variable | Description | Unit | Value labels / ranges |
| --- | --- | --- | --- |
| run_id | Sequencing run accession | n/a | free text / categorical (6870 unique values) |
| palm_id | Palmprint cluster identifier | n/a | free text / categorical (1703 unique values) |
| coverage | Read coverage | unitless | [0.45, 40683] |
| sOTU | sub-operational taxonomic unit identifier | n/a | free text / categorical (1165 unique values) |
| qseqid | Query sequence identifier | n/a | free text / categorical (35 unique values) |
| pident | Percent identity of the alignment | unitless | [39, 100] |
| evalue | Alignment expectation value | unitless | [0, 1] |
| sra_sequence | SRA hit sequence | n/a | free text / categorical (1748 unique values) |
| biosample_id | BioSample accession identifier | n/a | free text / categorical (5931 unique values) |
| scientific_name | Scientific name | n/a | free text / categorical (1009 unique values) |
| nickname | Nickname | n/a | free text / categorical (1165 unique values) |
| date | Collection or processing date | ISO-8601 | [2012-01-21, 2021-01-11] |
| lng | Longitude | degrees | [-180, 180] |
| lat | Latitude | degrees | [-90, 90] |
| file | Source file name or path | n/a | free text / categorical (35 unique values) |
* Missing data codes: *list code/symbol and definition*
  - NaN/blank parsed by pandas: found 87167 value(s)
* Specialized formats or other abbreviations used: datetime columns: date

# DATA-SPECIFIC INFORMATION FOR: [Supplementary_Table4.tsv]

* Number of variables: 164
* Number of cases/rows: 253
* Variable List: variable name, description, unit, and value labels/ranges
| Variable | Description | Unit | Value labels / ranges |
| --- | --- | --- | --- |
| run_accession | Sequencing run accession | n/a | free text / categorical (253 unique values) |
| study_accession | Study accession identifier | n/a | free text / categorical (69 unique values) |
| study_title | Study title | n/a | free text / categorical (69 unique values) |
| experiment_accession | Experiment accession identifier | n/a | free text / categorical (251 unique values) |
| experiment_title | Experiment title | n/a | free text / categorical (149 unique values) |
| experiment_desc | Experiment description | n/a | free text / categorical (149 unique values) |
| organism_taxid | NCBI taxonomic identifier for the organism | unitless | [2708, 2.69692e+06] |
| organism_name | Organism name | n/a | free text / categorical (51 unique values) |
| library_name | library name | n/a | free text / categorical (198 unique values) |
| library_strategy | Library strategy | n/a | {'RNA-Seq', 'WGA', 'WGS'} |
| library_source | Library source | n/a | {'GENOMIC', 'METAGENOMIC', 'METATRANSCRIPTOMIC', 'TRANSCRIPTOMIC', 'TRANSCRIPTOMIC SINGLE CELL'} |
| library_selection | Library selection method | n/a | free text / categorical (11 unique values) |
| library_layout | Library layout | n/a | {'PAIRED', 'SINGLE'} |
| sample_accession | Sample accession identifier | n/a | free text / categorical (247 unique values) |
| sample_title | Sample title | n/a | n/a |
| biosample | BioSample accession identifier | n/a | free text / categorical (247 unique values) |
| bioproject | BioProject accession identifier | n/a | free text / categorical (69 unique values) |
| instrument | Instrument | n/a | free text / categorical (11 unique values) |
| instrument_model | Instrument model | n/a | free text / categorical (11 unique values) |
| instrument_model_desc | Instrument model description | n/a | {'ILLUMINA'} |
| total_spots | Total number of spots in the dataset | n/a | [102819, 5.79302e+08] |
| total_size | Total size of the dataset | n/a | [3.29977e+07, 4.05363e+10] |
| run_total_spots | Total number of spots for the sequencing run | unitless | [102819, 1.70209e+08] |
| run_total_bases | Total number of bases for the sequencing run | unitless | [5.16151e+07, 3.65976e+10] |
| run_alias | Submitter-provided run alias | n/a | free text / categorical (250 unique values) |
| public_filename | Publicly available file name | n/a | free text / categorical (253 unique values) |
| public_size | Size of the public file | unitless | [1.32465e+07, 2.35906e+10] |
| public_date | Public release date | ISO-8601 | [2015-05-28, 2022-09-29] |
| public_md5 | MD5 checksum for the public file | n/a | free text / categorical (253 unique values) |
| public_version | Version number of the public file | unitless | [1, 2] |
| public_semantic_name | Public file format label | n/a | {'SRA Lite', 'SRA Normalized', 'cram', 'fastq'} |
| public_supertype | Public file distribution class | n/a | {'Original', 'Primary ETL'} |
| public_sratoolkit | Flag indicating whether the file is accessible via SRA Toolkit | unitless | {0, 1} |
| aws_url | Amazon Web Services URL | n/a | free text / categorical (253 unique values) |
| aws_free_egress | Amazon Web Services free-egress region or policy | n/a | {'-', 's3.us-east-1', 'worldwide'} |
| aws_access_type | Amazon Web Services access type | n/a | {'Use Cloud Data Delivery', 'anonymous', 'aws identity'} |
| public_url | public url | n/a | free text / categorical (253 unique values) |
| ncbi_url | NCBI download URL | n/a | free text / categorical (253 unique values) |
| ncbi_free_egress | NCBI free-egress policy | n/a | {'worldwide'} |
| ncbi_access_type | NCBI access type | n/a | {'anonymous'} |
| gcp_url | Google Cloud Storage URL | n/a | free text / categorical (253 unique values) |
| gcp_free_egress | Google Cloud free-egress region or policy | n/a | {'gs.us-east1'} |
| gcp_access_type | Google Cloud access type | n/a | {'gcp identity'} |
| experiment_alias | experiment alias | n/a | free text / categorical (19 unique values) |
| strain | microbial or eukaryotic strain name | n/a | {'C57BL/6', 'MdSGHV', 'fenpropathrin-resistant', 'fenpropathrin-susceptible'} |
| isolate | identification or description of the specific individual from which this sample was obtained | n/a | free text / categorical (33 unique values) |
| host | The natural (as opposed to laboratory) host to the organism from which the sample was obtained. Use the full taxonomi... | n/a | free text / categorical (16 unique values) |
| isolation_source | Describes the physical, environmental and/or local geographical source of the biological sample from which the sample... | n/a | free text / categorical (20 unique values) |
| collection_date | the date on which the sample was collected; date/time ranges are supported by providing two dates from among the supp... | ISO-8601 | [2008-07-06, 2019-01-01] |
| geo_loc_name | Geographical origin of the sample; use the appropriate name from this list https://www.insdc.org/submitting-standards... | n/a | free text / categorical (48 unique values) |
| sample_type | Sample type, such as cell culture, mixed culture, tissue sample, whole organism, single cell, metagenomic assembly | n/a | {'2786_HFKKJALXX_L5', 'cells', 'single cell'} |
| host_tissue_sampled | name of body site where the sample was obtained from, such as a specific organ or tissue, e.g., tongue, lung. For fou... | n/a | {'Gut'} |
| biosamplemodel | biosample model or type | n/a | {'Human', 'Invertebrate', 'MIGS/MIMS/MIMARKS.human-gut', 'MIGS/MIMS/MIMARKS.water', 'Metagenome or environmental', ..... |
| lat_lon | The geographical coordinates of the location where the sample was collected. Specify as degrees latitude and longitud... | degrees | formatted latitude/longitude pair |
| breed | breed name - chiefly used in domesticated animals or plants | n/a | {'Missing', 'bashibai sheep', 'citrus', 'not applicable', 'not collected'} |
| tissue | Type of tissue the sample was taken from. | n/a | free text / categorical (32 unique values) |
| gold ecosystem classification | gold ecosystem classification | n/a | {'Engineered | Bioreactor | Anaerobic | Digestate | Unclassified', 'Environmental | Terrestrial | Soil | Floodplain |... |
| cultivar | cultivar name - cultivated variety of plant | n/a | free text / categorical (14 unique values) |
| age | age at the time of sampling; relevant scale depends on species and study, e.g. could be seconds for amoebae or centur... | n/a | free text / categorical (16 unique values) |
| biomaterial_provider | name and address of the lab or PI, or a culture collection identifier | n/a | {'Carolina Biological Supply Company', 'Hui-Ling Liao', 'Mark Coggeshall, University of Missouri', 'Prof. Dr. med Geo... |
| collected_by | Name of persons or institute who collected the sample | n/a | {'Dr Brett Williams', 'Gretchen A Gerrish, James G Morin, Nicholai M Hensley, Trevor J Rivers, ...} (10 levels) |
| biological_replicate | biological replicate number; a biological replicate is a sample taken from a different individual or experimental uni... | n/a | [1, 3] |
| replicate | replicate | n/a | free text / categorical (12 unique values) |
| samp_mat_process | Processing applied to the sample during or after isolation | n/a | {'Total RNA was extracted from each individual using Trizol-Chloroform extractions', 'none'} |
| samp_size | Amount or size of sample (volume, mass or area) that was collected | n/a | {'45 individuals'} |
| ref_biomaterial | Primary publication or genome report | n/a | n/a |
| rel_to_oxygen | Is this organism an aerobe, anaerobe? Please note that aerobic and anaerobic are valid descriptors for microbial envi... | n/a | n/a |
| samp_collect_device | Method or device employed for collecting sample | n/a | n/a |
| ecotype | a population within a given species displaying genetically based, phenotypic traits that reflect adaptation to a loca... | n/a | {'B73', 'mising', 'not applicable', 'not collected'} |
| identified_by | name of the taxonomist who identified the specimen | n/a | {'Gretchen A Gerrish, James G Morin, Nicholai M Hensley, Trevor J Rivers, Todd H Oakley, ...} (6 levels) |
| sex | physical sex of sampled organism | n/a | {'female', 'male', 'not collected', 'not determined'} |
| dev_stage | Developmental stage at the time of sampling. | n/a | free text / categorical (14 unique values) |
| ena-first-public | ENA first public release date | n/a | {'2018-06-11 19:03:00', '2019-02-08 18:04:04', '2019-02-26 18:02:58', '2019-04-19 06:02:37'} |
| ena-last-update | ENA last update date | n/a | {'2018-06-05 17:33:43', '2018-06-22 15:33:29', '2018-12-17 15:58:17', '2019-02-21 10:18:26'} |
| external id | external id | n/a | {'SAMEA3923215', 'SAMEA4708865', 'SAMEA4753698', 'SAMEA5185551', 'SAMEA5365958'} |
| insdc center name | INSDC center name | n/a | {'Center for Biotechnology (CeBiTec) - Bielefeld University', 'EAWAG', 'Faculty of Medicine, Goethe-University Frankf... |
| insdc first public | INSDC first public release date | n/a | {'2018-03-30 19:02:14', '2018-06-11 19:03:00', '2019-02-08 18:04:04', '2019-02-26 18:02:58', '2019-04-19 06:02:37'} |
| insdc last update | INSDC last update date | n/a | {'2016-10-21 11:43:06', '2018-06-05 17:33:43', '2018-06-22 15:33:29', '2018-12-17 15:58:17', '2019-02-21 10:18:26'} |
| insdc status | INSDC status | n/a | {'public'} |
| submitter id | Submitter name | n/a | {'20730x94', '6bcb1c40-540d-11e6-9a3d-3c4a9275d6c8', 'E-MTAB-6943:IgG-ACM_1, 2, 3', ...} (7 levels) |
| common name | Common name | n/a | {'human', 'viral metagenome'} |
| sample name | sample name in source database | n/a | {'20730x94', '6bcb1c40-540d-11e6-9a3d-3c4a9275d6c8', 'E-MTAB-6943:IgG-ACM_1, 2, 3', ...} (7 levels) |
| scientific_name | Scientific name | n/a | {'Homo sapiens', 'soil metagenome', 'viral metagenome'} |
| phenotype | Phenotype of sampled organism. For Phenotypic quality Ontology (PATO) (v1.269) terms, please see http://bioportal.bio... | n/a | {'cyan gray', 'white', 'wild type phenotype'} |
| experimental factor: compound | experimental factor - compound | n/a | {'ryegrass silage'} |
| experimental factor: dose | experimental factor - dose | n\a | [1500, 1500] |
| ebi_url | EBI download URL | n/a | {'http://ftp.sra.ebi.ac.uk/vol1/run/ERR135/ERR1356733/RU122R2.fq.gz', 'http://ftp.sra.ebi.ac.uk/vol1/run/ERR266/ERR26... |
| ebi_free_egress | EBI free-egress policy | n/a | {'worldwide'} |
| ebi_access_type | EBI access type | n/a | {'anonymous'} |
| broker name | broker name | n/a | {'ArrayExpress'} |
| altitude | The altitude of the sample is the vertical distance between Earth's surface above Sea Level and the sampled position... | n/a | n/a |
| env_broad_scale | Add terms that identify the major environment type(s) where your sample was collected. Recommend subclasses of biome... | unitless | environment (broad scale) |
| env_local_scale | Add terms that identify environmental entities having causal influences upon the entity at time of sampling, multiple... | unitless | environment (local scale) |
| env_medium | Add terms that identify the material displaced by the entity at time of sampling. Recommend subclasses of environment... | unitless | environment (medium scale) |
| season | The season when sampling occurred. Any of the four periods into which the year is divided by the equinoxes and solsti... | n/a | n/a |
| type_of_nucleic_acid | Record type | n/a | n/a |
| arrayexpress-species | Taxonomic species | n/a | {'viral metagenome'} |
| sample description | Sample identifier | n/a | {'Total nucleic acid extraction from human feces'} |
| experimental factor: growth condition | experimental factor - growth condition | n/a | {'MCF-7 conditioned media'} |
| experimental factor: stimulus | experimental factor - stimulus | n/a | {'none'} |
| experimental factor: immunoprecipitate | experimental factor - immunoprecipitate | n/a | {'anti-IgG'} |
| cell type | Type of cell of the sample or from which the sample was obtained. | n/a | {'S2 cells', 'macrophage'} |
| developmental stage | Developmental stage at the time of sampling. | n/a | {'adult'} |
| genotype | observed genotype | n/a | {'PB1+ Pi54', 'Wildtype', 'wild type genotype'} |
| growth condition | growth condition | n/a | {'MCF-7 conditioned media'} |
| indivudal | individual | n/a | {'pool of donors 1 to 3'} |
| organism part | Type of tissue the sample was taken from. | n/a | {'blood'} |
| progenitor cell type | Record type | n/a | {'peripheral blood mononuclear cell'} |
| source_name | Sample source name or description | n/a | {'Cortical thymic epithelial cells', 'Leaf', 'Medullary thymic epithelial cells', 'S2 cells', 'Skin epithelial cells'... |
| treatment | treatment applied to the sample, such as chemical treatment, mechanical wounding, rehydration, etc. Please note that... | n/a | free text / categorical (79 unique values) |
| sample_name | sample name in source database | n/a | {'Sample_1086', 'Sample_1235'} |
| elev | The elevation of the sampling site as measured by the vertical distance from mean sea level. | n/a | elevation |
| env_feature | Add terms that identify environmental entities having causal influences upon the entity at time of sampling, multiple... | n/a | {'Amazon River', 'ENVO:2100002'} |
| depth | Depth is defined as the vertical distance below surface, e.g. for sediment or soil samples depth is measured from sed... | m in text | {'0.5 m', '10 m', '14 m', '33 m'} |
| env_material | Add terms that identify the material displaced by the entity at time of sampling. Recommend subclasses of environment... | n/a | {'ENVO:00002003', 'water'} |
| env_biome | Add terms that identify the major environment type(s) where your sample was collected. Recommend subclasses of biome... | n/a | {'ENVO:00009003', 'river'} |
| ena first public | ENA first public release date | n/a | {'2018-03-30'} |
| ena last update | ENA last update date | n/a | {'2016-10-21'} |
| insdc center alias | INSDC center alias | n/a | {'EAWAG'} |
| collection date | the date on which the sample was collected; date/time ranges are supported by providing two dates from among the supp... | ISO-8601 | [2014, 2014] |
| environment (biome) | Add terms that identify the major environment type(s) where your sample was collected. Recommend subclasses of biome... | n/a | {'wastewater'} |
| environment (feature) | Add terms that identify environmental entities having causal influences upon the entity at time of sampling, multiple... | n/a | {'wastewater treatment plant'} |
| environment (material) | Add terms that identify the material displaced by the entity at time of sampling. Recommend subclasses of environment... | n/a | {'activated sludge'} |
| geographic location (country and/or sea) | Geographical origin of the sample; use the appropriate name from this list https://www.insdc.org/submitting-standards... | n/a | {'Switzerland'} |
| geographic location (latitude) | Sampling latitude | degrees | [-90, 90] |
| geographic location (longitude) | Sampling longitude | degrees | [-180, 180] |
| investigation type | Nucleic Acid Sequence Report is the root element of all MIGS/MIMS compliant reports as standardized by Genomic Standa... | n/a | {'metagenome'} |
| project name | A concise name that describes the overall project, for example "Analysis of sequences collected from Antarctic soil" | n/a | {'What biological and operational factors determine the release of antibiotic resistance genes from WWTPs?'} |
| sequencing method | Sequencing method | n/a | {'Illumina HiSeq'} |
| wastewater/sludge environmental package | wastewater/sludge environmental package | n/a | {'wastewater/sludge'} |
| cell_line | Name of the cell line. | n/a | {'P493-4', 'P493-5', 'P493-6'} |
| cell_type | Type of cell of the sample or from which the sample was obtained. | n/a | {'B-cell'} |
| nucleic_acid | nucleic acid type | n/a | {'rna'} |
| sample_timepoint | Sample identifier | n/a | {'SF05'} |
| host_subject_id | a unique identifier by which each subject can be referred to, de-identified, e.g. #131 | n/a | [5.01762e+07, 5.01762e+07] |
| isolation source | Describes the physical, environmental and/or local geographical source of the biological sample from which the sample... | n/a | n\a |
| isol_growth_condt | PMID or url for isolation and growth condition specifications | unitless | n/a |
| propagation | phage: lytic/lysogenic/temperate/obligately lytic; plasmid: incompatibility group; eukaryote: asexual/sexual | n/a | n/a |
| number_of_cells | number of cells | unitless | [500000, 2e+06] |
| target_molecules | target molecules | n/a | {'mRNA'} |
| site_name | site name | n/a | {'Macapa North (MCPN)', 'Macapa South (MCPS)', 'Obidos (OB)', 'Tapajos Surface (TAPS)'} |
| chlorophyll | concentration of chlorophyll | n/a | {'0.87 ug/L', '1.67 ug/L', '1.99 ug/L', '3.23 ug/L', '3.82 ug/L'} |
| conduc | electrical conductivity of water | uS/cm in text | {'15.9 uS/cm', '55.3 uS/cm', '56.2 uS/cm', '56.4 uS/cm', '61.5 uS/cm'} |
| diss_inorg_carb | dissolved inorganic carbon concentration | umol C/kg in text | {'155 umol C/kg', '459 umol C/kg', '502 umol C/kg', '507 umol C/kg', '551 umol C/kg'} |
| bacterial_count | Count | n/a | {'3.58E06 cells/mL', '3.63E06 cells/mL', '3.76E06 cells/mL', '3.94E06 cells/mL', '3.96E06 cells/mL'} |
| samp_volume | sampling volume | L (liters) in text format | {'0.65 L', '0.7 L', '0.73 L', '0.76 L', '1.06 L', ...} (6 levels) |
| min_filter_size | minimum filter size | n/a | {'0.2 um', '2 um', '2.0 Âµm'} |
| max_filter_size | maximum filter size | n/a | {'156 Âµm', '2.0 um', '297 um'} |
| sequencing_method | Sequencing method | n/a | {'Illumina PE 150x150', 'Illumina Rapid Run'} |
| sequencing_machine | Sequencing machine | n/a | {'HiSeq 2500'} |
| std_norm_factor | std normalization factor | unitless | [4217, 43832] |
| specimen_voucher | Identifier for the physical specimen. Use format: "[<institution-code>:[<collection-code>:]]<specimen_id>", eg, "UAM:... | n/a | {'TH-08-6182'} |
| source_material_id | unique identifier assigned to a material sample used for extracting nucleic acids, and subsequent sequencing. The ide... | n/a | n/a |
| sample_number | Sample identifier | unitless | n/a |
| sample_source | Sample identifier | n/a | n/a |
| investigation_type | Nucleic Acid Sequence Report is the root element of all MIGS/MIMS compliant reports as standardized by Genomic Standa... | n/a | {'Non-selective Metatranscriptome'} |
| temp | temperature of the sample at time of sampling | n/a | {'AMBIENT'} |
| geographic_location | Geographical origin of the sample; use the appropriate name from this list https://www.insdc.org/submitting-standards... | n/a | {'not collected'} |
| ena_fastq_http | ena fastq http | n/a | free text / categorical (131 unique values) |
| ena_fastq_http_1 | ena fastq http 1 | n/a | free text / categorical (122 unique values) |
| ena_fastq_http_2 | ena fastq http 2 | n/a | free text / categorical (121 unique values) |
| ena_fastq_ftp | ena fastq ftp | n/a | free text / categorical (131 unique values) |
| ena_fastq_ftp_1 | ena fastq ftp 1 | n/a | free text / categorical (122 unique values) |
| ena_fastq_ftp_2 | ena fastq ftp 2 | n/a | free text / categorical (121 unique values) |
* Missing data codes: *list code/symbol and definition*
  - 'missing': found 12 value(s)
  - 'none': found 20 value(s)
  - 'not applicable': found 286 value(s)
  - 'not collected': found 42 value(s)
  - 'unknown': found 3 value(s)
  - NaN/blank parsed by pandas: found 27738 value(s)
* Specialized formats or other abbreviations used: datetime columns: public_date, ena-first-public, ena-last-update, insdc first public, insdc last update, ena first public, ena last update

# DATA-SPECIFIC INFORMATION FOR: [Supplementary_Table5.tsv]

* Number of variables: 54
* Number of cases/rows: 510
* Variable List: variable name, description, unit, and value labels/ranges
| Variable | Description | Unit | Value labels / ranges |
| --- | --- | --- | --- |
| contig_id | Contig identifier | n/a | free text / categorical (510 unique values) |
| sequence | Nucleotide sequence | n/a | free text / categorical (510 unique values) |
| corresponding_srr | corresponding SRA Run accession identifier | n/a | free text / categorical (131 unique values) |
| assembler | Assembler name | n/a | {'megahit', 'spades'} |
| cluster_membership | cluster membership identifier | n/a | free text / categorical (56 unique values) |
| known_or_potentially_novel_tobamovirus | flag indicating whether the contig is a known or potentially novel tobamovirus (1) or not (0) | unitless | {0, 1} |
| contig_length | Contig length | bp | [606, 7832] |
| cluster_representative | cluster representative flag | unitless | {0, 1} |
| orf1_complete | flag indicating whether the contig is complete (1) or not (0) | unitless | {0, 1} |
| orf1_partial | flag indicating whether the contig is partial (1) or not (0) | unitless | {0, 1} |
| orf1_length | Sequence length | bp | [0, 1135] |
| orf1_start | Start position | bp | [0, 1] |
| orf1_stop | Stop position | bp | [0, 1] |
| orf1_tree_representative | flag indicating whether the contig is a tree representative (1) or not (0) | unitless | {0, 1} |
| orf2_rdrp_complete | flag indicating whether the contig is complete (1) or not (0) | unitless | {0, 1} |
| orf2_rdrp_partial | flag indicating whether the contig is partial (1) or not (0) | unitless | {0, 1} |
| orf2_rdrp_length | Sequence length | bp | [0, 592] |
| orf2_start | Start position | bp | [0, 1] |
| orf2_stop | Stop position | bp | [0, 1] |
| orf2_tree_representative | flag indicating whether the contig is a tree representative (1) or not (0) | unitless | {0, 1} |
| orf3_mp_complete | flag indicating whether the contig is complete (1) or not (0) | unitless | {0, 1} |
| orf3_mp_partial | flag indicating whether the contig is partial (1) or not (0) | unitless | {0, 1} |
| orf3_mp_length | Sequence length | bp | [0, 306] |
| orf3_start | Start position | bp | [0, 1] |
| orf3_stop | Stop position | bp | [0, 1] |
| orf3_tree_representative | flag indicating whether the contig is a tree representative (1) or not (0) | unitless | {0, 1} |
| orf4_cp_complete | flag indicating whether the contig is complete (1) or not (0) | unitless | {0, 1} |
| orf4_cp_partial | flag indicating whether the contig is partial (1) or not (0) | unitless | {0, 1} |
| orf4_cp_length | Sequence length | bp | [0, 186] |
| orf4_start | Start position | bp | [0, 1] |
| orf4_stop | Stop position | bp | [0, 1] |
| orf4_tree_representative | flag indicating whether the contig is a tree representative (1) or not (0) | unitless | {0, 1} |
| notes | Free-text notes | n/a | {'missing orf3', 'missing stop codon in orf3', 'missing stop codon in orf4', 'missing stop codons in orf2 and orf3',... |
| ground_truth_category | Ground-truth category label | n/a | {'mas', 'oth1', 'oth2', 'tob1', 'tob2', ...} (6 levels) |
| model_prediction | Model prediction label | unitless | {0, 1} |
| model_prediction_probabiility | Model prediction probability | unitless | [0, 1] |
| first_diamond_blastx_hit_name | first diamond blastx hit name | n/a | free text / categorical (40 unique values) |
| first_diamond_blastx_hit_identity_percent | Percent identity of the top BLASTX hit | % | [2.8, 99.9] |
| first_diamond_blastx_hit_alignment_length_aa | Amino-acid alignment length of the top BLASTX hit | bp | [53, 2810] |
| first_blastx_hit_name_accession | Accession identifier | n/a | free text / categorical (237 unique values) |
| first_blastx_hit_identity_percent | Percent identity of the top BLASTX hit | % | [22, 100] |
| first_blastx_hit_alignment_length_aa | Amino-acid alignment length of the top BLASTX hit | bp | [50, 1889] |
| first_blastx_hit_e_value | E-value of the top BLASTX hit | unitless | [0, 1] |
| first_blastx_hit_bit_score | Bit score of the top BLASTX hit | unitless | [64.3, 3428] |
| first_blastx_hit_protein | First blastx hit protein result | n/a | free text / categorical (212 unique values) |
| study_accession | Study accession identifier | n/a | free text / categorical (57 unique values) |
| study_title | Study title | n/a | free text / categorical (57 unique values) |
| organism_name | Organism name | n/a | free text / categorical (44 unique values) |
| submitter | Submitter name | n/a | free text / categorical (46 unique values) |
| country | Geographical origin of the sample; use the appropriate name from this list https://www.insdc.org/submitting-standards... | n/a | free text / categorical (44 unique values) |
| collection_date | the date on which the sample was collected; date/time ranges are supported by providing two dates from among the supp... | ISO-8601 | [2008-07-06, 2019-01-01] |
| source_sample_category | Broad source sample category | n/a | {'Aquatic', 'Host-associated', 'Terrestrial'} |
| source_sample_subcategory | Source sample subcategory | n/a | {'Animal', 'Animal gut', 'Freshwater', 'Human', 'Human gut', ...} (9 levels) |
| genbank_accession_number | GenBank accession number | n/a | free text / categorical (63 unique values) |
* Missing data codes: *list code/symbol and definition*
  - 'missing': found 6 value(s)
  - NaN/blank parsed by pandas: found 11000 value(s)
* Specialized formats or other abbreviations used: none specified

# DATA-SPECIFIC INFORMATION FOR: [Supplementary_Table6.tsv]

* Number of variables: 16
* Number of cases/rows: 125
* Variable List: variable name, description, unit, and value labels/ranges
| Variable | Description | Unit | Value labels / ranges |
| --- | --- | --- | --- |
| type | group type | n/a | {'outgroup', 'tobamo'} |
| genus | Taxonomic genus | n/a | {'Furovirus', 'Goravirus', 'Hordeivirus', 'Pecluvirus', 'Pomovirus', ...} (7 levels) |
| species | Taxonomic species | n/a | free text / categorical (60 unique values) |
| virus_name | Virus or viral taxon | n/a | free text / categorical (88 unique values) |
| accession | Accession identifier | n/a | free text / categorical (125 unique values) |
| RNA | Ribonucleic acid | n/a | {'RNA1', 'RNA2', 'RNA3'} |
| training | Training-set membership flag (training 1, not in training 0) | unitless | {0, 1} |
| length | Sequence length | bp | [1799, 7226] |
| orf1_ref | Reference identifier for ORF1 | n/a | free text / categorical (88 unique values) |
| orf2_ref | Reference identifier for ORF2 | n/a | free text / categorical (88 unique values) |
| cp_ref | Reference identifier for coat protein | n/a | free text / categorical (87 unique values) |
| record_id | UNCLEAR - human review required | n/a | free text / categorical (125 unique values) |
| orf1_record_id | Record identifier for ORF1 | n/a | free text / categorical (88 unique values) |
| orf2_record_id | Record identifier for ORF2 | n/a | free text / categorical (88 unique values) |
| cp_record_id | Record identifier for coat protein | n/a | free text / categorical (87 unique values) |
| sampling_prob | Sampling probability used during training set construction | unitless | [0, 1] |
* Missing data codes: *list code/symbol and definition*
  - NaN/blank parsed by pandas: found 224 value(s)
* Specialized formats or other abbreviations used: none specified
