# DATA-SPECIFIC INFORMATION FOR: [Supplementary_Table3.tsv]
*repeat this section for each dataset, folder or file, as appropriate*

* Number of variables: 15
* Number of cases/rows: 136964
* Variable List: *list variable name(s), description(s), unit(s) and value labels as appropriate for each*
  - run_id: description=Sequencing run accession; unit=n/a; value labels/range=free text / categorical (6870 unique values); dtype=object
  - palm_id: description=Palmprint cluster identifier; unit=n/a; value labels/range=free text / categorical (1703 unique values); dtype=object
  - coverage: description=Read coverage; unit=x; value labels/range=[0.45, 40683]; dtype=float64; review_status=review_recommended
  - sOTU: description=UNCLEAR - human review required; unit=n/a; value labels/range=free text / categorical (1165 unique values); dtype=object; review_status=unclear
  - qseqid: description=Query sequence identifier; unit=n/a; value labels/range=free text / categorical (35 unique values); dtype=object
  - pident: description=Percent identity of the alignment; unit=unitless; value labels/range=[39, 100]; dtype=float64
  - evalue: description=Alignment expectation value; unit=unitless; value labels/range=[0, 1]; dtype=int64
  - sra_sequence: description=SRA accession; unit=n/a; value labels/range=free text / categorical (1748 unique values); dtype=object; review_status=review_recommended
  - biosample_id: description=BioSample accession identifier; unit=n/a; value labels/range=free text / categorical (5931 unique values); dtype=object
  - scientific_name: description=Scientific name; unit=n/a; value labels/range=free text / categorical (1009 unique values); dtype=object
  - nickname: description=Nickname; unit=n/a; value labels/range=free text / categorical (1165 unique values); dtype=object; review_status=review_recommended
  - date: description=Collection or processing date; unit=ISO-8601; value labels/range=[2012-01-21, 2021-01-11]; dtype=datetime64[ns]; review_status=review_recommended
  - lng: description=Longitude; unit=degrees; value labels/range=[-180, 180]; dtype=float64; review_status=review_recommended
  - lat: description=Latitude; unit=degrees; value labels/range=[-90, 90]; dtype=float64; review_status=review_recommended
  - file: description=Source file name or path; unit=n/a; value labels/range=free text / categorical (35 unique values); dtype=object; review_status=review_recommended
* Variable dictionary source: S3_serratus_variable_dictionary.review_copy.yaml (review_copy)
* Missing data codes: *list code/symbol and definition*
  - NaN/blank parsed by pandas: found 87167 value(s)
* Specialized formats or other abbreviations used: datetime columns: date

# DATA-SPECIFIC INFORMATION FOR: [Supplementary_Table4.tsv]
*repeat this section for each dataset, folder or file, as appropriate*

* Number of variables: 164
* Number of cases/rows: 253
* Variable List: *list variable name(s), description(s), unit(s) and value labels as appropriate for each*
  - run_accession: description=Sequencing run accession; unit=n/a; value labels/range=free text / categorical (253 unique values); dtype=object
  - study_accession: description=Study accession identifier; unit=n/a; value labels/range=free text / categorical (69 unique values); dtype=object
  - study_title: description=Study title; unit=n/a; value labels/range=free text / categorical (69 unique values); dtype=object
  - experiment_accession: description=Experiment accession identifier; unit=n/a; value labels/range=free text / categorical (251 unique values); dtype=object
  - experiment_title: description=Experiment title; unit=n/a; value labels/range=free text / categorical (149 unique values); dtype=object
  - experiment_desc: description=Experiment description; unit=n/a; value labels/range=free text / categorical (149 unique values); dtype=object
  - organism_taxid: description=NCBI taxonomic identifier for the organism; unit=unitless; value labels/range=[2708, 2.69692e+06]; dtype=int64
  - organism_name: description=Organism name; unit=n/a; value labels/range=free text / categorical (51 unique values); dtype=object
  - library_name: description=UNCLEAR - human review required; unit=n/a; value labels/range=free text / categorical (198 unique values); dtype=object; review_status=unclear
  - library_strategy: description=Library strategy; unit=n/a; value labels/range={'RNA-Seq', 'WGA', 'WGS'}; dtype=object
  - library_source: description=Library source; unit=n/a; value labels/range={'GENOMIC', 'METAGENOMIC', 'METATRANSCRIPTOMIC', 'TRANSCRIPTOMIC', 'TRANSCRIPTOMIC SINGLE CELL'}; dtype=object
  - library_selection: description=Library selection method; unit=n/a; value labels/range=free text / categorical (11 unique values); dtype=object
  - library_layout: description=Library layout; unit=n/a; value labels/range={'PAIRED', 'SINGLE'}; dtype=object
  - sample_accession: description=Sample accession identifier; unit=n/a; value labels/range=free text / categorical (247 unique values); dtype=object
  - sample_title: description=Sample title; unit=UNCLEAR - human review required; value labels/range=UNCLEAR - human review required; dtype=float64; review_status=unclear
  - biosample: description=BioSample accession identifier; unit=n/a; value labels/range=free text / categorical (247 unique values); dtype=object; review_status=review_recommended
  - bioproject: description=BioProject accession identifier; unit=n/a; value labels/range=free text / categorical (69 unique values); dtype=object
  - instrument: description=UNCLEAR - human review required; unit=n/a; value labels/range=free text / categorical (11 unique values); dtype=object; review_status=unclear
  - instrument_model: description=Instrument model; unit=n/a; value labels/range=free text / categorical (11 unique values); dtype=object
  - instrument_model_desc: description=Instrument model description; unit=n/a; value labels/range={'ILLUMINA'}; dtype=object
  - total_spots: description=UNCLEAR - human review required; unit=UNCLEAR - human review required; value labels/range=[102819, 5.79302e+08]; dtype=int64; review_status=unclear
  - total_size: description=UNCLEAR - human review required; unit=UNCLEAR - human review required; value labels/range=[3.29977e+07, 4.05363e+10]; dtype=int64; review_status=unclear
  - run_total_spots: description=Total number of spots for the sequencing run; unit=unitless; value labels/range=[102819, 1.70209e+08]; dtype=int64
  - run_total_bases: description=Total number of bases for the sequencing run; unit=unitless; value labels/range=[5.16151e+07, 3.65976e+10]; dtype=int64
  - run_alias: description=Submitter-provided run alias; unit=n/a; value labels/range=free text / categorical (250 unique values); dtype=object
  - public_filename: description=Publicly available file name; unit=n/a; value labels/range=free text / categorical (253 unique values); dtype=object
  - public_size: description=Size of the public file; unit=unitless; value labels/range=[1.32465e+07, 2.35906e+10]; dtype=int64
  - public_date: description=Public release date; unit=ISO-8601; value labels/range=[2015-05-28, 2022-09-29]; dtype=datetime64[ns]
  - public_md5: description=MD5 checksum for the public file; unit=n/a; value labels/range=free text / categorical (253 unique values); dtype=object
  - public_version: description=Version number of the public file; unit=unitless; value labels/range=[1, 2]; dtype=int64
  - public_semantic_name: description=Public file format label; unit=n/a; value labels/range={'SRA Lite', 'SRA Normalized', 'cram', 'fastq'}; dtype=object
  - public_supertype: description=Public file distribution class; unit=n/a; value labels/range={'Original', 'Primary ETL'}; dtype=object
  - public_sratoolkit: description=Flag indicating whether the file is accessible via SRA Toolkit; unit=unitless; value labels/range={0, 1}; dtype=int64
  - aws_url: description=Amazon Web Services URL; unit=n/a; value labels/range=free text / categorical (253 unique values); dtype=object
  - aws_free_egress: description=Amazon Web Services free-egress region or policy; unit=n/a; value labels/range={'-', 's3.us-east-1', 'worldwide'}; dtype=object
  - aws_access_type: description=Amazon Web Services access type; unit=n/a; value labels/range={'Use Cloud Data Delivery', 'anonymous', 'aws identity'}; dtype=object
  - public_url: description=UNCLEAR - human review required; unit=n/a; value labels/range=free text / categorical (253 unique values); dtype=object; review_status=unclear
  - ncbi_url: description=NCBI download URL; unit=n/a; value labels/range=free text / categorical (253 unique values); dtype=object
  - ncbi_free_egress: description=NCBI free-egress policy; unit=n/a; value labels/range={'worldwide'}; dtype=object
  - ncbi_access_type: description=NCBI access type; unit=n/a; value labels/range={'anonymous'}; dtype=object
  - gcp_url: description=Google Cloud Storage URL; unit=n/a; value labels/range=free text / categorical (253 unique values); dtype=object
  - gcp_free_egress: description=Google Cloud free-egress region or policy; unit=n/a; value labels/range={'gs.us-east1'}; dtype=object
  - gcp_access_type: description=Google Cloud access type; unit=n/a; value labels/range={'gcp identity'}; dtype=object
  - experiment_alias: description=UNCLEAR - human review required; unit=n/a; value labels/range=free text / categorical (19 unique values); dtype=object; review_status=unclear
  - strain: description=UNCLEAR - human review required; unit=n/a; value labels/range={'C57BL/6', 'MdSGHV', 'fenpropathrin-resistant', 'fenpropathrin-susceptible'}; dtype=object; review_status=unclear
  - isolate: description=Isolate identifier or isolate name; unit=n/a; value labels/range=free text / categorical (33 unique values); dtype=object
  - host: description=Host organism; unit=n/a; value labels/range=free text / categorical (16 unique values); dtype=object; review_status=review_recommended
  - isolation_source: description=Material or environment from which the sample was isolated; unit=n/a; value labels/range=free text / categorical (20 unique values); dtype=object
  - collection_date: description=Collection date; unit=ISO-8601; value labels/range=[2012-01-01, 2017-01-01]; dtype=object
  - geo_loc_name: description=Geographic location name; unit=n/a; value labels/range=free text / categorical (48 unique values); dtype=object
  - sample_type: description=Sample identifier; unit=n/a; value labels/range={'2786_HFKKJALXX_L5', 'cells', 'single cell'}; dtype=object; review_status=review_recommended
  - host_tissue_sampled: description=Host organism; unit=n/a; value labels/range={'Gut'}; dtype=object; review_status=review_recommended
  - biosamplemodel: description=UNCLEAR - human review required; unit=n/a; value labels/range={'Human', 'Invertebrate', 'MIGS/MIMS/MIMARKS.human-gut', 'MIGS/MIMS/MIMARKS.water', 'Metagenome or environmental', 'Microbe, viral or environmental', 'Model organism or animal', 'Plant', 'Viral'}; dtype=object; review_status=unclear
  - lat_lon: description=Latitude/longitude coordinate pair; unit=degrees; value labels/range=formatted latitude/longitude pair; dtype=object
  - breed: description=UNCLEAR - human review required; unit=n/a; value labels/range={'Missing', 'bashibai sheep', 'citrus', 'not applicable', 'not collected'}; dtype=object; review_status=unclear
  - tissue: description=UNCLEAR - human review required; unit=n/a; value labels/range=free text / categorical (32 unique values); dtype=object; review_status=unclear
  - gold ecosystem classification: description=UNCLEAR - human review required; unit=n/a; value labels/range={'Engineered | Bioreactor | Anaerobic | Digestate | Unclassified', 'Environmental | Terrestrial | Soil | Floodplain | Unclassified', 'Environmental | Terrestrial | Soil | Riparian zone | Unclassified', 'Host-associated | Mammals | Digestive system | Stomach | Rumen'}; dtype=object; review_status=unclear
  - cultivar: description=UNCLEAR - human review required; unit=n/a; value labels/range=free text / categorical (14 unique values); dtype=object; review_status=unclear
  - age: description=UNCLEAR - human review required; unit=n/a; value labels/range=free text / categorical (16 unique values); dtype=object; review_status=unclear
  - biomaterial_provider: description=UNCLEAR - human review required; unit=n/a; value labels/range={'Carolina Biological Supply Company', 'Hui-Ling Liao', 'Mark Coggeshall, University of Missouri', 'Prof. Dr. med Georg Bornkamm, HelmholtzZentrum Munich', 'Sugar Research Australia'}; dtype=object; review_status=unclear
  - collected_by: description=UNCLEAR - human review required; unit=n/a; value labels/range={'Dr Brett Williams', 'Gretchen A Gerrish, James G Morin, Nicholai M Hensley, Trevor J Rivers, Todd H Oakley, Emily A Ellis', 'NOTRE DAME', 'QAAFI, University of Queensland'}; dtype=object; review_status=unclear
  - biological_replicate: description=UNCLEAR - human review required; unit=UNCLEAR - human review required; value labels/range=[1, 3]; dtype=float64; review_status=unclear
  - replicate: description=UNCLEAR - human review required; unit=n/a; value labels/range=free text / categorical (12 unique values); dtype=object; review_status=unclear
  - samp_mat_process: description=UNCLEAR - human review required; unit=n/a; value labels/range={'Total RNA was extracted from each individual using Trizol-Chloroform extractions', 'none'}; dtype=object; review_status=unclear
  - samp_size: description=UNCLEAR - human review required; unit=n/a; value labels/range={'45 individuals'}; dtype=object; review_status=unclear
  - ref_biomaterial: description=UNCLEAR - human review required; unit=UNCLEAR - human review required; value labels/range=UNCLEAR - human review required; dtype=float64; review_status=unclear
  - rel_to_oxygen: description=UNCLEAR - human review required; unit=UNCLEAR - human review required; value labels/range=UNCLEAR - human review required; dtype=float64; review_status=unclear
  - samp_collect_device: description=UNCLEAR - human review required; unit=UNCLEAR - human review required; value labels/range=UNCLEAR - human review required; dtype=float64; review_status=unclear
  - ecotype: description=UNCLEAR - human review required; unit=n/a; value labels/range={'B73', 'mising', 'not applicable', 'not collected'}; dtype=object; review_status=unclear
  - identified_by: description=UNCLEAR - human review required; unit=n/a; value labels/range={'Gretchen A Gerrish, James G Morin, Nicholai M Hensley, Trevor J Rivers, Todd H Oakley, Emily A Ellis'}; dtype=object; review_status=unclear
  - sex: description=UNCLEAR - human review required; unit=n/a; value labels/range={'female', 'male', 'not collected', 'not determined'}; dtype=object; review_status=unclear
  - dev_stage: description=UNCLEAR - human review required; unit=n/a; value labels/range=free text / categorical (14 unique values); dtype=object; review_status=unclear
  - ena-first-public: description=ENA first public release date; unit=n/a; value labels/range={'2018-06-11 19:03:00', '2019-02-08 18:04:04', '2019-02-26 18:02:58', '2019-04-19 06:02:37'}; dtype=datetime64[ns]
  - ena-last-update: description=ENA last update date; unit=n/a; value labels/range={'2018-06-05 17:33:43', '2018-06-22 15:33:29', '2018-12-17 15:58:17', '2019-02-21 10:18:26'}; dtype=datetime64[ns]
  - external id: description=UNCLEAR - human review required; unit=n/a; value labels/range={'SAMEA3923215', 'SAMEA4708865', 'SAMEA4753698', 'SAMEA5185551', 'SAMEA5365958'}; dtype=object; review_status=unclear
  - insdc center name: description=UNCLEAR - human review required; unit=n/a; value labels/range={'Center for Biotechnology (CeBiTec) - Bielefeld University', 'EAWAG', 'Faculty of Medicine, Goethe-University Frankfurt', 'SC', 'UNIVERSITY OF EDINBURGH'}; dtype=object; review_status=unclear
  - insdc first public: description=INSDC first public release date; unit=n/a; value labels/range={'2018-03-30 19:02:14', '2018-06-11 19:03:00', '2019-02-08 18:04:04', '2019-02-26 18:02:58', '2019-04-19 06:02:37'}; dtype=datetime64[ns]
  - insdc last update: description=INSDC last update date; unit=n/a; value labels/range={'2016-10-21 11:43:06', '2018-06-05 17:33:43', '2018-06-22 15:33:29', '2018-12-17 15:58:17', '2019-02-21 10:18:26'}; dtype=datetime64[ns]
  - insdc status: description=UNCLEAR - human review required; unit=n/a; value labels/range={'public'}; dtype=object; review_status=unclear
  - submitter id: description=Submitter name; unit=n/a; value labels/range={'20730x94', '6bcb1c40-540d-11e6-9a3d-3c4a9275d6c8', 'E-MTAB-6943:IgG-ACM_1,2,3', 'E-MTAB-7533:Sample 16', 'U121'}; dtype=object; review_status=review_recommended
  - common name: description=Common name; unit=n/a; value labels/range={'human', 'viral metagenome'}; dtype=object
  - sample name: description=Sample name; unit=n/a; value labels/range={'20730x94', '6bcb1c40-540d-11e6-9a3d-3c4a9275d6c8', 'E-MTAB-6943:IgG-ACM_1,2,3', 'E-MTAB-7533:Sample 16', 'U121'}; dtype=object
  - scientific_name: description=Scientific name; unit=n/a; value labels/range={'Homo sapiens', 'soil metagenome', 'viral metagenome'}; dtype=object
  - phenotype: description=UNCLEAR - human review required; unit=n/a; value labels/range={'cyan gray', 'white', 'wild type phenotype'}; dtype=object; review_status=unclear
  - experimental factor: compound: description=UNCLEAR - human review required; unit=n/a; value labels/range={'ryegrass silage'}; dtype=object; review_status=unclear
  - experimental factor: dose: description=UNCLEAR - human review required; unit=UNCLEAR - human review required; value labels/range=[1500, 1500]; dtype=float64; review_status=unclear
  - ebi_url: description=EBI download URL; unit=n/a; value labels/range={'http://ftp.sra.ebi.ac.uk/vol1/run/ERR135/ERR1356733/RU122R2.fq.gz', 'http://ftp.sra.ebi.ac.uk/vol1/run/ERR266/ERR2660556/10_S2_R1_001.fastq.gz', 'http://ftp.sra.ebi.ac.uk/vol1/run/ERR273/ERR2737479/20730_1#94.cram', 'http://ftp.sra.ebi.ac.uk/vol1/run/ERR301/ERR3010915/YS-16_CCGTCC_R2.fastq.gz', 'http://ftp.sra.ebi.ac.uk/vol1/run/ERR317/ERR3179625/20730_1#94.cram'}; dtype=object
  - ebi_free_egress: description=EBI free-egress policy; unit=n/a; value labels/range={'worldwide'}; dtype=object
  - ebi_access_type: description=EBI access type; unit=n/a; value labels/range={'anonymous'}; dtype=object
  - broker name: description=UNCLEAR - human review required; unit=n/a; value labels/range={'ArrayExpress'}; dtype=object; review_status=unclear
  - altitude: description=UNCLEAR - human review required; unit=UNCLEAR - human review required; value labels/range=UNCLEAR - human review required; dtype=float64; review_status=unclear
  - env_broad_scale: description=UNCLEAR - human review required; unit=UNCLEAR - human review required; value labels/range=UNCLEAR - human review required; dtype=float64; review_status=unclear
  - env_local_scale: description=UNCLEAR - human review required; unit=UNCLEAR - human review required; value labels/range=UNCLEAR - human review required; dtype=float64; review_status=unclear
  - env_medium: description=UNCLEAR - human review required; unit=UNCLEAR - human review required; value labels/range=UNCLEAR - human review required; dtype=float64; review_status=unclear
  - season: description=UNCLEAR - human review required; unit=UNCLEAR - human review required; value labels/range=UNCLEAR - human review required; dtype=float64; review_status=unclear
  - type_of_nucleic_acid: description=Record type; unit=UNCLEAR - human review required; value labels/range=UNCLEAR - human review required; dtype=float64; review_status=unclear
  - arrayexpress-species: description=Taxonomic species; unit=n/a; value labels/range={'viral metagenome'}; dtype=object; review_status=review_recommended
  - sample description: description=Sample identifier; unit=n/a; value labels/range={'Total nucleic acid extraction from human feces'}; dtype=object; review_status=review_recommended
  - experimental factor: growth condition: description=UNCLEAR - human review required; unit=n/a; value labels/range={'MCF-7 conditioned media'}; dtype=object; review_status=unclear
  - experimental factor: stimulus: description=UNCLEAR - human review required; unit=n/a; value labels/range={'none'}; dtype=object; review_status=unclear
  - experimental factor: immunoprecipitate: description=UNCLEAR - human review required; unit=n/a; value labels/range={'anti-IgG'}; dtype=object; review_status=unclear
  - cell type: description=Cell type; unit=n/a; value labels/range={'S2 cells', 'macrophage'}; dtype=object
  - developmental stage: description=UNCLEAR - human review required; unit=n/a; value labels/range={'adult'}; dtype=object; review_status=unclear
  - genotype: description=UNCLEAR - human review required; unit=n/a; value labels/range={'PB1+ Pi54', 'Wildtype', 'wild type genotype'}; dtype=object; review_status=unclear
  - growth condition: description=UNCLEAR - human review required; unit=n/a; value labels/range={'MCF-7 conditioned media'}; dtype=object; review_status=unclear
  - indivudal: description=UNCLEAR - human review required; unit=n/a; value labels/range={'pool of donors 1 to 3'}; dtype=object; review_status=unclear
  - organism part: description=UNCLEAR - human review required; unit=n/a; value labels/range={'blood'}; dtype=object; review_status=unclear
  - progenitor cell type: description=Record type; unit=n/a; value labels/range={'peripheral blood mononuclear cell'}; dtype=object; review_status=review_recommended
  - source_name: description=UNCLEAR - human review required; unit=n/a; value labels/range={'Cortical thymic epithelial cells', 'Leaf', 'Medullary thymic epithelial cells', 'S2 cells', 'Skin epithelial cells', 'acute myeloid leukemia', 'stage 1: one week before water deficit', 'whole body'}; dtype=object; review_status=unclear
  - treatment: description=UNCLEAR - human review required; unit=n/a; value labels/range=free text / categorical (79 unique values); dtype=object; review_status=unclear
  - sample_name: description=Sample name; unit=n/a; value labels/range={'Sample_1086', 'Sample_1235'}; dtype=object
  - elev: description=UNCLEAR - human review required; unit=UNCLEAR - human review required; value labels/range=UNCLEAR - human review required; dtype=float64; review_status=unclear
  - env_feature: description=UNCLEAR - human review required; unit=n/a; value labels/range={'Amazon River', 'ENVO:2100002'}; dtype=object; review_status=unclear
  - depth: description=Read depth; unit=n/a; value labels/range={'0.5 m', '10 m', '14 m', '33 m'}; dtype=object; review_status=review_recommended
  - env_material: description=UNCLEAR - human review required; unit=n/a; value labels/range={'ENVO:00002003', 'water'}; dtype=object; review_status=unclear
  - env_biome: description=UNCLEAR - human review required; unit=n/a; value labels/range={'ENVO:00009003', 'river'}; dtype=object; review_status=unclear
  - ena first public: description=ENA first public release date; unit=n/a; value labels/range={'2018-03-30'}; dtype=datetime64[ns]
  - ena last update: description=ENA last update date; unit=n/a; value labels/range={'2016-10-21'}; dtype=datetime64[ns]
  - insdc center alias: description=UNCLEAR - human review required; unit=n/a; value labels/range={'EAWAG'}; dtype=object; review_status=unclear
  - collection date: description=Collection date; unit=ISO-8601; value labels/range=[2014, 2014]; dtype=float64
  - environment (biome): description=UNCLEAR - human review required; unit=n/a; value labels/range={'wastewater'}; dtype=object; review_status=unclear
  - environment (feature): description=UNCLEAR - human review required; unit=n/a; value labels/range={'wastewater treatment plant'}; dtype=object; review_status=unclear
  - environment (material): description=UNCLEAR - human review required; unit=n/a; value labels/range={'activated sludge'}; dtype=object; review_status=unclear
  - geographic location (country and/or sea): description=Geographic location: country and/or sea; unit=n/a; value labels/range={'Switzerland'}; dtype=object
  - geographic location (latitude): description=Sampling latitude; unit=degrees; value labels/range=[-90, 90]; dtype=float64
  - geographic location (longitude): description=Sampling longitude; unit=degrees; value labels/range=[-180, 180]; dtype=float64
  - investigation type: description=Record type; unit=n/a; value labels/range={'metagenome'}; dtype=object; review_status=review_recommended
  - project name: description=UNCLEAR - human review required; unit=n/a; value labels/range={'What biological and operational factors determine the release of antibiotic resistance genes from WWTPs?'}; dtype=object; review_status=unclear
  - sequencing method: description=Sequencing method; unit=n/a; value labels/range={'Illumina HiSeq'}; dtype=object
  - wastewater/sludge environmental package: description=UNCLEAR - human review required; unit=n/a; value labels/range={'wastewater/sludge'}; dtype=object; review_status=unclear
  - cell_line: description=UNCLEAR - human review required; unit=n/a; value labels/range={'P493-4', 'P493-5', 'P493-6'}; dtype=object; review_status=unclear
  - cell_type: description=Cell type; unit=n/a; value labels/range={'B-cell'}; dtype=object
  - nucleic_acid: description=UNCLEAR - human review required; unit=n/a; value labels/range={'rna'}; dtype=object; review_status=unclear
  - sample_timepoint: description=Sample identifier; unit=n/a; value labels/range={'SF05'}; dtype=object; review_status=review_recommended
  - host_subject_id: description=Host subject identifier; unit=unitless; value labels/range=[5.01762e+07, 5.01762e+07]; dtype=float64
  - isolation source: description=Material or environment from which the sample was isolated; unit=UNCLEAR - human review required; value labels/range=UNCLEAR - human review required; dtype=float64; review_status=unclear
  - isol_growth_condt: description=UNCLEAR - human review required; unit=UNCLEAR - human review required; value labels/range=UNCLEAR - human review required; dtype=float64; review_status=unclear
  - propagation: description=UNCLEAR - human review required; unit=UNCLEAR - human review required; value labels/range=UNCLEAR - human review required; dtype=float64; review_status=unclear
  - number_of_cells: description=UNCLEAR - human review required; unit=UNCLEAR - human review required; value labels/range=[500000, 2e+06]; dtype=float64; review_status=unclear
  - target_molecules: description=UNCLEAR - human review required; unit=n/a; value labels/range={'mRNA'}; dtype=object; review_status=unclear
  - site_name: description=UNCLEAR - human review required; unit=n/a; value labels/range={'Macapa North (MCPN)', 'Macapa South (MCPS)', 'Obidos (OB)', 'Tapajos Surface (TAPS)'}; dtype=object; review_status=unclear
  - chlorophyll: description=UNCLEAR - human review required; unit=n/a; value labels/range={'0.87 ug/L', '1.67 ug/L', '1.99 ug/L', '3.23 ug/L', '3.82 ug/L'}; dtype=object; review_status=unclear
  - conduc: description=UNCLEAR - human review required; unit=n/a; value labels/range={'15.9 uS/cm', '55.3 uS/cm', '56.2 uS/cm', '56.4 uS/cm', '61.5 uS/cm'}; dtype=object; review_status=unclear
  - diss_inorg_carb: description=UNCLEAR - human review required; unit=n/a; value labels/range={'155 umol C/kg', '459 umol C/kg', '502 umol C/kg', '507 umol C/kg', '551 umol C/kg'}; dtype=object; review_status=unclear
  - bacterial_count: description=Count; unit=n/a; value labels/range={'3.58E06 cells/mL', '3.63E06 cells/mL', '3.76E06 cells/mL', '3.94E06 cells/mL', '3.96E06 cells/mL'}; dtype=object; review_status=review_recommended
  - samp_volume: description=UNCLEAR - human review required; unit=n/a; value labels/range={'0.65 L', '0.7 L', '0.73 L', '0.76 L', '1.06 L', '1.59 L'}; dtype=object; review_status=unclear
  - min_filter_size: description=UNCLEAR - human review required; unit=n/a; value labels/range={'0.2 um', '2 um', '2.0 Âµm'}; dtype=object; review_status=unclear
  - max_filter_size: description=UNCLEAR - human review required; unit=n/a; value labels/range={'156 Âµm', '2.0 um', '297 um'}; dtype=object; review_status=unclear
  - sequencing_method: description=Sequencing method; unit=n/a; value labels/range={'Illumina PE 150x150', 'Illumina Rapid Run'}; dtype=object
  - sequencing_machine: description=Sequencing machine; unit=n/a; value labels/range={'HiSeq 2500'}; dtype=object
  - std_norm_factor: description=UNCLEAR - human review required; unit=UNCLEAR - human review required; value labels/range=[4217, 43832]; dtype=float64; review_status=unclear
  - specimen_voucher: description=UNCLEAR - human review required; unit=n/a; value labels/range={'TH-08-6182'}; dtype=object; review_status=unclear
  - source_material_id: description=UNCLEAR - human review required; unit=UNCLEAR - human review required; value labels/range=UNCLEAR - human review required; dtype=float64; review_status=unclear
  - sample_number: description=Sample identifier; unit=UNCLEAR - human review required; value labels/range=UNCLEAR - human review required; dtype=float64; review_status=unclear
  - sample_source: description=Sample identifier; unit=UNCLEAR - human review required; value labels/range=UNCLEAR - human review required; dtype=float64; review_status=unclear
  - investigation_type: description=Record type; unit=n/a; value labels/range={'Non-selective Metatranscriptome'}; dtype=object; review_status=review_recommended
  - temp: description=UNCLEAR - human review required; unit=n/a; value labels/range={'AMBIENT'}; dtype=object; review_status=unclear
  - geographic_location: description=Sampling location; unit=n/a; value labels/range={'not collected'}; dtype=object; review_status=review_recommended
  - ena_fastq_http: description=UNCLEAR - human review required; unit=n/a; value labels/range=free text / categorical (131 unique values); dtype=object; review_status=unclear
  - ena_fastq_http_1: description=UNCLEAR - human review required; unit=n/a; value labels/range=free text / categorical (122 unique values); dtype=object; review_status=unclear
  - ena_fastq_http_2: description=UNCLEAR - human review required; unit=n/a; value labels/range=free text / categorical (121 unique values); dtype=object; review_status=unclear
  - ena_fastq_ftp: description=UNCLEAR - human review required; unit=n/a; value labels/range=free text / categorical (131 unique values); dtype=object; review_status=unclear
  - ena_fastq_ftp_1: description=UNCLEAR - human review required; unit=n/a; value labels/range=free text / categorical (122 unique values); dtype=object; review_status=unclear
  - ena_fastq_ftp_2: description=UNCLEAR - human review required; unit=n/a; value labels/range=free text / categorical (121 unique values); dtype=object; review_status=unclear
* Variable dictionary source: S4_253SRR_metadata_variable_dictionary.review_copy.yaml (review_copy)
* Missing data codes: *list code/symbol and definition*
  - 'none': found 20 value(s)
  - 'Unknown': found 1 value(s)
  - 'unknown': found 2 value(s)
  - NaN/blank parsed by pandas: found 27738 value(s)
* Specialized formats or other abbreviations used: datetime columns: public_date, ena-first-public, ena-last-update, insdc first public, insdc last update, ena first public, ena last update

# DATA-SPECIFIC INFORMATION FOR: [Supplementary_Table5.tsv]
*repeat this section for each dataset, folder or file, as appropriate*

* Number of variables: 54
* Number of cases/rows: 510
* Variable List: *list variable name(s), description(s), unit(s) and value labels as appropriate for each*
  - contig_id: description=Contig identifier; unit=n/a; value labels/range=free text / categorical (510 unique values); dtype=object; review_status=review_recommended
  - sequence: description=Nucleotide sequence; unit=n/a; value labels/range=free text / categorical (510 unique values); dtype=object; review_status=review_recommended
  - corresponding_srr: description=UNCLEAR - human review required; unit=n/a; value labels/range=free text / categorical (131 unique values); dtype=object; review_status=unclear
  - assembler: description=Assembler name; unit=n/a; value labels/range={'megahit', 'spades'}; dtype=object; review_status=review_recommended
  - cluster_membership: description=UNCLEAR - human review required; unit=n/a; value labels/range=free text / categorical (56 unique values); dtype=object; review_status=unclear
  - known_or_potentially_novel_tobamovirus: description=UNCLEAR - human review required; unit=unitless; value labels/range={0, 1}; dtype=bool; review_status=unclear
  - contig_length: description=Contig length; unit=bp; value labels/range=[606, 7832]; dtype=int64
  - cluster_representative: description=UNCLEAR - human review required; unit=unitless; value labels/range={0, 1}; dtype=float64; review_status=unclear
  - orf1_complete: description=UNCLEAR - human review required; unit=unitless; value labels/range={0, 1}; dtype=float64; review_status=unclear
  - orf1_partial: description=UNCLEAR - human review required; unit=unitless; value labels/range={0, 1}; dtype=float64; review_status=unclear
  - orf1_length: description=Sequence length; unit=bp; value labels/range=[0, 1135]; dtype=float64; review_status=review_recommended
  - orf1_start: description=Start position; unit=bp; value labels/range=[0, 1]; dtype=float64; review_status=review_recommended
  - orf1_stop: description=Stop position; unit=bp; value labels/range=[0, 1]; dtype=float64; review_status=review_recommended
  - orf1_tree_representative: description=UNCLEAR - human review required; unit=unitless; value labels/range={0, 1}; dtype=float64; review_status=unclear
  - orf2_rdrp_complete: description=UNCLEAR - human review required; unit=unitless; value labels/range={0, 1}; dtype=float64; review_status=unclear
  - orf2_rdrp_partial: description=UNCLEAR - human review required; unit=unitless; value labels/range={0, 1}; dtype=float64; review_status=unclear
  - orf2_rdrp_length: description=Sequence length; unit=bp; value labels/range=[0, 592]; dtype=float64; review_status=review_recommended
  - orf2_start: description=Start position; unit=bp; value labels/range=[0, 1]; dtype=float64; review_status=review_recommended
  - orf2_stop: description=Stop position; unit=bp; value labels/range=[0, 1]; dtype=float64; review_status=review_recommended
  - orf2_tree_representative: description=UNCLEAR - human review required; unit=unitless; value labels/range={0, 1}; dtype=float64; review_status=unclear
  - orf3_mp_complete: description=UNCLEAR - human review required; unit=unitless; value labels/range={0, 1}; dtype=float64; review_status=unclear
  - orf3_mp_partial: description=UNCLEAR - human review required; unit=unitless; value labels/range={0, 1}; dtype=float64; review_status=unclear
  - orf3_mp_length: description=Sequence length; unit=bp; value labels/range=[0, 306]; dtype=float64; review_status=review_recommended
  - orf3_start: description=Start position; unit=bp; value labels/range=[0, 1]; dtype=float64; review_status=review_recommended
  - orf3_stop: description=Stop position; unit=bp; value labels/range=[0, 1]; dtype=float64; review_status=review_recommended
  - orf3_tree_representative: description=UNCLEAR - human review required; unit=unitless; value labels/range={0, 1}; dtype=float64; review_status=unclear
  - orf4_cp_complete: description=UNCLEAR - human review required; unit=unitless; value labels/range={0, 1}; dtype=float64; review_status=unclear
  - orf4_cp_partial: description=UNCLEAR - human review required; unit=unitless; value labels/range={0, 1}; dtype=float64; review_status=unclear
  - orf4_cp_length: description=Sequence length; unit=bp; value labels/range=[0, 186]; dtype=float64; review_status=review_recommended
  - orf4_start: description=Start position; unit=bp; value labels/range=[0, 1]; dtype=float64; review_status=review_recommended
  - orf4_stop: description=Stop position; unit=bp; value labels/range=[0, 1]; dtype=float64; review_status=review_recommended
  - orf4_tree_representative: description=UNCLEAR - human review required; unit=unitless; value labels/range={0, 1}; dtype=float64; review_status=unclear
  - notes: description=Free-text notes; unit=n/a; value labels/range={'missing orf3', 'missing stop codon in orf3', 'missing stop codon in orf4', 'missing stop codons in orf2 and orf3', 'orf1 disrupted by premature stop codon', 'orf2 disrupted by premature stop codon', 'two start codons in orf3', 'two start codons in orf3 and orf4', 'two start codons in orf4'}; dtype=object; review_status=review_recommended
  - ground_truth_category: description=Ground-truth category label; unit=n/a; value labels/range={'mas', 'oth1', 'oth2', 'tob1', 'tob2', 'tob3'}; dtype=object
  - model_prediction: description=Model prediction label; unit=unitless; value labels/range={0, 1}; dtype=int64
  - model_prediction_probabiility: description=Model prediction probability; unit=unitless; value labels/range=[0, 1]; dtype=float64
  - first_diamond_blastx_hit_name: description=UNCLEAR - human review required; unit=n/a; value labels/range=free text / categorical (40 unique values); dtype=object; review_status=unclear
  - first_diamond_blastx_hit_identity_percent: description=Percent identity of the top BLASTX hit; unit=%; value labels/range=[2.8, 99.9]; dtype=float64
  - first_diamond_blastx_hit_alignment_length_aa: description=Amino-acid alignment length of the top BLASTX hit; unit=bp; value labels/range=[53, 2810]; dtype=int64
  - first_blastx_hit_name_accession: description=Accession identifier; unit=n/a; value labels/range=free text / categorical (237 unique values); dtype=object; review_status=review_recommended
  - first_blastx_hit_identity_percent: description=Percent identity of the top BLASTX hit; unit=%; value labels/range=[22, 100]; dtype=float64
  - first_blastx_hit_alignment_length_aa: description=Amino-acid alignment length of the top BLASTX hit; unit=bp; value labels/range=[50, 1889]; dtype=int64
  - first_blastx_hit_e_value: description=E-value of the top BLASTX hit; unit=unitless; value labels/range=[0, 1]; dtype=float64
  - first_blastx_hit_bit_score: description=Bit score of the top BLASTX hit; unit=unitless; value labels/range=[64.3, 3428]; dtype=float64
  - first_blastx_hit_protein: description=UNCLEAR - human review required; unit=n/a; value labels/range=free text / categorical (212 unique values); dtype=object; review_status=unclear
  - study_accession: description=Study accession identifier; unit=n/a; value labels/range=free text / categorical (57 unique values); dtype=object
  - study_title: description=Study title; unit=n/a; value labels/range=free text / categorical (57 unique values); dtype=object
  - organism_name: description=Organism name; unit=n/a; value labels/range=free text / categorical (44 unique values); dtype=object
  - submitter: description=Submitter name; unit=n/a; value labels/range=free text / categorical (46 unique values); dtype=object; review_status=review_recommended
  - country: description=Country of origin; unit=n/a; value labels/range=free text / categorical (44 unique values); dtype=object; review_status=review_recommended
  - collection_date: description=Collection date; unit=ISO-8601; value labels/range=[2012-01-01, 2017-01-01]; dtype=object
  - source_sample_category: description=Broad source sample category; unit=n/a; value labels/range={'Aquatic', 'Host-associated', 'Terrestrial'}; dtype=object
  - source_sample_subcategory: description=Source sample subcategory; unit=n/a; value labels/range={'Animal', 'Animal gut', 'Freshwater', 'Human', 'Human gut', 'Other', 'Plant', 'Soil', 'Wastewater'}; dtype=object
  - genbank_accession_number: description=GenBank accession number; unit=n/a; value labels/range=free text / categorical (63 unique values); dtype=object
* Variable dictionary source: S5_510discovered_contigs_variable_dictionary.review_copy.yaml (review_copy)
* Missing data codes: *list code/symbol and definition*
  - NaN/blank parsed by pandas: found 11000 value(s)
* Specialized formats or other abbreviations used: none specified

# DATA-SPECIFIC INFORMATION FOR: [Supplementary_Table6.tsv]
*repeat this section for each dataset, folder or file, as appropriate*

* Number of variables: 16
* Number of cases/rows: 125
* Variable List: *list variable name(s), description(s), unit(s) and value labels as appropriate for each*
  - type: description=Record type; unit=n/a; value labels/range={'outgroup', 'tobamo'}; dtype=object; review_status=review_recommended
  - genus: description=Taxonomic genus; unit=n/a; value labels/range={'Furovirus', 'Goravirus', 'Hordeivirus', 'Pecluvirus', 'Pomovirus', 'Tobamovirus', 'Tobravirus'}; dtype=object; review_status=review_recommended
  - species: description=Taxonomic species; unit=n/a; value labels/range=free text / categorical (60 unique values); dtype=object; review_status=review_recommended
  - virus_name: description=Virus or viral taxon; unit=n/a; value labels/range=free text / categorical (88 unique values); dtype=object; review_status=review_recommended
  - accession: description=Accession identifier; unit=n/a; value labels/range=free text / categorical (125 unique values); dtype=object; review_status=review_recommended
  - RNA: description=Reference genome segment or RNA component; unit=n/a; value labels/range={'RNA1', 'RNA2', 'RNA3'}; dtype=object
  - training: description=Training-set membership flag; unit=unitless; value labels/range={0, 1}; dtype=int64; review_status=review_recommended
  - length: description=Sequence length; unit=bp; value labels/range=[1799, 7226]; dtype=int64; review_status=review_recommended
  - orf1_ref: description=Reference identifier for ORF1; unit=n/a; value labels/range=free text / categorical (88 unique values); dtype=object
  - orf2_ref: description=Reference identifier for ORF2; unit=n/a; value labels/range=free text / categorical (88 unique values); dtype=object
  - cp_ref: description=Reference identifier for coat protein; unit=n/a; value labels/range=free text / categorical (87 unique values); dtype=object
  - record_id: description=UNCLEAR - human review required; unit=n/a; value labels/range=free text / categorical (125 unique values); dtype=object; review_status=unclear
  - orf1_record_id: description=Record identifier for ORF1; unit=n/a; value labels/range=free text / categorical (88 unique values); dtype=object
  - orf2_record_id: description=Record identifier for ORF2; unit=n/a; value labels/range=free text / categorical (88 unique values); dtype=object
  - cp_record_id: description=Record identifier for coat protein; unit=n/a; value labels/range=free text / categorical (87 unique values); dtype=object
  - sampling_prob: description=Sampling probability used during reference set construction; unit=unitless; value labels/range=[0, 1]; dtype=float64
* Variable dictionary source: S6_reference_database_variable_dictionary.review_copy.yaml (review_copy)
* Missing data codes: *list code/symbol and definition*
  - NaN/blank parsed by pandas: found 224 value(s)
* Specialized formats or other abbreviations used: none specified
