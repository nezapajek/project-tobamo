# DATA-SPECIFIC INFORMATION FOR: [Supplementary_Table3.tsv]
*repeat this section for each dataset, folder or file, as appropriate*

* Number of variables: 15
* Number of cases/rows: 136964
* Variable List: *list variable name(s), description(s), unit(s) and value labels as appropriate for each*
  - run_id: description=Sequencing run accession; unit=n/a; value labels/range=free text / categorical (6870 unique values); dtype=object
  - palm_id: description=Palmprint cluster identifier; unit=n/a; value labels/range=free text / categorical (1703 unique values); dtype=object
  - coverage: description=Read coverage; unit=unitless; value labels/range=[0.45, 40683]; dtype=float64; review_status=review completed
  - sOTU: description=sub-operational taxonomic unit identifier; unit=n/a; value labels/range=free text / categorical (1165 unique values); dtype=object; review_status=manual_review_completed
  - qseqid: description=Query sequence identifier; unit=n/a; value labels/range=free text / categorical (35 unique values); dtype=object
  - pident: description=Percent identity of the alignment; unit=unitless; value labels/range=[39, 100]; dtype=float64
  - evalue: description=Alignment expectation value; unit=unitless; value labels/range=[0, 1]; dtype=int64
  - sra_sequence: description=SRA hit sequence; unit=n/a; value labels/range=free text / categorical (1748 unique values); dtype=object; review_status=review completed
  - biosample_id: description=BioSample accession identifier; unit=n/a; value labels/range=free text / categorical (5931 unique values); dtype=object
  - scientific_name: description=Scientific name; unit=n/a; value labels/range=free text / categorical (1009 unique values); dtype=object
  - nickname: description=Nickname; unit=n/a; value labels/range=free text / categorical (1165 unique values); dtype=object; review_status=review completed
  - date: description=Collection or processing date; unit=ISO-8601; value labels/range=[2012-01-21, 2021-01-11]; dtype=datetime64[ns]; review_status=review completed
  - lng: description=Longitude; unit=degrees; value labels/range=[-180, 180]; dtype=float64; review_status=review completed
  - lat: description=Latitude; unit=degrees; value labels/range=[-90, 90]; dtype=float64; review_status=review completed
  - file: description=Source file name or path; unit=n/a; value labels/range=free text / categorical (35 unique values); dtype=object; review_status=review completed
* Variable dictionary source: S3_serratus_variable_dictionary.reviewed.yaml (review_copy)
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
  - library_name: description=library name; unit=n/a; value labels/range=free text / categorical (198 unique values); dtype=object; review_status=review completed
  - library_strategy: description=Library strategy; unit=n/a; value labels/range={'RNA-Seq', 'WGA', 'WGS'}; dtype=object
  - library_source: description=Library source; unit=n/a; value labels/range={'GENOMIC', 'METAGENOMIC', 'METATRANSCRIPTOMIC', 'TRANSCRIPTOMIC', 'TRANSCRIPTOMIC SINGLE CELL'}; dtype=object
  - library_selection: description=Library selection method; unit=n/a; value labels/range=free text / categorical (11 unique values); dtype=object
  - library_layout: description=Library layout; unit=n/a; value labels/range={'PAIRED', 'SINGLE'}; dtype=object
  - sample_accession: description=Sample accession identifier; unit=n/a; value labels/range=free text / categorical (247 unique values); dtype=object
  - sample_title: description=Sample title; unit=n/a; value labels/range=n/a; dtype=float64; review_status=review completed
  - biosample: description=BioSample accession identifier; unit=n/a; value labels/range=free text / categorical (247 unique values); dtype=object; review_status=review_recommended
  - bioproject: description=BioProject accession identifier; unit=n/a; value labels/range=free text / categorical (69 unique values); dtype=object
  - instrument: description=Instrument; unit=n/a; value labels/range=free text / categorical (11 unique values); dtype=object; review_status=review completed
  - instrument_model: description=Instrument model; unit=n/a; value labels/range=free text / categorical (11 unique values); dtype=object
  - instrument_model_desc: description=Instrument model description; unit=n/a; value labels/range={'ILLUMINA'}; dtype=object
  - total_spots: description=Total number of spots in the dataset; unit=n/a; value labels/range=[102819, 5.79302e+08]; dtype=int64; review_status=review completed
  - total_size: description=Total size of the dataset; unit=n/a; value labels/range=[3.29977e+07, 4.05363e+10]; dtype=int64; review_status=review completed
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
  - public_url: description=public url; unit=n/a; value labels/range=free text / categorical (253 unique values); dtype=object; review_status=review completed
  - ncbi_url: description=NCBI download URL; unit=n/a; value labels/range=free text / categorical (253 unique values); dtype=object
  - ncbi_free_egress: description=NCBI free-egress policy; unit=n/a; value labels/range={'worldwide'}; dtype=object
  - ncbi_access_type: description=NCBI access type; unit=n/a; value labels/range={'anonymous'}; dtype=object
  - gcp_url: description=Google Cloud Storage URL; unit=n/a; value labels/range=free text / categorical (253 unique values); dtype=object
  - gcp_free_egress: description=Google Cloud free-egress region or policy; unit=n/a; value labels/range={'gs.us-east1'}; dtype=object
  - gcp_access_type: description=Google Cloud access type; unit=n/a; value labels/range={'gcp identity'}; dtype=object
  - experiment_alias: description=experiment alias; unit=n/a; value labels/range=free text / categorical (19 unique values); dtype=object; review_status=review completed
  - strain: description=microbial or eukaryotic strain name; unit=n/a; value labels/range={'C57BL/6', 'MdSGHV', 'fenpropathrin-resistant', 'fenpropathrin-susceptible'}; dtype=object
  - isolate: description=identification or description of the specific individual from which this sample was obtained; unit=n/a; value labels/range=free text / categorical (33 unique values); dtype=object
  - host: description=The natural (as opposed to laboratory) host to the organism from which the sample was obtained. Use the full taxonomic name, eg, "Homo sapiens".; unit=n/a; value labels/range=free text / categorical (16 unique values); dtype=object
  - isolation_source: description=Describes the physical, environmental and/or local geographical source of the biological sample from which the sample was derived.; unit=n/a; value labels/range=free text / categorical (20 unique values); dtype=object
  - collection_date: description=the date on which the sample was collected; date/time ranges are supported by providing two dates from among the supported value formats, delimited by a forward-slash character; collection times are supported by adding "T", then the hour and minute after the date, and must be in Coordinated Universal Time (UTC), otherwise known as "Zulu Time" (Z); supported formats include "DD-Mmm-YYYY", "Mmm-YYYY", "YYYY" or ISO 8601 standard "YYYY-mm-dd", "YYYY-mm", "YYYY-mm-ddThh:mm:ss"; e.g., 30-Oct-1990, Oct-1990, 1990, 1990-10-30, 1990-10,  21-Oct-1952/15-Feb-1953, 2015-10-11T17:53:03Z; valid non-ISO dates will be automatically transformed to ISO format; unit=ISO-8601; value labels/range=[2008-07-06, 2019-01-01]; dtype=object
  - geo_loc_name: description=Geographical origin of the sample; use the appropriate name from this list https://www.insdc.org/submitting-standards/geo_loc_name-qualifier-vocabulary/. Use a colon to separate the country or ocean from more detailed information about the location, eg "Canada: Vancouver" or "Germany: halfway down Zugspitze, Alps"; unit=n/a; value labels/range=free text / categorical (48 unique values); dtype=object
  - sample_type: description=Sample type, such as cell culture, mixed culture, tissue sample, whole organism, single cell, metagenomic assembly; unit=n/a; value labels/range={'2786_HFKKJALXX_L5', 'cells', 'single cell'}; dtype=object
  - host_tissue_sampled: description=name of body site where the sample was obtained from, such as a specific organ or tissue, e.g., tongue, lung. For foundational model of anatomy ontology (fma) (v 4.11.0) or Uber-anatomy ontology (UBERON) (v releases/2014-06-15) terms, please see http://purl.bioontology.org/ontology/FMA or http://purl.bioontology.org/ontology/UBERON; unit=n/a; value labels/range={'Gut'}; dtype=object
  - biosamplemodel: description=biosample model or type; unit=n/a; value labels/range={'Human', 'Invertebrate', 'MIGS/MIMS/MIMARKS.human-gut', 'MIGS/MIMS/MIMARKS.water', 'Metagenome or environmental', 'Microbe, viral or environmental', 'Model organism or animal', 'Plant', 'Viral'}; dtype=object; review_status=review completed
  - lat_lon: description=The geographical coordinates of the location where the sample was collected. Specify as degrees latitude and longitude in format "d[d.dddd] N|S d[dd.dddd] W|E", eg, 38.98 N 77.11 W; unit=degrees; value labels/range=formatted latitude/longitude pair; dtype=object
  - breed: description=breed name - chiefly used in domesticated animals or plants; unit=n/a; value labels/range={'Missing', 'bashibai sheep', 'citrus', 'not applicable', 'not collected'}; dtype=object
  - tissue: description=Type of tissue the sample was taken from.; unit=n/a; value labels/range=free text / categorical (32 unique values); dtype=object
  - gold ecosystem classification: description=gold ecosystem classification; unit=n/a; value labels/range={'Engineered | Bioreactor | Anaerobic | Digestate | Unclassified', 'Environmental | Terrestrial | Soil | Floodplain | Unclassified', 'Environmental | Terrestrial | Soil | Riparian zone | Unclassified', 'Host-associated | Mammals | Digestive system | Stomach | Rumen'}; dtype=object; review_status=review completed
  - cultivar: description=cultivar name - cultivated variety of plant; unit=n/a; value labels/range=free text / categorical (14 unique values); dtype=object
  - age: description=age at the time of sampling; relevant scale depends on species and study, e.g. could be seconds for amoebae or centuries for trees; unit=n/a; value labels/range=free text / categorical (16 unique values); dtype=object
  - biomaterial_provider: description=name and address of the lab or PI, or a culture collection identifier; unit=n/a; value labels/range={'Carolina Biological Supply Company', 'Hui-Ling Liao', 'Mark Coggeshall, University of Missouri', 'Prof. Dr. med Georg Bornkamm, HelmholtzZentrum Munich', 'Sugar Research Australia'}; dtype=object
  - collected_by: description=Name of persons or institute who collected the sample; unit=n/a; value labels/range={'Dr Brett Williams', 'Gretchen A Gerrish, James G Morin, Nicholai M Hensley, Trevor J Rivers, Todd H Oakley, Emily A Ellis', 'NOTRE DAME', 'QAAFI, University of Queensland'}; dtype=object
  - biological_replicate: description=biological replicate number; a biological replicate is a sample taken from a different individual or experimental unit than other samples in the dataset, but which is treated identically to other samples in the dataset; biological replicates are used to capture biological variation and to support statistical inference about the population from which the samples were drawn; unit=n/a; value labels/range=[1, 3]; dtype=float64; review_status=review completed
  - replicate: description=replicate; unit=n/a; value labels/range=free text / categorical (12 unique values); dtype=object; review_status=review completed
  - samp_mat_process: description=Processing applied to the sample during or after isolation; unit=n/a; value labels/range={'Total RNA was extracted from each individual using Trizol-Chloroform extractions', 'none'}; dtype=object
  - samp_size: description=Amount or size of sample (volume, mass or area) that was collected; unit=n/a; value labels/range={'45 individuals'}; dtype=object
  - ref_biomaterial: description=Primary publication or genome report; unit=n/a; value labels/range=n/a; dtype=float64; review_status=review completed
  - rel_to_oxygen: description=Is this organism an aerobe, anaerobe? Please note that aerobic and anaerobic are valid descriptors for microbial environments, eg, aerobe, anaerobe, facultative, microaerophilic, microanaerobe, obligate aerobe, obligate anaerobe, missing, not applicable, not collected, not provided, restricted access; unit=n/a; value labels/range=n/a; dtype=float64; review_status=review completed
  - samp_collect_device: description=Method or device employed for collecting sample; unit=n/a; value labels/range=n/a; dtype=float64; review_status=review completed
  - ecotype: description=a population within a given species displaying genetically based, phenotypic traits that reflect adaptation to a local habitat, e.g., Columbia; unit=n/a; value labels/range={'B73', 'mising', 'not applicable', 'not collected'}; dtype=object
  - identified_by: description=name of the taxonomist who identified the specimen; unit=n/a; value labels/range={'Gretchen A Gerrish, James G Morin, Nicholai M Hensley, Trevor J Rivers, Todd H Oakley, Emily A Ellis'}; dtype=object
  - sex: description=physical sex of sampled organism; unit=n/a; value labels/range={'female', 'male', 'not collected', 'not determined'}; dtype=object
  - dev_stage: description=Developmental stage at the time of sampling.; unit=n/a; value labels/range=free text / categorical (14 unique values); dtype=object
  - ena-first-public: description=ENA first public release date; unit=n/a; value labels/range={'2018-06-11 19:03:00', '2019-02-08 18:04:04', '2019-02-26 18:02:58', '2019-04-19 06:02:37'}; dtype=datetime64[ns]
  - ena-last-update: description=ENA last update date; unit=n/a; value labels/range={'2018-06-05 17:33:43', '2018-06-22 15:33:29', '2018-12-17 15:58:17', '2019-02-21 10:18:26'}; dtype=datetime64[ns]
  - external id: description=external id; unit=n/a; value labels/range={'SAMEA3923215', 'SAMEA4708865', 'SAMEA4753698', 'SAMEA5185551', 'SAMEA5365958'}; dtype=object; review_status=review completed
  - insdc center name: description=INSDC center name; unit=n/a; value labels/range={'Center for Biotechnology (CeBiTec) - Bielefeld University', 'EAWAG', 'Faculty of Medicine, Goethe-University Frankfurt', 'SC', 'UNIVERSITY OF EDINBURGH'}; dtype=object; review_status=review completed
  - insdc first public: description=INSDC first public release date; unit=n/a; value labels/range={'2018-03-30 19:02:14', '2018-06-11 19:03:00', '2019-02-08 18:04:04', '2019-02-26 18:02:58', '2019-04-19 06:02:37'}; dtype=datetime64[ns]
  - insdc last update: description=INSDC last update date; unit=n/a; value labels/range={'2016-10-21 11:43:06', '2018-06-05 17:33:43', '2018-06-22 15:33:29', '2018-12-17 15:58:17', '2019-02-21 10:18:26'}; dtype=datetime64[ns]
  - insdc status: description=INSDC status; unit=n/a; value labels/range={'public'}; dtype=object; review_status=review completed
  - submitter id: description=Submitter name; unit=n/a; value labels/range={'20730x94', '6bcb1c40-540d-11e6-9a3d-3c4a9275d6c8', 'E-MTAB-6943:IgG-ACM_1,2,3', 'E-MTAB-7533:Sample 16', 'U121'}; dtype=object; review_status=review_recommended
  - common name: description=Common name; unit=n/a; value labels/range={'human', 'viral metagenome'}; dtype=object
  - sample name: description=sample name in source database; unit=n/a; value labels/range={'20730x94', '6bcb1c40-540d-11e6-9a3d-3c4a9275d6c8', 'E-MTAB-6943:IgG-ACM_1,2,3', 'E-MTAB-7533:Sample 16', 'U121'}; dtype=object
  - scientific_name: description=Scientific name; unit=n/a; value labels/range={'Homo sapiens', 'soil metagenome', 'viral metagenome'}; dtype=object
  - phenotype: description=Phenotype of sampled organism. For Phenotypic quality Ontology (PATO) (v1.269) terms, please see http://bioportal.bioontology.org/visualize/44601; unit=n/a; value labels/range={'cyan gray', 'white', 'wild type phenotype'}; dtype=object
  - experimental factor: compound: description=experimental factor - compound; unit=n/a; value labels/range={'ryegrass silage'}; dtype=object; review_status=review completed
  - experimental factor: dose: description=experimental factor - dose; unit=n\a; value labels/range=[1500, 1500]; dtype=float64; review_status=review completed
  - ebi_url: description=EBI download URL; unit=n/a; value labels/range={'http://ftp.sra.ebi.ac.uk/vol1/run/ERR135/ERR1356733/RU122R2.fq.gz', 'http://ftp.sra.ebi.ac.uk/vol1/run/ERR266/ERR2660556/10_S2_R1_001.fastq.gz', 'http://ftp.sra.ebi.ac.uk/vol1/run/ERR273/ERR2737479/20730_1#94.cram', 'http://ftp.sra.ebi.ac.uk/vol1/run/ERR301/ERR3010915/YS-16_CCGTCC_R2.fastq.gz', 'http://ftp.sra.ebi.ac.uk/vol1/run/ERR317/ERR3179625/20730_1#94.cram'}; dtype=object
  - ebi_free_egress: description=EBI free-egress policy; unit=n/a; value labels/range={'worldwide'}; dtype=object
  - ebi_access_type: description=EBI access type; unit=n/a; value labels/range={'anonymous'}; dtype=object
  - broker name: description=broker name; unit=n/a; value labels/range={'ArrayExpress'}; dtype=object; review_status=review completed
  - altitude: description=The altitude of the sample is the vertical distance between Earth's surface above Sea Level and the sampled position in the air.; unit=n/a; value labels/range=n/a; dtype=float64; review_status=review completed
  - env_broad_scale: description=Add terms that identify the major environment type(s) where your sample was collected. Recommend subclasses of biome [ENVO:00000428]. Multiple terms can be separated by one or more pipes e.g.:  mangrove biome [ENVO:01000181]|estuarine biome [ENVO:01000020]
; unit=unitless; value labels/range=environment (broad scale); dtype=float64; review_status=review completed
  - env_local_scale: description=Add terms that identify environmental entities having causal influences upon the entity at time of sampling, multiple terms can be separated by pipes, e.g.:  shoreline [ENVO:00000486]|intertidal zone [ENVO:00000316]
; unit=unitless; value labels/range=environment (local scale); dtype=float64; review_status=review completed
  - env_medium: description=Add terms that identify the material displaced by the entity at time of sampling. Recommend subclasses of environmental material [ENVO:00010483]. Multiple terms can be separated by pipes e.g.: estuarine water [ENVO:01000301]|estuarine mud [ENVO:00002160]
; unit=unitless; value labels/range=environment (medium scale); dtype=float64; review_status=review completed
  - season: description=The season when sampling occurred. Any of the four periods into which the year is divided by the equinoxes and solstices. This field accepts terms listed under season (http://purl.obolibrary.org/obo/NCIT_C94729); unit=n/a; value labels/range=n/a; dtype=float64; review_status=review completed
  - type_of_nucleic_acid: description=Record type; unit=n/a; value labels/range=n/a; dtype=float64; review_status=review completed
  - arrayexpress-species: description=Taxonomic species; unit=n/a; value labels/range={'viral metagenome'}; dtype=object; review_status=review_recommended
  - sample description: description=Sample identifier; unit=n/a; value labels/range={'Total nucleic acid extraction from human feces'}; dtype=object; review_status=review_recommended
  - experimental factor: growth condition: description=experimental factor - growth condition; unit=n/a; value labels/range={'MCF-7 conditioned media'}; dtype=object; review_status=review completed
  - experimental factor: stimulus: description=experimental factor - stimulus; unit=n/a; value labels/range={'none'}; dtype=object; review_status=review completed
  - experimental factor: immunoprecipitate: description=experimental factor - immunoprecipitate; unit=n/a; value labels/range={'anti-IgG'}; dtype=object; review_status=review completed
  - cell type: description=Type of cell of the sample or from which the sample was obtained.; unit=n/a; value labels/range={'S2 cells', 'macrophage'}; dtype=object
  - developmental stage: description=Developmental stage at the time of sampling.; unit=n/a; value labels/range={'adult'}; dtype=object
  - genotype: description=observed genotype; unit=n/a; value labels/range={'PB1+ Pi54', 'Wildtype', 'wild type genotype'}; dtype=object
  - growth condition: description=growth condition; unit=n/a; value labels/range={'MCF-7 conditioned media'}; dtype=object; review_status=review completed
  - indivudal: description=individual; unit=n/a; value labels/range={'pool of donors 1 to 3'}; dtype=object; review_status=review completed
  - organism part: description=Type of tissue the sample was taken from.; unit=n/a; value labels/range={'blood'}; dtype=object
  - progenitor cell type: description=Record type; unit=n/a; value labels/range={'peripheral blood mononuclear cell'}; dtype=object; review_status=review_recommended
  - source_name: description=Sample source name or description; unit=n/a; value labels/range={'Cortical thymic epithelial cells', 'Leaf', 'Medullary thymic epithelial cells', 'S2 cells', 'Skin epithelial cells', 'acute myeloid leukemia', 'stage 1: one week before water deficit', 'whole body'}; dtype=object
  - treatment: description=treatment applied to the sample, such as chemical treatment, mechanical wounding, rehydration, etc. Please note that allowed values include any treatment applied to the sample, including but not limited to chemical treatment, mechanical wounding, rehydration, etc.; unit=n/a; value labels/range=free text / categorical (79 unique values); dtype=object; review_status=review completed
  - sample_name: description=sample name in source database; unit=n/a; value labels/range={'Sample_1086', 'Sample_1235'}; dtype=object
  - elev: description=The elevation of the sampling site as measured by the vertical distance from mean sea level.; unit=n/a; value labels/range=elevation; dtype=float64; review_status=review completed
  - env_feature: description=Add terms that identify environmental entities having causal influences upon the entity at time of sampling, multiple terms can be separated by pipes, e.g.:  shoreline [ENVO:00000486]|intertidal zone [ENVO:00000316]
; unit=n/a; value labels/range={'Amazon River', 'ENVO:2100002'}; dtype=object
  - depth: description=Depth is defined as the vertical distance below surface, e.g. for sediment or soil samples depth is measured from sediment or soil surface, respectivly. Depth can be reported as an interval for subsurface samples.; unit=m in text; value labels/range={'0.5 m', '10 m', '14 m', '33 m'}; dtype=object
  - env_material: description=Add terms that identify the material displaced by the entity at time of sampling. Recommend subclasses of environmental material [ENVO:00010483]. Multiple terms can be separated by pipes e.g.: estuarine water [ENVO:01000301]|estuarine mud [ENVO:00002160]
; unit=n/a; value labels/range={'ENVO:00002003', 'water'}; dtype=object
  - env_biome: description=Add terms that identify the major environment type(s) where your sample was collected. Recommend subclasses of biome [ENVO:00000428]. Multiple terms can be separated by one or more pipes e.g.:  mangrove biome [ENVO:01000181]|estuarine biome [ENVO:01000020]
; unit=n/a; value labels/range={'ENVO:00009003', 'river'}; dtype=object
  - ena first public: description=ENA first public release date; unit=n/a; value labels/range={'2018-03-30'}; dtype=datetime64[ns]
  - ena last update: description=ENA last update date; unit=n/a; value labels/range={'2016-10-21'}; dtype=datetime64[ns]
  - insdc center alias: description=INSDC center alias; unit=n/a; value labels/range={'EAWAG'}; dtype=object; review_status=review completed
  - collection date: description=the date on which the sample was collected; date/time ranges are supported by providing two dates from among the supported value formats, delimited by a forward-slash character; collection times are supported by adding "T", then the hour and minute after the date, and must be in Coordinated Universal Time (UTC), otherwise known as "Zulu Time" (Z); supported formats include "DD-Mmm-YYYY", "Mmm-YYYY", "YYYY" or ISO 8601 standard "YYYY-mm-dd", "YYYY-mm", "YYYY-mm-ddThh:mm:ss"; e.g., 30-Oct-1990, Oct-1990, 1990, 1990-10-30, 1990-10,  21-Oct-1952/15-Feb-1953, 2015-10-11T17:53:03Z; valid non-ISO dates will be automatically transformed to ISO format; unit=ISO-8601; value labels/range=[2014, 2014]; dtype=float64
  - environment (biome): description=Add terms that identify the major environment type(s) where your sample was collected. Recommend subclasses of biome [ENVO:00000428]. Multiple terms can be separated by one or more pipes e.g.:  mangrove biome [ENVO:01000181]|estuarine biome [ENVO:01000020]
; unit=n/a; value labels/range={'wastewater'}; dtype=object
  - environment (feature): description=Add terms that identify environmental entities having causal influences upon the entity at time of sampling, multiple terms can be separated by pipes, e.g.:  shoreline [ENVO:00000486]|intertidal zone [ENVO:00000316]
; unit=n/a; value labels/range={'wastewater treatment plant'}; dtype=object
  - environment (material): description=Add terms that identify the material displaced by the entity at time of sampling. Recommend subclasses of environmental material [ENVO:00010483]. Multiple terms can be separated by pipes e.g.: estuarine water [ENVO:01000301]|estuarine mud [ENVO:00002160]
; unit=n/a; value labels/range={'activated sludge'}; dtype=object
  - geographic location (country and/or sea): description=Geographical origin of the sample; use the appropriate name from this list https://www.insdc.org/submitting-standards/geo_loc_name-qualifier-vocabulary/. Use a colon to separate the country or ocean from more detailed information about the location, eg "Canada: Vancouver" or "Germany: halfway down Zugspitze, Alps"; unit=n/a; value labels/range={'Switzerland'}; dtype=object
  - geographic location (latitude): description=Sampling latitude; unit=degrees; value labels/range=[-90, 90]; dtype=float64
  - geographic location (longitude): description=Sampling longitude; unit=degrees; value labels/range=[-180, 180]; dtype=float64
  - investigation type: description=Nucleic Acid Sequence Report is the root element of all MIGS/MIMS compliant reports as standardized by Genomic Standards Consortium. This field is either eukaryote,bacteria,virus,plasmid,organelle, metagenome, miens-survey or miens-culture; unit=n/a; value labels/range={'metagenome'}; dtype=object
  - project name: description=A concise name that describes the overall project, for example "Analysis of sequences collected from Antarctic soil"; unit=n/a; value labels/range={'What biological and operational factors determine the release of antibiotic resistance genes from WWTPs?'}; dtype=object
  - sequencing method: description=Sequencing method; unit=n/a; value labels/range={'Illumina HiSeq'}; dtype=object
  - wastewater/sludge environmental package: description=wastewater/sludge environmental package; unit=n/a; value labels/range={'wastewater/sludge'}; dtype=object; review_status=review completed
  - cell_line: description=Name of the cell line.; unit=n/a; value labels/range={'P493-4', 'P493-5', 'P493-6'}; dtype=object
  - cell_type: description=Type of cell of the sample or from which the sample was obtained.; unit=n/a; value labels/range={'B-cell'}; dtype=object
  - nucleic_acid: description=nucleic acid type; unit=n/a; value labels/range={'rna'}; dtype=object; review_status=review completed
  - sample_timepoint: description=Sample identifier; unit=n/a; value labels/range={'SF05'}; dtype=object; review_status=review_recommended
  - host_subject_id: description=a unique identifier by which each subject can be referred to, de-identified, e.g. #131; unit=n/a; value labels/range=[5.01762e+07, 5.01762e+07]; dtype=float64
  - isolation source: description=Describes the physical, environmental and/or local geographical source of the biological sample from which the sample was derived.; unit=n/a; value labels/range=n\a; dtype=float64; review_status=review completed
  - isol_growth_condt: description=PMID or url for isolation and growth condition specifications; unit=unitless; value labels/range=n/a; dtype=float64; review_status=review completed
  - propagation: description=phage: lytic/lysogenic/temperate/obligately lytic; plasmid: incompatibility group; eukaryote: asexual/sexual; unit=n/a; value labels/range=n/a; dtype=float64; review_status=review completed
  - number_of_cells: description=number of cells; unit=unitless; value labels/range=[500000, 2e+06]; dtype=float64; review_status=review completed
  - target_molecules: description=target molecules; unit=n/a; value labels/range={'mRNA'}; dtype=object; review_status=review completed
  - site_name: description=site name; unit=n/a; value labels/range={'Macapa North (MCPN)', 'Macapa South (MCPS)', 'Obidos (OB)', 'Tapajos Surface (TAPS)'}; dtype=object; review_status=review completed
  - chlorophyll: description=concentration of chlorophyll; unit=n/a; value labels/range={'0.87 ug/L', '1.67 ug/L', '1.99 ug/L', '3.23 ug/L', '3.82 ug/L'}; dtype=object
  - conduc: description=electrical conductivity of water; unit=uS/cm in text; value labels/range={'15.9 uS/cm', '55.3 uS/cm', '56.2 uS/cm', '56.4 uS/cm', '61.5 uS/cm'}; dtype=object
  - diss_inorg_carb: description=dissolved inorganic carbon concentration; unit=umol C/kg in text; value labels/range={'155 umol C/kg', '459 umol C/kg', '502 umol C/kg', '507 umol C/kg', '551 umol C/kg'}; dtype=object
  - bacterial_count: description=Count; unit=n/a; value labels/range={'3.58E06 cells/mL', '3.63E06 cells/mL', '3.76E06 cells/mL', '3.94E06 cells/mL', '3.96E06 cells/mL'}; dtype=object; review_status=review_recommended
  - samp_volume: description=sampling volume; unit=L (liters) in text format; value labels/range={'0.65 L', '0.7 L', '0.73 L', '0.76 L', '1.06 L', '1.59 L'}; dtype=object; review_status=review completed
  - min_filter_size: description=minimum filter size; unit=n/a; value labels/range={'0.2 um', '2 um', '2.0 Âµm'}; dtype=object; review_status=review completed
  - max_filter_size: description=maximum filter size; unit=n/a; value labels/range={'156 Âµm', '2.0 um', '297 um'}; dtype=object; review_status=review completed
  - sequencing_method: description=Sequencing method; unit=n/a; value labels/range={'Illumina PE 150x150', 'Illumina Rapid Run'}; dtype=object
  - sequencing_machine: description=Sequencing machine; unit=n/a; value labels/range={'HiSeq 2500'}; dtype=object
  - std_norm_factor: description=std normalization factor; unit=unitless; value labels/range=[4217, 43832]; dtype=float64; review_status=review completed
  - specimen_voucher: description=Identifier for the physical specimen. Use format: "[<institution-code>:[<collection-code>:]]<specimen_id>", eg, "UAM:Mamm:52179". Intended as a reference to the physical specimen that remains after it was analyzed. If the specimen was destroyed in the process of analysis, electronic images (e-vouchers) are an adequate substitute for a physical voucher specimen. Ideally the specimens will be deposited in a curated museum, herbarium, or frozen tissue collection, but often they will remain in a personal or laboratory collection for some time before they are deposited in a curated collection. There are three forms of specimen_voucher qualifiers. If the text of the qualifier includes one or more colons it is a 'structured voucher'. Structured vouchers include institution-codes (and optional collection-codes) taken from a controlled vocabulary maintained by the INSDC that denotes the museum or herbarium collection where the specimen resides, please visit: http://www.insdc.org/controlled-vocabulary-specimenvoucher-qualifier.; unit=n/a; value labels/range={'TH-08-6182'}; dtype=object
  - source_material_id: description=unique identifier assigned to a material sample used for extracting nucleic acids, and subsequent sequencing. The identifier can refer either to the original material collected or to any derived sub-samples.; unit=n/a; value labels/range=n/a; dtype=float64; review_status=review completed
  - sample_number: description=Sample identifier; unit=unitless; value labels/range=n/a; dtype=float64; review_status=review completed
  - sample_source: description=Sample identifier; unit=n/a; value labels/range=n/a; dtype=float64; review_status=review completed
  - investigation_type: description=Nucleic Acid Sequence Report is the root element of all MIGS/MIMS compliant reports as standardized by Genomic Standards Consortium. This field is either eukaryote,bacteria,virus,plasmid,organelle, metagenome, miens-survey or miens-culture; unit=n/a; value labels/range={'Non-selective Metatranscriptome'}; dtype=object
  - temp: description=temperature of the sample at time of sampling; unit=n/a; value labels/range={'AMBIENT'}; dtype=object
  - geographic_location: description=Geographical origin of the sample; use the appropriate name from this list https://www.insdc.org/submitting-standards/geo_loc_name-qualifier-vocabulary/. Use a colon to separate the country or ocean from more detailed information about the location, eg "Canada: Vancouver" or "Germany: halfway down Zugspitze, Alps"; unit=n/a; value labels/range={'not collected'}; dtype=object
  - ena_fastq_http: description=ena fastq http; unit=n/a; value labels/range=free text / categorical (131 unique values); dtype=object; review_status=review completed
  - ena_fastq_http_1: description=ena fastq http 1; unit=n/a; value labels/range=free text / categorical (122 unique values); dtype=object; review_status=review completed
  - ena_fastq_http_2: description=ena fastq http 2; unit=n/a; value labels/range=free text / categorical (121 unique values); dtype=object; review_status=review completed
  - ena_fastq_ftp: description=ena fastq ftp; unit=n/a; value labels/range=free text / categorical (131 unique values); dtype=object; review_status=review completed
  - ena_fastq_ftp_1: description=ena fastq ftp 1; unit=n/a; value labels/range=free text / categorical (122 unique values); dtype=object; review_status=review completed
  - ena_fastq_ftp_2: description=ena fastq ftp 2; unit=n/a; value labels/range=free text / categorical (121 unique values); dtype=object; review_status=review completed
* Variable dictionary source: S4_253SRR_metadata_variable_dictionary.reviewed.yaml (review_copy)
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
  - corresponding_srr: description=corresponding SRA Run accession identifier; unit=n/a; value labels/range=free text / categorical (131 unique values); dtype=object; review_status=review completed
  - assembler: description=Assembler name; unit=n/a; value labels/range={'megahit', 'spades'}; dtype=object; review_status=review_recommended
  - cluster_membership: description=cluster membership identifier; unit=n/a; value labels/range=free text / categorical (56 unique values); dtype=object; review_status=review completed
  - known_or_potentially_novel_tobamovirus: description=flag indicating whether the contig is a known or potentially novel tobamovirus (1) or not (0); unit=unitless; value labels/range={0, 1}; dtype=bool; review_status=review completed
  - contig_length: description=Contig length; unit=bp; value labels/range=[606, 7832]; dtype=int64
  - cluster_representative: description=cluster representative flag; unit=unitless; value labels/range={0, 1}; dtype=float64; review_status=review completed
  - orf1_complete: description=flag indicating whether the contig is complete (1) or not (0); unit=unitless; value labels/range={0, 1}; dtype=float64; review_status=review completed
  - orf1_partial: description=flag indicating whether the contig is partial (1) or not (0); unit=unitless; value labels/range={0, 1}; dtype=float64; review_status=review completed
  - orf1_length: description=Sequence length; unit=bp; value labels/range=[0, 1135]; dtype=float64; review_status=review_recommended
  - orf1_start: description=Start position; unit=bp; value labels/range=[0, 1]; dtype=float64; review_status=review_recommended
  - orf1_stop: description=Stop position; unit=bp; value labels/range=[0, 1]; dtype=float64; review_status=review_recommended
  - orf1_tree_representative: description=flag indicating whether the contig is a tree representative (1) or not (0); unit=unitless; value labels/range={0, 1}; dtype=float64; review_status=review completed
  - orf2_rdrp_complete: description=flag indicating whether the contig is complete (1) or not (0); unit=unitless; value labels/range={0, 1}; dtype=float64; review_status=review completed
  - orf2_rdrp_partial: description=flag indicating whether the contig is partial (1) or not (0); unit=unitless; value labels/range={0, 1}; dtype=float64; review_status=review completed
  - orf2_rdrp_length: description=Sequence length; unit=bp; value labels/range=[0, 592]; dtype=float64; review_status=review_recommended
  - orf2_start: description=Start position; unit=bp; value labels/range=[0, 1]; dtype=float64; review_status=review_recommended
  - orf2_stop: description=Stop position; unit=bp; value labels/range=[0, 1]; dtype=float64; review_status=review_recommended
  - orf2_tree_representative: description=flag indicating whether the contig is a tree representative (1) or not (0); unit=unitless; value labels/range={0, 1}; dtype=float64; review_status=review completed
  - orf3_mp_complete: description=flag indicating whether the contig is complete (1) or not (0); unit=unitless; value labels/range={0, 1}; dtype=float64; review_status=review completed
  - orf3_mp_partial: description=flag indicating whether the contig is partial (1) or not (0); unit=unitless; value labels/range={0, 1}; dtype=float64; review_status=review completed
  - orf3_mp_length: description=Sequence length; unit=bp; value labels/range=[0, 306]; dtype=float64; review_status=review_recommended
  - orf3_start: description=Start position; unit=bp; value labels/range=[0, 1]; dtype=float64; review_status=review_recommended
  - orf3_stop: description=Stop position; unit=bp; value labels/range=[0, 1]; dtype=float64; review_status=review_recommended
  - orf3_tree_representative: description=flag indicating whether the contig is a tree representative (1) or not (0); unit=unitless; value labels/range={0, 1}; dtype=float64; review_status=review completed
  - orf4_cp_complete: description=flag indicating whether the contig is complete (1) or not (0); unit=unitless; value labels/range={0, 1}; dtype=float64; review_status=review completed
  - orf4_cp_partial: description=flag indicating whether the contig is partial (1) or not (0); unit=unitless; value labels/range={0, 1}; dtype=float64; review_status=review completed
  - orf4_cp_length: description=Sequence length; unit=bp; value labels/range=[0, 186]; dtype=float64; review_status=review_recommended
  - orf4_start: description=Start position; unit=bp; value labels/range=[0, 1]; dtype=float64; review_status=review_recommended
  - orf4_stop: description=Stop position; unit=bp; value labels/range=[0, 1]; dtype=float64; review_status=review_recommended
  - orf4_tree_representative: description=flag indicating whether the contig is a tree representative (1) or not (0); unit=unitless; value labels/range={0, 1}; dtype=float64; review_status=review completed
  - notes: description=Free-text notes; unit=n/a; value labels/range={'missing orf3', 'missing stop codon in orf3', 'missing stop codon in orf4', 'missing stop codons in orf2 and orf3', 'orf1 disrupted by premature stop codon', 'orf2 disrupted by premature stop codon', 'two start codons in orf3', 'two start codons in orf3 and orf4', 'two start codons in orf4'}; dtype=object; review_status=review_recommended
  - ground_truth_category: description=Ground-truth category label; unit=n/a; value labels/range={'mas', 'oth1', 'oth2', 'tob1', 'tob2', 'tob3'}; dtype=object
  - model_prediction: description=Model prediction label; unit=unitless; value labels/range={0, 1}; dtype=int64
  - model_prediction_probabiility: description=Model prediction probability; unit=unitless; value labels/range=[0, 1]; dtype=float64
  - first_diamond_blastx_hit_name: description=first diamond blastx hit name; unit=n/a; value labels/range=free text / categorical (40 unique values); dtype=object; review_status=review completed
  - first_diamond_blastx_hit_identity_percent: description=Percent identity of the top BLASTX hit; unit=%; value labels/range=[2.8, 99.9]; dtype=float64
  - first_diamond_blastx_hit_alignment_length_aa: description=Amino-acid alignment length of the top BLASTX hit; unit=bp; value labels/range=[53, 2810]; dtype=int64
  - first_blastx_hit_name_accession: description=Accession identifier; unit=n/a; value labels/range=free text / categorical (237 unique values); dtype=object; review_status=review_recommended
  - first_blastx_hit_identity_percent: description=Percent identity of the top BLASTX hit; unit=%; value labels/range=[22, 100]; dtype=float64
  - first_blastx_hit_alignment_length_aa: description=Amino-acid alignment length of the top BLASTX hit; unit=bp; value labels/range=[50, 1889]; dtype=int64
  - first_blastx_hit_e_value: description=E-value of the top BLASTX hit; unit=unitless; value labels/range=[0, 1]; dtype=float64
  - first_blastx_hit_bit_score: description=Bit score of the top BLASTX hit; unit=unitless; value labels/range=[64.3, 3428]; dtype=float64
  - first_blastx_hit_protein: description=First blastx hit protein result; unit=n/a; value labels/range=free text / categorical (212 unique values); dtype=object; review_status=review completed
  - study_accession: description=Study accession identifier; unit=n/a; value labels/range=free text / categorical (57 unique values); dtype=object
  - study_title: description=Study title; unit=n/a; value labels/range=free text / categorical (57 unique values); dtype=object
  - organism_name: description=Organism name; unit=n/a; value labels/range=free text / categorical (44 unique values); dtype=object
  - submitter: description=Submitter name; unit=n/a; value labels/range=free text / categorical (46 unique values); dtype=object; review_status=review_recommended
  - country: description=Geographical origin of the sample; use the appropriate name from this list https://www.insdc.org/submitting-standards/geo_loc_name-qualifier-vocabulary/. Use a colon to separate the country or ocean from more detailed information about the location, eg "Canada: Vancouver" or "Germany: halfway down Zugspitze, Alps"; unit=n/a; value labels/range=free text / categorical (44 unique values); dtype=object
  - collection_date: description=the date on which the sample was collected; date/time ranges are supported by providing two dates from among the supported value formats, delimited by a forward-slash character; collection times are supported by adding "T", then the hour and minute after the date, and must be in Coordinated Universal Time (UTC), otherwise known as "Zulu Time" (Z); supported formats include "DD-Mmm-YYYY", "Mmm-YYYY", "YYYY" or ISO 8601 standard "YYYY-mm-dd", "YYYY-mm", "YYYY-mm-ddThh:mm:ss"; e.g., 30-Oct-1990, Oct-1990, 1990, 1990-10-30, 1990-10,  21-Oct-1952/15-Feb-1953, 2015-10-11T17:53:03Z; valid non-ISO dates will be automatically transformed to ISO format; unit=ISO-8601; value labels/range=[2008-07-06, 2019-01-01]; dtype=object
  - source_sample_category: description=Broad source sample category; unit=n/a; value labels/range={'Aquatic', 'Host-associated', 'Terrestrial'}; dtype=object
  - source_sample_subcategory: description=Source sample subcategory; unit=n/a; value labels/range={'Animal', 'Animal gut', 'Freshwater', 'Human', 'Human gut', 'Other', 'Plant', 'Soil', 'Wastewater'}; dtype=object
  - genbank_accession_number: description=GenBank accession number; unit=n/a; value labels/range=free text / categorical (63 unique values); dtype=object
* Variable dictionary source: S5_510discovered_contigs_variable_dictionary.reviewed.yaml (review_copy)
* Missing data codes: *list code/symbol and definition*
  - NaN/blank parsed by pandas: found 11000 value(s)
* Specialized formats or other abbreviations used: none specified

# DATA-SPECIFIC INFORMATION FOR: [Supplementary_Table6.tsv]
*repeat this section for each dataset, folder or file, as appropriate*

* Number of variables: 16
* Number of cases/rows: 125
* Variable List: *list variable name(s), description(s), unit(s) and value labels as appropriate for each*
  - type: description=group type; unit=n/a; value labels/range={'outgroup', 'tobamo'}; dtype=object; review_status=review_recommended
  - genus: description=Taxonomic genus; unit=n/a; value labels/range={'Furovirus', 'Goravirus', 'Hordeivirus', 'Pecluvirus', 'Pomovirus', 'Tobamovirus', 'Tobravirus'}; dtype=object; review_status=review_recommended
  - species: description=Taxonomic species; unit=n/a; value labels/range=free text / categorical (60 unique values); dtype=object; review_status=review_recommended
  - virus_name: description=Virus or viral taxon; unit=n/a; value labels/range=free text / categorical (88 unique values); dtype=object; review_status=review_recommended
  - accession: description=Accession identifier; unit=n/a; value labels/range=free text / categorical (125 unique values); dtype=object; review_status=review_recommended
  - RNA: description=Ribonucleic acid; unit=n/a; value labels/range={'RNA1', 'RNA2', 'RNA3'}; dtype=object
  - training: description=Training-set membership flag (training 1, not in training 0); unit=unitless; value labels/range={0, 1}; dtype=int64; review_status=review_recommended
  - length: description=Sequence length; unit=bp; value labels/range=[1799, 7226]; dtype=int64; review_status=review_recommended
  - orf1_ref: description=Reference identifier for ORF1; unit=n/a; value labels/range=free text / categorical (88 unique values); dtype=object
  - orf2_ref: description=Reference identifier for ORF2; unit=n/a; value labels/range=free text / categorical (88 unique values); dtype=object
  - cp_ref: description=Reference identifier for coat protein; unit=n/a; value labels/range=free text / categorical (87 unique values); dtype=object
  - record_id: description=UNCLEAR - human review required; unit=n/a; value labels/range=free text / categorical (125 unique values); dtype=object; review_status=unclear
  - orf1_record_id: description=Record identifier for ORF1; unit=n/a; value labels/range=free text / categorical (88 unique values); dtype=object
  - orf2_record_id: description=Record identifier for ORF2; unit=n/a; value labels/range=free text / categorical (88 unique values); dtype=object
  - cp_record_id: description=Record identifier for coat protein; unit=n/a; value labels/range=free text / categorical (87 unique values); dtype=object
  - sampling_prob: description=Sampling probability used during training set construction; unit=unitless; value labels/range=[0, 1]; dtype=float64
* Variable dictionary source: S6_reference_database_variable_dictionary.reviewed.yaml (review_copy)
* Missing data codes: *list code/symbol and definition*
  - NaN/blank parsed by pandas: found 224 value(s)
* Specialized formats or other abbreviations used: none specified
