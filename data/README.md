# Data

No data are included in this repository. All seven cohorts are distributed by their own
providers; none may be redistributed.

Request the **harmonised** files from the Gateway to Global Aging Data (https://g2aging.org/)
after registering with each cohort. The analysis expects the following layout under
`data/` (override with `AGING_DB_ROOT`):

```
data/
├── 1. CHARLS  中国/
│   └── CHARLS_中国/
│       ├── Working_data/charls.dta
│       └── Raw_data/Harmonized_CHARLS_D/H_CHARLS_D_Data.dta
├── 2. ELSA 英国/
│   ├── Working_data/elsa.dta
│   └── Raw_data/Harmonized ELSA/h_elsa_g3.dta
├── 3. HRS  美国/
│   └── HRS_美国/Working_data/hrs.dta
├── 4. KLoSA 韩国/
│   ├── KLoSA_韩国/Working_data/klosa.dta
│   └── KLoSA_韩国/Raw_data/H_KLoSA/H_KLoSA_e2.dta
├── 5. LASI 印度/
│   ├── LASI_印度/Working_data/lasi.dta
│   └── LASI_印度/Raw_data/Harmonized LASI (A.3)/H_LASI_a3.dta
├── 6. MHAS 墨西哥/
│   ├── MHAS_墨西哥/Working_data/mhas.dta
│   └── MHAS_墨西哥/Raw_data/Harmonized MHAS File/H_MHAS_c2.dta
└── 7.SHARE 欧洲/
    ├── SHARE_欧洲/Working_data/share.dta
    └── SHARE_欧洲/Raw_data/Harmonized SHARE/H_SHARE_f2.dta
```

The CHARLS wave 6 validation additionally reads the raw 2021–2023 health-status module
(`da011_15_`), expected under `data_raw/` (override with `AGING_DB_RAW_ROOT`):

```
data_raw/CHARLS/2021/Health_Status_and_Functioning.dta
```

If your files sit elsewhere, either mirror this tree or set the environment variables before
running any script:

```bash
export AGING_DB_ROOT="/path/to/harmonised"
export AGING_DB_RAW_ROOT="/path/to/charls_w6"
export OOPE02_OUT="/path/to/outputs"
```

Cohort acknowledgement and citation requirements are set by each provider; cite the cohorts as
required by their data-use agreements when publishing results based on this code.
