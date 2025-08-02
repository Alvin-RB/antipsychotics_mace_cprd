# Risk of major adverse cardiovascular events with aripiprazole versus olanzapine, quetiapine, and risperidone in severe mental illness: a target trial emulation.

<i>Status: Under review</i>

> Richards-Belle A, Launders N, Hardoon S, Man KKC, Davies NM, Bramon E, Hayes JF, Osborn DPJ. Risk of major adverse cardiovascular events with aripiprazole versus olanzapine, quetiapine, and risperidone in severe mental illness: a target trial emulation.

This study aimed to emulate a trial on the following research question: In adults diagnosed with SMI who are prescribed a new antipsychotic in primary care, what is the effect of aripiprazole, as compared to olanzapine, quetiapine, and risperidone, on the long-term risk of major adverse cardiovascular events? The data source is [Clinical Practice Research Datalink (CPRD)](https://www.cprd.com/). To support transparency and open science, this repository hosts code and materials which support the manuscript - including:

#### Code lists

1. SMI - [CPRD Aurum](https://github.com/Alvin-RB/antipsychotics_descriptive_study_cprd/blob/main/Aurum_SMI_codelist_21032024.txt) and [CPRD GOLD](https://github.com/Alvin-RB/antipsychotics_descriptive_study_cprd/blob/main/GOLD_SMI_codelist_21032024.txt)
   - used to define SMI
   
2. Antipsychotics - [CPRD Aurum](https://github.com/Alvin-RB/antipsychotics_descriptive_study_cprd/blob/main/antipsychotics_AURUM_250324.txt) and [CPRD GOLD](https://github.com/Alvin-RB/antipsychotics_descriptive_study_cprd/blob/main/antipsychotics_GOLD_250324.txt)
    - used to identify prescriptions of antipsychotics

#### R Scripts

1. [Cohort derivation](https://github.com/Alvin-RB/antipsychotics_mace_cprd/blob/main/R%20Scripts/1.0%20AP%20and%20MACE%20-%20Cohort%20derivation_github.R)
   
2. [Multiple imputation](https://github.com/Alvin-RB/antipsychotics_mace_cprd/blob/main/R%20Scripts/2.0%20AP%20and%20MACE%20-%20Multiple%20imputation_github.R)
   
3. [Overlap weights](https://github.com/Alvin-RB/antipsychotics_mace_cprd/blob/main/R%20Scripts/3.0%20AP%20and%20MACE%20-%20Overlap%20weights_github.R)

4. [Descriptive analysis](https://github.com/Alvin-RB/antipsychotics_mace_cprd/blob/main/R%20Scripts/4.0%20AP%20and%20MACE%20-%20Descriptives_github.R)

5. [Cox models](https://github.com/Alvin-RB/antipsychotics_mace_cprd/blob/main/R%20Scripts/5.0%20AP%20and%20MACE%20-%20Cox%20regression%20hazard%20ratios_github.R)

6. [Pooled logistic regression](https://github.com/Alvin-RB/antipsychotics_mace_cprd/blob/main/R%20Scripts/6.0%20AP%20and%20MACE%20-%20Overlap%20weights%20and%20PLR%20-%20primary%20analysis_github.R)

7. [Inverse probability of censoring weighting - sensitivity analysis](https://github.com/Alvin-RB/antipsychotics_mace_cprd/blob/main/R%20Scripts/7.0%20AP%20and%20MACE%20-%20Sensitivity_IPCW_timevarying_github.R)

#### Contact

If you would like any further information, then please [contact me](https://github.com/Alvin-RB).
