# Chemogenomic-DTI-Prediction-Methods

## Table of Contents

* [About This Repository](#about-this-repository)
* [Alternative Online Source Codes for DTI Prediction Algorithms](#alternative-online-source-codes-for-dti-prediction-algorithms)
* [Web Servers](#web-servers)
* [Databases](#databases)
* [Datasets](#datasets)
* [Surveys](#surveys)

----

## About This Repository
This repository has been initially established to hold the source code used for conducting the experiments explained in the following paper:

**[Computational Prediction of Drug-Target Interactions using Chemogenomic Approaches: An Empirical Survey [2018]](https://doi.org/10.1093/bib/bby002)**  
Briefings in Bioinformatics  
*Ali Ezzat, Min Wu, Xiao-Li Li, Chee-Keong Kwoh*

This repository consists of two folders:
* **"1. Main":** contains the source code that was used to generate the cross validation results displayed in the main text. The prediction methods that were compared in the main text belong to the category of *similarity-based methods*.
* **"2. Feature-based Methods":** contains the source code that was used to generate the results displayed in the supplementary text accompanying the publication. The prediction methods involved belong to the category of *feature-based methods*.

The source code is entirely written in the MATLAB programming language. In the future, it is our intention to regularly add source codes for more DTI prediction algorithms, especially those that were not released in their respective papers. Stay tuned!

----

## Alternative Online Source Codes for DTI Prediction Algorithms
The list below consists of other links containing source codes that have been made available online along with their corresponding publications (sorted in chronological order): 

|Algorithm|Source Code|Publication|
|---------|:---------:|----------:|
|Jacob et al.       |[Source Code](https://academic.oup.com/bioinformatics/article/24/19/2149/247731)
|SVM-based BLM      |[Source Code](http://members.cbio.mines-paristech.fr/~yyamanishi/bipartitelocal/)  |[Publication](https://academic.oup.com/bioinformatics/article/25/18/2397/197654)|
|RLS-avg & RLS-kron |[Source Code](http://cs.ru.nl/~tvanlaarhoven/drugtarget2011/)                      |[Publication](https://academic.oup.com/bioinformatics/article/27/21/3036/216840)|
|KBMF2K             |[Source Code](http://users.ics.aalto.fi/gonen/bioinfo12.php)                       |[Publication](https://academic.oup.com/bioinformatics/article/28/18/2304/241817)|
|DT-Hybrid          |[Source Code](https://alpha.dmi.unict.it/dtweb/dthybrid.php)                       |[Publication](https://academic.oup.com/bioinformatics/article/29/16/2004/199066)|
|RLS-WNN            |[Source Code](http://cs.ru.nl/~tvanlaarhoven/drugtarget2013/)                      |[Publication](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0066952)|
|MHL1SVM            |[Source Code](https://sites.google.com/site/interactminhash/)                      |[Publication](https://bmcsystbiol.biomedcentral.com/articles/10.1186/1752-0509-7-S6-S3)|
|Nanni et al.       |[Source Code](https://www.dropbox.com/s/gv37ujz3f93h5ye/ToolProteinDrug.rar)       |[Publication](https://www.sciencedirect.com/science/article/pii/S0022519314003452)|
|semEP              |[Source Code](https://github.com/gpalma/semep)                                     |[Publication](https://link.springer.com/chapter/10.1007%2F978-3-319-11964-9_9)|
|PSL                |[Source Code](https://github.com/shobeir/fakhraei_tcbb2014)                        |[Publication](https://ieeexplore.ieee.org/document/6817596/)|
|PUCPI              |[Source Code](http://admis.fudan.edu.cn/projects/pucpi.html)                       |[Publication](https://ieeexplore.ieee.org/document/7471459/)|
|GIFT               |[Source Code](http://bioinfo.au.tsinghua.edu.cn/software/GIFT/)                    |[Publication](https://academic.oup.com/bioinformatics/article/31/15/2523/188618)|
|KronRLS-MKL        |[Source Code](http://www.cin.ufpe.br/~acan/kronrlsmkl/)                            |[Publication](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-0890-3)|
|Coelho et al.      |[Source Code](http://bioinformatics.ua.pt/software/dtipred/)                       |[Publication](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005219)|
|RLS-KF             |[Source Code](https://github.com/minghao2016/RLS-KF)                               |[Publication](https://www.sciencedirect.com/science/article/pii/S0003267016300630)|
|NRLMF              |[Source Code](https://github.com/stephenliu0423/PyDTI)                             |[Publication](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004760)|
|DTINet             |[Source Code](https://github.com/luoyunan/DTINet)                                  |[Publication](https://www.nature.com/articles/s41467-017-00680-8)|
|DeepWalk           |[Source Code](https://github.com/zongnansu1982/drug-target-prediction)             |[Publication](https://academic.oup.com/bioinformatics/article-abstract/33/15/2337/3738543?redirectedFrom=fulltext)|
|DeepDTIs           |[Source Code](https://github.com/Bjoux2/DeepDTIs_DBN)                              |[Publication](https://pubs.acs.org/doi/abs/10.1021/acs.jproteome.6b00618)|
|DNILMF             |[Source Code](https://github.com/minghao2016/DNILMF)                               |[Publication](https://www.nature.com/articles/srep40376)|
|GRMF               |[Source Code](http://www.comp.nus.edu.sg/~lixl/GRMF/index.html)                    |[Publication](https://doi.org/10.1109/TCBB.2016.2530062)|
|pairwiseMKL        |[Source Code](https://github.com/aalto-ics-kepaco/pairwiseMKL)                     |[Publication](https://doi.org/10.1093/bioinformatics/bty277)|


----

## Web Servers
Following is a list of online DTI prediction servers:

* [BalestraWeb](http://balestra.csb.pitt.edu/)
* [D3TPredictor](https://www.d3pharma.com/d3tpredictor)
* [DASPfind](http://www.cbrc.kaust.edu.sa/daspfind/)
* [DINIES](http://www.genome.jp/tools/dinies/)
* [DT-Web](https://alpha.dmi.unict.it/dtweb/)
* [iDTI-ESBoost](http://farshidrayhan.pythonanywhere.com/iDTI-ESBoost/)
* [iGPCR-Drug](http://www.jci-bioinfo.cn/iGPCR-Drug/)
* [SuperPred](http://prediction.charite.de/)
* [SwissTargetPrediction](http://www.swisstargetprediction.ch/)

----

## Databases
Below is a list of databases that have previously been used in efforts pertaining to drug-target interaction prediction and, more generally, drug discovery. This list was compiled  with the help of [this paper](https://academic.oup.com/bib/article-abstract/17/4/696/2240330).

* [ASDCD: Antifungal Synergistic Drug Combination Database](http://asdcd.amss.ac.cn/)
* [BindingDB](www.bindingdb.org/)
* [BioLip](https://zhanglab.ccmb.med.umich.edu/BioLiP/)
* [CancerDR](http://crdd.osdd.net/raghava/cancerdr/)
* [canSAR](http://cansar.icr.ac.uk/)
* [ChemBank](http://chembank.broadinstitute.org/)
* [ChEMBL](https://www.ebi.ac.uk/chembl/)
* [DCDB: Drug Combination Database](http://www.cls.zju.edu.cn/dcdb/)
* [DrugBank](http://drugbank.ca/)
* [Drug Target Commons (DTC)](https://drugtargetcommons.fimm.fi/)
* [FAERS: FDA Adverse Event Reporting System](https://www.fda.gov/Drugs/GuidanceComplianceRegulatoryInformation/Surveillance/AdverseDrugEffects/default.htm)
* [Integrity](integrity.thomson-pharma.com)
* [IUPHAR/BPS Guide to Pharmacology](http://www.guidetopharmacology.org/)
* [JAPIC: Japan Pharmaceutical Information Center](http://www.japic.or.jp/)
* [KEGG: Kyoto Encyclopedia of Genes and Genomes](http://www.kegg.jp/)
* [MATADOR: Manually Annotated Targets and Drugs Online Resource](http://matador.embl.de/)
* [NRDTD: ncRNA Drug Targets Database](http://chengroup.cumt.edu.cn/NRDTD/)
* [Open Targets](https://www.targetvalidation.org/)
* [PDSP: Psychoactive Drug Screening Program](https://pdsp.unc.edu/)
* [SIDER](http://sideeffects.embl.de/)
* [STITCH: Search Tool for Interactions of Chemicals](http://stitch.embl.de/)
* [SuperTarget](http://insilico.charite.de/supertarget/) ([alternative link](http://bioinf-apache.charite.de/supertarget/))
* [TDR (Tropical Disease Research) Targets](http://tdrtargets.org/)
* [TTD: Therapeutic Target Database](https://db.idrblab.org/ttd/)
* [ZINC](http://zinc.docking.org/)

For more data sources that are used in drug discovery in general (and not just DTI prediction), please refer to *Table 2* of the following paper:

**[Toward better drug repositioning: prioritizing and integrating existing methods into efficient pipelines [2014]](https://doi.org/10.1016/j.drudis.2013.11.005)**  
Drug Discovery Today  
*Guangxu Jin, Stephen T.C. Wong*

----

## Datasets

Listed below are datasets that have been compiled by other researchers and used in DTI prediction efforts.

|Dataset|Link|Publication|
|-------|:---------:|----------:|
|Yamanishi et al., 2008 |[Download Link](http://web.kuicr.kyoto-u.ac.jp/supp/yoshi/drugtarget/)                               |[Publication](https://doi.org/10.1093/bioinformatics/btn162) |
|Yamanishi et al., 2010 |[Download Link](http://members.cbio.mines-paristech.fr/~yyamanishi/pharmaco/)                        |[Publication](https://doi.org/10.1093/bioinformatics/btq176) |
|Tabei et al., 2012     |[Download Link](http://members.cbio.mines-paristech.fr/~yyamanishi/l1binary/)                        |[Publication](https://doi.org/10.1093/bioinformatics/bts412) |
|Tabei et al., 2013     |[Download Link](https://sites.google.com/site/interactminhash/)                                      |[Publication](https://doi.org/10.1186/1752-0509-7-S6-S3)     |
|Xiao et al., 2013a     |[Download Link](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0072234)            |[Publication](https://doi.org/10.1371/journal.pone.0072234)  |
|Xiao et al., 2013b     |[Download Link](http://www.sciencedirect.com/science/article/pii/S0022519313003871)                  |[Publication](https://doi.org/10.1016/j.jtbi.2013.08.013)    |
|Min et al., 2013       |[Download Link](https://www.hindawi.com/journals/bmri/2013/701317/sup/)                              |[Publication](http://dx.doi.org/10.1155/2013/701317)         |
|Fan et al., 2014       |[Download Link](http://www.mdpi.com/1422-0067/15/3/4915#supplementary)                               |[Publication](https://doi.org/10.3390/ijms15034915)          |
|Nanni et al., 2014     |[Download Link](https://www.dropbox.com/s/gv37ujz3f93h5ye/ToolProteinDrug.rar)                       |[Publication](https://doi.org/10.1016/j.jtbi.2014.06.008)    |
|Ezzat et al., 2016     |[Download Link](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1377-y#Sec2) |[Publication](https://doi.org/10.1186/s12859-016-1377-y)     |
|Cheng et al., 2016     |[Download Link](http://admis.fudan.edu.cn/projects/pucpi.html)                                       |[Publication](https://doi.org/10.1109/TCBB.2016.2570211)     |
|Nascimento et al., 2016|[Download Link](http://www.cin.ufpe.br/~acan/kronrlsmkl/)                                            |[Publication](https://doi.org/10.1186/s12859-016-0890-3)     |
|Coelho et al., 2016    |[Download Link](http://bioinformatics.ua.pt/software/dtipred/)                                       |[Publication](https://doi.org/10.1371/journal.pcbi.1005219)  |
|Li at al., 2016        |[Download Link](https://academic.oup.com/bioinformatics/article/32/7/1057/1743938)                   |[Publication](https://academic.oup.com/bioinformatics/article/32/7/1057/1743938) |
|Zong et al., 2017      |[Download Link](https://github.com/zongnansu1982/drug-target-prediction)                             |[Publication](https://doi.org/10.1093/bioinformatics/btx160) |
|Wen at al., 2017       |[Download Link](https://github.com/Bjoux2/DeepDTIs_DBN)                                              |[Publication](https://doi.org/10.1021/acs.jproteome.6b00618) |


----

## Surveys

Following are a list of surveys on DTI prediction that have been published over the years. Note that some of these surveys contain alternative source codes for the DTI prediction algorithms being surveyed.

**[Similarity-based machine learning methods for predicting drug-target interactions - a brief review [2013]](https://doi.org/10.1093/bib/bbt056)**  
Briefings in Bioinformatics  
*Hao Ding, Ichigaku Takigawa, Hiroshi Mamitsuka, Shanfeng Zhu*

**[Toward more realistic drug–target interaction predictions [2014]](https://doi.org/10.1093/bib/bbu010)**  
Briefings in Bioinformatics  
*Tapio Pahikkala, Antti Airola, Sami Pietilä, Sushil Shakyawar, Agnieszka Szwajda, Jing Tang, Tero Aittokallio*

**[Drug-target interaction prediction via chemogenomic space - learning-based methods [2014]](https://doi.org/10.1517/17425255.2014.950222)**  
Expert Opinion on Drug Metabolism & Toxicology  
*Zaynab Mousavian, Ali Masoudi-Nejad*

**[Drug–target interaction prediction: databases, web servers and computational models  [2016]](https://doi.org/10.1093/bib/bbv066)**  
Briefings in Bioinformatics  
*Xing Chen, Chenggang Clarence Yan, Xiaotian Zhang, Xu Zhang, Feng Dai, Jian Yin, Yongdong Zhang*

**[Large-Scale Prediction of Drug-Target Interaction: a Data-Centric Review [2017]](https://doi.org/10.1208/s12248-017-0092-6)**  
The AAPS Journal, Springer  
*Tiejun Cheng, Ming Hao, Takako Takeda, Stephen H. Bryant, Yanli Wang*

**[Open-source chemogenomic data-driven algorithms for predicting drug–target interactions [2018]](https://doi.org/10.1093/bib/bby010)**  
Briefings in Bioinformatics  
*Ming Hao, Stephen H. Bryant, Yanli Wang*

**[Computational Prediction of Drug-Target Interactions using Chemogenomic Approaches: An Empirical Survey [2018]](https://doi.org/10.1093/bib/bby002)**  
Briefings in Bioinformatics  
*Ali Ezzat, Min Wu, Xiao-Li Li, Chee-Keong Kwoh*

