# ============================================================================
# references.R — Master reference list for the UKB MS proteomics manuscript
# Keys are short alphanumeric IDs; values are formatted citation strings.
# Citation order (number) is determined by first call to cite() in the document,
# NOT by order here. Edit/add/delete freely; numbers auto-update on next knit.
# ============================================================================

REFERENCES <- list(

  # ── Introduction refs (expected order 1–14) ──────────────────────────────
  walton2020      = "Walton C, et al. Rising prevalence of multiple sclerosis worldwide: insights from the Atlas of MS, third edition. *Mult Scler* **26**, 1816–1821 (2020).",
  gbd2019         = "GBD 2016 Multiple Sclerosis Collaborators. Global, regional, and national burden of multiple sclerosis 1990–2016. *Lancet Neurol* **18**, 269–285 (2019).",
  demeo2025       = "De Meo E, et al. Disability worsening phenotypes in multiple sclerosis and impact of disease-modifying treatments. *Neurology* **105**, e214408 (2025).",
  wijnands2017    = "Wijnands JMA, et al. Health-care use before a first demyelinating event suggestive of a multiple sclerosis prodrome. *Lancet Neurol* **16**, 445–451 (2017).",
  ruizalguero2025 = "Ruiz-Algueró M, et al. Health care use before multiple sclerosis symptom onset. *JAMA Netw Open* **8**, e2524635 (2025).",
  marrie2022      = "Marrie RA, et al. From the prodromal stage of multiple sclerosis to disease prevention. *Nat Rev Neurol* **18**, 559–572 (2022).",
  abdelhak2026    = "Abdelhak A, et al. Myelin injury precedes axonal injury and symptomatic onset in multiple sclerosis. *Nat Med* **32**, 362–368 (2026).",
  bjornevik2020   = "Bjornevik K, et al. Serum neurofilament light chain levels in patients with presymptomatic multiple sclerosis. *JAMA Neurol* **77**, 58–64 (2020).",
  disanto2018     = "Disanto G, et al. Prodromal symptoms of multiple sclerosis in primary care. *Ann Neurol* **83**, 1162–1173 (2018).",
  khalil2024      = "Khalil M, et al. Neurofilaments as biomarkers in neurological disorders: towards clinical application. *Nat Rev Neurol* **20**, 269–287 (2024).",
  goris2022       = "Goris A, et al. Genetics of multiple sclerosis: lessons from polygenicity. *Lancet Neurol* **21**, 830–842 (2022).",
  bjornevik2022   = "Bjornevik K, et al. Longitudinal analysis reveals high prevalence of Epstein-Barr virus associated with multiple sclerosis. *Science* **375**, 296–301 (2022).",
  olsson2017      = "Olsson T, Barcellos LF, Alfredsson L. Interactions between genetic, lifestyle and environmental risk factors for multiple sclerosis. *Nat Rev Neurol* **13**, 25–36 (2017).",
  gross2024       = "Gross CC, et al. Multiple sclerosis endophenotypes identified by high-dimensional blood signatures are associated with distinct disease trajectories. *Sci Transl Med* **16**, eade8560 (2024).",

  # ── Discussion refs (expected order 15–44) ────────────────────────────────
  kim2003         = "Kim JY, et al. The role of ErbB2 signaling in the onset of terminal differentiation of oligodendrocytes in vivo. *J Neurosci* **23**, 5561–5571 (2003).",
  kataria2021     = "Kataria H, et al. Neuregulin-1 beta 1 is implicated in pathogenesis of multiple sclerosis. *Brain* **144**, 162–185 (2021).",
  kataria2019     = "Kataria H, Alizadeh A, Karimi-Abdolrezaee S. Neuregulin-1/ErbB network: an emerging modulator of nervous system injury and repair. *Prog Neurobiol* **180**, 101643 (2019).",
  indari2021      = "Indari O, Chandramohanadas R, Jha HC. Epstein-Barr virus infection modulates blood-brain barrier cells. *Pathog Dis* **79**, ftaa080 (2021).",
  tsai2013        = "Tsai SC, et al. Autocrine CCL3 and CCL4 induced by the oncoprotein LMP1 promote Epstein-Barr virus-triggered B cell proliferation. *J Virol* **87**, 9041–9052 (2013).",
  etingin1993     = "Etingin OR, Silverstein RL, Hajjar DP. Von Willebrand factor mediates platelet adhesion to virally infected endothelial cells. *Proc Natl Acad Sci USA* **90**, 5153–5156 (1993).",
  balashov1999    = "Balashov KE, et al. CCR5+ and CXCR3+ T cells are increased in multiple sclerosis and their ligands MIP-1α and IP-10 are expressed in demyelinating brain lesions. *Proc Natl Acad Sci USA* **96**, 6873–6878 (1999).",
  puthenparampil2020 = "Puthenparampil M, et al. Wide cytokine analysis in cerebrospinal fluid at diagnosis identified CCL-3 as a possible prognostic factor for multiple sclerosis. *Front Immunol* **11**, 174 (2020).",
  ziliotto2018    = "Ziliotto N, et al. Hemostasis biomarkers in multiple sclerosis. *Eur J Neurol* **25**, 1169–1176 (2018).",
  zierfuss2024    = "Zierfuss B, Larochelle C, Prat A. Blood-brain barrier dysfunction in multiple sclerosis: causes, consequences, and potential effects of therapies. *Lancet Neurol* **23**, 95–109 (2024).",
  dorrier2021     = "Dorrier CE, et al. CNS fibroblasts form a fibrotic scar in response to immune cell infiltration. *Nat Neurosci* **24**, 234–244 (2021).",
  pikor2015       = "Pikor NB, et al. Integration of Th17- and lymphotoxin-derived signals initiates meningeal-resident stromal cell remodeling to propagate neuroinflammation. *Immunity* **43**, 1160–1173 (2015).",
  cheng2011       = "Cheng H, et al. Crystal structure of leukocyte Ig-like receptor LILRB4 (ILT3/LIR-5/CD85k). *J Biol Chem* **286**, 18013–18025 (2011).",
  katz2007        = "Katz HR. Inhibition of pathologic inflammation by leukocyte Ig-like receptor B4 and related inhibitory receptors. *Immunol Rev* **217**, 222–230 (2007).",
  hou2024         = "Hou J, et al. Antibody-mediated targeting of human microglial leukocyte Ig-like receptor B4 attenuates amyloid pathology in a mouse model. *Sci Transl Med* **16**, eadj9052 (2024).",
  imsgc2007       = "International Multiple Sclerosis Genetics Consortium. Risk alleles for multiple sclerosis identified by a genomewide study. *N Engl J Med* **357**, 851–862 (2007).",
  maier2009       = "Maier LM, et al. Soluble IL-2RA levels in multiple sclerosis subjects and the effect of soluble IL-2RA on immune responses. *J Immunol* **182**, 1541–1547 (2009).",
  hartmann2014    = "Hartmann FJ, et al. Multiple sclerosis-associated IL2RA polymorphism controls GM-CSF production in human TH cells. *Nat Commun* **5**, 5056 (2014).",
  kosa2025        = "Kosa P, et al. Longitudinal proteomic profiling of cerebrospinal fluid in untreated multiple sclerosis defines evolving disease biology. *Nat Commun* **16**, 8891 (2025).",
  perwieniec2021  = "Perwieniec J, et al. Gender-related differences in prodromal multiple sclerosis characteristics: a 7-year observation study. *J Clin Med* **10**, 3821 (2021).",
  kreft2026       = "Kreft KL, et al. Genetic subtypes associated with multiple sclerosis severity and response to treatment. *J Neurol Neurosurg Psychiatry* **97**, 413–421 (2026).",
  okuda2023       = "Okuda DT, et al. Dimethyl fumarate delays multiple sclerosis in radiologically isolated syndrome. *Ann Neurol* **93**, 604–614 (2023).",
  lebrunfreney2023 = "Lebrun-Frénay C, et al. Teriflunomide and time to clinical multiple sclerosis in patients with radiologically isolated syndrome: the TERIS randomized clinical trial. *JAMA Neurol* **80**, 1080–1088 (2023).",
  herold2019      = "Herold KC, et al. An anti-CD3 antibody, teplizumab, in relatives at risk for type 1 diabetes. *N Engl J Med* **381**, 603–613 (2019).",
  cope2024        = "Cope AP, et al. Abatacept in individuals at high risk of rheumatoid arthritis (APIPPRA): a randomised, double-blind, multicentre, parallel, placebo-controlled, phase 2b clinical trial. *Lancet* **403**, 838–849 (2024).",
  fry2017         = "Fry A, et al. Comparison of sociodemographic and health-related characteristics of UK Biobank participants with those of the general population. *Am J Epidemiol* **186**, 1026–1034 (2017).",
  nova2024        = "Nova A, Di Caprio G, Bernardinelli L, Fazia T. Genetic and early life factors influence on time-to-multiple sclerosis diagnosis: a UK Biobank study. *Mult Scler* **30**, 994–1003 (2024).",
  schoeler2023    = "Schoeler T, et al. Participation bias in the UK Biobank distorts genetic associations and downstream analyses. *Nat Hum Behav* **7**, 1216–1227 (2023).",
  rahimov2025     = "Rahimov F, et al. Common diseases in clinical cohorts: not always what they seem. *N Engl J Med* **393**, 1589–1598 (2025).",
  rooney2025      = "Rooney MR, et al. Correlations within and between highly multiplexed proteomic assays of human plasma. *Clin Chem* **71**, hvaf030 (2025)."
)
