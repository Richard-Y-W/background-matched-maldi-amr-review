# Master Paper Journal And Display Strategy

This note records the editorial choices behind `master_research_paper.md`. It is not part of the manuscript.

## Target Interpretation

The user asked about Journal of Clinical Microbiology and "NPJ Biology." There is no single journal called NPJ Biology in the Nature Portfolio npj series. I interpreted the closest relevant target as `npj Systems Biology and Applications`, because the paper is a computational and systems-level analysis of clinical bacterial phenotypes, model behavior, co-resistance networks, and lineage-associated spectral structure.

## Journal Fit

### Journal of Clinical Microbiology

Best fit if the paper is framed as a clinical-laboratory validation guardrail for MALDI-TOF AMR prediction. The strongest JCM pitch is:

> Before MALDI-TOF AMR models are used for clinical decisions, external AUC should be decomposed by co-resistance background to determine whether prediction survives within clinically observable resistant-population contexts.

This route should keep the language practical: model output, AST panel, external site, audit decision, deployment risk. The paper should not over-emphasize method novelty at the expense of clinical microbiology value.

### npj Systems Biology and Applications

Best fit if the paper is framed as a systems-biology/methods contribution:

> MALDI-TOF spectra, AST labels, co-resistance ecology, and bacterial lineage form a coupled system. The background-matched audit quantifies which part of AMR prediction survives after controlling one observable layer of that system.

This route can foreground cross-resistance networks, public WGS-linked lineage evidence, and model-agnostic evaluation. It should still be concise and avoid becoming a software paper.

## Display Rule Applied

Both journal styles favor a clean main story and supplemental support. The main manuscript should stand on its own, while large tables, secondary model checks, stress tests, and implementation details belong in supplementary material.

## Revision Response To External Critique

The Overleaf manuscript treats the manuscript as a background-sensitivity audit, not a proof of causal confounding. The title and major claims avoid saying that the model definitively learned the resistance mechanism or that background caused the raw AUC. The central claim is narrower: raw MALDI-TOF AMR performance can be sensitive to resistant-population background, and the audit quantifies how much focal-drug ranking survives matched co-resistance backgrounds.

The closest existing method is Weis et al.'s hierarchical stratification. The revised Introduction and Discussion now compare it directly with this audit: hierarchical stratification is spectrum-level and split-level, whereas this audit is prediction-table-level and post hoc. That makes it complementary and easier for an independent auditor to run when only isolate-level predictions and AST labels are available.

The paper now explicitly frames the audit as a domain-specific extension to TRIPOD+AI/PROBAST-style prediction-model reporting. It also foregrounds Yu et al. as the WGS precedent: they showed population-structure confounding for WGS AMR prediction; this manuscript tests and audits the analogous risk for MALDI-TOF AMR prediction.

The ciprofloxacin/norfloxacin circularity concern was checked empirically in `scripts/revision_audit_checks.py`. The critique's expected pattern was not observed: valid ciprofloxacin strata at A-2018, DRIAMS-C, and DRIAMS-D were not mostly norfloxacin-susceptible; norfloxacin was unknown in those valid strata. This weakens a fluoroquinolone-pair interpretation, but it strengthens the primary cipro result in another way: the interpretable cipro rows are effectively cross-class tests of ciprofloxacin discrimination within beta-lactam and cephalosporin co-resistance backgrounds. The manuscript now presents this as a positive reinterpretation in the Results and Discussion, while still noting that norfloxacin missingness prevents a fully controlled fluoroquinolone-pair analysis.

The n >= 3 threshold is now supported by descriptive sensitivity analyses at n >= 5 and n >= 10. The ciprofloxacin versus amoxicillin-clavulanic acid contrast persists at stricter thresholds, but the manuscript describes DRIAMS-D ciprofloxacin as weak retained signal rather than strong retention. DRIAMS-B ciprofloxacin is treated as statistically uninformative because it has one valid stratum and only 25 matched isolates.

The critique's recommendation to add `Staphylococcus aureus`/oxacillin has now been addressed with the CNN/Mega Sa/Oxa panel export and background audit. The result should be framed as a focused second-organism generality check, not as broad multi-organism validation: A-2018 and DRIAMS-C retain interpretable within-background signal, while DRIAMS-B is sparse and DRIAMS-D has no usable oxacillin focal rows in the export. Weis LR official-panel oxacillin remains a published-workflow compatibility check rather than the primary Sa/Oxa evidence.

## Recommended Main Display Items

1. **Figure 1: Audit framework and metric logic.**
   - Use `manuscript/figures/figure_1_framework.pdf`.
   - Purpose: define raw AUC, matched AUC, background-centered AUC, pairwise within-background accuracy, and matched retention.

2. **Figure 2: Three-way decomposition.**
   - Use `manuscript/figures/figure_6_three_way_decomposition.pdf`.
   - Purpose: make the central decomposition visually immediate: raw MALDI AUC, no-spectrum co-resistance-only AUC, and background-centered MALDI AUC for `E. coli`/ciprofloxacin, `E. coli`/amoxicillin-clavulanic acid, and `S. aureus`/oxacillin.
   - This should be the centerpiece because it directly answers the reviewer question: how much apparent MALDI performance is already predictable from the AST background alone?

3. **Figure 3: Primary DRIAMS background audit.**
   - Use `manuscript/figures/figure_2_primary_background_audit.pdf`.
   - Purpose: show ciprofloxacin retention and amoxicillin-clavulanic acid collapse without burying the reader in all drugs.

4. **Figure 4: Public WGS-linked MALDI support.**
   - Use `manuscript/figures/figure_5_public_wgs_proteomic_support.pdf`.
   - Purpose: demonstrate that MALDI spectra strongly encode ST131 lineage and that resistance-associated discriminative bins are enriched for ST131 biomarker mass neighborhoods.
   - Keep the claim narrow: this is independent biological plausibility, not proof that DRIAMS cipro predictions are ST131-driven.

5. **Table 1: Primary interpretable audit rows.**
   - Use a shortened version of `outputs/final_framework_outputs/table_1_primary_background_matched_audit.csv`.
   - Include the `E. coli` ciprofloxacin and amoxicillin-clavulanic acid rows plus explicit caution labels for low-support DRIAMS-B rows. Move the full drug-site table to supplement.

## Move To Supplement

- Model-family comparison: important, but secondary to the primary claim. Move `manuscript/figures/figure_3_model_family_replication.pdf` and Table 2 to Supplementary Figure/Table unless the target journal asks for a larger main display set.
- Full cross-resistance network: move `manuscript/figures/figure_4_cross_resistance_network.pdf` to supplement if Figure 2 already establishes the co-resistance-only shortcut. The Results text can cite the strongest phi edges without making the heatmap a main display item.
- Falsification controls: important guardrail, but too technical for main flow. Move `manuscript/figures/figure_6_falsification_controls.pdf` and Table 15 to supplement; mention the headline in Results and Discussion.
- `S. aureus`/oxacillin detail figure: keep `manuscript/figures/figure_6_saureus_oxacillin_audit.pdf` in supplement because the main three-way decomposition already includes the second-organism check.
- Weis/Borgwardt official LR parity and six-drug stress test: supplement only. It strengthens reproducibility and audit portability, but the main paper should not look like it is mainly a Weis replication paper.
- MARISMa external stress test: supplement only. It is useful as a boundary condition and failed transfer case, not part of the main positive argument.
- Deployment decision flow and framework comparison: supplement or final discussion schematic only if the target journal wants translational guidance. Otherwise, state the decision rules in text.
- Calibration, temporal reliability, and readiness tables: supplement.

## Claim Guardrails

- Do claim that raw MALDI-TOF AMR performance can contain resistant-population background signal.
- Do claim that background-matched evaluation distinguishes retained within-background signal from background-sensitive transfer.
- Do claim that public WGS-linked MALDI data make lineage-associated spectral shortcuts biologically plausible.
- Do not claim direct ST131 detection inside DRIAMS without DRIAMS-linked WGS/MLST.
- Do not claim protein identity for mass-bin overlaps without MS/MS or LC-MS/MS.
- Do not claim exact Weis replication beyond the official 8-row LR parity panel.
- Do not present MARISMa as successful external validation.

## Online Sources Consulted

Primary journal and editorial sources checked on 2026-05-12:

- Nature Portfolio, `npj Systems Biology and Applications`, For authors: https://www.nature.com/npjsba/for-authors-and-referees/submisions
  - Initial submissions do not need special formatting if the study is suitable for editorial assessment and peer review.
  - Figures may be inserted near the relevant text for reviewer readability.
- Nature Portfolio, `npj Systems Biology and Applications`, Aims and Scope: https://www.nature.com/npjsba/aims
  - The journal considers computational and mathematical approaches to complex biological systems, including network biology and disease applications.
- Nature formatting guide: https://www.nature.com/nature/for-authors/formatting-guide
  - Nature-style Articles use a concise summary paragraph, short subheadings, and a limited main display set.
  - Typical biological/clinical Articles with 5-6 modest display items are around 4,300 words, so the master manuscript was kept near that scale.
- Nature Portfolio reporting standards for `npj Systems Biology and Applications`: https://www.nature.com/npjsba/editorial-policies/reporting-standards
  - Data, materials, code, and protocols must be made available or restrictions must be disclosed.
  - Data availability statements should identify primary and referenced datasets with access information and accession identifiers where relevant.
- ICMJE manuscript preparation recommendations: https://www.icmje.org/recommendations/browse/manuscript-preparation/preparing-for-submission.html
  - Original biomedical research usually follows IMRAD structure.
  - References should prioritize original sources, and published articles should cite persistent identifiers for datasets.
- ASM announcement that ASM journals accept initial submissions in any format: https://asm.org/press-releases/2018/asm-journals-now-accepting-submissions-in-any-form
  - This supports focusing the first submission on clarity and scientific structure rather than premature formatting.
- ASM scientific writing guidance: https://asm.org/articles/2019/november/get-your-scientific-paper-off-the-ground
  - The paper was rebuilt from results outline to figures, title/abstract, Results, Introduction, then Discussion.
- ASM figure guidance: https://asm.org/articles/2019/april/use-clear-figures-to-tell-a-story-with-your-data
  - The main display plan prioritizes a small number of figures that each carry one conclusion, with detailed audit tables in supplement.
- Journal of Clinical Microbiology archived instructions: https://pmc.ncbi.nlm.nih.gov/articles/PMC2620845/
  - JCM scope is clinical microbiology, diagnosis, and epidemiology of infection.
  - Supplemental material is appropriate for large or complex datasets, but main conclusions should be supported by the manuscript without requiring the supplement.
  - The current journals.asm.org pages were not fully accessible through the crawler, so current ASM-wide author resources were paired with this archived JCM-specific instruction page.
