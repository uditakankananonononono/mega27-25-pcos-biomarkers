# Source provenance and rights documentation gaps (pcos child, 2026-10-08 audit record)

This file records what this repository's own documents say about source hashes and terms. It is a documentation record, not a license verdict, not an integrity certificate, and not a clearance of any recorded rights hold. A matching claim below means two recorded hash strings are equal at a named field in pinned repository documents. It does not mean any assay bytes were re-downloaded or verified. No rights record was found for these sources, so reuse and redistribution are NOT cleared by this file.

Audited child commit: `e94e6d244c6c0812c7442ab354eb588d4ff54cbf`  
Compared shared-parent repo: `mega27-25-biomarkers-underserved-diseases` at `36ff87995f4fe0e1d08f9e2b2cd885d341b88937`

Rights: No source-specific rights check found in inspected disease-child documentation. Hash matches are not legal clearance or byte verification.

## 1. Result-file hash records in this child (9)

Hash values are shown as 12-hex display prefixes only; the full value is in the cited file and field. Classes are kept separate on purpose: an expression/assay source hash, a GEO series-matrix hash (metadata that may include expression), a reference annotation, and a published artifact are different kinds of record.

- `results/pcos_p26_result.json` `/matrix_sha256` hash `93335008afe2...` (GSE267287); class: expression_or_assay_source_hash_record; equal to shared-parent field
- `results/pcos_p26_result.json` `/count_sha256` hash `38caa96acb45...` (GSE267287); class: expression_or_assay_source_hash_record; equal to shared-parent field
- `results/pcos_p8.json` `/matrix_sha256` hash `00d0b39c1087...` (GSE293353); class: expression_or_assay_source_hash_record; equal to shared-parent field
- `results/pcos_p8.json` `/sample_info_sha256` hash `71bbef8f05e0...` (GSE293353); class: sample_or_reference_metadata_not_assay; equal to shared-parent field
- `results/pcos_p46_result.json` `/data_sha256` hash `c260e9e1a446...` (GSE271363); class: expression_or_assay_source_hash_record; equal to shared-parent field
- `results/pcos_p46_result.json` `/meta_sha256` hash `fc3d45b68187...` (GSE271363); class: sample_or_reference_metadata_not_assay; equal to shared-parent field
- `results/pcos_p27_result.json` `/sha256/GSE294074_series_matrix.txt.gz` hash `575fef37c3ef...` (GSE294074); class: GEO_series_matrix_file_includes_metadata_and_may_include_expression; equal to shared-parent field
- `results/pcos_p27_result.json` `/sha256/GSE294074_RAW.tar` hash `aced9e17fed2...` (GSE294074); class: source_hash_record_target_needs_context; equal to shared-parent field
- `results/pcos_p27_result.json` `/sha256/GSE294074_annotation.xls.gz` hash `f2ab5c472dfb...` (GSE294074); class: sample_or_reference_metadata_not_assay; equal to shared-parent field

## 4. Metadata crosswalk tables (not assay bytes)

These per-sample tables carry hash columns describing GEO source-response metadata. They are not blanket expression-matrix integrity.

- `projects/pcos/sources/GSE5090_used_sample_crosswalk.csv`: columns ['sha256'], 17 records; file identical to a shared-parent file: True; scope: GEO sample/series source-response metadata; not blanket expression matrix integrity
- `projects/pcos/sources/GSE271363_used_sample_crosswalk.csv`: columns ['sha256'], 16 records; file identical to a shared-parent file: False; scope: GEO sample/series source-response metadata; not blanket expression matrix integrity
- `projects/pcos/sources/gse43322/GSE43322_used_sample_crosswalk.csv`: columns ['source_sha256', 'matrix_sha256'], 31 records; file identical to a shared-parent file: False; scope: GEO sample/series source-response metadata; not blanket expression matrix integrity

## Coverage boundary

- Source accessions listed in the child manifest: 20.
- Of those, 16 have no mapped result-file payload hash in this audit: GSE106724, GSE124226, GSE137684, GSE155489, GSE262735, GSE277906, GSE304677, GSE34526, GSE43264, GSE43266, GSE43322, GSE5090, GSE54248, GSE54250, GSE6798, GSE80432.
- 0 hash fields inherited from other diseases' records are not counted toward this child.
- Records absent from child manifest can still have result-file hashes. Neither presence nor absence proves assay acquisition by this scout. Other-disease copied hashes never count toward this child.
