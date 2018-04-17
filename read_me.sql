SELECT si.id, si.spectrum_id AS sid,
pep1.seq_mods AS pep1, pep2.seq_mods AS pep2, pep1.link_site AS linkpos1, pep2.link_site AS linkpos2,
si.charge_state AS charge,
sp.frag_tol AS fragTolerance,
si.pass_threshold, si.rank, si.ions, si.scores, pep1.crosslinker_modmass AS modmass1, pep2.crosslinker_modmass AS crosslinker_modmass,
pep1_ev.decoy AS is_decoy1,
pep2_ev.decoy AS is_decoy2,
pep1_ev.protein AS protein1,
pep2_ev.protein AS protein2,
sp.peak_list_file_name AS file,
sp.scan_id AS scan_id,
si.spectrum_id AS peakList_id
FROM spectrum_identifications AS si
LEFT JOIN spectra AS sp ON (si.spectrum_id = sp.id)
LEFT JOIN peptides AS pep1 ON (si.pep1_id = pep1.id)
LEFT JOIN
  (SELECT peptide_ref, group_concat(DISTINCT protein_accession) AS protein, group_concat(DISTINCT is_decoy) AS decoy
  FROM peptide_evidences GROUP BY peptide_ref)
AS pep1_ev ON (si.pep1_id = pep1_ev.peptide_ref)
LEFT JOIN peptides AS pep2 ON (si.pep2_id = pep2.id)
LEFT JOIN
  (SELECT peptide_ref, group_concat(DISTINCT protein_accession) AS protein, group_concat(DISTINCT is_decoy) AS decoy
  FROM peptide_evidences GROUP BY peptide_ref)
AS pep2_ev ON (si.pep2_id = pep2_ev.peptide_ref)
