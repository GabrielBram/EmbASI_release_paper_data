This dataset contains the input file, generation, HPC submission scripts, and output files for the Manuscript:

A General Pythonic Framework for DFT-in-DFT and WF-in-DFT Embedding

by Gabriel Bramley, Pavel Stishenko, Volker Blum, and Andrew Logsdail

The datasets contain calculation performed with the EmbASI QM/QM embedding wrapper (https://github.com/tamm-cci/EmbASI)
interfaced to the FHI-aims all-electron, full-potential software package.

The directory structure and the corresponding Sections in the manuscript are specified below:

Section 4.1: EMBASI_RELEASE_PAPER_DATA/Pentanol_embedding/

Each subdirectory contains calculations performed with the EmbASI QM/QM embedding wrapper for the dissocation
of a pentanol dimer over a range of distances (1.26-4.36 \AA). Simulations are performed with the NAO-VCC-4Z
and NAO-VCC-5Z basis sets of FHI-aims, and final energies are evaluated using the two-point extrapolation
scheme detailed in the manuscript. Each directory contains a .csv file with the final energies, timings,
and relavent intermediate values required for post-processing calculations.

The following calculations are performed using the PBE, PBE0, and RPA levels of theory for the full pentanol
dimer/pentanol monomer units, and are equivalent to a normal wavefunction/DFT calculation:

Pentanol_RPA/
Pentanol_PBE/
Pentanol_PBE0/

The following calculations are performed using the QM/QM embedding via EmbASI at the following levels of theory:
PBE-in-PBE, PBE0-in-PBE, and RPA@PBE-in-PBE. Projection is performed with the Huzinaga level-shift operator
and charge localisation/partitioning is performed with the Pipek-Mizey iterative localisation scheme. Unless
otherwise stated, the -OH moieties are evaluated at the high-level of theory, while the rest of the pentanol
dimer/monomer is evaluated at the stated low-level of theory. The directories for these tests are listed below:

Pentanol_PBEinPBE_SPADE_huzinaga/
Pentanol_PBE0inPBE_SPADE_huzinaga/
Pentanol_RPA@PBEinPBE_qmloc_huzinaga/

Additional calculations were performed for the RPA@PBE-in-PBE level of theory to determine the effect of the
size of the QM region on accuracy with respect to the full RPA@PBE calculation. The number of -[CH2]- units
included in the final evaluation are denoted as "nnX":

Pentanol_RPA@PBEinPBE_qmloc_huzinaga/nn0_BASIS/
Pentanol_RPA@PBEinPBE_qmloc_huzinaga/nn1_BASIS/
Pentanol_RPA@PBEinPBE_qmloc_huzinaga/nn2_BASIS/
Pentanol_RPA@PBEinPBE_qmloc_huzinaga/nn3_BASIS/
Pentanol_RPA@PBEinPBE_qmloc_huzinaga/nn4_BASIS/


Section 4.2: EMBASI_RELEASE_PAPER_DATA/Acid_dimer_embedding

Each subdirectory contains calculations performed with the EmbASI QM/QM embedding wrapper for the dissocation
of an organic acid dimer with various levels of basis truncation (1.0 |e| - 10^{-6} |e|) . Simulations are
performed with the "tight" NAO  basis sets of FHI-aims. Each directory contains a .csv file with the final
energies, timings, and relavent intermediate values required for post-processing calculations.

The following calculations are performed using the PBE, PBE0, and RPA levels of theory for the full organic acid
dimer/monomer units, and are equivalent to a normal wavefunction/DFT calculation:

C22COOH_dimer_full_mol_tight/trunc-None_tight_PBE_PBE_None/
C22COOH_dimer_full_mol_tight/trunc-None_tight_PBE0_PBE0_None/
C22COOH_dimer_full_mol_tight/trunc-None_tight_PBE_PBE_RPA/

The following calculations are performed using the QM/QM embedding via EmbASI at the following levels of theory:
PBE-in-PBE, PBE0-in-PBE, and RPA@PBE-in-PBE. Projection is performed with the Huzinaga level-shift operator
and charge localisation/partitioning is performed SPADE paritioning scheme. Unless otherwise stated, the -COOH
moieties are evaluated at the high-level of theory, while the rest of the organic acid dimer dimer/monomer is
evaluated at the stated low-level of theory. The charge parameter used to perform truncation is denoted as:
"trunc-<trunc>_<basis>_<high-level>_<low-level>_<post-HF-method>". The directories for these tests are listed below:

Acid_dimer_embedding/C22COOH_dimer_truncation_schuz_tight/trunc-<trunc>_tight_PBE_PBE_None/
Acid_dimer_embedding/C22COOH_dimer_truncation_schuz_tight/trunc-<trunc>_tight_PBE0_PBE_None/
Acid_dimer_embedding/C22COOH_dimer_truncation_schuz_tight/trunc-<trunc>_tight_PBE_PBE_RPA/
