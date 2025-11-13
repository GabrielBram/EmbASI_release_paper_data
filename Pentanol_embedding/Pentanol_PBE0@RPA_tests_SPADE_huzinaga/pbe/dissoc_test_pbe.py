import os
from asiembedding.embedding import ProjectionEmbedding, StandardDFT
from ase.calculators.aims import Aims
from ase.build import molecule
from ase.visualize import view
from ase.data.s22 import s22, s26, create_s22_system
from ase.io import read
import numpy as np
from csv_writer import write_csv_data_row, run_single_fragment_full_molecule, run_binding_energy_test_full_molecule

os.environ['ASI_LIB_PATH'] = "/scratch/c.sacgb4/CCI-GEP_tests/Software/FHIaims/_build_avx2_lib_2018/libaims.240826.scalapack.mpi.so"

calc_ll = Aims(xc='PBE',
    output_level="full",
    relativistic='atomic_zora scalar',
    occupation_type="gaussian 0.01",
    mixer="pulay",
    n_max_pulay=10,
    KS_method="parallel",
    RI_method="LVL",
    collect_eigenvectors = True,
    postprocess_anyway = True,
    density_update_method='density_matrix', # for DM export
    lmo_init_guess="canon",
    compensate_multipole_errors=True,
    atomic_solver_xc="PBE",
#    lmo_pm_charge_metric="mulliken",
    compute_kinetic=True,
    basis_threshold=1e-5,
    override_illconditioning=True,
    )

calc_hl = Aims(xc='PBE',
    relativistic='atomic_zora scalar',
    occupation_type="gaussian 0.01",
    mixer="pulay",
    KS_method="parallel",
    RI_method="LVL",
    n_max_pulay=10,
    collect_eigenvectors=True,
    postprocess_anyway = True,
    density_update_method='density_matrix', # for DM export
    lmo_pm_charge_metric="mulliken",
    compensate_multipole_errors=True,
    atomic_solver_xc="PBE",
    compute_kinetic=True,
    basis_threshold=1e-5,
    override_illconditioning=True,
  )

pentanol_dimer = read('pentanol.xyz')
pentanol = pentanol_dimer[0:18]

os.environ['AIMS_SPECIES_DIR'] = "/scratch/c.sacgb4/CCI-GEP_tests/Software/FHIaims/species_defaults/NAO-VCC-nZ/NAO-VCC-5Z/"
run_single_fragment_full_molecule(pentanol, hl_calc=calc_hl, ll_calc=calc_ll, hl_xc="PBE", ll_xc="PBE", post_hf_method=None, basis="NAO-VCC-5Z")
run_binding_energy_test_full_molecule(pentanol_dimer, hl_calc=calc_hl, ll_calc=calc_ll,hl_xc="PBE", ll_xc="PBE", post_hf_method=None, basis="NAO-VCC-5Z")

os.environ['AIMS_SPECIES_DIR'] = "/scratch/c.sacgb4/CCI-GEP_tests/Software/FHIaims/species_defaults/NAO-VCC-nZ/NAO-VCC-4Z/"
run_single_fragment_full_molecule(pentanol, hl_calc=calc_hl, ll_calc=calc_ll, hl_xc="PBE", ll_xc="PBE", post_hf_method=None, basis="NAO-VCC-4Z")
run_binding_energy_test_full_molecule(pentanol_dimer, hl_calc=calc_hl, ll_calc=calc_ll,hl_xc="PBE", ll_xc="PBE", post_hf_method=None, basis="NAO-VCC-4Z")

os.environ['AIMS_SPECIES_DIR'] = "/scratch/c.sacgb4/CCI-GEP_tests/Software/FHIaims/species_defaults/NAO-VCC-nZ/NAO-VCC-3Z/"
run_single_fragment_full_molecule(pentanol, hl_calc=calc_hl, ll_calc=calc_ll, hl_xc="PBE", ll_xc="PBE", post_hf_method=None, basis="NAO-VCC-3Z")
run_binding_energy_test_full_molecule(pentanol_dimer, hl_calc=calc_hl, ll_calc=calc_ll,hl_xc="PBE", ll_xc="PBE", post_hf_method=None, basis="NAO-VCC-3Z")

