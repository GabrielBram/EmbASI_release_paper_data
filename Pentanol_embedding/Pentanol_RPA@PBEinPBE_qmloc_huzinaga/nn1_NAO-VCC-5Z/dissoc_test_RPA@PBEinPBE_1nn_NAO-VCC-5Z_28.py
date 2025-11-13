import os
from embasi.embedding import ProjectionEmbedding
from ase.calculators.aims import Aims, AimsProfile 
from ase.build import molecule
from ase.visualize import view
from ase.data.s22 import s22, s26, create_s22_system
from ase.io import read
import numpy as np
from csv_writer import write_csv_data_row, run_single_fragment, get_nn_index

os.environ['ASI_LIB_PATH'] = "/work/e05/e05/gabram/FHIaims/_build-gnu-lib/libaims.250403.scalapack.mpi.so"
os.environ['ASE_AIMS_COMMAND'] = "/work/e05/e05/gabram/FHIaims/_build-gnu-lib/libaims.250403.scalapack.mpi.so"

calc_ll = Aims(profile=AimsProfile(command=os.environ['ASE_AIMS_COMMAND']),
    xc='PBE0',
    output_level="full",
    relativistic='atomic_zora scalar',
    occupation_type="gaussian 0.01",
    mixer="pulay",
    n_max_pulay=10,
    KS_method="parallel",
    RI_method="LVL",
    collect_eigenvectors = True,
    density_update_method='density_matrix', # for DM export
    atomic_solver_xc="PBE",
    basis_threshold=1e-5,
    override_illconditioning=True,
    override_initial_charge_check=True,
    use_2d_corr=False
    )

calc_hl = Aims(profile=AimsProfile(command=os.environ['ASE_AIMS_COMMAND']),
    xc='PBE0',
    relativistic='atomic_zora scalar',
    occupation_type="gaussian 0.01",
    mixer="pulay",
    KS_method="parallel",
    RI_method="LVL",
    n_max_pulay=10,
    collect_eigenvectors=True,
    density_update_method='density_matrix', # for DM export
    lmo_pm_charge_metric="mulliken",
    atomic_solver_xc="PBE",
    basis_threshold=1e-5,
    override_illconditioning=True,
    override_initial_charge_check=True,
    frequency_points=50,
    use_2d_corr=False
  )

pentanol_dimer = read('pentanol_molecules.xyz', index=28)
region2 = get_nn_index(pentanol_dimer, 5, 1) + get_nn_index(pentanol_dimer, 23, 1)
embed_mask = [2]*len(pentanol_dimer)
for atom_idx in region2:
    embed_mask[atom_idx] = 1
os.environ['AIMS_SPECIES_DIR'] = "/work/e05/e05/gabram/FHIaims/species_defaults/NAO-VCC-nZ//NAO-VCC-5Z"
run_single_fragment(pentanol_dimer, embed_mask=embed_mask, hl_calc=calc_hl, ll_calc=calc_ll, hl_xc="PBE", ll_xc="PBE", post_hf_method="RPA", basis="NAO-VCC-5Z", run_dir="nn1_NAO-VCC-5Z")
