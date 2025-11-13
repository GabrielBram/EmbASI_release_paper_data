import os
from embasi.embedding import ProjectionEmbedding, StandardDFT
from ase.calculators.aims import Aims, AimsProfile
from ase.build import molecule
from ase.visualize import view
from ase.data.s22 import s22, s26, create_s22_system
import numpy as np

def get_nn_index(atoms, root, nn):

    from carmm.analyse.neighbours import neighbours

    neighb_list, junk = neighbours(atoms, [root], nn)

    backbone_list = [at for at in neighb_list if atoms[at].symbol in ["C","O"]]
    H_list = [at.index for at in atoms if at.symbol == "H" ]

    nn_idx_list = []
    for bb_idx in backbone_list:
        bb_neighb, junk = neighbours(atoms, [bb_idx], 1)

        H_local_list = list(set(bb_neighb).intersection(set(H_list)))
        nn_idx_list = list(set(nn_idx_list).union(set(H_local_list)))

    nn_idx_list = list(set(nn_idx_list).union(set(backbone_list)))

    return nn_idx_list

def write_csv_data_row(filename, datadict):
    import os
    import csv
    from mpi4py import MPI

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    if (rank == 0):
        field_names = [ key for key in datadict.keys() ]
        file_exist = os.path.exists(filename)

        with open(filename, 'a') as csvfile:

            writer = csv.DictWriter(csvfile, fieldnames = field_names)

            if not file_exist:
                writer.writeheader()

            writer.writerow(datadict)

os.environ['ASE_AIMS_COMMAND'] = "/work/e05/e05/gabram/FHIaims/_build-gnu_lib_21102025/libaims.251014.scalapack.mpi.so"

calc_ll = Aims(xc='PBE',
    profile=AimsProfile(command=os.environ['ASE_AIMS_COMMAND']),
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
    override_initial_charge_check=True,
    )

calc_hl = Aims(xc='PBE',
    profile=AimsProfile(command=os.environ['ASE_AIMS_COMMAND']),
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
    override_initial_charge_check=True,
  )

def run_single_fragment_full_molecule(atoms, hl_calc=calc_hl, ll_calc=calc_ll, hl_xc="pbe", ll_xc="pbe", post_hf_method=None, basis="N/A", run_dir="./Embasi_calc", frag_name="Default"):

    data_dict_list = []

    hl_calc.parameters["xc"] = hl_xc
    ll_calc.parameters["xc"] = ll_xc

    if post_hf_method is not None:
        ll_calc.parameters["total_energy_method"] = post_hf_method

    dist = atoms.info["dist"]
    run_dir = os.path.join(run_dir, str(round(dist, ndigits=2)))

    Projection = StandardDFT(atoms, embed_mask=None, calc_base_ll=ll_calc, calc_base_hl=hl_calc, run_dir=run_dir)
    Projection.run()

    data_dict = {}
    data_dict["Frag Name"] = frag_name
    data_dict["Trunc"] = "N/A"
    data_dict["LL XC"] = ll_calc.parameters["xc"]
    data_dict["HL XC"] = hl_calc.parameters["xc"]
    if post_hf_method is not None:
        data_dict["Post-HF Method"] = post_hf_method
    else:
        data_dict["Post-HF Method"] = "N/A"
    data_dict["Basis"] = basis
    data_dict["Total DFT Energy"] = Projection.AB_LL.total_energy
    data_dict["Total Time"] = Projection.time_tot

    write_csv_data_row("Values.csv", data_dict)


def run_single_fragment(atoms, embed_mask=None, hl_calc=calc_hl, ll_calc=calc_ll, hl_xc="pbe", ll_xc="pbe", post_hf_method=None, basis="N/A", trunc=None, run_dir="./EmbASI_calc", frag_name="Default"):

    hl_calc.parameters["xc"] = hl_xc
    ll_calc.parameters["xc"] = ll_xc

    if post_hf_method is not None:
        hl_calc.parameters["total_energy_method"] = "rpa"

    if trunc is None:
        run_dir = os.path.join(run_dir, "NoTrunc")
        Projection = ProjectionEmbedding(atoms, embed_mask=embed_mask, calc_base_ll=ll_calc, calc_base_hl=hl_calc, mu_val=1.e+6, gc=True, parallel=True, projection="huzinaga", run_dir=run_dir, localisation="SPADE")
    else:
        run_dir = os.path.join(run_dir, str(trunc))
        Projection = ProjectionEmbedding(atoms, embed_mask=embed_mask, calc_base_ll=ll_calc, calc_base_hl=hl_calc, mu_val=1.e+6, gc=True, parallel=True, projection="huzinaga", run_dir=run_dir, localisation="SPADE", truncate_basis_thresh=trunc)

    Projection.run()

    data_dict={}
    data_dict["Frag Name"] = frag_name
    if trunc is None:
        data_dict["Trunc"] = "N/A"
    else:
        data_dict["Trunc"] = trunc
    data_dict["LL XC"] = ll_calc.parameters["xc"]
    data_dict["HL XC"] = hl_calc.parameters["xc"]
    data_dict["HL Atoms"] = embed_mask.count(1)
    if post_hf_method is not None:
        data_dict["Post-HF Method"] = post_hf_method
    else:
        data_dict["Post-HF Method"] = "N/A"
    data_dict["Basis"] = basis
    data_dict["Untrunc Nbasis"] = Projection.basis_info.full_natoms
    data_dict["Trunc Nbasis"] = Projection.basis_info.trunc_natoms
    data_dict["Untrunc Natoms"] = Projection.basis_info.full_nbasis
    data_dict["Trunc Natoms"] = Projection.basis_info.trunc_nbasis
    data_dict["Total DFT Energy (A-in-B)"] = Projection.DFT_AinB_total_energy
    data_dict["AB LL Energy"] = Projection.subsys_AB_lowlvl_totalen
    data_dict["A HL Energy"] = Projection.subsys_A_highlvl_totalen
    data_dict["A LL Energy"] = Projection.subsys_A_lowlvl_totalen
    data_dict["PB Correction"] = Projection.PB_corr
    data_dict["AB_LL Time"] = Projection.time_ab_lowlevel
    data_dict["A_LL Time"] = Projection.time_a_lowlevel
    data_dict["A_HL Time"] = Projection.time_a_highlevel

    write_csv_data_row("Values.csv", data_dict)
