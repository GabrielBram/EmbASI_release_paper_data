import os
from asiembedding.embedding import ProjectionEmbedding, StandardDFT
from ase.calculators.aims import Aims
from ase.build import molecule
from ase.visualize import view
from ase.data.s22 import s22, s26, create_s22_system
import numpy as np

def get_nn_index(atoms, root, nn):

    from carmm.analyse.neighbours import neighbours

    neighb_list, junk = neighbours(atoms, [root], nn)

    backbone_list = [at for at in neighb_list if atoms[at].symbol in ["C","O"]]
    H_list = [at.index for at in pentanol if at.symbol == "H" ]

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

def run_single_fragment_full_molecule(atoms, hl_calc=calc_hl, ll_calc=calc_ll, hl_xc="pbe", ll_xc="pbe", post_hf_method=None, basis="N/A"):

    data_dict_list = []

    hl_calc.parameters["xc"] = hl_xc
    ll_calc.parameters["xc"] = ll_xc

    if post_hf_method is not None:
        ll_calc.parameters["total_energy_method"] = post_hf_method
    dist = 0.

    Projection = StandardDFT(atoms, embed_mask=None, calc_base_ll=ll_calc, calc_base_hl=hl_calc)
    Projection.run()

    data_dict = {}
    data_dict["LL XC"] = ll_calc.parameters["xc"]
    data_dict["HL XC"] = hl_calc.parameters["xc"]
    if post_hf_method is not None:
        data_dict["Post-HF Method"] = post_hf_method
    else:
        data_dict["Post-HF Method"] = "N/A"
    data_dict["Basis"] = basis
    data_dict["OH-H Distance"] = 0.
    data_dict["Total DFT Energy"] = Projection.AB_LL.total_energy
    data_dict["Total Time"] = Projection.time_tot

    write_csv_data_row("Values.csv", data_dict)
    del Projection
    import gc
    gc.collect()

def run_binding_energy_test_full_molecule(atoms, hl_calc=calc_hl, ll_calc=calc_ll, hl_xc="pbe", ll_xc="pbe", post_hf_method=None, basis="N/A"):

    hl_calc.parameters["xc"] = hl_xc
    ll_calc.parameters["xc"] = ll_xc

    if post_hf_method is not None:
        ll_calc.parameters["total_energy_method"] = post_hf_method

    atoms_old = atoms.copy()
    for dist_offset in np.arange(-0.5,1.5,0.11):
        H_bond_vec = atoms.positions[5] - atoms.positions[35]
        H_bond_vec_norm = H_bond_vec/np.linalg.norm(H_bond_vec)

        atoms.positions[:18] = atoms.positions[:18] + ( dist_offset * H_bond_vec_norm )
        dist = np.linalg.norm(atoms.positions[35] - atoms.positions[5])

        Projection = StandardDFT(atoms, embed_mask=None, calc_base_ll=ll_calc, calc_base_hl=hl_calc)
        Projection.run()

        data_dict = {}
        data_dict["LL XC"] = ll_calc.parameters["xc"]
        data_dict["HL XC"] = hl_calc.parameters["xc"]
        if post_hf_method is not None:
            data_dict["Post-HF Method"] = post_hf_method
        else:
            data_dict["Post-HF Method"] = "N/A"
        data_dict["Basis"] = basis
        data_dict["OH-H Distance"] = dist
        data_dict["Total DFT Energy"] = Projection.AB_LL.total_energy
        data_dict["Total Time"] = Projection.time_tot

        atoms.positions = atoms_old.positions
        write_csv_data_row("Values.csv", data_dict)
        del Projection
        import gc
        gc.collect()


def run_binding_energy_test(atoms, embed_mask=None, frag_charge=0, hl_calc=calc_hl, ll_calc=calc_ll, hl_xc="pbe", ll_xc="pbe", post_hf_method=None, basis="N/A"):

    hl_calc.parameters["xc"] = hl_xc
    ll_calc.parameters["xc"] = ll_xc

    if post_hf_method is not None:
        hl_calc.parameters["total_energy_method"] = "rpa"

    atoms_old = atoms.copy()
    for dist_offset in np.arange(-0.5,1.5,0.11):
        H_bond_vec = atoms.positions[5] - atoms.positions[35]
        H_bond_vec_norm = H_bond_vec/np.linalg.norm(H_bond_vec)

        atoms.positions[:18] = atoms.positions[:18] + ( dist_offset * H_bond_vec_norm )
        dist = np.linalg.norm(atoms.positions[5] - atoms.positions[35])

        data_dict = run_single_fragment(atoms, embed_mask=embed_mask, frag_charge=frag_charge, hl_calc=hl_calc, ll_calc=ll_calc, hl_xc=hl_xc, ll_xc=ll_xc, post_hf_method=post_hf_method, basis=basis)
        data_dict["OH-H Distance"] = dist

        write_csv_data_row("Values.csv", data_dict)

        atoms.positions = atoms_old.positions

def run_single_fragment(atoms, embed_mask=None, frag_charge=0, hl_calc=calc_hl, ll_calc=calc_ll, hl_xc="pbe", ll_xc="pbe", post_hf_method=None, basis="N/A"):

    hl_calc.parameters["xc"] = hl_xc
    ll_calc.parameters["xc"] = ll_xc

    dist = 0.

    if post_hf_method is not None:
        hl_calc.parameters["total_energy_method"] = "rpa"

    Projection = ProjectionEmbedding(atoms, embed_mask=embed_mask, calc_base_ll=ll_calc, calc_base_hl=hl_calc, frag_charge=frag_charge, mu_val=1.e+6)
    Projection.run()

    data_dict={}
    data_dict["LL XC"] = ll_calc.parameters["xc"]
    data_dict["HL XC"] = hl_calc.parameters["xc"]
    if post_hf_method is not None:
        data_dict["Post-HF Method"] = post_hf_method
    else:
        data_dict["Post-HF Method"] = "N/A"
    data_dict["Basis"] = basis
    data_dict["OH-H Distance"] = dist
    data_dict["Total DFT Energy (A-in-B)"] = Projection.DFT_AinB_total_energy
    data_dict["AB LL Energy"] = Projection.AB_LL_PP.ev_corr_total_energy
    data_dict["A HL Energy"] = Projection.A_HL_PP.ev_corr_total_energy
    data_dict["A LL Energy"] = Projection.A_LL.ev_corr_total_energy
    data_dict["PB Correction"] = Projection.PB_corr
    data_dict["AB_LL Time"] = Projection.time_ab_lowlevel
    data_dict["A_LL Time"] = Projection.time_a_lowlevel
    data_dict["A_HL Time"] = Projection.time_a_highlevel
    data_dict["A_LL_PP Time"] = Projection.time_a_lowlevel_pp
    data_dict["A_HL_PP Time"] = Projection.time_a_highlevel_pp
    data_dict["AB_LL_PP Time"] = Projection.time_ab_lowlevel_pp

    del Projection
    import gc
    gc.collect()

    return data_dict

def run_single_fragment_and_print(atoms, embed_mask=None, frag_charge=0, hl_calc=calc_hl, ll_calc=calc_ll, hl_xc="pbe", ll_xc="pbe", post_hf_method=None, basis="N/A"):

    hl_calc.parameters["xc"] = hl_xc
    ll_calc.parameters["xc"] = ll_xc

    dist = 0.

    if post_hf_method is not None:
        hl_calc.parameters["total_energy_method"] = "rpa"

    Projection = ProjectionEmbedding(atoms, embed_mask=embed_mask, calc_base_ll=ll_calc, calc_base_hl=hl_calc, frag_charge=frag_charge, mu_val=1.e+6)
    Projection.run()

    data_dict={}
    data_dict["LL XC"] = ll_calc.parameters["xc"]
    data_dict["HL XC"] = hl_calc.parameters["xc"]
    if post_hf_method is not None:
        data_dict["Post-HF Method"] = post_hf_method
    else:
        data_dict["Post-HF Method"] = "N/A"
    data_dict["Basis"] = basis
    data_dict["OH-H Distance"] = dist
    data_dict["Total DFT Energy (A-in-B)"] = Projection.DFT_AinB_total_energy
    data_dict["AB LL Energy"] = Projection.AB_LL_PP.ev_corr_total_energy
    data_dict["A HL Energy"] = Projection.A_HL_PP.ev_corr_total_energy
    data_dict["A LL Energy"] = Projection.A_LL.ev_corr_total_energy
    data_dict["PB Correction"] = Projection.PB_corr
    data_dict["AB_LL Time"] = Projection.time_ab_lowlevel
    data_dict["A_LL Time"] = Projection.time_a_lowlevel
    data_dict["A_HL Time"] = Projection.time_a_highlevel
    data_dict["A_LL_PP Time"] = Projection.time_a_lowlevel_pp
    data_dict["A_HL_PP Time"] = Projection.time_a_highlevel_pp
    data_dict["AB_LL_PP Time"] = Projection.time_ab_lowlevel_pp

    del Projection
    import gc
    gc.collect()

    write_csv_data_row("Values.csv", data_dict)
