import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def get_germline_sequence(ab):
    '''
    Get full germline sequence. Note this includes some differences to 9114 and
    6261 that are ignored in the paper as they lie away from the paratope.
    Germlines are copied from figure in paper and extracted from screenshot (some Qs->O need correcting)
    https://brandfolder.com/workbench/extract-text-from-image

    :param ab: str with "9114" or "6261"
    :returns: str of germline H seq
    '''

    # the two germlines differ mainly in the H3
    if ab == "9114":
        germline = "QVQLVQSGAEVKKPGSSVKVSCKASGGTFSSYAISWVRQAPGQGLEWMGGIIPIFGTANYAQKFQGRVTITADKSTSTAYMELSSLRSEDTAVYYCARHGNYYYYYGMDVWGQGTTVTVSS"
    elif ab == "6261":
        germline = "QVQLVQSGAEVKKPGSSVKVSCKASGGTFSSYAISWVRQAPGQGLEWMGGIIPIFGTANYAQKFQGRVTITADKSTSTAYMELSSLRSEDTAVYYCAKHMGYQLRETMDVWGQGTTVTVSS"
    else:
        raise ValueError("ab must be either '9114' or 6261")
    
    return germline


def get_heavy_sequence(ab):
    '''
    Get Fv sequence of ab heavy chain

    :param ab: str with "9114" or "6261"
    :returns: str of H seq
    '''

    if ab == "9114":
        # taken from pdb 4fqi
        H_seq = "QVQLVQSGAEVKKPGSSVKVSCKSSGGTSNNYAISWVRQAPGQGLDWMGGISPIFGSTAYAQKFQGRVTISADIFSNTAYMELNSLTSEDTAVYFCARHGNYYYYSGMDVWGQGTTVTVSS"
    elif ab == "6261":
        # taken from pdb 3gbn
        H_seq = "EVQLVESGAEVKKPGSSVKVSCKASGGPFRSYAISWVRQAPGQGPEWMGGIIPIFGTTKYAPKFQGRVTITADDFAGTVYMELSSLRSEDTAMYYCAKHMGYQVRETMDVWGKGTTVTVSS"
    else:
        raise ValueError("ab must be either '9114' or 6261")
    
    return H_seq


def get_anarci_numbering():
    '''
    Hard coding to avoid anarci env dependency. Recreate using:
        raw_numbering = anarci.number(<seq>)
        [(str(res[0]) + res[1]).strip() for res, AA in raw_numbering[0] if AA != "-"]
    Note - this is the same for noth 6261 and 9114

    :returns: list of imgt numbering
    '''
    return ['1', '2', '3', '4', '5', '6', '7', '8', '9', '11', '12', '13', '14', '15', '16', '17',
            '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '35', '36',
            '37', '38', '39', '40', '41', '42', '43', '44', '45', '46', '47', '48', '49', '50', '51',
            '52', '53', '54', '55', '56', '57', '58', '59', '62', '63', '64', '65', '66', '67', '68',
            '69', '70', '71', '72', '74', '75', '76', '77', '78', '79', '80', '81', '82', '83', '84',
            '85', '86', '87', '88', '89', '90', '91', '92', '93', '94', '95', '96', '97', '98', '99',
            '100', '101', '102', '103', '104', '105', '106', '107', '108', '109', '110', '111',
            '112A','112', '113', '114', '115', '116', '117', '118', '119', '120', '121', '122',
            '123', '124', '125', '126', '127', '128']


def get_positions_mutated_from_germline(ab, ignore_same_residues_as_paper=True):
    '''
    Get imgt positions that are mutated from germline for given ab

    :param ab: str with "9114" or "6261"
    :param ignore_same_residues_as_paper: bool, ignore non-paratope mutations
    :returns: list of imgt numbers
    '''

    numbering = get_anarci_numbering()
    germline = get_germline_sequence(ab)
    somatic = get_heavy_sequence(ab)

    positions_mutated = []
    for aa_germline, aa_somatic, imgt in zip(germline, somatic, numbering):
        if aa_germline != aa_somatic:
            positions_mutated.append(imgt)

    positions_to_ignore = ["1", "6", "25", "50", "51", "101", "120"]
    if ignore_same_residues_as_paper:
        for position in positions_to_ignore:
            try:
                positions_mutated.remove(position)
            except ValueError:  # not the same in both abs
                pass

    return positions_mutated


def get_mutated_germline_vs_somatic_residues(ab, ignore_same_residues_as_paper=True):
    '''
    Get AAs that differ between germline and somatic

    :param ab: str with "9114" or "6261"
    :param ignore_same_residues_as_paper: bool, ignore non-paratope mutations
    :returns: two lists - first is AAs in germline, second is AAs in somatic ab
    '''

    numbering = get_anarci_numbering()
    germline = get_germline_sequence(ab)
    somatic = get_heavy_sequence(ab)
    
    positions_mutated = \
        get_positions_mutated_from_germline(ab, ignore_same_residues_as_paper=ignore_same_residues_as_paper)

    germline_AAs = [aa for position, aa in zip(numbering, germline) if position in positions_mutated]
    somatic_AAs = [aa for position, aa in zip(numbering, somatic) if position in positions_mutated]

    return germline_AAs, somatic_AAs


def load_raw_affinity_data(ab):
    '''
    Load raw -log(Kd) data for given ab

    :param ab: str with "9114" or "6261"
    :returns: df with affinity info
    '''
    if ab == "9114":
        date = "20210427"
        extra = ""
    elif ab == "6261":
        date = "20210323"
        extra = f"_{ab}"
    else:
        raise ValueError("ab must be either '9114' or 6261")
    
    # taken from bnab-landscapes/Figures/Figure1/20210329 - Figure 1_6261.ipynb
    df = pd.DataFrame(pd.read_csv(f'CR{ab}/Kd_meanbin/kd_processed/{date}{extra}_HA_unadj_fil_merg.csv',
                                  delimiter=',',
                                  dtype={'variant': str}))
    
    return df


def get_columns_of_interest(ab):
    '''
    9114 and 6261 bind different Ags and have different column names

    :param ab: str with "9114" or "6261"
    :returns: list of present column names of interest
    '''
    if ab == "9114":
        cols = ["variant", "h1_mean", "h3_mean", "fluB_mean", "som_mut"]
    elif ab == "6261":
        cols = ["variant", "h1_mean", "h9_mean", "som_mut"]
    else:
        raise ValueError("ab must be either '9114' or 6261")

    return cols


def binary_to_aa_str(binary, ab, ignore_same_residues_as_paper=True):
    '''
    Convert binary 10010... etc to matching amino acids

    :param binary: str of 1s and 0s that indicate germline vs somatic residues
    :param ab: str with "9114" or "6261"
    :param ignore_same_residues_as_paper: bool, ignore non-paratope mutations
    :returns: str of single letter amino acids for mutations only
    '''

    germ, som = get_mutated_germline_vs_somatic_residues(ab, ignore_same_residues_as_paper=True)

    assert ignore_same_residues_as_paper is True, "Must ignore residues to match length of binary"
    assert len(binary) == len(som), "lenght of binary str and list of somatic AAs must be the same"
   

    aa_str = ""
    for bit, germ_aa, som_aa in zip([b for b in binary], germ, som):
        if bit == "0":
            aa_str += germ_aa
        elif bit == "1":
            aa_str += som_aa
        else:
            raise ValueError("binary should contain string ones and zeros only")
        
    return aa_str


def binary_to_full_Hseq(binary, ab, ignore_same_residues_as_paper=True):
    '''
    Convert binary 10010... etc to matching amino acids and expand to whole H seq

    :param binary: str of 1s and 0s that indicate germline vs somatic residues
    :param ab: str with "9114" or "6261"
    :param ignore_same_residues_as_paper: bool, ignore non-paratope mutations
    :returns: str of single letter amino acids for whole H seq
    '''
    _, som = get_mutated_germline_vs_somatic_residues(ab, ignore_same_residues_as_paper=True)
    mutated_positions = get_positions_mutated_from_germline(ab, ignore_same_residues_as_paper=True)

    # note - ignoring any somatic mutations away from paratope by taking germline as base seq
    full_germline = get_germline_sequence(ab)
    full_numbering = get_anarci_numbering()

    assert ignore_same_residues_as_paper is True, "Must ignore residues to match length of binary"
    assert len(binary) == len(som), "lenght of binary str and list of somatic AAs must be the same"

    full_H_seq = ""
    for imgt, germ_aa in zip(full_numbering, full_germline):
        if imgt in mutated_positions:
            mut_idx = mutated_positions.index(imgt)
            if binary[mut_idx] == "1":
                full_H_seq += som[mut_idx]
            else:
                full_H_seq += germ_aa
        else:
            full_H_seq += germ_aa

    return full_H_seq


def load_processed_affinity_data(ab):
    '''
    Create df with only columns of interest, nicely ordered and named

    :param ab: str with "9114" or "6261"
    :returns: pretty df with affinity info
    '''
    raw_df = load_raw_affinity_data(ab)
    processed_df = raw_df[get_columns_of_interest(ab)]
    processed_df = processed_df.rename(columns={"variant": "binary_id"})

    processed_df["str_aa_id"] = processed_df.apply(lambda row: binary_to_aa_str(row["binary_id"], ab), axis=1)
    processed_df["H_seq"] = processed_df.apply(lambda row: binary_to_full_Hseq(row["binary_id"], ab), axis=1)

    max_mutations = 16 if ab == "9114" else 11
    processed_df["edit_distance"] = processed_df.apply(lambda row: max_mutations-row["som_mut"], axis=1)

    unordered_cols = list(processed_df.columns)
    ordered_cols = ["binary_id", "str_aa_id", "edit_distance", "H_seq"]
    new_cols = ordered_cols + list(set(unordered_cols)-set(ordered_cols))
    new_cols.remove("som_mut")

    return processed_df[new_cols]


def plot_affinity_by_edit_distance(ab, processed_df, target):
    '''
    Plot affinity of ab against a target, similar to paper

    :param processed_df: nicely processed df with edit distance
    :param target: antigen name e.g. "h9"
    '''

    sns.boxplot(processed_df, x="edit_distance", y=f"{target}_mean")
    plt.title(f"{ab} affinity against {target}")
    plt.xlabel(f"edit distance from {ab}")
    plt.ylabel("-log(Kd)")
    plt.show()
