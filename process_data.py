import os
import pandas as pd


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
