import os
import pandas as pd


def get_germline_sequence():
    '''
    Get full germline sequence. Note this includes some differences to 9114 and 6261
    that are ignored as they lie away from the paratope

    :returns: str of germline H seq
    '''
    return "QVQLVQSGAEVKKPGSSVKVSCKASGGTFSSYAISWVRQAPGQGLEWMGGIIPIFGTANYAQKFQGRVTITADKSTSTAYMELSSLRSEDTAVYYCARHGNYYYYYGMDVWGQGTTVTVSS"


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
