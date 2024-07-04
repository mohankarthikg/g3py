import numpy as np
from scipy.interpolate import RegularGridInterpolator as RGI


class sig:

    def __init__(self):

        pass

    def create_SpatialEnergypdf(self, data, bins=[30, 100], range=None):

        bg_Hist_dec_log10Ne, bg_binedges = np.histogramdd(
            np.transpose(
                np.array(
                    [np.sin(np.radians(data2["EvDec1"])), np.log10(data2["NKGSize"])]
                )
            ),
            bins=bins,
            range=range,
            density=True,
        )

        bg_Hist_dec_log10Ne[bg_Hist_dec_log10Ne == 0] = np.finfo(float).eps

        logbg_spt_en_pdf = RGI(
            ((bgsindec_cnt), bglog10Ne_bins_cnt),
            np.log(bg_Hist_dec_log10Ne),
            method="linear",
            bounds_error=False,
            fill_value=None,
        )

        with open("bgspt_en_pdf.pkl", "wb") as f:

            pickle.dump(
                [
                    logbg_spt_en_pdf,
                    bg_binedges[0][0],
                    bg_binedges[0][-1],
                    bg_binedges[2][0],
                    bg_binedges[2][-1],
                ],
                f,
            )