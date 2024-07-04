import numpy as np
from scipy.interpolate import RegularGridInterpolator as RGI


class bg:

    def __init__(self):

        pass

    def create_Spatialpdf(self, data, bins=50):

        hist, histbins = np.histogram(
            np.sin(np.radians(data2["EvDec1"])),
            cumulative=False,
            density=True,
            bins=bins,
        )

        log_bckgpdf = scipy.interpolate.InterpolatedUnivariateSpline(
            (histbins[1:] + histbins[:-1]) / 2, np.log(hist), k=1
        )

        with open("bg_spt_pdf.pkl", "wb") as file:
            pickle.dump([log_bckgpdf, histbins[0], histbins[-1]], file)

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


def get_Spatialpdf():

    with open("bg_spt_pdf.pkl", "rb") as file:
        log_bckgpdf, histbins_min, histbins_max = pickle.load(file)
    return log_bckgpdf, histbins_min, histbins_max


def get_SpatialEnergypdf():

    with open("bgspt_en_pdf.pkl", "rb") as file:
        (
            logbg_spt_en_pdf,
            bg_binedges_min0,
            bg_binedges_max0,
            bg_binedges_min1,
            bg_binedges_max1,
        ) = pickle.load(file)
    return (
        logbg_spt_en_pdf,
        bg_binedges_min0,
        bg_binedges_max0,
        bg_binedges_min1,
        bg_binedges_max1,
    )
