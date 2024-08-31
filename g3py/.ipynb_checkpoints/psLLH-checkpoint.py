from iminuit import cost, Minuit
import numpy as np
from .funcs import spaceAngle as sA
from . import CoordTransform as CT
import pickle
from . import IRF
import pandas as pd


class ps:
    """
    An point source Likelihood obj. which can generate TS distribution for
    bg. only hypothesis and sig. + bg. hypothesis.
    ...

    Attributes
    ----------

    Srcdec : float
        Declination of point src. in deg.

    Srcra : float
        RightAscension of point src. in deg.

    data : dict
        Dict. consisting of numpy arrays of evt. info. atleast EVDec, EvRA, NKGSize

    psrcsim : dict
        Dict. consisting of numpy arrays of simulated gamma evt. info. atleast EVDec, EvRA, NKGSize.

    n: int
        n is the average of poisson distrubtion of source events injected after drawing from psrcsim

    inj: bool
        If TRUE point source simulation will be injected before calculating TS.

    Methods
    -------

    """

    def __init__(self, Srcdec, SrcRA, data, psrcdata=None, n=0, spec_index = 2.5,inj = False):
        """
        Constructs all the necessary attributes for the point source Likelihood obj.

        Parameters
        ----------
        Srcdec : float
            Declination of point src. in deg.

        Srcra : float
            RightAscension of point src. in deg.

        data : dict
            Dict. consisting of numpy arrays of evt. info. atleast EVDec, EvRA, NKGSize

        psrcsim : dict
            Dict. consisting of numpy arrays of simulated gamma evt. info. atleast EVDec, EvRA, NKGSize.

        n: int
            n is the average of poisson distrubtion of source events injected after drawing from psrcsim
        
        inj: bool
            If TRUE point source simulation will be injected before calculating TS.
        """

        self.Srcdec = Srcdec
        self.SrcRA = SrcRA

        (
            self.bckg_spt_pdf,
            self.bckg_spt_cdf,
            self.bckg_spt_inv_cdf,
            self.sindec_min,
            self.sindec_max,
        ) = IRF.spt_pdf()
        (
            self.bckg_spt_en_pdf2,
            self.sig_en_pdf2,
            self.sindec_min2,
            self.sindec_max2,
            self.log10Ne_min2,
            self.log10Ne_max2,
        ) = IRF.spt_en_pdf2()

        self.angres = (
            sA(
                data["ThetaOdd1"],
                data["PhiOdd1"],
                data["ThetaEven1"],
                data["PhiEven1"],
                rad="False",
            )
            / 2
        )
        self.angres = np.radians(self.angres)

        mask1 = self.angres == self.angres
        mask2 = np.log10(data["NKGSize"]) < self.log10Ne_max2

        mask = (mask1) * (mask2)
        self.angres = self.angres[mask]

        self.data = {}
        
        self.data["evsindec"] = np.sin(np.radians(data["EvDec1"][mask]))
        self.data["evra"] = np.radians(data["EvRa1"][mask])
        self.data["evlog10Ne"] = np.log10(data["NKGSize"][mask])

        if psrcdata is not None:

            self.psrc = {}
            self.psrc["evsindec"] = np.sin(np.radians(psrcdata["EvDec1"]))
            self.psrc["evra"] = np.radians(psrcdata["EvRa1"])
            self.psrc["evlog10Ne"] = np.log10(psrcdata["NKGSize"])
            self.psrc["PrimaryEnergy"] = psrcdata["PrimaryEnergy"]


        self.inj = inj
        self.n = n
        self.spec_index = spec_index   

        # self.gam = 2.2

    def ClassicPdf(self, x, f):

        Srcdec = self.Srcdec
        SrcRA = self.SrcRA
        sindec, RA = x
        angres = np.radians(0.83)

        spcAng = sA(
            np.arcsin(sindec), RA, np.radians(Srcdec), np.radians(SrcRA), rad=True
        )

        s = (
            (1 / (2 * np.pi * (angres**2)))
            * np.exp(-(spcAng**2) / (2 * angres**2))
            * (f)
        )

        b = (1 / (2 * np.pi)) * np.exp(self.bckg_spt_pdf(sindec)) * (1 - f)

        return s + b

    def SpatialEnergyPdf_gamma(self, x, f, gam):

        Srcdec = self.Srcdec
        SrcRA = self.SrcRA

        evlog10Ne, evsindec, evra, evangres = x

        evangres = np.radians(0.83)  # evangres

        spcAng = sA(
            np.arcsin(evsindec), evra, np.radians(Srcdec), np.radians(SrcRA), rad=True
        )

        s = (
            (1 / (2 * np.pi * (evangres**2)))
            * np.exp(-(spcAng**2) / (2 * evangres**2))
            * np.exp(self.sig_en_pdf2((np.sin(np.radians(Srcdec)), evlog10Ne, gam)))
            * f
        )

        b = (
            (1 / (2 * np.pi))
            * np.exp(self.bckg_spt_en_pdf2((evsindec, evlog10Ne)))
            * (1 - f)
        )

        return s + b

    def TSClassicpdf(self, seed):
        """

        Calculates Test Statistic -2*ln( LLH(0)/LLH( ns, gamma) ).

        Parameters
        ----------

        seed: int
            seed for the scrambling of right ascension of background events (data)

        Returns
        -------

        TS: float

        NLL: float
            minimum value of Negative Log Likelihood

        f: float
            optimum value of source strength - ns, lies in [-1,1]

        valid: bool
            True / 1 if minimum of NLL satisfies iminuit valid minimum conditions.
        seed: int
            same as the seed set above
        Ntot: int
            total number of events

        method: int
            method used to minimize.
            1 - 'Minuit.migrad'
            2 - 'Powell'
            3 - 'Nelder-Mead'
            4 - 'L-BFGS-B'
            5 - 'SLSQP'
            6 - 'Minuit.scan'

            tries to minimize the NLL in this order if a method fails to minimize.
        """

        np.random.seed(seed)

        
        
        if self.inj:

            psrc = self.inject(seed, self.psrc)

            sindec = np.concatenate([self.data["evsindec"], psrc["evsindec"]])
            ra = np.concatenate([np.random.permutation(self.data["evra"]),  psrc["evra"]])
            
        else:
            sindec = self.data["evsindec"]
            ra = np.random.permutation(self.data["evra"])

        Ntot = sindec.shape[0]
            

        NLL = cost.UnbinnedNLL((sindec, ra), self.ClassicPdf)

        m = Minuit(NLL, f=np.random.uniform(0, 0.1))
        m.limits["f"] = (-0.0001, 1)
        m.tol = 10**-5
        m.precision = 2**-100
        m.scan(ncall=10)

        m.migrad()

        if (NLL(0) - NLL(m.values["f"])) > 0 and m.valid:
            return (
                (NLL(0) - NLL(m.values["f"])),
                NLL(0),
                NLL(1),
                NLL(m.values["f"]),
                m.values["f"],
                m.valid,
                seed,
                Ntot,
                1
            )  # , m, NLL

        m.scipy(method="Powell")

        if (NLL(0) - NLL(m.values["f"])) > 0 and m.valid:
            return (
                (NLL(0) - NLL(m.values["f"])),
                NLL(0),
                NLL(1),
                NLL(m.values["f"]),
                m.values["f"],
                m.valid,
                seed,
                Ntot,
                2
            )  # , m, NLL

        m.scipy(method="Nelder-Mead")
        if (NLL(0) - NLL(m.values["f"])) > 0 and m.valid:
            return (
                (NLL(0) - NLL(m.values["f"])),
                NLL(0),
                NLL(1),
                NLL(m.values["f"]),
                m.values["f"],
                m.valid,
                seed,
                Ntot,
                3
            )  # , m, NLL

        m.scipy(method="L-BFGS-B")

        if (NLL(0) - NLL(m.values["f"])) > 0 and m.valid:
            return (
                (NLL(0) - NLL(m.values["f"])),
                NLL(0),
                NLL(1),
                NLL(m.values["f"]),
                m.values["f"],
                m.valid,
                seed,
                Ntot,
                4
            )  # , m, NLL

        m.scipy(method="SLSQP")

        if (NLL(0) - NLL(m.values["f"])) > 0 and m.valid:
            return (
                (NLL(0) - NLL(m.values["f"])),
                NLL(0),
                NLL(1),
                NLL(m.values["f"]),
                m.values["f"],
                m.valid,
                seed,
                Ntot,
                5
            )  # , m, NLL

        m.scan(ncall=10)

        return (
            (NLL(0) - NLL(m.values["f"])),
            NLL(0),
            NLL(1),
            NLL(m.values["f"]),
            m.values["f"],
            m.valid,
            seed,
            Ntot,
            6
        )  # , m, NLL

    def TSSpatialEnergyPdf(self, seed):
        """

        Calculates Test Statistic -2*ln( LLH(0)/LLH( ns, gamma) ).

        Parameters
        ----------

        seed: int
            seed for the scrambling of right ascension of background events (data)

        Returns
        -------

        TS: float

        NLL: float
            minimum value of Negative Log Likelihood

        f: float
            optimum value of source strength - ns, lies in [-1,1]

        gam:
            optimum value of spectral index - gamma, lies in [0,4]

        valid: bool
            True / 1 if minimum of NLL satisfies iminuit valid minimum conditions.
        seed: int
            same as the seed set above
        Ntot: int
            total number of events

        method: int
            method used to minimize.
            1 - 'Minuit.migrad'
            2 - 'Powell'
            3 - 'Nelder-Mead'
            4 - 'L-BFGS-B'
            5 - 'SLSQP'
            6 - 'Minuit.scan'

            tries to minimize the NLL in this order if a method fails to minimize.
        """

        np.random.seed(seed)

        if self.inj:

            psrc = self.inject(seed, self.psrc)
            
            sindec = np.concatenate([self.data["evsindec"], psrc["evsindec"]])
            ra = np.concatenate([np.random.permutation(self.data["evra"]),  psrc["evra"]])
            log10Ne = np.concatenate([self.data["evlog10Ne"],  psrc["evlog10Ne"]])
            evangres = np.radians(.83)*np.ones_like(sindec) #self.angres
            
        else:
            sindec = self.data["evsindec"]
            ra = np.random.permutation(self.data["evra"])
            log10Ne = self.data["evlog10Ne"]
            evangres = np.radians(.83)*np.ones_like(sindec)

        Ntot = sindec.shape[0]

        NLL = cost.UnbinnedNLL(
            (log10Ne, sindec, ra, evangres), self.SpatialEnergyPdf_gamma
        )

        m = Minuit(NLL, f=np.random.uniform(0, 0.1), gam=np.random.uniform(1.5, 3.5))
        m.limits["f"] = (-0.0001, 1)
        m.limits["gam"] = (1, 4)
        m.tol = 10**-5
        m.precision = 2**-100
        m.scan(ncall=10)

        m.migrad()

        if (NLL(0, 2) - NLL(m.values["f"], m.values["gam"])) > 0 and m.valid:
            return (
                (NLL(0, 2) - NLL(m.values["f"], m.values["gam"])),
                NLL(m.values["f"], m.values["gam"]),
                m.values["f"],
                m.values["gam"],
                m.valid,
                seed,
                Ntot,
                1,
            )  # , m, NLL

        m.scipy(method="Powell")

        if (NLL(0, 2) - NLL(m.values["f"], m.values["gam"])) > 0 and m.valid:
            return (
                (NLL(0, 2) - NLL(m.values["f"], m.values["gam"])),
                NLL(m.values["f"], m.values["gam"]),
                m.values["f"],
                m.values["gam"],
                m.valid,
                seed,
                Ntot,
                2,
            )  # , m, NLL

        m.scipy(method="Nelder-Mead")
        if (NLL(0, 2) - NLL(m.values["f"], m.values["gam"])) > 0 and m.valid:
            return (
                (NLL(0, 2) - NLL(m.values["f"], m.values["gam"])),
                NLL(m.values["f"], m.values["gam"]),
                m.values["f"],
                m.values["gam"],
                m.valid,
                seed,
                Ntot,
                3,
            )  # , m, NLL

        m.scipy(method="L-BFGS-B")

        if (NLL(0, 2) - NLL(m.values["f"], m.values["gam"])) > 0 and m.valid:
            return (
                (NLL(0, 2) - NLL(m.values["f"], m.values["gam"])),
                NLL(m.values["f"], m.values["gam"]),
                m.values["f"],
                m.values["gam"],
                m.valid,
                seed,
                Ntot,
                4,
            )  # , m, NLL

        m.scipy(method="SLSQP")

        if (NLL(0, 2) - NLL(m.values["f"], m.values["gam"])) > 0 and m.valid:
            return (
                (NLL(0, 2) - NLL(m.values["f"], m.values["gam"])),
                NLL(m.values["f"], m.values["gam"]),
                m.values["f"],
                m.values["gam"],
                m.valid,
                seed,
                Ntot,
                5,
            )  # , m, NLL

        m.scan(ncall=10)

        return (
            (NLL(0, 2) - NLL(m.values["f"], m.values["gam"])),
            NLL(m.values["f"], m.values["gam"]),
            m.values["f"],
            m.values["gam"],
            m.valid,
            seed,
            Ntot,
            6,
        )  # , m, NLL

    def inject(self, seed, psrcdata):
        """

        Takes in Monte-Carlo generated point source events and randomly choose events
        following poisson distribution with average number of events = n.

        Parameters
        ----------

        data : dict
            contains numpy arrays of declination, right ascension, NKGSize etc.

        n : int
            average of poisson distribution

        Returns
        -------

        out: dict
            contains numpy arrays of same info. as input (data).
        """
        #np.random.seed(seed)
        NSrcev = np.random.poisson(self.n)
        #psrcdatadf = pd.DataFrame(psrcdata)
        
        p = psrcdata['PrimaryEnergy']**(2.5 - self.spec_index)/np.sum( psrcdata['PrimaryEnergy']**(2.5 - self.spec_index))
        Nevs = range(len(psrcdata[list(psrcdata.keys())[0]]))
        
        index = np.random.choice(Nevs, NSrcev, p = p) 

        out = {}
        for key in psrcdata:
            out[key] = psrcdata[key][index]

        #out =  psrcdatadf.sample(NSrcev, weights= psrcdatadf['PrimaryEnergy']**(2.5 - self.spec_index))

        return out
