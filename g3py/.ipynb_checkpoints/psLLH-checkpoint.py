from iminuit import cost, Minuit
import numpy as np
from funcs import spaceAngle as sA
import CoordTransform as CT
from astropy.coordinates import SkyCoord
import astropy.units as u
import pickle
import IRF


class ps:
    def __init__(self, Srcdec, SrcRA, data):

        self.Srcdec = Srcdec
        self.SrcRA = SrcRA

        with open(f"tesrSrcAltAz_{self.Srcdec}_{self.SrcRA}.pkl", "rb") as file:
            self.SrcAz, self.SrcAlt = pickle.load(file)
        # self.SrcAz, self.SrcAlt =  CT.RADecToAltAz(np.degrees(self.SrcRA), np.degrees(self.Srcdec) ,\
        #                                            data['EvDate'], data['EvTime1'], data['EvTime2'])
        # with open('tesrSrcAltAz.pkl', 'wb') as file:
        #     pickle.dump([self.SrcAz, self.SrcAlt],file)

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

        mask1 = self.angres > 0
        mask2 = self.SrcAlt > 45
        mask = mask1  # *(mask2)
        self.angres = self.angres[mask]

        self.SrcAz = np.radians(self.SrcAz)[mask]
        self.SrcAlt = np.radians(self.SrcAlt)[mask]

        self.Srcsinalt = np.sin(self.SrcAlt)

        (
            self.bckg_spt_pdf,
            self.bckg_spt_cdf,
            self.bckg_spt_inv_cdf,
            self.sindec_min,
            self.sindec_max,
        ) = IRF.spt_pdf()

        (
            self.bckg_spt_pdf2,
            self.bckg_spt_pdf3,
            self.bckg_spt_pdf4,
            self.alt_min,
            self.alt_max,
            self.sinalt_min,
            self.sinalt_max,
        ) = IRF.spt_pdf2()

        (
            self.bckg_spt_en_pdf,
            self.sig_en_pdf,
            self.log10Ne_min,
            self.log10Ne_max,
            self.alt_min2,
            self.alt_max2,
            self.bg_norm,
        ) = IRF.spt_en_pdf()

        self.bckg_phi_pdf, self.phi_min, self.phi_max = IRF.spt_phi_pdf()

        self.evsindec = np.sin(np.radians(data["EvDec1"][mask]))
        self.evra = np.radians(data["EvRa1"][mask])

        self.evsinAlt = np.sin(np.radians(90 - data["Theta1"][mask]))
        self.evlog10Ne = np.log10(data["NKGSize"][mask])
        self.evAz = np.radians(data["Phi1"][mask])

        self.gam = 2.2

    @staticmethod
    def AltAzcache(Srcdec, SrcRA):
        SrcAz, SrcAlt = CT.RADecToAltAz(
            SrcRA, Srcdec, data["EvDate"], data["EvTime1"], data["EvTime2"]
        )
        with open(f"tesrSrcAltAz_{Srcdec}_{SrcRA}.pkl", "wb") as file:
            pickle.dump([SrcAz, SrcAlt], file)

    def ClassicPdf(self, x, f):

        Srcdec = self.Srcdec
        SrcRA = self.SrcRA
        sindec, RA = x
        # angres = np.radians(.83)

        # srcCord = SkyCoord(ra = SrcRA * u.rad, dec = Srcdec * u.rad, frame='icrs')
        # coords = SkyCoord(ra = RA * u.rad, dec = np.arcsin(sindec) * u.rad, frame='icrs')
        # spcAng = srcCord.separation(coords).radian

        spcAng = sA(np.arcsin(sindec), RA, Srcdec, SrcRA, rad=True)

        s = (
            (1 / (2 * np.pi * (self.angres**2)))
            * np.exp(-(spcAng**2) / (2 * self.angres**2))
            * (f)
        )

        b = (1 / (2 * np.pi)) * np.exp(self.bckg_spt_pdf(sindec)) * (1 - f)

        return s + b

    def SpatialEnergyPdf(self, x, f):

        Srcsinalt = self.Srcsinalt
        SrcAz = self.SrcAz

        evlog10Ne, evalt, evaz, Srcsinalt, SrcAz, angres = x

        spcAng = sA(evalt, evaz, np.arcsin(Srcsinalt), SrcAz, rad=True)

        s = (
            (1 / (2 * np.pi * (angres**2)))
            * np.exp(-(spcAng**2) / (2 * angres**2))
            * self.sig_en_pdf((self.gam, evlog10Ne, np.sin(evalt)))
            * f
        )

        b = (
            np.exp(self.bckg_phi_pdf(evaz))
            * self.bckg_spt_en_pdf((evlog10Ne, np.sin(evalt)))
            * (1 - f)
        )

        return s + b

    def SpatialEnergyPdf2(self, x, f):

        Srcsinalt = self.Srcsinalt
        SrcAz = self.SrcAz

        evlog10Ne, evalt, evaz, Srcsinalt, SrcAz, angres = x

        spcAng = sA(evalt, evaz, np.arcsin(Srcsinalt), SrcAz, rad=True)

        s = (
            (1 / (2 * np.pi * (angres**2)))
            * np.exp(-(spcAng**2) / (2 * angres**2))
            * self.sig_en_pdf((self.gam, evlog10Ne, np.sin(evalt)))
            * f
        )

        b = (
            np.exp(self.bckg_phi_pdf(evaz))
            * self.bckg_spt_en_pdf((evlog10Ne, np.sin(evalt)))
            * (1 - f)
        )

        return s + b

    def LocalSpatialPdf(self, x, f):

        evalt, evaz, Srcsinalt, SrcAz, angres = x

        spcAng = sA(evalt, evaz, np.arcsin(Srcsinalt), SrcAz, rad=True)
        s = (
            (1 / (2 * np.pi * (angres**2)))
            * np.exp(-(spcAng**2) / (2 * angres**2))
            * (np.degrees(np.arcsin(Srcsinalt)) > 45)
            * f
        )

        # b = (1/(2*np.pi))*np.exp(self.bckg_spt_pdf2(evalt))*(1-f)/np.cos(evalt)
        # b = (1/(2*np.pi))*self.bckg_spt_pdf4(np.sin(evalt))*(1-f)
        b = (
            np.exp(self.bckg_phi_pdf(evaz))
            * self.bckg_spt_pdf4(np.sin(evalt))
            * (1 - f)
        )

        return s + b

    def fun(self, i):

        np.random.seed(i)

        # Ntot = 1337857
        # sindec = np.random.uniform(self.sindec_min,self.sindec_max,2*Ntot)
        # sindec = sindec[ np.random.uniform(0,1.3,2*Ntot) < np.exp(self.bckg_pdf(sindec)) ]
        # ra = np.random.uniform(0,2*np.pi, sindec.shape[0])

        sindec = self.evsindec
        ra = np.random.permutation(self.evra)

        NLL = cost.UnbinnedNLL((sindec, ra), self.ClassicPdf)

        m = Minuit(NLL, f=np.random.uniform(-0.0002, 0.1))
        m.limits["f"] = (0, 1)
        m.tol = 10**-5
        m.errordef = 0.5
        m.precision = 2**-90
        m.print_level = 0
        m.scan(ncall=10)
        m.migrad()

        return (
            (NLL(0) - NLL(m.values["f"])),
            NLL(0),
            NLL(1),
            NLL(m.values["f"]),
            m.values["f"],
            m.valid,
            i,
        )

    def help(self, NLL, tol, method, init):

        m = Minuit(NLL, f=init)
        m.limits["f"] = (0, 1)  # (-.0001,1)
        m.tol = tol
        m.errordef = 0.5
        m.precision = 2**-100

        m.scipy(method=method)

        return m

    def fun2(self, i):

        np.random.seed(i)

        Ntot = self.Srcsinalt.shape[0]
        evlog10Ne = self.evlog10Ne
        # evlog10Ne = np.random.uniform(self.log10Ne_min, self.log10Ne_max, Ntot)
        evalt = np.arcsin(self.evsinAlt)

        # evlog10Ne = np.random.uniform(self.log10Ne_min, self.log10Ne_max, 3*Ntot)
        # evalt = np.random.uniform(self.alt_min2,self.alt_max2,3*Ntot)

        # mask = np.random.uniform(0,3,3*Ntot) <  np.exp(self.bckg_spt_en_pdf((evlog10Ne, evalt)))
        # evalt = evalt[mask]
        # evlog10Ne = evlog10Ne[mask]

        evaz = np.random.permutation(self.evAz)

        # evaz = np.random.uniform(0,2*np.pi, evalt.shape[0])
        # angres = np.radians(.83)*np.ones(evalt.shape[0])
        NLL = cost.UnbinnedNLL(
            (evlog10Ne, evalt, evaz, self.Srcsinalt, self.SrcAz, self.angres),
            self.SpatialEnergyPdf,
        )
        # NLL = cost.UnbinnedNLL((evlog10Ne, evalt, evaz, .85*np.ones(evalt.shape[0]), np.radians(180)*np.ones(evalt.shape[0]),\
        #                        np.radians(.83)*np.ones(evalt.shape[0])), self.SpatialEnergyPdf)

        m = self.help(NLL, tol=10**-5, method="Powell", init=0)

        if (NLL(0) - NLL(m.values["f"])) > 0 and m.valid:
            return (
                (NLL(0) - NLL(m.values["f"])),
                NLL(0),
                NLL(1),
                NLL(m.values["f"]),
                m.values["f"],
                m.valid,
                i,
                Ntot,
                m,
                NLL,
            )

        m = self.help(NLL, tol=10**-5, method="Powell", init=np.random.uniform(0, 0.1))

        if (NLL(0) - NLL(m.values["f"])) > 0 and m.valid:
            return (
                (NLL(0) - NLL(m.values["f"])),
                NLL(0),
                NLL(1),
                NLL(m.values["f"]),
                m.values["f"],
                m.valid,
                i,
                Ntot,
                m,
                NLL,
            )

        m = self.help(NLL, tol=10**-4, method="Powell", init=np.random.uniform(0, 0.1))

        if (NLL(0) - NLL(m.values["f"])) > 0 and m.valid:
            return (
                (NLL(0) - NLL(m.values["f"])),
                NLL(0),
                NLL(1),
                NLL(m.values["f"]),
                m.values["f"],
                m.valid,
                i,
                Ntot,
                m,
                NLL,
            )

        m = self.help(NLL, tol=10**-5, method="Nelder-Mead", init=0)

        if (NLL(0) - NLL(m.values["f"])) > 0 and m.valid:
            return (
                (NLL(0) - NLL(m.values["f"])),
                NLL(0),
                NLL(1),
                NLL(m.values["f"]),
                m.values["f"],
                m.valid,
                i,
                Ntot,
                m,
                NLL,
            )

        m = self.help(
            NLL, tol=10**-4, method="Nelder-Mead", init=np.random.uniform(0, 0.1)
        )

        if (NLL(0) - NLL(m.values["f"])) > 0 and m.valid:
            return (
                (NLL(0) - NLL(m.values["f"])),
                NLL(0),
                NLL(1),
                NLL(m.values["f"]),
                m.values["f"],
                m.valid,
                i,
                Ntot,
                m,
                NLL,
            )

        m.scan(ncall=10)
        return (
            (NLL(0) - NLL(m.values["f"])),
            NLL(0),
            NLL(1),
            NLL(m.values["f"]),
            m.values["f"],
            m.valid,
            i,
            Ntot,
            m,
            NLL,
        )

    def fun3(self, i):

        np.random.seed(i)

        Ntot = self.Srcsinalt.shape[0]  # 1337857

        evalt = np.arcsin(self.evsinAlt)

        # evalt = np.random.uniform(self.alt_min,self.alt_max,3*Ntot)
        # evalt = evalt[ np.random.uniform(0,9.5,3*Ntot) <  np.exp(self.bckg_spt_pdf2(evalt)) ]
        # evaz = np.random.uniform(0,2*np.pi, evalt.shape[0])

        evaz = np.random.permutation(self.evAz)
        # NLL = cost.UnbinnedNLL((evalt, evaz, .9*np.ones(evalt.shape[0]), np.radians(180)*np.ones(evalt.shape[0]),\
        #                      np.radians(.83)*np.ones(evalt.shape[0])), self.LocalSpatialPdf)

        NLL = cost.UnbinnedNLL(
            (evalt, evaz, self.Srcsinalt, self.SrcAz, self.angres), self.LocalSpatialPdf
        )
        # NLL = cost.UnbinnedNLL((evsinalt, evaz, self.Srcsinalt[:evsinalt.shape[0]], self.SrcAz[:evsinalt.shape[0]],\
        #                      self.angres[:evsinalt.shape[0]] ), self.LocalSpatialPdf)

        m = self.help(NLL, tol=10**-5, method="Powell", init=0)

        if (NLL(0) - NLL(m.values["f"])) > 0 and m.valid:
            return (
                (NLL(0) - NLL(m.values["f"])),
                NLL(0),
                NLL(1),
                NLL(m.values["f"]),
                m.values["f"],
                m.valid,
                i,
                Ntot,
                m,
                NLL,
            )

        m = self.help(NLL, tol=10**-5, method="Powell", init=np.random.uniform(0, 0.1))

        if (NLL(0) - NLL(m.values["f"])) > 0 and m.valid:
            return (
                (NLL(0) - NLL(m.values["f"])),
                NLL(0),
                NLL(1),
                NLL(m.values["f"]),
                m.values["f"],
                m.valid,
                i,
                Ntot,
                m,
                NLL,
            )

        m = self.help(NLL, tol=10**-5, method="Nelder-Mead", init=0)

        if (NLL(0) - NLL(m.values["f"])) > 0 and m.valid:
            return (
                (NLL(0) - NLL(m.values["f"])),
                NLL(0),
                NLL(1),
                NLL(m.values["f"]),
                m.values["f"],
                m.valid,
                i,
                Ntot,
                m,
                NLL,
            )

        # m = self.help(NLL, tol = 10**-5, method='COBYLA', init = 0 )

        # if (NLL(0)- NLL(m.values['f']))>0 and m.valid:
        #     return (NLL(0)- NLL(m.values['f'])), NLL(0), NLL(1), NLL(m.values['f']), m.values['f'], m.valid, i,  Ntot, m, NLL

        m = self.help(NLL, tol=10**-4, method="Powell", init=np.random.uniform(0, 0.1))

        if (NLL(0) - NLL(m.values["f"])) > 0 and m.valid:
            return (
                (NLL(0) - NLL(m.values["f"])),
                NLL(0),
                NLL(1),
                NLL(m.values["f"]),
                m.values["f"],
                m.valid,
                i,
                Ntot,
                m,
                NLL,
            )

        m = self.help(
            NLL, tol=10**-4, method="Nelder-Mead", init=np.random.uniform(0, 0.1)
        )

        if (NLL(0) - NLL(m.values["f"])) > 0 and m.valid:
            return (
                (NLL(0) - NLL(m.values["f"])),
                NLL(0),
                NLL(1),
                NLL(m.values["f"]),
                m.values["f"],
                m.valid,
                i,
                Ntot,
                m,
                NLL,
            )

        m.scan(ncall=50)
        return (
            (NLL(0) - NLL(m.values["f"])),
            NLL(0),
            NLL(1),
            NLL(m.values["f"]),
            m.values["f"],
            m.valid,
            i,
            Ntot,
            m,
            NLL,
        )

    def NLL(self, i):

        np.random.seed(3)

        Ntot = self.Srcsinalt.shape[0]  # 1337857
        evlog10Ne = self.evlog10Ne
        evalt = np.arcsin(self.evsinAlt)

        evaz = np.random.uniform(0, 2 * np.pi, evalt.shape[0])

        NLL_ = cost.UnbinnedNLL(
            (evlog10Ne, evalt, evaz, self.Srcsinalt, self.SrcAz, self.angres),
            self.SpatialEnergyPdf,
        )

        return NLL_(i)
