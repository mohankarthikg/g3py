from Sim import *


class PSrcSim(Sim):

    def __init__(self):

        super().__init__()

        self._MuonVars = g3dir.MC_PSim_MuonVars

        self._ShowerVars = g3dir.MC_PSim_ShowerVars
        self._ShowerPTVars = g3dir.MC_PSim_ShowerPTVars

        self.cut = {
            "NKGFitFlag": [0, 3],
            "Theta1": [0, 45],
            "Age": [0.2, 1.8],
            "nMuonLarge": [-1, 1],
        }

        self._Cutvar = None
        self._infiducial = None
        self._distcut = None
        self._loadedvar = None

    def loadfiles(self):

        files1 = {}
        files2 = {}
        files3 = {}

        CA = glob.glob("/boson/users/mohan/PSim/RECO/CORSIKA_ANALYSIS_00????.root")
        CA = sorted(CA)
        PT = glob.glob("/boson/users/mohan/PSim/RECO/CORSIKA_ANALYSIS_PT_00????.root")
        PT = sorted(PT)
        MU = glob.glob("/boson/users/mohan/PSim/RECO/MuTrack00????.root")
        MU = sorted(MU)

        self._MuonFiles = {key: "mupar" for key in MU}
        self._ShowerFiles = {key: "cor" for key in CA}
        self._ShowerPTFiles = {key: "pt" for key in PT}

    def loadCutvar(
        self,
        Cutvar=[
            "nMuonLarge",
            "nMuonSmall",
            "Age",
            "NKGFitFlag",
            "NKGSize",
            "NKGX",
            "NKGY",
            "Theta1",
            "Phi1",
        ],
        cache=True,
    ):

        return super().loadCutvar(Cutvar=Cutvar, cache=cache)

    def loadvar(
        self,
        loadvar=[
            "nMuonLarge",
            "nMuonSmall",
            "Age",
            "NKGFitFlag",
            "NKGSize",
            "NKGX",
            "NKGY",
            "Theta1",
            "Phi1",
            "ThetaOdd1",
            "ThetaEven1",
            "Theta6",
            "Phi6",
            "ThetaOdd6",
            "ThetaEven6",
            "Density",
            "DetNo",
            "nHitDetectors",
        ],
        cache=True,
    ):

        return super().loadvar(loadvar=loadvar, cache=cache)

    def getvar(
        self,
        var=[
            "nMuonLarge",
            "nMuonSmall",
            "Age",
            "NKGFitFlag",
            "NKGSize",
            "NKGX",
            "NKGY",
            "Theta1",
            "Phi1",
            "ThetaOdd1",
            "ThetaEven1",
            "Theta6",
            "Phi6",
            "ThetaOdd6",
            "ThetaEven6",
            "Density",
            "DetNo",
            "nHitDetectors",
        ],
        cache=True,
    ):

        return super().getvar(var=var, cache=cache)


