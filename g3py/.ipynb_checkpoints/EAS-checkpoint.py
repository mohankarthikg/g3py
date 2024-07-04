from Sim import *


class EAS(Sim):

    def __init__(self):

        super().__init__()

        self._MuonVars = g3dir.Data_MuonVars
        self._ShowerVars = g3dir.Data_ShowerVars
        self._ShowerPTVars = g3dir.Data_ShowerPTVars

        self.cut = {
            "NKGFitFlag": [0, 3],
            "Theta1": [0, 45],
            "Age": [0.2, 1.8],
            "MuStatus": [0, 2],
            "nMuonLarge": [-1, 1],
        }

        self._Cutvar = None
        self._infiducial = None
        self._distcut = None
        self._loadedvar = None

    def loadfiles(
        self,
        date=[
            "20220110",
            "20220210",
            "20220310",
            "20220410",
            "20220510",
            "20220610",
            "20220710",
            "20220810",
            "20220910",
            "20221010",
            "20221110",
            "20221210",
        ],
    ):

        if isinstance(date, str):
            dates = [date]

        else:
            dates = date

        files1 = {}
        files2 = {}
        files3 = {}
        for date in dates:

            year = date[:4]
            self._muonDir = (
                "/boson/users/common/FINAL_EAS/" + str(year) + "/" + "MUPAR/"
            )
            self._showerDir = (
                "/boson/users/common/FINAL_EAS/" + str(year) + "/" + "SHPAR/STAGE2/"
            )
            self._showerPTDir = "/boson/users/mohan/DataPT/"

            if os.path.exists(f"{self._muonDir}mupar{date}.root") and os.path.exists(
                f"{self._showerDir}sh{date}.root"
            ):

                files1[f"{self._muonDir}mupar{date}.root"] = "mupar"
                files2[f"{self._showerDir}sh{date}.root"] = "sh"
                files3[f"{self._showerPTDir}pt{date}.root"] = "pt"

        self._MuonFiles = files1
        self._ShowerFiles = files2
        self._ShowerPTFiles = files3

    def loadCutvar(
        self,
        Cutvar=[
            "nMuonLarge",
            "nMuonSmall",
            "Age",
            "NKGFitFlag",
            "MuStatus",
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
            "MuStatus",
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
            "Density",
            "DetNo",
            "nHitDetectors",
        ],
        cache=True,
    ):

        return super().getvar(var=var, cache=cache)
