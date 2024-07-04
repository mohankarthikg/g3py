import os
import numpy as np
import uproot as up
from functools import reduce
import pickle
from . import CoordTransform as CT
from .fiducial import check_within_fiducial
from .funcs import distance_based_cut
from . import g3dir


class Sim:
    """
    An obj which can get/handle Monte Carlo data after
    applying basic cuts on Monte Carlo data.

    ...

    Attributes
    ----------
    name : str
        first name of the person
    surname : str
        family name of the person
    age : int
        age of the person

    Methods
    -------
    info(additional=""):
        Prints the person's name and age.
    """

    def __init__(self):

        self.fiducial = True
        self.cut = {"NKGFitFlag": [0, 3], "Theta1": [0, 45], "Age": [0.2, 1.8]}
        self.distcut = True

        self._MuonFiles = None
        self._ShowerFiles = None
        self._ShowerPTFiles = None

        self._MuonVars = g3dir.MC_MuonVars
        self._ShowerVars = g3dir.MC_ShowerVars
        self._ShowerPTVars = g3dir.MC_ShowerPTVars

        # self._DetNo, self._detX , self._detY = self._getDetInfo()

        self.loadfiles()

        self._Cutvar = None
        self._infiducial = None
        self._distcut = None
        self._loadedvar = None

    def loadfiles(self, primary=["gamma"], ebins=[10, 11, 12, 13, 14, 15]):

        files1 = {}
        files2 = {}
        files3 = {}

        filepathsmuon = g3dir.MC_pathsmuon
        filepathsshower = g3dir.MC_pathsshower
        filepathsshowerPT = g3dir.MC_pathsshowerPT

        for prim in primary:

            if prim == "gamma":

                for binno in ebins:
                    for i in range(1, 201):
                        num = str(i).zfill(3)
                        if os.path.exists(
                            f"{filepathsmuon[prim]}bin_{binno}/MuTrack000{num}.root"
                        ):
                            files1[
                                f"{filepathsmuon[prim]}bin_{binno}/MuTrack000{num}.root"
                            ] = "mupar"
                        else:
                            print(
                                f"{filepathsmuon[prim]}bin_{binno}/MuTrack000{num}.root does not exist."
                            )

            else:
                for binno in ebins:
                    if os.path.exists(f"{filepathsmuon[prim]}MuTrack{binno}.root"):
                        files1[f"{filepathsmuon[prim]}MuTrack{binno}.root"] = "mupar"
                    else:
                        print(
                            f"File {filepathsmuon[prim]}MuTrack{binno}.root doesnt exist."
                        )

        for prim in primary:

            if prim == "gamma":

                for binno in ebins:
                    for i in range(1, 201):
                        num = str(i).zfill(3)
                        if os.path.exists(
                            f"{filepathsshower[prim]}bin_{binno}/CORSIKA_ANALYSIS_000{num}.root"
                        ):
                            files2[
                                f"{filepathsshower[prim]}bin_{binno}/CORSIKA_ANALYSIS_000{num}.root"
                            ] = "data"
                        else:
                            print(
                                f"{filepathsshower[prim]}bin_{binno}/CORSIKA_ANALYSIS_000{num}.root doesnt exist"
                            )

            else:
                for binno in ebins:
                    if os.path.exists(
                        f"{filepathsshower[prim]}CORSIKA_ANALYSIS_{binno}.root"
                    ):
                        files2[
                            f"{filepathsshower[prim]}CORSIKA_ANALYSIS_{binno}.root"
                        ] = "data"
                    else:
                        print(
                            f"File {filepathsshower[prim]}CORSIKA_ANALYSIS_{binno}.root doesnt exist."
                        )

        for prim in primary:

            if prim == "gamma":

                for binno in ebins:
                    for i in range(1, 201):
                        num = str(i).zfill(3)
                        if os.path.exists(
                            f"{filepathsshowerPT[prim]}bin_{binno}/CORSIKA_ANALYSIS_PT_000{num}.root"
                        ):
                            files3[
                                f"{filepathsshowerPT[prim]}bin_{binno}/CORSIKA_ANALYSIS_PT_000{num}.root"
                            ] = "pt"
                        else:
                            print(
                                f"{filepathsshowerPT[prim]}bin_{binno}/CORSIKA_ANALYSIS_PT_000{num}.root doesnt exist"
                            )

            else:
                for binno in ebins:
                    if os.path.exists(
                        f"{filepathsshowerPT[prim]}CORSIKA_ANALYSIS_PT_{binno}.root"
                    ):
                        files3[
                            f"{filepathsshowerPT[prim]}CORSIKA_ANALYSIS_PT_{binno}.root"
                        ] = "pt"
                    else:
                        print(
                            f"File {filepathsshowerPT[prim]}CORSIKA_ANALYSIS_PT_{binno}.root doesnt exist."
                        )

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
            "NKGSize",
            "NKGX",
            "NKGY",
            "Theta1",
            "Phi1",
        ],
        cache=True,
    ):

        if cache == True and self._Cutvar is not None:
            return self._Cutvar

        else:

            mucutvar2 = list(np.intersect1d(Cutvar, self._MuonVars))
            shcutvar2 = list(np.intersect1d(Cutvar, self._ShowerVars))
            shPTcutvar2 = list(np.intersect1d(Cutvar, self._ShowerPTVars))

            muCutvar, shCutvar, shPTCutvar = {}, {}, {}

            if len(mucutvar2) != 0:
                muCutvar = up.concatenate(
                    self._MuonFiles, mucutvar2, allow_missing=True, library="np"
                )
            if len(shcutvar2) != 0:
                shCutvar = up.concatenate(
                    self._ShowerFiles, shcutvar2, allow_missing=True, library="np"
                )
            if len(shPTcutvar2) != 0:
                shPTCutvar = up.concatenate(
                    self._ShowerPTFiles, shPTcutvar2, allow_missing=True, library="np"
                )

            self._Cutvar = muCutvar, shCutvar, shPTCutvar

        return muCutvar, shCutvar, shPTCutvar

    def cutconditions(self, muCutvar, shCutvar, shPTCutvar, cache=False):

        if self.fiducial == True or self.distcut == True:
            CoreX = shCutvar["NKGX"]
            CoreY = shCutvar["NKGY"]

        if self.distcut == True:
            log10NKGSize = np.log10(shCutvar["NKGSize"])

        if self.fiducial == True:
            if cache == True and self._infiducial is not None:
                infiducial = self._infiducial

            else:
                infiducial = check_within_fiducial(CoreX, CoreY)
                self._infiducial = infiducial

        if self.distcut == True:

            if cache == True and self._distcut is not None:
                indist_based_cut = self._distcut
            else:
                indist_based_cut = distance_based_cut(CoreX, CoreY, log10NKGSize)
                self._distcut = indist_based_cut

        temp = []

        for key in self.cut:

            if key in self._MuonVars:
                low = muCutvar[key] > self.cut[key][0]
                high = muCutvar[key] < self.cut[key][1]
                temp.append(np.all([low, high], axis=0))

            elif key in self._ShowerVars:
                low = shCutvar[key] > self.cut[key][0]
                high = shCutvar[key] < self.cut[key][1]
                temp.append(np.all([low, high], axis=0))

            elif key in self._ShowerPTVars:
                low = shPTCutvar[key] > self.cut[key][0]
                high = shPTCutvar[key] < self.cut[key][1]
                temp.append(np.all([low, high], axis=0))

            if self.fiducial == True:
                temp.append(infiducial)
            if self.distcut == True:
                temp.append(indist_based_cut)
        if temp != []:
            condition = np.all(temp, axis=0)
        else:
            return None
        return condition

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
            "PrimaryEnergy",
            "CoreX",
            "CoreY",
            "Density",
            "DetNo",
            "nHitDetectors",
        ],
        cache=True,
    ):

        if cache == True and self._loadedvar is not None:

            return self._loadedvar

        else:
            datamu = {}
            datash = {}
            datashPT = {}

            muvars = list(np.intersect1d(loadvar, self._MuonVars))
            shvars = list(np.intersect1d(loadvar, self._ShowerVars))
            shPTvars = list(np.intersect1d(loadvar, self._ShowerPTVars))

            if len(muvars) != 0:
                datamu = up.concatenate(
                    self._MuonFiles, muvars, allow_missing=True, library="np"
                )

            if len(shvars) != 0:
                datash = up.concatenate(
                    self._ShowerFiles, shvars, allow_missing=True, library="np"
                )

            if len(shPTvars) != 0:
                datashPT = up.concatenate(
                    self._ShowerPTFiles, shPTvars, allow_missing=True, library="np"
                )

            self._loadedvar = datamu, datash, datashPT
            return datamu, datash, datashPT

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
            "PrimaryEnergy",
            "CoreX",
            "CoreY",
            "Density",
            "DetNo",
            "nHitDetectors",
        ],
        cache=True,
    ):

        union = reduce(
            np.union1d, (self._MuonVars, self._ShowerVars, self._ShowerPTVars)
        )

        if not np.all(np.in1d(var, union)) or var == []:
            print(f"Some variables in {var} doesnt exist.")
            return {}

        elif len(self._MuonFiles) == 0 or len(self._MuonFiles) == 0:

            out = {}
            for key in var:
                out[key] = []

            return out

        else:

            var2 = []
            for key in self.cut:
                var2.append(key)

            if self.fiducial == True or self.distcut == True:
                var2.append("NKGX")
                var2.append("NKGY")

            if self.distcut == True:
                var2.append("NKGSize")

            muCutvar, shCutvar, shPTCutvar = self.loadCutvar(var2, cache=cache)

            try:
                condition = self.cutconditions(
                    muCutvar, shCutvar, shPTCutvar, cache=cache
                )

            except:
                out = {}
                for key in var:
                    out[key] = []
                return out

            datamu, datash, datashPT = self.loadvar(loadvar=var, cache=cache)

            out = {}

            for data in [datamu, datash, datashPT]:
                for key in data:
                    if condition is not None:
                        out[key] = data[key][condition]

                    else:
                        out[key] = data[key]

            return out
