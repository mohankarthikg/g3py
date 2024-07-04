import os
import pickle
import numpy as np
from . import CoordTransform as CT
from .EAS import *
from .funcs import generate_dates
from .g3dir import disk

class FullEASdata:

    def __init__(self):
        self.var = [
            "EvDate",
            "EvTime1",
            "EvTime2",
            "NKGSize",
            "Theta6",
            "Phi6",
            "Theta1",
            "Phi1",
            "ThetaOdd1",
            "PhiOdd1",
            "ThetaEven1",
            "PhiEven1",
            "ThetaOdd6",
            "PhiOdd6",
            "ThetaEven6",
            "PhiEven6",
        ]

    def _help(self, date="20140102"):

        # print(date)
        a = EAS()
        a.loadfiles(date)
        a.cut = {
            "NKGFitFlag": [0, 3],
            "Theta1": [0, 45],
            "Age": [0.2, 1.8],
            "MuStatus": [0, 2],
        }
        a.distcut = True

        # fp = open(f'{date}.pkl','wb')
        # pickle.dump(a.getvar(self.var),fp)

        return a.getvar(self.var)

    def load(
        self,
        var=[
            "EvDate",
            "EvTime1",
            "EvTime2",
            "NKGSize",
            "Theta6",
            "Phi6",
            "Theta1",
            "Phi1",
            "ThetaOdd1",
            "PhiOdd1",
            "ThetaEven1",
            "PhiEven1",
            "ThetaOdd6",
            "PhiOdd6",
            "ThetaEven6",
            "PhiEven6",
        ],
        year=2013,
        cache=True,
    ):

        if cache == True and os.path.exists(
            f"/{disk}/users/mohan/G3PY/CacheDir/FullEASdata{year}.pkl"
        ):
            file = open(
                f"/{disk}/users/mohan/G3PY/CacheDir/FullEASdata{year}.pkl", "rb"
            )
            cacheout, cachevar = pickle.load(file)
            if np.all(np.in1d(var, cachevar)):
                return cacheout

        else:

            self.var = var
            from multiprocessing import Pool, cpu_count

            p = Pool(cpu_count())

            dates = generate_dates(f"{year}0101", f"{year}1231")

            outdicts = p.map(self._help, dates)

            p.close()
            p.join()

            output = {}
            for key in self.var:
                output[key] = []

            for dicts in outdicts:
                for key in self.var:
                    output[key].extend(dicts[key])

            file = open(
                f"/{disk}/users/mohan/G3PY/CacheDir/FullEASdata{year}.pkl", "wb"
            )
            pickle.dump([output, self.var], file)

            return output

    def merge(
        self,
        startyear=2013,
        endyear=2022,
        var=[
            "EvDate",
            "EvTime1",
            "EvTime2",
            "NKGSize",
            "Theta6",
            "Phi6",
            "Theta1",
            "Phi1",
            "ThetaOdd1",
            "PhiOdd1",
            "ThetaEven1",
            "PhiEven1",
            "ThetaOdd6",
            "PhiOdd6",
            "ThetaEven6",
            "PhiEven6",
        ],
    ):

        output = {}
        for key in var:
            output[key] = []

        for year in range(startyear, endyear + 1):
            file = open(
                f"/{disk}/users/mohan/G3PY/CacheDir/FullEASdata{year}.pkl", "rb"
            )

            out, var = pickle.load(file)
            for key in var:
                output[key].extend(out[key])

        for key in var:
            output[key] = np.array(output[key])

        # output['EvRa6'], output['EvDec6'] = CT.ZenPhiToRADec(output['Theta6'],\
        #                                                   output['Phi6'], output['EvDate'],
        #                                                   output['EvTime1'], output['EvTime2'] )

        output["EvRa1"], output["EvDec1"] = CT.ZenPhiToRADec(
            output["Theta1"],
            output["Phi1"],
            output["EvDate"],
            output["EvTime1"],
            output["EvTime2"],
        )

        with open(
            f"/{disk}/users/mohan/G3PY/CacheDir/FullEASdata_RaDec.pkl", "wb"
        ) as f:
            pickle.dump([output, var], f)

    def out(
        self,
        var=[
            "EvDate",
            "EvTime1",
            "EvTime2",
            "NKGSize",
            "Theta6",
            "Phi6",
            "Theta1",
            "Phi1",
            "ThetaOdd1",
            "PhiOdd1",
            "ThetaEven1",
            "PhiEven1",
            "ThetaOdd6",
            "PhiOdd6",
            "ThetaEven6",
            "PhiEven6",
        ],
        cache=True,
    ):

        if cache == True and os.path.exists(
            f"/{disk}/users/mohan/G3PY/CacheDir/FullEASdata_RaDec.pkl"
        ):

            file = open(
                f"/{disk}/users/mohan/G3PY/CacheDir/FullEASdata_RaDec.pkl", "rb"
            )

            cacheout, cachevar = pickle.load(file)

            if np.all(np.in1d(var, cachevar)):
                return cacheout
            else:
                print("var requested are not in cache")

        else:
            print("Cache not available. Create Cache")
