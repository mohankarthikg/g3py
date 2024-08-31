from crflux.models import HillasGaisser2012


def H4Awts(primary, E, org_spec_ind=2.5):

    model = HillasGaisser2012(model="H4a")
    corsikaId = {
        "proton": 14,
        "helium": 402,
        "nitrogen": 1206,
        "aluminium": 2814,
        "iron": 5426,
    }

    weights = model.nucleus_flux(corsikaId[primary], E) / E ** (-org_spec_ind)

    return weights


def gammawt(E, org_spec_ind=2):

    model = HillasGaisser2012(model="H4a")
    corsikaId = {
        "proton": 14,
        "helium": 402,
        "nitrogen": 1206,
        "aluminium": 2814,
        "iron": 5426,
    }

    weights = model.total_flux(E) / E ** (-org_spec_ind)

    return weights
