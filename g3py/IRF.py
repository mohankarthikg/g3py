import pickle

def spt_pdf():
    with open("/boson/users/mohan/G3PY/bg_spt_pdf.pkl", 'rb') as file:
        bckg_spt_pdf, bckg_spt_cdf, bckg_spt_inv_cdf, sindec_min, sindec_max = pickle.load(file)
    
    return bckg_spt_pdf, bckg_spt_cdf, bckg_spt_inv_cdf, sindec_min, sindec_max
    
def spt_en_pdf2():
    with open("/boson/users/mohan/G3PY/spt_en_pdf2.pkl", "rb") as f:
        
        bg_spt_en_pdf, sig_en_pdf, sindecmin, sindecmax, log10Nemin, log10Nemax = pickle.load(f)
    return bg_spt_en_pdf, sig_en_pdf, sindecmin, sindecmax, log10Nemin, log10Nemax