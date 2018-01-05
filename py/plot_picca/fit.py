### Python lib
import scipy as sp
import scipy.stats
import copy
import matplotlib.pyplot as plt
import h5py
import os.path

from . import utils


raw_dic_class = {
    "path" : "",
    "name" : "",
}

class Fit:

    def __init__(self,dic=None):

        if (dic is None):
            dic = copy.deepcopy(raw_dic_class)
    
        self._name     = dic["name"]
        self._param    = None
        self._fitAtrrs = None
        self._da       = None
        self.read_fit_results(dic["path"],dic["name"])
    
        return
    
    def read_fit_results(self,path,name):
    
        path = os.path.expandvars(path)
        f = h5py.File(path,"r")
        
        ### Parameters
        self._param    = {}
        self._fitAtrrs = {}
        self._fit      = {}
        lst_forFit = ["list of free pars","list of fixed pars","fval","ndata","npar","cov["]
        for el in list(f["best fit"].attrs):
            if any( ell in el for ell in lst_forFit):
                if str(el)=="list of free pars" or str(el)=="list of fixed pars":
                    self._fitAtrrs[str(el)]=[ ell.decode('UTF-8') for ell in f["best fit"].attrs[el]]
                else:
                    self._fitAtrrs[str(el)]=f["best fit"].attrs[el]
                continue
            else:
                tmp = {}
                tmp["value"] = f["best fit"].attrs[el][0]
                tmp["error"] = f["best fit"].attrs[el][1]
                tmp["name"]  = str(el)
                self._param[str(el)] = tmp
        
        ### Set errors to zero for unfitted param
        for el in self._fitAtrrs["list of fixed pars"]:
            self._param[el]["error"] = 0.
        
        ### Best fit
        self._fitAtrrs["own_chi2"] = f[name].attrs["chi2"]
        self._da = f[name]["fit"].value
        
        ### Set proba
        self._fitAtrrs["proba"] = 1.-sp.stats.chi2.cdf(self._fitAtrrs["fval"],self._fitAtrrs["ndata"]-self._fitAtrrs["npar"])
        
        f.close()
    
        return
    def print_fitted_par(self,lst=None,coeffBias=1.,header=True):
    
        if lst is None:
            lst = self._fitAtrrs["list of free pars"]
    
        if header:
            ### 
            to_print0  = " || " + "".ljust(30) + " || "
            to_print1  = " || " + "".ljust(30) + " || "
            for p in lst:
                to_print0 += self._param[p]["name"].ljust(20)
                to_print0 += " || "
                to_print1 += "".ljust(20)
                to_print1 += " || "
            to_print0 += "".ljust(30) + " || "
            to_print1 += "".ljust(30) + " || "
            print(to_print0)
            print(to_print1)
        ###
        to_print  = " || " + self._name.ljust(30) + " || "

        for p in lst:
            val = self._param[p]["value"]
            err = self._param[p]["error"]
            if len(p)>len("bias") and p[:4]=="bias":
                val *= coeffBias
                err *= coeffBias
            val = utils.format_number_with_precision(val,err)
            err = utils.format_number_with_precision(err,err)
            to_print += (val + " +/- " + err).ljust(20)
            to_print += " || "
            
        val = self._fitAtrrs["fval"]
        err = 0.1
        val = utils.format_number_with_precision(val,err)
        nbBin   = str(self._fitAtrrs["ndata"])
        nbParam = str(self._fitAtrrs["npar"])
        proba = self._fitAtrrs["proba"]
        proba = utils.format_number_with_precision(proba,proba)
        s       = val + ' / (' + nbBin + '-' + nbParam + '),  p = ' + proba
        to_print  += s.ljust(30) + " || "
            
        print(to_print)
    
        return














