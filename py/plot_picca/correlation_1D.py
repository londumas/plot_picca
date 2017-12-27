### Python lib
import fitsio
import scipy as sp
import copy
import matplotlib.pyplot as plt

from . import utils

raw_dic_class = {
    "correlation"             : "",
    "f1"                      : "",
    "f2"                      : "",
    "path"                    : "",
    "title"                   : "",
}

class Correlation1D:

    def __init__(self,dic=None):

        if (dic is None):
            dic = copy.deepcopy(raw_dic_class)

        ### info from dic
        self._correlation = dic['correlation']
        self._f1 = dic["f1"]
        self._f2 = dic["f2"]
        self._title = dic["title"]

        self._llmin = None
        self._llmax = None
        self._dll   = None
        self._n1d   = None
        self._var   = None
        self._cor   = None
        self._mat   = None
        self.read_from_do_cor(dic["path"])

        return

    def read_from_do_cor(self,path):

        vac = fitsio.FITS(path)
        head  = vac[1].read_header()
        self._llmin = head['LLMIN']
        self._llmax = head['LLMAX']
        self._dll   = head['DLL']
        self._n1d   = int((self._llmax-self._llmin)/self._dll+1)

        ### All Matrix
        self._mat = {}
        self._mat["DA"] = vac[2]['DA'][:]
        self._mat["WE"] = vac[2]['WE'][:]
        self._mat["NB"] = vac[2]['NB'][:]

        ### Variance
        self._var = sp.zeros( (self._n1d,4) )
        self._var[:,0] = 10.**( sp.arange(self._n1d)*self._dll+self._llmin )
        self._var[:,1] = sp.diag(self._mat["DA"])
        self._var[:,2] = sp.diag(self._mat["WE"])
        self._var[:,3] = sp.diag(self._mat["NB"])

        ### Correlation
        self._cor = sp.zeros( (self._n1d,4) )
        self._cor[:,0] = 10.**( sp.arange(self._n1d)*self._dll )

        inUp = False
        upperTri = sp.triu(self._mat["NB"])
        if (upperTri>0.).sum()==0:
            inUp=True

        norm=1.
        if sp.trace(self._mat["NB"],offset=0)>0:
            norm = sp.sum(sp.diag(self._mat["DA"],k=0)*sp.diag(self._mat["WE"],k=0))/sp.trace(self._mat["WE"],offset=0)

        for i in range(self._n1d):
            d = i
            if inUp:
                d = 1+i-self._n1d
            tda = sp.diag(self._mat["DA"],k=d)
            twe = sp.diag(self._mat["WE"],k=d)
            tnb = sp.trace(self._mat["NB"],offset=d)
            if tnb>0:
                self._cor[i,1] = sp.sum(tda*twe)/sp.sum(twe)/norm
                self._cor[i,2] = sp.sum(twe)/norm
                self._cor[i,3] = tnb

        vac.close()
    
        return
    def plot_var(self):

        x = self._var[:,0]
        y = self._var[:,1]
        w = (self._var[:,2]>0.) & (self._var[:,3]>10)
        x = x[w]
        y = y[w]
        if x.size==0: return
        plt.plot(x,y,linewidth=4)
        plt.xlabel(r'$\lambda_{\mathrm{Obs.}} \, \mathrm{\AA{}}$',fontsize=30)
        plt.ylabel(r'$\xi^{ff,1D}$',fontsize=30)
        plt.grid()
        plt.show()

        return
    def plot_cor(self):

        x = self._cor[:,0]
        y = self._cor[:,1]
        w = (self._cor[:,2]>0.) & (self._cor[:,3]>10)
        x = x[w]
        y = y[w]
        if x.size==0: return
        plt.plot(x,y,linewidth=4)
        plt.xlabel(r'$\lambda_{1}/\lambda_{2}$',fontsize=30)
        plt.ylabel(r'$\mathrm{normalized} \, \xi^{ff,1D}(\lambda_{1},\lambda_{2})$',fontsize=30)
        plt.grid()
        plt.show()

        return
    def plot_mat(self):

        ###
        for k in range(5):
            x = sp.arange(sp.diag(self._mat["DA"],k=k).size)
            y = sp.diag(self._mat["DA"],k=k)
            w = (sp.diag(self._mat["WE"],k=k)>0.) & (sp.diag(self._mat["NB"],k=k)>=2.)
            x = x[w]
            y = y[w]
            plt.plot(x,y,linewidth=4,alpha=0.7)
        plt.grid()
        plt.show()

        return










