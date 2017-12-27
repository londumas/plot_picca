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
    "o1"                      : "",
    "o2"                      : "",
    "l1"                      : "",
    "l2"                      : "",
    "path"                    : "",
    "title"                   : "",
}

class Correlation3D:

    def __init__(self,dic=None):

        if (dic is None):
            dic = copy.deepcopy(raw_dic_class)

        ### info from dic
        self._correlation = dic['correlation']
        self._f1 = dic["f1"]
        self._f2 = dic["f2"]
        self._o1 = dic["o1"]
        self._o2 = dic["o2"]
        self._l1 = dic["l1"]
        self._l2 = dic["l2"]
        self._title = dic["title"]

        ### bin size (only square)
        self._binSize = None

        ### Grid
        self._nt = None
        self._np = None
        self._rp = None
        self._rt = None
        self._r  = None
        self._z  = None
        self._nb = None
        self._we = None
        self._co = None
        self._dm = None

        ### Correlation
        self._da = None
        self._er = None


        self.read_from_do_cor(dic["path"])

        return

    def __add__(self,other):

        self._rp = self._rp*self._we
        self._rt = self._rt*self._we
        self._z  = self._z*self._we
        self._da = self._da*self._we
        
        self._nb += other._nb
        self._we += other._we
        self._rp += other._rp*other._we
        self._rt += other._rt*other._we
        self._z  += other._z*other._we
        self._da += other._da*other._we

        cut = (self._we > 0.)
        self._rp[cut] /= self._we[cut]
        self._rt[cut] /= self._we[cut]
        self._z[cut]  /= self._we[cut]
        self._da[cut] /= self._we[cut]
        
        self._r  = sp.sqrt(self._rp**2. + self._rt**2.)
        
        return self

    def multiply(self,scalar):
        self._da *= scalar

        return

    def read_from_do_cor(self,path):

        vac = fitsio.FITS(path)

        ### bin size (only square)
        head = vac[1].read_header()
        self._binSize = head['RTMAX'] / head['NT']
    
        ### Grid
        self._nt = head['NT']
        self._np = head['NP']
        self._rt_min = 0.
        self._rt_max = head['RTMAX']
        try:
            self._rp_min = head['RPMIN']
        except:
            self._rp_min = -head['RPMAX']
        self._rp_max = head['RPMAX']
        self._rp = vac[1]['RP'][:]
        self._rt = vac[1]['RT'][:]
        self._r = sp.sqrt(self._rp**2. + self._rt**2.)
        self._z  = vac[1]['Z'][:]
        self._nb = vac[1]['NB'][:]
    
        ### Correlation
        we  = vac[2]['WE'][:]
        da  = vac[2]['DA'][:]
        self._we = we.sum(axis=0)
        cut = (self._we>0.)
        self._da       = (da*we).sum(axis=0)
        self._da[cut] /= self._we[cut]
        
        vac.close()

        return
    def read_from_export(self,path):
    
        vac = fitsio.FITS(path)
        self._rp = vac[1]['RP'][:]
        self._rt = vac[1]['RT'][:]
        self._r  = sp.sqrt(self._rp**2. + self._rt**2.)
        self._z  = vac[1]['Z'][:]
        self._da = vac[1]['DA'][:]
        self._co = vac[1]['CO'][:]
        self._dm = vac[1]['DM'][:]
        self._nb = vac[1]['NB'][:]
        
        self._er = sp.copy(sp.diag(self._co))
        cut = (self._er>0.)
        self._er[cut] = sp.sqrt(self._er[cut])
        
        vac.close()
        
        return
    def get_mean_redshift(self,rmin=80.,rmax=120.):
        
        cut = (self._r>rmin) & (self._r<rmax)
        mz = sp.sum( self._z[cut]*self._we[cut] )/sp.sum( self._we[cut] )

        return mz
    def print_bias(self,omega_M_0=0.315,rmin=80.,rmax=120.):

        z = self.get_mean_redshift(rmin=80.,rmax=120.)
        val = round(z, 4)
        print("z        = ",val)
        val = round(utils.croom(z), 4)
        print("b        = ",val)
        val = round(utils.growthStructure(z, omega_M_0), 4)
        print("g        = ",val)
        val = round(utils.growthStructureSimple(z, omega_M_0), 4)
        print("g_simple = ",val)
        val = round(utils.growthRateStructure(z, omega_M_0), 4)
        print("f        = ",val)
        val = round(utils.growthRateStructure(z, omega_M_0)/utils.croom(z), 4)
        print("beta     = ",val)

        return
    def plot_2d(self,x_power=0):

        origin='lower'
        extent=[self._rt_min, self._rt_max, self._rp_min, self._rp_max]
        if (self._correlation=='o_f' or self._correlation=='f_f2'):
            origin='upper'
            extent=[self._rt_min, self._rt_max, self._rp_max, self._rp_min]
    

        xxx = utils.convert1DTo2D(self._r,self._np,self._nt)
        yyy = utils.convert1DTo2D(self._da,self._np,self._nt)
        wee = utils.convert1DTo2D(self._we,self._np,self._nt)
    
        cut = (wee<=0.)
        if (xxx[cut].size==xxx.size):
            return
        yyy[ cut ] = float('nan')
        coef = sp.power(xxx,x_power)
    
        fig = plt.figure()
        ax = fig.add_subplot(111)
        #ax.set_xticks([ i for i in sp.arange(self._minX2D-50., self._maxX2D+50., 50.) ])
        #ax.set_yticks([ i for i in sp.arange(self._minY2D-50., self._maxY2D+50., 50.) ])

        plt.imshow(coef*yyy, origin=origin,extent=extent, interpolation='nearest')
        cbar = plt.colorbar()
    
        if (x_power==0):
            cbar.set_label(r'$\xi(\, r_{\parallel},r_{\perp} \,)$',size=40)
        if (x_power==1):
            cbar.set_label(r'$|r|.\xi(\, r_{\parallel},r_{\perp} \,)$',size=40)
        if (x_power==2):
            cbar.set_label(r'$|r|^{2}.\xi(\, r_{\parallel},r_{\perp} \,)$',size=40)


        plt.xlabel(r'$r_{\perp} \, [h^{-1} \, \rm{Mpc}]$', fontsize=40)
        plt.ylabel(r'$r_{\parallel} \, [h^{-1} \, \rm{Mpc}]$', fontsize=40)
        plt.grid(True)
        cbar.formatter.set_powerlimits((0, 0))
        cbar.update_ticks()

        plt.show()

        return
    def plot_1d(self,x_power=0, other=[]):

        list_corr = [self] + other

        for el in list_corr:
            r_max = sp.amax(el._r)
            r_min = 0.
            size = int((r_max-r_min)/el._binSize)

            xxx = sp.zeros(size)
            yyy = sp.zeros(size)
            wee = sp.zeros(size)

            we   = el._we
            wer  = el._we*el._r
            weda = el._we*el._da

            for i in sp.arange(size):
                cut = (el._r>i*el._binSize) & (el._r<(i+1)*el._binSize)
                xxx[i] += sp.sum(wer[cut])
                yyy[i] += sp.sum(weda[cut])
                wee[i] += sp.sum(we[cut])
            cut = (wee>0.)
            xxx[cut] /= wee[cut]
            yyy[cut] /= wee[cut]

            coef = sp.power(xxx,x_power)
            plt.errorbar(xxx,coef*yyy,marker='o',label=r'$'+el._title+'$')


        plt.xlabel(r'$r \, [\mathrm{h^{-1}\,Mpc}]$',fontsize=30)
        if (x_power==0):
            plt.ylabel(r'$\xi$',fontsize=30)
        if (x_power==1):
            plt.ylabel(r'$r \cdot \xi \, [\mathrm{h^{-1}\,Mpc}]$',fontsize=30)
        if (x_power==2):
            plt.ylabel(r'$r^{2} \cdot \xi \, [(\mathrm{h^{-1}\,Mpc})^{2}]$',fontsize=30)
        plt.legend(fontsize=20, numpoints=1,ncol=2, loc=1)
        plt.grid()
        plt.show()

        return
    def plot_slice_2d(self,sliceX=None,sliceY=None, other=[]):

        list_corr = [self] + other

        for el in list_corr:
            if (sliceX is not None):
                cut = (el._rt>el._binSize*sliceX) & (el._rt<el._binSize*(sliceX+1))
                xxx = el._rp[cut]
            if (sliceY is not None):
                cut = (el._rp>el._binSize*sliceY) & (el._rp<el._binSize*(sliceY+1))
                xxx = el._rt[cut]
            yyy = el._da[cut]
            if not el._er is None:
                yer = el._er[cut]
                
            if el._er is None:
                plt.errorbar(xxx,yyy,linewidth=4,label=r'$'+el._title+'$')
            else:
                plt.errorbar(xxx,yyy,yerr=yer,linewidth=4,label=r'$'+el._title+'$')

        if (sliceX is not None):
            plt.title(r"$0<r_{\perp}<4 \, [\mathrm{h^{-1}\,Mpc}]$",fontsize=30)
            plt.xlabel(r'$r_{\parallel} \, [\mathrm{h^{-1}\,Mpc}]$',fontsize=30)
        if (sliceY is not None):
            plt.xlabel(r'$r_{\perp} \, [\mathrm{h^{-1}\,Mpc}]$',fontsize=30)
        plt.ylabel(r'$\xi$',fontsize=30)
        plt.legend(fontsize=20, numpoints=1,ncol=2, loc=1)
        plt.grid()
        plt.show()

        return
    def plot_cov(self):

        cov = self._co

        if False:
            ###
            plt.plot(sp.diag(cov))
            plt.grid()
            plt.show()
            ###
            plt.plot(sp.diag(cov)*self._nb)
            plt.grid()
            plt.show()
            ###
            plt.plot(sp.diag(cov)*self._rt)
            plt.grid()
            plt.show()
        if True:
            cor = utils.getCorrelationMatrix(cov) 
            ###
            #plt.imshow(cor, interpolation='nearest')
            #plt.show()
            ###
            for i in range(3): 
                mcor = sp.asarray( [ sp.mean(sp.diag(cor,k=i+self._nt*k)) for k in sp.arange(self._np) ]  )
                plt.plot(sp.arange(mcor.size)*self._binSize,mcor,linewidth=2)
            plt.grid()
            plt.show()

        return






























