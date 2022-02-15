import pylab as pl
import matplotlib
from matplotlib import ticker
import numpy as np
from numpy import log10 as log
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rc('font', family='serif', size=10)
matplotlib.rcParams['xtick.direction'] = 'in'
matplotlib.rcParams['ytick.direction'] = 'in'
matplotlib.rcParams['xtick.top'] = True
matplotlib.rcParams['ytick.right'] = True
fig_width_half = 3.4
from astropy import units as un, constants as cons
from astropy.cosmology import Planck15 as cosmo
import HaloPotential_old as Halo
import cooling_flow as CF
import WiersmaCooling as Cool

figDir = '/home/jonathan/Dropbox/jonathanmain/CGM/rapidCoolingCGM/figures/'
fig_size_half = [3.39375300954753, 2.09745470932262]
slantlineprops = {'color': '.5', 'headlength': 2, 'headwidth': 1.5, 'width': 0.1}

Z2Zsun = 1/3.
z = 1.
cooling = Cool.Wiersma_Cooling(Z2Zsun,z)
m=0 # power-law index
Rvir = 200.*un.kpc # potential.rvir() #virial radius
R_phi0 = 100*Rvir # radius where Potential is zero (for calculating Bernoulli Parameter)
R_min = 0.1*un.kpc           #inner radius of supersonic part of solution
R_max = 9.*Rvir              #outer radius of integration
max_step = 0.1               #resolution of solution in ln(r)
R_circ = 10.*un.kpc   #circularization radius
v0 = 3*un.km/un.s     #radial velocity at circularization radius
X = 0.7
tHubble = cosmo.age(z)
fb = cosmo.Ob0 / cosmo.Om0

class Res:
    pass
def spin(lMhalo,z,rvir,vc):
    Mhalo = 10**lMhalo*un.Msun
    potential = Halo.DK14_NFW(lMhalo,0.*un.Msun,z)
    rho_bs = fb * potential.rho(potential._Rs*un.kpc)
    nHs = (X * rho_bs / cons.m_p).to('cm**-3')
    tcools = cooling.tcool(potential.Tvir(), nHs).to('Gyr')
    rcool = np.interp(tHubble,tcools,potential._Rs)
    Mcool = potential.enclosedMass(rcool)
    if Mcool > Mhalo: Mcool=0.99*Mhalo
    
    spin0 = 0.035
    mu = 10**-0.6+1  #avg value in Bullock+01
    b = -mu * np.log(1-mu**-1)-1 #eqn. 10 in Bullock+01
    j0 = 2**0.5* spin0 * vc *rvir / b
    J = j0 *b
    j_to_j0 = Mcool/Mhalo  / (mu - Mcool / Mhalo)    
    if j_to_j0/b > 1: spin = spin0
    else: spin = j_to_j0/b * spin0
    
    print('lMhalo=%.1f mu=%.2f rfrac=%.2f Mfrac=%.2f b=%.2f j/j0=%.2f spin=%.4f'%(lMhalo,mu,rcool/rvir,Mcool/Mhalo,b,j_to_j0,spin))
    return spin
def Mdot_baryoncomplete(lMhalo,z,rvir):    
    potential = Halo.DK14_NFW(lMhalo,0.*un.Msun,z)
    rho_bs = fb * potential.rho(potential._Rs*un.kpc)
    nHs = (X * rho_bs / cons.m_p).to('cm**-3')
    tcools = cooling.tcool(potential.Tvir(), nHs).to('Gyr')
    rcool = np.interp(tHubble,tcools,potential._Rs)
    if rcool < rvir:
        rho_rcool = np.interp(tHubble,tcools,rho_bs.to('g/cm**3').value)*un.g/un.cm**3
        Mdot = (4*np.pi*rcool**3 / tHubble * rho_rcool).to('Msun/yr')
        print('lMhalo=%.1f Rcool=%d kpc, rho(R_cool)=%.2g, Mdot=%.2f Msun/yr'%(lMhalo,rcool.value,rho_rcool.value,Mdot.value))
    else:
        Mdot = (fb * 10**lMhalo*un.Msun / tHubble).to('Msun/yr')
        print('lMhalo=%.1f Mdot=%.2f Msun/yr'%(lMhalo,Mdot.value))
    return Mdot
def honeInonMdot(Mdot,PLpotentialAM,R_sonic_min,R_sonic_max,tol=0.05):
    _Mdot = 1e5*un.Msun/un.yr
    while abs(log(Mdot/_Mdot))>tol:
        R_sonic = (R_sonic_max*R_sonic_min)**0.5            
        res = CF.shoot_from_sonic_point(PLpotentialAM,cooling,R_sonic,R_max=500*un.kpc,R_min=R_min,max_step=max_step,pr=False)
        _Mdot = res.Mdot
        if _Mdot > Mdot: R_sonic_max = R_sonic
        else: R_sonic_min = R_sonic
        print('R_sonic=%.1f Mdot=%.2f goal M_dot=%.2f'%(R_sonic.value,_Mdot.value,Mdot.value))
    return res
def calcSolutions(vcs,fhot=1):    
    ress = []
    for ivc,vc in enumerate(vcs):
        print(vc)
        res  = Res()
        res.vc = vc
        res.M_halo = 4e11*un.Msun * (vc/(100.*un.km/un.s))**3.
        lMhalo = log(res.M_halo.value)        
        
        res.Rvir = 200.*un.kpc * (vc/(100.*un.km/un.s))
        res.Tvir = 290e3*un.K * (vc/(100.*un.km/un.s))**2.
        res.Mdot = Mdot_baryoncomplete(lMhalo, z, res.Rvir)*fhot**2
        print(res.Mdot)
        res.Rcirc2Rvir = 0.05#max(10.*un.kpc / res.Rvir,0.05) #2**0.5*spin(lMhalo, z, res.Rvir, res.vc)
        if vc>120*un.km/un.s:
            res.PLpotential = Halo.PowerLaw(m,res.vc,Rvir,R_phi0) #for shoot from Rcirc
            res.res = CF.shoot_from_R_circ(res.PLpotential,cooling,res.Rcirc2Rvir*res.Rvir,res.Mdot,
                                               R_max,v0,max_step=max_step,pr=True,epsilon=0.03,T_high=res.Tvir)
        else:
            res.PLpotentialAM = Halo.PowerLaw_with_AngularMomentum(m,res.vc,Rvir,R_phi0,res.Rcirc2Rvir*res.Rvir) #for shoot from sonic point
            if vc>50*un.km/un.s:                
                res.res = honeInonMdot(res.Mdot,res.PLpotentialAM,res.Rcirc2Rvir*res.Rvir,res.Rvir)
            else:
                res.res = CF.shoot_from_sonic_point(res.PLpotentialAM,cooling,R_sonic=150*un.kpc,R_max=500*un.kpc,R_min=R_min,max_step=max_step,pr=False)
        ress.append(res)
    #     ress.append( )    
    goods = [res for res in ress if res.res!=None]
    print([res.M_halo for res in goods])
    return goods
def plotSolutions(goods):
    fs=12; alpha=0.1
    for ifig,resClass in enumerate(goods):
        res = resClass.res
        f=pl.figure(figsize=(4.5,5.8))
        pl.subplots_adjust(hspace=0.15,top=0.92,left=0.25,right=0.85)
        for iPanel in range(3):
            c = 'k'
            ax = pl.subplot(3,1,iPanel+1)
            if iPanel==0: ys = res.Ts()
            if iPanel==1: ys = res.Ms()
            if iPanel==2: ys = res.Ks()
            pl.plot(res.Rs()/resClass.Rvir, log(ys.value),c=c,ls='-',lw=1)
            R_circ = resClass.Rcirc2Rvir
            pl.axvline(R_circ,c='.5',lw=0.5,ls='--')
            xls = 0.025,1.
            pl.xlim(*xls)
    
            pl.semilogx()
            ax.yaxis.set_major_locator(ticker.MultipleLocator(1))
    #         ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.5))
            pl.yticks(range(-5,8),[r'$10^{-5}$',r'$10^{-4}$',r'$10^{-3}$',r'$0.01$',r'$0.1$',r'$1$',r'$10$',r'$100$',r'$10^{3}$',r'$10^{4}$',r'$10^{5}$',r'$10^{6}$',r'$10^{7}$'])
    #         ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x,pos: r'$%d$'%x))  
            if iPanel==2:
                pl.xlabel(r'$r/R_{\rm vir}$',fontsize=13)        
                #ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x,pos: (r'$%d$'%x,'')[x in (4,6,7,8,9,40,60,70,80,90)]))    
            else:
                ax.xaxis.set_major_formatter(ticker.NullFormatter())
            if res.R_sonic()!=None:
                x0 = res.R_sonic().value[0]/resClass.Rvir.value
                pl.fill_between(np.array([R_circ,x0]),np.array([-10,-10]),np.array([10,10]),facecolor='b',alpha=alpha,zorder=-100)
            else: 
                x0 = R_circ        
            if x0<xls[1]:
                pl.fill_between(np.array([x0,xls[1]]),np.array([-10,-10]),np.array([10,10]),facecolor='r',alpha=alpha,zorder=-100)
    
            if iPanel==0:
                pl.text(-0.25,0.5,'temperature\n[K]',fontsize=fs,ha='center',va='center',rotation=90,transform=ax.transAxes)
                coeff,index = ('%.2g'%resClass.M_halo.value)[0],('%.2g'%resClass.M_halo.value)[-2:]
                pl.text(0.95,0.075,r'$%s\cdot 10^{%s}{\rm M}_\odot,\ z=0,\ 0.3{\rm Z}_\odot$'%(coeff,index),
                        ha='right',fontsize=fs-1,transform=ax.transAxes)
                pl.ylim(3.8,7.2) 
                pl.text(5/resClass.Rvir.value,7.35,'rotational\nsupport',color='r',ha='center',fontsize=fs)
                pl.annotate(r'',(10/resClass.Rvir.value,6.8),(8/resClass.Rvir.value,7.4),arrowprops=slantlineprops)
                pl.axhline(log(resClass.Tvir.value),c='k',lw=0.5,ls='--')
                pl.text(210/resClass.Rvir.value,log(resClass.Tvir.value)-0.1,r'$T_{\rm vir}$',va='center',fontsize=fs)
    
                if res.R_sonic()!=None:
                    pl.text((x0*R_circ)**0.5,7.35,'free\nfall',color='b',ha='center',fontsize=fs)
                if x0<xls[1]:
                    pl.text((x0*xls[1])**0.5,7.35,'cooling\nflow',color='r',ha='center',fontsize=fs)
            if iPanel==1: 
                pl.ylabel('Mach\n',fontsize=fs)
                pl.ylim(-2,1.5)
                pl.axhline(0.,c='k',lw=0.5,ls='--')  
                ax2 = pl.twinx()
                pl.ylabel(r'$\sim t_{\rm cool}/t_{\rm ff}$',fontsize=fs)
                pl.yticks(range(-5,7),[r'$10^{-5}$',r'$10^{-4}$',r'$10^{-3}$',r'$0.01$',r'$0.1$',r'$1$',r'$10$',r'$100$',r'$10^{3}$',r'$10^{4}$',r'$10^{5}$',r'$10^{6}$'])
                pl.ylim(2,-1.5)                      
                ax2.tick_params(axis='both',labelsize=13)
            #if iPanel==2:
                #pl.ylabel('density\n[cm$^{-3}$]',fontsize=fs)
                #pl.ylim(-5.7,-1.7)  
                #rs = 10.**np.arange(log(xls[0]),log(xls[1]),.01)*un.kpc
                #potential = Halo.DK14_NFW(log(resClass.M_halo.value),0.*un.Msun,z)                
                #pl.plot(rs,log(0.7*potential.dk14.rho_b(rs).value/cons.m_p.to('g').value),c='k',lw=0.5,ls='--')
                #pl.text(215,-3.8,r'$f_{\rm b}\rho_{\rm DM}$',fontsize=fs)
                #pl.annotate(r'',(130,-4.5),(215,-3.8),arrowprops=slantlineprops)
            if iPanel==2:
                pl.ylabel('entropy\n[keV cm$^{-2}$]',fontsize=fs)
                pl.ylim(-1.5,2.5)
            ax.xaxis.set_major_locator(ticker.LogLocator(subs=range(1,10)))
            ax.tick_params(axis='both',labelsize=13)
        pl.savefig(figDir+'steady_state_presentation_F%05d.pdf'%ifig)    


def plotSolutions2(goods,nPanels=2,letter='H',fhot=1,suffix='pdf'):
    fs=12; alpha=0.1
    for ifig,resClass in enumerate(goods):
        res = resClass.res
        f=pl.figure(figsize=(4.5,(5.8,5)[nPanels==2]))
        pl.subplots_adjust(hspace=0.15,top=0.8,left=0.25,right=0.85)
        for iPanel in range(nPanels):
            c = 'k'
            ax = pl.subplot(nPanels,1,iPanel+1)
            if iPanel==0: ys = res.Ts()
            if iPanel==1: ys = res.Ms()
            if iPanel==2: ys = res.Ks()
            pl.plot(res.Rs(), log(ys.value),c=c,ls='-',lw=1)
            R_circ = resClass.Rcirc2Rvir*resClass.Rvir.value
            pl.axvline(R_circ,c='.5',lw=0.5,ls='--')
            xls = 0.0,100.
            pl.xlim(*xls)
    
            #pl.semilogx()
            ax.yaxis.set_major_locator(ticker.MultipleLocator(1))
            ax.yaxis.set_minor_locator(ticker.FixedLocator(log(np.concatenate([np.arange(2,10)*0.01,
                                                       np.arange(2,10)*0.1,
                                                       np.arange(2,10)*1,
                                                       np.arange(2,10)*10,
                                                       np.arange(2,10)*100,
                                                       np.arange(2,10)*1e3,
                                                       np.arange(2,10)*1e4,
                                                       np.arange(2,10)*1e5,
                                                       np.arange(2,10)*1e6,
                                                       np.arange(2,10)*1e7]))))
            pl.yticks(range(-5,8),[r'$10^{-5}$',r'$10^{-4}$',r'$10^{-3}$',r'$0.01$',r'$0.1$',r'$1$',r'$10$',r'$100$',r'$10^{3}$',r'$10^{4}$',r'$10^{5}$',r'$10^{6}$',r'$10^{7}$'])
    #         ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x,pos: r'$%d$'%x))  
            if ax.is_last_row():
                pl.xlabel(r'radius [kpc]',fontsize=13)        
                #ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x,pos: (r'$%d$'%x,'')[x in (4,6,7,8,9,40,60,70,80,90)]))    
            else:
                ax.xaxis.set_major_formatter(ticker.NullFormatter())
            if res.R_sonic()!=None:
                x0 = res.R_sonic().value[0]
                pl.fill_between(np.array([R_circ,x0]),np.array([-10,-10]),np.array([10,10]),facecolor='b',alpha=alpha,zorder=-100)
            else: 
                x0 = R_circ        
            if x0<xls[1]:
                pl.fill_between(np.array([x0,xls[1]]),np.array([-10,-10]),np.array([10,10]),facecolor='r',alpha=alpha,zorder=-100)
    
            if iPanel==0:                                
                coeff,index = ('%.2g'%resClass.M_halo.value)[0],('%.2g'%resClass.M_halo.value)[-2:]
                if letter!='I':
                    pl.text(-0.25,0.5,'temperature\n[K]',fontsize=fs,ha='center',va='center',rotation=90,transform=ax.transAxes)
                    pl.text(0.95,0.075,r'$f_{\rm CGM}=1,\ Z={\rm Z}_{\rm MZR}$',
                        ha='right',fontsize=fs-1,transform=ax.transAxes)
                    pl.ylim(4,7.2) 
                    pl.text(0,7.35,r'$R_{\rm circ}$',color='r',ha='center',fontsize=fs)
                    pl.annotate(r'',(R_circ,6.8),(6.5,7.3),arrowprops=slantlineprops)
                else:
                    pl.text(-0.25,0.5,r'$T\ [{\rm K}]$',fontsize=fs,ha='center',va='center',rotation=90,transform=ax.transAxes)
                    pl.text(0.95,0.075,r'$f_{\rm hot}=%.1f,\ Z={\rm Z}_\odot/3$'%fhot,
                        ha='right',fontsize=fs-1,transform=ax.transAxes)
                    pl.ylim(5,6.5) 
                    pl.text(0,6.65,r'$R_{\rm circ}$',color='r',ha='center',fontsize=fs)
                    pl.annotate(r'',(R_circ,6.3),(6.5,6.6),arrowprops=slantlineprops)
                
                pl.axhline(log(resClass.Tvir.value),c='k',lw=0.5,ls='--')
                pl.text(102,log(resClass.Tvir.value)-0.1,r'$T_{\rm vir}$',va='center',fontsize=fs)
    
                if res.R_sonic()!=None:
                    pl.text((x0+R_circ)/2.,7.35,"free\nfall",color='b',ha='center',fontsize=fs)
                if x0<xls[1] and letter!='I':
                    pl.text(70,7.35,"cooling\nflow",color='r',ha='center',fontsize=fs)
            if iPanel==1:                 
                if letter!='I': 
                    pl.ylim(-1.5,1.)
                    pl.ylabel('inflow Mach number\n',fontsize=fs)
                else: 
                    pl.ylim(-2,0.3)
                    pl.ylabel('inflow Mach number\n',fontsize=fs)                    
                pl.axhline(0.,c='k',lw=0.5,ls='--')  
                ax2 = pl.twinx()
                pl.ylabel(r'$\approx t_{\rm cool}/t_{\rm ff}$',fontsize=fs)
                pl.yticks(range(-5,7),[r'$10^{-5}$',r'$10^{-4}$',r'$10^{-3}$',r'$0.01$',r'$0.1$',r'$1$',r'$10$',r'$100$',r'$10^{3}$',r'$10^{4}$',r'$10^{5}$',r'$10^{6}$'])
                if letter!='I':  pl.ylim(1.5,-1.)                      
                else: pl.ylim(2,-0.3)
                ax2.tick_params(axis='both',labelsize=13)
            if iPanel==2:
                if letter!='I':
                    pl.ylabel(r"'entropy' $(\propto P\rho^{-\gamma})$"+"\n[keV cm$^2$]",fontsize=fs)
                else:
                    pl.ylabel(r"$K\ [{\rm keV}\ {\rm cm}^2$]",fontsize=fs)
                    
                pl.ylim(-0.,2.)
            
            #ax.xaxis.set_major_locator(ticker.LogLocator(subs=range(1,10)))
            ax.tick_params(axis='both',labelsize=13)
        pl.text(0.475,0.935,r'solution for ',transform=f.transFigure,fontsize=14, ha='right')
        pl.text(0.475,0.935,r'$z=0, \log\ M_{\rm h}=$',transform=f.transFigure,fontsize=16, ha='left')
        pl.text(0.9,0.935,r'$%.1f$'%log(resClass.M_halo.value),transform=f.transFigure,fontsize=16, ha='right',color='red')
        
        if suffix=='pdf': pl.savefig(figDir+'steady_state_presentation_%c%05d.pdf'%(letter,ifig))
        if suffix=='png': pl.savefig(figDir+'steady_state_presentation_%c%05d.png'%(letter,ifig),dpi=600)
