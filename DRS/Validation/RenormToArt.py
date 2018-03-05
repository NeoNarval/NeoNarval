import numpy as np
from scipy import interpolate
import pyfits
import matplotlib.pyplot as plt



def Renorm(lbd,flx):
#
#--- input is:
#          lbd..... wavelength grid
#          flx..... initial (unproperly normalized) flux values
#
#--- returns renormalized profile [see FitCont procedure]
#
    lbdo=lbd
    flxo=flx

    for i in np.arange(10):
        lbd,flx=FitCont(lbd,flx)

    flin=interpolate.interp1d(lbd,flx,kind='linear',bounds_error=False,fill_value=1.0)
    coeff=np.polyfit(lbd,flx,3)
    #s=interpolate.InterpolatedUnivariateSpline(lbd,flx)
    
    fcp=flin(lbdo)
    p = np.poly1d(coeff)
    
    plt.subplot(1,2,1)
    plt.plot(lbdo,flxo)
    plt.plot(lbdo,fcp)
    plt.plot(lbdo,p(lbdo))
    
    plt.subplot(1,2,2)
    plt.plot(lbdo,flxo/p(lbdo))
    #plt.plot(lbdo,s(lbdo))
    plt.show()
    
    return flxo/fcp

def FitCont(lbd,flx):
#
#--- input is:
#          lbd..... wavelength grid
#          flx..... initial (unproperly normalized) flux values
#
#--- returns estimate of continuum level (to be reinterpolated onto original grid for correction)
#
    pp=np.polyfit(lbd,flx,8)
    fcp=pp[0]*lbd*lbd*lbd*lbd*lbd*lbd*lbd*lbd + pp[1]*lbd*lbd*lbd*lbd*lbd*lbd*lbd + pp[2]*lbd*lbd*lbd*lbd*lbd*lbd + pp[3]*lbd*lbd*lbd*lbd*lbd + pp[4]*lbd*lbd*lbd*lbd + pp[5]*lbd*lbd*lbd + pp[6]*lbd*lbd + pp[7]*lbd + pp[8]

    dif=flx-fcp
    mmf=np.mean(dif)
    ect=np.std(dif)
    crit= (dif > mmf-0.5*ect) & (dif < (mmf+3*ect) )

    lbd=lbd[crit]
    flx=flx[crit]

    return lbd,flx



file_Lune='luna_calibered.fits'
l=pyfits.open(file_Lune)
a=l[1].data
l.close()

Orderlimit=np.concatenate(([0],np.where(np.diff(a['wavelength_lane1'])<0)[0],[len(a['wavelength_lane1'])-1]))

for order in np.arange(1,30):
	Oprof=a['intensity_lane1'][Orderlimit[order-1]+1:Orderlimit[order]]
	Owv=a['wavelength_lane1'][Orderlimit[order-1]+1:Orderlimit[order]]
	out=Renorm(Owv,Oprof)

