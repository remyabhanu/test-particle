# Imports
import spacepy.time as spt
import spacepy.omni as om
import numpy as np
import tsygFort
import weimFort
import datetime as dt
from datetime import timedelta
from scipy.interpolate import interp1d
### Debug
class avedrift:
    #pcy2B = 0.0
    # Input parameters
    # 1- Date and time of the injection
    #date = dt.datetime(2014,9,27)
    #--------------------------------------------------------------------------------------#
    #imod = 1        # magnetic field model: 0=dipole, 1=T96_01, 2=T01_01, 3=T04_s
    #ipot = 1        # electric field model: 1=Weimer, 2=CIMI 
    #itrace = -1      # 1=forward, -1=backward tracing
    #tsec0 = 3600       # initial time in second
    #tend = 0     # end time in second
    #ro = 6          # initial radial distance in RE
    #mlt = 6         # initial mlt in hour
    #ekev0 = 60      # initial energy in keV
    #Kin = 0         # K index, 0=equatorially mirroring
    #ion= 1          # 1=ion, -1=electron
    #tf =60         # time resolution in seconds in updating omni parameters
    def __init__(self,date,ro=None,mlt=None,ekev0=None,tsec0=None,tend=None,sinA=1,imod = 1,ipot = 1,itrace = 1,ion= 1,tf =60):
        self.pcy2B = 0.0
        self.K = 0.0
        self.date =date
        self.sinA = sinA
        self.ro = ro
        self.mlt = mlt
        self.ekev0 = ekev0
        self.tsec = tsec0
        self.tend = tend
        self.imod = imod
        self.ipot = ipot
        self.itrace = itrace
        self.ion = ion
        self.tf = tf
    #---------------------------------------------------------------------------
    def mapping(self,imod,parmod,xeq,yeq,zeq,lati,mlti,switch):
    #-----------------------------------------------------------------------------
    # Routine mapping field line from equator to earth surface and vice versa
    # input:  imod,parmod
    #         if switch=1,  xeq,yeq,zeq
    #         if switch=-1, lati,mlti
    # output:  
    #         if switch=1,  lati,mlti
    #         if switch=-1, xeq,yeq,zeq,iout
            exname = ['zeroB','T96_01','T01_01','T04_s']  # The name should be the same than
            inname = ['DIP_08','IGRF_GSW_08'] # in routine rhand_08. !!Check capital letters!!
            dsmax=1.33
            err=0.00001
            Lmax = 1000
            pi=np.pi
            rlim=15.
            r0=1.
            iopt=1                          # dummy variable in TRACE
            ps=0.                           # force ps = 0
            tsygFort.geopack1.sps=np.sin(ps)
            tsygFort.geopack1.cps=np.cos(ps)
            if (switch ==1):        # mapping from equator to earth surface
                DIR=-1.
                xf,yf,zf,xx,yy,zz,npf = tsygFort.trace_08(xeq,yeq,zeq,DIR,dsmax,err,rlim,r0,\
                    iopt,parmod,exname[imod],'DIP_08',Lmax) 
                rf=np.sqrt(xf*xf+yf*yf+zf*zf)
                lati=np.arcsin(zf/rf)            # latitude in radian at ionosphere
                mlti=np.arctan2(yf,xf)+pi        # local time in radian at ionosphere
                return lati,mlti
            else :                       # mapping from Earth surface to equator
                DIR=1.
                rr=1.0001
                xi=-rr*np.cos(lati)*np.cos(mlti)
                yi=-rr*np.cos(lati)*np.sin(mlti)
                zi=rr*np.sin(lati)   
                xf,yf,zf,xx,yy,zz,npf,iout=tsygFort.trace1_08(xi,yi,zi,DIR,dsmax,err,rlim,r0,\
                                    iopt,parmod,exname[imod],'DIP_08',Lmax)   
                xeq=xf
                yeq=yf
                zeq=zf
                return xf,yf,zf,iout
    def TsyParmod(self,date_t,imod):
        # Get SW data from omni
        omni_data = om.get_omni(date_t,dbase='qd5min')
        # Solar wind speed
        vxgse = omni_data['Vsw']
        #vxgse = 400.
        vygse = 0.
        vzgse = 0.
        # Magnetic activity [SW pressure (nPa), Dst, ByIMF, BzIMF]
        parmod = np.zeros(10)
        parmod[0:4] = [omni_data['Pdyn'], omni_data['Dst'], omni_data['ByIMF'], omni_data['BzIMF']]
        if imod == 2:
            parmod[4:6]=omni_data['G'][0][0:2]
        elif imod == 3:
            parmod[4:10]=omni_data['W']
        nsw = omni_data['Den_P']
        #nsw = 5
        return parmod,nsw,vxgse,vygse,vzgse
    def Tsy_w_dip(self,imod,parmod,xeq,yeq,zeq,psi=0,setpsi=1):
        #psi tilt specify by user
        #If setpsi > 0 use psi value, if setpsi < 0 use psi as recalc_08 
        # Wrapper to the Tsyganenko model plus dipole
        exname = ['zeroB','T96_01','T01_01','T04_s']  # The name should be the same than
        inname = ['DIP_08','IGRF_GSW_08'] # in routine rhand_08. !!Check capital letters!!
        bx,by,bz,Bmag = tsygFort.tsymod(xeq,yeq,zeq,1,parmod,psi,setpsi,exname[imod],'DIP_08')
        return bx,by,bz,Bmag 
    #-----------------------------------------------------------------------------
    def Vdrift(self,imod,ipot,parmod,vswi,xnswi,xme,x0,Eo,ion,dra,sinA,veldip=0):
    #-----------------------------------------------------------------------------
    # Routine calculates bounce averaged drift velocity across field lines 
    # labelled by the ionospheric foot points: 
    # vl = vel(1) = v along latitude in rad/s
    # vp = vel(2) = v along mlt in rad/s
    # input: imod,parmod,vswi,xnswi,xme,x0,ekev,ion,dra
    # output: vel
        pi=np.pi
        pi2=2.*pi
        m_one=-1
        UseAL=False
        ALindex=10.          # arbitrary value
        Tilt=0.              # dipole tilt angle in degree, force it to zero
        re=6.371e6           # earth's radius (m)
        Eo2=Eo*Eo
        # Find 4 adjacent points of x0
        dra1=0.25*dra
        Xadj=np.zeros([4,2])
        Xadj[0,0]=x0[0]-dra1
        Xadj[0,1]=x0[1]
        Xadj[1,0]=x0[0]+dra1
        Xadj[1,1]=x0[1]
        Xadj[2,0]=x0[0]
        Xadj[2,1]=x0[1]-dra1
        Xadj[3,0]=x0[0]
        Xadj[3,1]=x0[1]+dra1
        #Find kinetic energy at Xadj
        EeV = np.zeros(4)
        for i in range(4):
            #### Include here calculation of pitch angle (or B at the mirror point)
            #### Using 2nd Invariant K(Bm)
            #ipdb.set_trace()
            if sinA == 1:
                _sinA = sinA
                xeq,yeq,zeq,iout = self.mapping(imod,parmod,0,0,0,Xadj[i,0],Xadj[i,1],m_one)
                bx,by,bz,_Bmag = self.Tsy_w_dip(imod,parmod,xeq,yeq,zeq)
                #print("B at zsm=0: %6.2f" %(_Bmag))
            else:
                _Alpha,_Bmag = self.AlphaOfK(self.K,Xadj[i,0],Xadj[i,1],imod,parmod)
                _sinA = np.sin(_Alpha)
                #ipdb.set_trace()
                #print("B at zsm=0: %6.2f" %(_Bmag))
            pc_sq=self.pcy2B*_Bmag/_sinA/_sinA #pcy2B is constant, first invariant 
            _ekev=np.sqrt(pc_sq+Eo2)-Eo  
            EeV[i]=1000.*_ekev
        #ipdb.set_trace()
        # Find potentials at Xadj
        if (ipot==1):   # Setup Weimer potential
            angle=np.arctan2(parmod[2],parmod[3])*180./pi   # deg from northward toward +Y 
            Bt=np.sqrt(parmod[2]*parmod[2]+parmod[3]*parmod[3]) # IMF mag in Y-Z plane
            weimFort.setmodel(angle,Bt,Tilt,vswi,xnswi,ALindex,UseAL)
        potent=np.zeros(4)
        for i in range(4):
            if (Xadj[i,1] <= 0.): Xadj[i,1]=Xadj[i,1]+pi2   
            if (Xadj[i,1] >= pi2): Xadj[i,1]=Xadj[i,1]-pi2    
            if (ipot==1):                # Weimer potential
                gLat=Xadj[i,0]*180./pi         # invariant lat. in degree
                gMLT=Xadj[i,1]*12./pi          # convert mlt from radians to hour 0-24.
                BnLat=weimFort.boundarylat(gMLT)        # boundary latitude
                if (gLat >= BnLat): 
                    potent[i]=weimFort.epotval(gLat,gMLT)*1000.  # potent(i) in V
            else:        # use cimi potential
                print('TODO: setup CIMI potential. Please use only Weimer')
        # Calculate sum of kinetic and potential energy (normalized by q)
        Esum=np.zeros(4)
        for i in range(4):
            Esum[i]=EeV[i]+ion*potent[i]
        # Calculate the drift velocity
        #ipdb.set_trace()
        cor=2.*pi/86400.                          # corotation speed in rad/s
        dra2=dra1*2.
        ksai=xme*np.sin(2.*x0[0])/re
        vel=np.zeros(2)
        veldipi_lat = veldip*np.cos(x0[1]) # dipolarization velocity
        veldipi_long = veldip*np.sin(x0[1])
        vel[0]=-(Esum[3]-Esum[2])/dra2/ksai/ion - veldipi_lat          # V along latitude
        vel[1]=(Esum[1]-Esum[0])/dra2/ksai/ion + cor + veldipi_long    # V along mlt
        #print(veldipi_lat,veldipi_long)
        return vel
    def calc_EeV(): ## Not finished. This function should be used to paralelize the first for loop in Vdrift
        xeq,yeq,zeq,iout = self.mapping(imod,parmod,0,0,0,Xadj[i,0],Xadj[i,1],m_one)
        bx,by,bz,_Bmag = self.Tsy_w_dip(imod,parmod,xeq,yeq,zeq)
        #### Include here calculation of pitch angle (or B at themirror point)
        #### Using 2nd Invariant K(Bm)
        #ipdb.set_trace()
        if sinA == 1:
            _sinA = sinA
        else:
            _sinA = np.sin(self.AlphaOfK(self.K,Xadj[i,0],Xadj[i,1],imod,parmod))
        pc_sq=self.pcy2B*_Bmag/_sinA/_sinA #pcy2B is constant, first invariant 
        _ekev=np.sqrt(pc_sq+Eo2)-Eo  
        EeV[i]=1000.*_ekev
    #-----------------------------------------------------------------------------
    def rk4(self,imod,ipot,parmod,vswi,xnswi,ion,h,Eo,xme,x0,sinA,veldip=0):
    #-----------------------------------------------------------------------------
    #   *** FOURTH-ORDER RUNGE-KUTTA ***
    #   Solve xend = x0 + fn*h                ! fn is the unit vector of ff
    #   input:  Vdrift,imod,ipot,parmod,vswi,xnswi,ion,h,Eo,xme,x0      
    #   output: xend,dt
        nd=2
        pi2=2.*np.pi
        xend = np.zeros(2)
        xwrk = np.zeros([4,2])
        ff = np.zeros(2)
        ff = self.Vdrift(imod,ipot,parmod,vswi,xnswi,xme,x0,Eo,ion,h,sinA,veldip=veldip)
        fmag0=np.sqrt(ff[0]*ff[0]+ff[1]*ff[1])
        for i in range(nd):
            xwrk[0,i]=h*ff[i]/fmag0
            xend[i]=x0[i]+xwrk[0,i]/2.
        ff = self.Vdrift(imod,ipot,parmod,vswi,xnswi,xme,xend,Eo,ion,h,sinA,veldip=veldip)
        fmag1=np.sqrt(ff[0]*ff[0]+ff[1]*ff[1])
        for i in range(nd):
            xwrk[1,i]=h*ff[i]/fmag1
            xend[i]=x0[i]+xwrk[1,i]/2.
        ff = self.Vdrift(imod,ipot,parmod,vswi,xnswi,xme,xend,Eo,ion,h,sinA,veldip=veldip)
        fmag2=np.sqrt(ff[0]*ff[0]+ff[1]*ff[1])
        for i in range(nd):
            xwrk[2,i]=h*ff[i]/fmag2
            xend[i]=x0[i]+xwrk[2,i]
        ff = self.Vdrift(imod,ipot,parmod,vswi,xnswi,xme,xend,Eo,ion,h,sinA,veldip=veldip)
        fmag3=np.sqrt(ff[0]*ff[0]+ff[1]*ff[1])
        for i in range(nd):
            xwrk[3,i]=h*ff[i]/fmag3
            xend[i]=x0[i]+(xwrk[0,i]+2.*xwrk[1,i]+2.*xwrk[2,i]+xwrk[3,i])/6.
        if (xend[1]<0.): xend[1]=xend[1]+pi2
        if (xend[1]>pi2): xend[1]=xend[1]-pi2
        fave=(fmag0+2.*fmag1+2.*fmag2+fmag3)/6.
        dt=h/fave
    #    xend[0]=x0[0]+0.001
    #    print(xend)
        return xend,dt  
    def trace_drift(self):
        # Input parameters
        # 1- Date and time of the injection
        #-------------------------------------------------------------------------------#
        _date = self.date
        _imod = self.imod        # magnetic field model: 0=dipole, 1=T96_01, 2=T01_1, 3=T04_s
        _ipot = self.ipot        # electric field model: 1=Weimer, 2=CIMI 
        _itrace = self.itrace      # 1=forward, -1=backward tracing
        _tsec0 = self.tsec       # initial time in second
        _tend = self.tend     # end time in second
        _ro = self.ro          # initial radial distance in RE
        _mlt = self.mlt         # initial mlt in hour
        _ekev0 = self.ekev0      # initial energy in keV
        _Kin = 0         # K index, 0=equatorially mirroring (reserved, not used now)
        _ion= self.ion          # 1=ion, -1=electron
        _tf = self.tf         # time resolution in seconds in updating omni parameters
        #-------------------------------------------------------------------------------#
        ### Main trace_drift subroutine 
        pi=np.pi
        dra=_itrace*0.5*pi/180.           # drift path tracing step size in radian
        iprint=1                         # print result every iprint steps
        re=6.371e6                       # earth's radius (m)
        if (_ion == 1): Eo=9.38269e5      # proton rest energy
        if (_ion == -1): Eo=511.          # electron rest energy
        Eo2=Eo*Eo
        pc_sq=(_ekev0+Eo)**2-Eo2         # (pc)^2 in keV^2
        latmax=72.*pi/180.              # max allowable latitude at the ionosphere
        reqmax=10.                      # max allowable equatorial distance
        istepmax=3000                   # max allowable istep
        rc=1.
        #Find dipole moment
        iyear = _date.year
        dts = spt.Ticktock(_date)
        iday = dts.DOY
        tsec0_dt = _date+timedelta(seconds=_tsec0)
        parmod,nsw,vxgse,vygse,vzgse = self.TsyParmod(_date+timedelta(seconds=_tsec0),_imod)# vy and vz are 0 to ensure GSM coordinate system
        tsygFort.recalc_08(iyear,iday,tsec0_dt.hour,tsec0_dt.minute,tsec0_dt.second,vxgse,vygse,vzgse)
        GGG=tsygFort.geopack2.g
        HHH=tsygFort.geopack2.h
        DIPMOM=np.sqrt(GGG[1]**2+GGG[2]**2+HHH[2]**2)   # DIPMOM in (nT RE^3)
        xme=DIPMOM*re**3*1.e-9                       # dipole moment in (T m^3)
        #xme = np.exp(36.58626097)
        # Find the initial ionospheric projection (along field line) point
        phi=_mlt*pi/12.               
        xeq=-_ro*np.cos(phi)              # +x --> Sun
        yeq=-_ro*np.sin(phi)
        zeq=0.
        lati,mlti = self.mapping(_imod,parmod,xeq,yeq,zeq,0,0,1)
        x0 = np.zeros(2)
        x0[0]=lati                 # latitude in radian
        x0[1]=mlti                 # local time in radian

        # find Bmag, sinA at the initial position and the invariant pcy2B
        bx,by,bz,Bmag=self.Tsy_w_dip(_imod,parmod,xeq,yeq,zeq)
        _sinA=self.sinA
        #ipdb.set_trace()
        self.pcy2B=pc_sq*_sinA*_sinA/Bmag  # First invariant, Bmag = B at equator
        self.K = self.KofAlpha(np.arcsin(_sinA),lati,mlti,_imod,parmod) # Second invariant K
        # Print the initial condition
        pitchA=np.arcsin(_sinA)*180./pi
        Lshell=rc/(np.cos(lati))**2
        mlti_d=mlti*12./pi              # mlt at ionosphere in hour
        print('%i %i %i %i %i %4.2f ! iyear,iday,imod,ipot,ion,tf' %(iyear,iday,_imod,_ipot,_ion,_tf))
        print('  tsec             Dst     Lshell     mlti        ro      mlto      ekeV   PA    Vsw   Pdyn')
        print('%6.2f %7.0f %9.3f %9.3f  %9.3f %6.3f  %6.1f %6.2f' %(_tsec0,parmod[1],Lshell,mlti_d,_ro,_mlt,_ekev0,pitchA))
        # ipdb.set_trace()
        # Start drift path tracing
        tsec=_tsec0
        istep=0
        iexit=0
        lastprint=0
        iexit=0
        vswi=vxgse
        xnswi = nsw
        Dst0 = parmod[1]
        xeq0=0;yeq0=0;zeq0=0
        ##
        roarr=np.zeros(istepmax)
        mltoarr=np.zeros(istepmax)
        tsecarr=np.zeros(istepmax)
        Eoarr=np.zeros(istepmax)
        ## Start debugging 
        while (iexit != 1):
            #ipdb.set_trace() 
            xend,dtsec=self.rk4(_imod,_ipot,parmod,vswi,xnswi,_ion,dra,Eo,xme,x0,_sinA)
            if (_sinA!=1):
                _Alpha,_Bmag = self.AlphaOfK(self.K,xend[0],xend[1],_imod,parmod)
                _sinA = np.sin(_Alpha)
            # update istep, tsec, Bmag
            istep=istep+1
            tsec=tsec+dtsec
            xeq,yeq,zeq,iout=self.mapping(_imod,parmod,0,0,0,xend[0],xend[1],-1)
            req=np.sqrt(xeq*xeq+yeq*yeq)
            # test for exit or continue
            if (istep==istepmax): iexit=1
            if (iout==1): iexit=1            
            if (xend[0]>=latmax): 
                iexit=1
                print('Warning: xend[0]>=latmax')
            if (req>=reqmax): iexit=1
            if ((_itrace==1) and (tsec>=_tend)): iexit=1            
            if ((_itrace==-1) and (tsec<=_tend)): iexit=1
            # Reasign values if iexit.eq.1
            if (iexit==1):
                #ipdb.set_trace()
                tsec=tsec-dtsec
                xend=x0
                parmod[1]=Dst0
                xeq=xeq0
                yeq=yeq0
                zeq=zeq0
                req=np.sqrt(xeq*xeq+yeq*yeq)
                if (tsec>=tprint): lastprint=1
            # print result every iprint step and at the end of the trace
            if ((np.mod(istep,iprint)==0) or (lastprint==1)):
                bx,by,bz,Bmag = self.Tsy_w_dip(_imod,parmod,xeq,yeq,zeq)
                ipdb.set_trace()
                pc_sq=self.pcy2B*Bmag/_sinA/_sinA
                ekev=np.sqrt(pc_sq+Eo2)-Eo
                pitchA=np.arcsin(_sinA)*180./pi
                Lshell=rc/(np.cos(xend[0]))**2
                mlti_d=xend[1]*12./pi          # mlt at ionosphere in hour
                phi=np.arctan2(yeq,xeq)
                mlteq=phi*12./pi+12.
                tprint=tsec
                print('%6.2f %7.0f %9.3f %9.3f  %9.3f %6.3f  %6.1f %6.2f %7.2f %6.4f'\
                      %(tsec,parmod[1],Lshell,mlti_d,req,mlteq,ekev,pitchA,vswi,parmod[0]))
                roarr[istep-1]=req ## variables for ploting
                mltoarr[istep-1]=mlteq
                tsecarr[istep-1]=tsec
                Eoarr[istep-1]=ekev
            # Save values 
            Dst0=parmod[1]
            xeq0=xeq
            yeq0=yeq
            zeq0=zeq
            # Update x0 and Tsyganenko parmod
            x0=xend
            if istep==1: dstep = int(_tf/dtsec) # choose an aproximate dt step for changing SW parameters
            if np.mod(istep,dstep)==0: 
                parmod,xnswi,vswi,vygse,vzgse = self.TsyParmod(_date+timedelta(seconds=tsec),_imod)
        return tsecarr,mltoarr,roarr,Eoarr,istep
    def KofAlpha(self,Alpha,lati,mlti,imod,parmod): 
        # lati and mlti are the footpoint of the FL at iono and Alpha is PA at magnetic equator
        exname = ['zeroB','T96_01','T01_01','T04_s']  # The name should be the same than
        inname = ['DIP_08','IGRF_GSW_08'] # in routine rhand_08. !!Check capital letters!!
        dsmax=.1
        err=0.0000001
        Lmax = 1000
        pi=np.pi
        rlim=15.
        r0=1.
        iopt=1                          # dummy variable in TRACE
        #ps=0.                           # force ps = 0
        #tsygFort.geopack1.sps=np.sin(ps)
        #tsygFort.geopack1.cps=np.cos(ps)
        DIR=1.
        rr=1.0001
        xi=-rr*np.cos(lati)*np.cos(mlti)
        yi=-rr*np.cos(lati)*np.sin(mlti)
        zi=rr*np.sin(lati)   # Mapping from North to South Hemisphere
        xf,yf,zf,xx,yy,zz,npf=tsygFort.trace_08(xi,yi,zi,DIR,dsmax,err,rlim,r0,\
                                    iopt,parmod,exname[imod],'DIP_08',Lmax)        
        Bmag = np.zeros(npf)
        ds = np.zeros(npf)
        for i in range(npf):
            bx,by,bz,Bmag[i] = self.Tsy_w_dip(imod,parmod,xx[i],yy[i],zz[i])
            ds[i] = np.sqrt((xx[i+1]-xx[i])**2+(yy[i+1]-yy[i])**2+(zz[i+1]-zz[i])**2)
        idx_eq = Bmag.argmin() 
        Beq = Bmag[idx_eq] 
        _sinA = np.sin(Alpha)
        Bm = Beq/(_sinA**2)
        idx_m1 = (np.abs(Bmag[0:idx_eq]-Bm)).argmin()
        idx_m2 = (np.abs(Bmag[idx_eq:-1]-Bm)).argmin() + idx_eq
        ### Do integration K ###
        I = np.sqrt(np.abs(Bm - Bmag[idx_m1:idx_m2+1]))
        #K  = np.dot(I,ds[idx_m1:idx_m2+1])
        S = np.cumsum(ds)
        K = np.trapz(I,S[idx_m1:idx_m2+1])
        #ipdb.set_trace()
        return K
    def AlphaOfK(self,Kt,lati,mlti,imod,parmod): # lati and mlti are the footpoint of the FL at iono
        exname = ['zeroB','T96_01','T01_01','T04_s']  # The name should be the same than
        inname = ['DIP_08','IGRF_GSW_08'] # in routine rhand_08. !!Check capital letters!!
        dsmax=.1
        err=0.0000001
        Lmax = 1000
        pi=np.pi
        rlim=15.
        r0=1.
        iopt=1                          # dummy variable in TRACE
        #ps=0.                           # force ps = 0
        #tsygFort.geopack1.sps=np.sin(ps)
        #tsygFort.geopack1.cps=np.cos(ps)
        DIR=1.
        rr=1.0001
        xi=-rr*np.cos(lati)*np.cos(mlti)
        yi=-rr*np.cos(lati)*np.sin(mlti)
        zi=rr*np.sin(lati)   # Mapping from North to South Hemisphere
        xf,yf,zf,xx,yy,zz,npf=tsygFort.trace_08(xi,yi,zi,DIR,dsmax,err,rlim,r0,\
                                    iopt,parmod,exname[imod],'DIP_08',Lmax)        
        Bmag = np.zeros(npf)
        ds = np.zeros(npf)
        for i in range(npf):
            bx,by,bz,Bmag[i] = self.Tsy_w_dip(imod,parmod,xx[i],yy[i],zz[i])
            ds[i] = np.sqrt((xx[i+1]-xx[i])**2+(yy[i+1]-yy[i])**2+(zz[i+1]-zz[i])**2)
        idx_eq = Bmag.argmin() 
        Beq = Bmag[idx_eq] 
        #if np.abs(zz[idx_eq])>0.1:
        #    print ("Warning: Min B out of equator at z=%6.2f, Bmin=%6.2f" %(zz[idx_eq],Beq))
            #ipdb.set_trace()
        Alpha = pi/180*np.array([20,25,30,35,40,45,50,55,60,70,80])
        _sinA = np.sin(Alpha)
        Bm = Beq/(_sinA**2)
        K = np.zeros(len(Alpha))
        for i in range(len(Alpha)): # look for Bm index
            idx_m1 = (np.abs(Bmag[0:idx_eq]-Bm[i])).argmin()
            idx_m2 = (np.abs(Bmag[idx_eq:-1]-Bm[i])).argmin() + idx_eq
            ### Do integration K ###
            I = np.sqrt(np.abs(Bm[i] - Bmag[idx_m1:idx_m2+1]))
            #K[i]  = np.dot(I,ds[idx_m1:idx_m2+1])
            S = np.cumsum(ds)
            K[i] = np.trapz(I,S[idx_m1:idx_m2+1])
        #ipdb.set_trace()
        f = interp1d(K, Alpha, kind='cubic',fill_value="extrapolate")
        alpha = f(Kt)
        #ipdb.set_trace()
        return alpha,Beq
    def trace_drift_shell(self,sm_loc,t_loc,Ech,sinA_loc,sim_dt,verbose=True,Vdip_inj=None):
        # This function is an improve version of trace_drift, it accept coordinates in SM as a vector
        # And the time as datetime parameter
        # Input parameters
        # 1- Date and time of the injection
        #-------------------------------------------------------------------------------#
        #_date = self.date
        _imod = self.imod        # magnetic field model: 0=dipole, 1=T96_01, 2=T01_1, 3=T04_s
        _ipot = self.ipot        # electric field model: 1=Weimer, 2=CIMI 
        _itrace = self.itrace      # 1=forward, -1=backward tracing
        #_tsec0 = self.tsec       # initial time in second
        #_tend = self.tend     # end time in second
        #_ro = self.ro          # initial radial distance in RE
        #_mlt = self.mlt         # initial mlt in hour
        _ekev0 = Ech      # initial energy in keV
        _ion= self.ion          # 1=ion, -1=electron
        _tf = self.tf         # time resolution in seconds in updating omni parameters
        #-------------------------------------------------------------------------------#
        ### Main trace_drift subroutine 
        pi=np.pi
        dra=_itrace*0.5*pi/180.           # drift path tracing step size in radian
        iprint=1                         # print result every iprint steps
        re=6.371e6                       # Earth's radius (m)
        if (_ion == 1): Eo=9.38269e5      # proton rest energy
        if (_ion == -1): Eo=511.          # electron rest energy
        Eo2=Eo*Eo
        pc_sq=(_ekev0+Eo)**2-Eo2         # (pc)^2 in keV^2
        latmax=72.*pi/180.              # max allowable latitude at the ionosphere
        reqmax=10.                      # max allowable equatorial distance
        istepmax=3000                   # max allowable istep
        rc=1.
        #Find dipole moment
        iyear = t_loc.year
        dts = spt.Ticktock(t_loc)
        iday = dts.DOY
        tsec0_dt = t_loc
        parmod,nsw,vxgse,vygse,vzgse = self.TsyParmod(tsec0_dt,_imod)# vy and vz are 0 to ensure GSM coordinate system
        tsygFort.recalc_08(iyear,iday,tsec0_dt.hour,tsec0_dt.minute,tsec0_dt.second,vxgse,vygse,vzgse)
        GGG=tsygFort.geopack2.g
        HHH=tsygFort.geopack2.h
        DIPMOM=np.sqrt(GGG[1]**2+GGG[2]**2+HHH[2]**2)   # DIPMOM in (nT RE^3)
        xme=DIPMOM*re**3*1.e-9                       # dipole moment in (T m^3)
        #xme = np.exp(36.58626097)
        
        # Find the initial ionospheric projection (along field line) point
        #phi=_mlt*pi/12.               
        #xeq=-_ro*np.cos(phi)              # +x --> Sun
        #yeq=-_ro*np.sin(phi)
        #zeq=0.
        lati,mlti = self.mapping(_imod,parmod,sm_loc[0],sm_loc[1],sm_loc[2],0,0,1)
        x0 = np.zeros(2)
        x0[0]=lati                 # latitude in radian
        x0[1]=mlti                 # local time in radian
        xeq,yeq,zeq,iout=self.mapping(_imod,parmod,0,0,0,x0[0],x0[1],-1)
        # find Bmag, sinA at the initial position and the invariant pcy2B
        bx,by,bz,Bmag=self.Tsy_w_dip(_imod,parmod,xeq,yeq,zeq) # find B at magnetic equator
        bx,by,bz,Bloc=self.Tsy_w_dip(_imod,parmod,sm_loc[0],sm_loc[1],sm_loc[2]) # find B at local position
        _sinA=np.sqrt(Bmag/Bloc*sinA_loc*sinA_loc) # Pitch angle at the magnetic equator
        if _sinA > 1:
            _sinA = 1
            print('Warning: B at the equator > B loc, possible orbit bifurcation. Temp solution make sin(PA)=1')
            #ipdb.set_trace()
        if (np.arcsin(_sinA) >= 80/180*pi): # if PA > 80 we consider as 90 degrees.
            print('Warning: real PA is %4.2f, as it is > 80 we assume 90' %(np.arcsin(_sinA)*180/pi) )
            _sinA = 1
            self.pcy2B=pc_sq*_sinA*_sinA/Bmag  # First invariant, Bmag = B at equator
            self.K = 0.
        else:
            self.pcy2B=pc_sq*_sinA*_sinA/Bmag  # First invariant, Bmag = B at equator
            self.K = self.KofAlpha(np.arcsin(_sinA),lati,mlti,_imod,parmod) # Second invariant K
        # Print the initial condition
        #ipdb.set_trace()
        pitchA=np.arcsin(_sinA)*180./pi
        Lshell=rc/(np.cos(lati))**2
        mlti_d=mlti*12./pi              # mlt at ionosphere in hour
        ro_loc=np.sqrt(sm_loc[0]*sm_loc[0]+sm_loc[1]*sm_loc[1]+sm_loc[2]*sm_loc[2])
        phi=np.arctan2(sm_loc[1],sm_loc[0])
        mlt_loc=phi*12./pi+12.
        print('Info: the first line show Ro local and mlt local, i.e at the pos of the S/C and not at the equator')
        print('%i %i %i %i %i %4.2f ! iyear,iday,imod,ipot,ion,tf' %(iyear,iday,_imod,_ipot,_ion,_tf))
        print('  tsec             Dst    Lshell    mlti       ro     mlto     ekeV    PA   Vsw   Pdyn')
        print('%s %7.0f %9.3f %9.3f  %9.3f %6.3f  %6.1f %6.2f' %(tsec0_dt.isoformat()\
                                                                 ,parmod[1],Lshell,mlti_d,ro_loc,mlt_loc,_ekev0,pitchA))
        #ipdb.set_trace()
        # Start drift path tracing
        epoch = dt.datetime(1970,1,1)
        _tend=(tsec0_dt-epoch).total_seconds() + _itrace*sim_dt # UTC seconds since epoch time
        tsec=(tsec0_dt-epoch).total_seconds() # UTC seconds since epoch time
        istep=0
        iexit=0
        lastprint=0
        iexit=0
        vswi=vxgse
        xnswi = nsw
        Dst0 = parmod[1]
        xeq0=0;yeq0=0;zeq0=0
        ##
        roarr=np.zeros(istepmax)
        mltoarr=np.zeros(istepmax)
        #tsecarr=np.zeros(istepmax)
        tsecarr = []
        Eoarr = np.zeros(istepmax)
        pitchAarr = np.zeros(istepmax)
        if Vdip_inj!=None:
            #ipdb.set_trace()
            vdip = self.Vdip_i(Vdip_inj,_imod,parmod) ## Added May 25, simulate injection
        else:
            vdip = 0
        ## Start debugging 
        while (iexit != 1):
            #ipdb.set_trace() 
            xend,dtsec=self.rk4(_imod,_ipot,parmod,vswi,xnswi,_ion,dra,Eo,xme,x0,_sinA,veldip=vdip)
            if (np.arcsin(_sinA) < 80/180*pi):
                _Alpha,_Bmag = self.AlphaOfK(self.K,xend[0],xend[1],_imod,parmod)
                _sinA = np.sin(_Alpha)
            else: _sinA = 1
            # update istep, tsec, Bmag
            istep=istep+1
            tsec=tsec+dtsec
            xeq,yeq,zeq,iout=self.mapping(_imod,parmod,0,0,0,xend[0],xend[1],-1)
            req=np.sqrt(xeq*xeq+yeq*yeq)
            # test for exit or continue
            if (istep==istepmax): iexit=1
            if (iout==1): iexit=1            
            if (xend[0]>=latmax): 
                iexit=1
                print('Warning: xend[0]>=latmax')
            if (req>=reqmax): iexit=1
            if ((_itrace==1) and (tsec>=_tend)): iexit=1            
            if ((_itrace==-1) and (tsec<=_tend)): iexit=1
            # Reasign values if iexit.eq.1
            if (iexit==1):
                #ipdb.set_trace()
                tsec=tsec-dtsec
                xend=x0
                parmod[1]=Dst0
                xeq=xeq0
                yeq=yeq0
                zeq=zeq0
                req=np.sqrt(xeq*xeq+yeq*yeq)
                if (tsec>=tprint): lastprint=1
            # print result every iprint step and at the end of the trace
            if ((np.mod(istep,iprint)==0) or (lastprint==1)):
                bx,by,bz,Bmag = self.Tsy_w_dip(_imod,parmod,xeq,yeq,zeq)
                #ipdb.set_trace()
                pc_sq=self.pcy2B*Bmag/_sinA/_sinA
                ekev=np.sqrt(pc_sq+Eo2)-Eo
                pitchA=np.arcsin(_sinA)*180./pi
                Lshell=rc/(np.cos(xend[0]))**2
                mlti_d=xend[1]*12./pi ## mlt at ionosphere in hour
                phi=np.arctan2(yeq,xeq)
                mlteq=phi*12./pi+12.
                tprint=tsec
                roarr[istep-1]=req ## variables for ploting
                mltoarr[istep-1]=mlteq
                #tsecarr[istep-1]= tsec
                tsecarr.append(epoch + timedelta(seconds=tsec))
                Eoarr[istep-1]=ekev
                pitchAarr[istep-1] = pitchA
                if verbose == True:
                    print('%s %7.0f %9.3f %9.3f  %9.3f %6.3f  %6.1f %6.2f %7.2f %6.4f'\
                          %(tsecarr[istep-1].isoformat(),parmod[1],Lshell,mlti_d,req,mlteq,ekev,pitchA,vswi,parmod[0]))
            # Save values 
            Dst0=parmod[1]
            xeq0=xeq
            yeq0=yeq
            zeq0=zeq
            # Update x0 and Tsyganenko parmod
            x0=xend
            if istep==1: dstep = int(_tf/dtsec) # choose an aproximate dt step for changing SW parameters
            if np.mod(istep,dstep)==0: 
                parmod,xnswi,vswi,vygse,vzgse = self.TsyParmod(tsecarr[istep-1],_imod)
        return tsecarr,mltoarr,roarr,Eoarr,pitchAarr,istep
    def Vdip_i(self,V,imod,parmod):
        lati1,mlti1 = self.mapping(imod,parmod,8,0,0,0,0,1) # equatorial position at x,y,z = 8,0,0
        lati2,mlti2 = self.mapping(imod,parmod,7,0,0,0,0,1) # equatorial position at x,y,z = 7,0,0
        if V==0:
            vel_rad=0
        else:
            delt = 1*6371/V
            vel_rad = abs(lati1-lati2)/delt
        return vel_rad
