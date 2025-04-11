#!/bin/python3

#import os
#import sys
import struct
from typing import Optional, Dict, List, Tuple, Union
from datetime import datetime, timedelta
import numpy as np
#import xarray as xr

def _decompress(xwork):
    xdat = []
    for value in xwork:
        if value > 0.0:
            xdat.append(value)
        else:
            xdat.extend([0.0]*int(-value))
    return xdat

def _skip_n_lines(f,n):
    for _ in range(n):
        rmark= struct.unpack("i", f.read(4))[0] 
        f.seek(rmark,1)
        f.seek(4    ,1)


class CalpuffOutput:
    def __init__(self):

        # File metadata
        self.filepath: str = ""
        self.file_type: str = ""                    # 'CONC', 'DFLX', 'WFLX', 'VISB'. 'RHO', 'T2D', 'FOG'
        self.model_version: str = ""                # > 7.2.1
        self.title: str = ""

        #Internal vars: (used to manipulate the data)
        self.NCOM: int = 0                          # Number of commented lines at header section.
        self.is_compressed: bool = 0                # Is main data compressed?

        # Date-time specifications
        self.start_date: datetime = None            # first reported date-time "YYYY-MM-DD HH:MM:SS"
        self.end_date: datetime = None              # last reported date-time  "YYYY-MM-DD HH:MM:SS"
        self.time_step: timedelta = None            # time step in seconds.
        self.run_length: int = 0                    # number of time steps on run.
        self.timezone: timedelta = timedelta(0)     # timezone

        # Grid specifications
        self.grid: Dict = {
            'nx': 0,   'ny': 0,                     # Grid size
            'x0': 0.0, 'y0': 0.0,                   # Grid origin (south-west corner)
            'dx': 0.0, 'dy': 0.0,                   # Grid spacing
        }

        # Proj specifications
        self.proj: Dict = {
            'crs': 'UTM',                           # Coordinate system
            'zone': 0,                              # UTM zone 
            'hemis': '',                            # UTM zone 
            'datum': 'WGS84'                        # Datum
        }

        # Species information
        self.nsp: int = 0                           # Number of species reported
        self.species: List[str] = []                # List of species with their properties

        # Source information
        self.nsrctype: int = 0                      # Number of source types: POINT, LINE, ...
        self.nsrcbytype: List[int] = []             # Number of sources of each type

        # Receptors information
        self.has_gridded_receptors: bool = 0        # Is there a grid of receptors?
        self.nrec: int = 0                          # Number of discrete receptors
        self.x_r: Optional[np.ndarray]
        self.y_r: Optional[np.ndarray]
        self.z_r: Optional[np.ndarray]                       #ground elevation
        self.h_r: Optional[np.ndarray]                       #above ground elevation
        self.igrp: List[int]
        self.ihill_r: List[int]
        self.nctrec: int = 0                        # Number of complex-terrain receptors

        self.rgroups: List[int]

        ## The actual data - using xarray for labeled multi-dimensional data
        #self.data: Optional[xr.Dataset] = None

    def info(self):

        idate=self.start_date.strftime("%Y %j %H:%M:%S")
        edate=self.end_date.strftime("%Y %j %H:%M:%S")

        print(f"Calpuff Output File \"{self.filepath}\":")
        if (self.is_compressed): print(f"(Compressed)")
        print(f"{self.title}")
        print(f"   Start Date: {idate}.")
        print(f"   End   Date: {edate}.")
        print(f"   Time-Step:  {self.time_step}.")    
        print(f"   Run Lenght: {self.run_length}.") #Avg. Time: {self.time_ave}. 
        print(f"")
        print(f"   Grid parameters:")
        print(f"      Grid size:  {self.grid['nx']} x {self.grid['ny']}.")
        print(f"      Cell size:  {self.grid['dx']} x {self.grid['dy']} m.")
        print(f"      Xmin, Ymin: ({self.grid['x0']}, {self.grid['y0']})")
        print(f"")
        print(f"   Proj. parameters:\n       {self.proj['datum']}, {self.proj['crs']}, {self.proj['zone']} {self.proj['hemis']}")
        print(f"")
        print(f"   Species: ({self.nsp})")
        print(f"      Species names: {self.species}")
        print(f"      Species units: {self.units}")
        print(f"")
        print(f"   Receptors:")
        print(f"      Number of gridded         rec.: {self.grid['nx']} x {self.grid['ny']}")
        print(f"      Number of discrete        rec.: {self.nrec}")
        print(f"      Number of complex terrain rec.: {self.nctrec}")
        print(f"")
        print(f"   Sources:")
        print(f"      Number of source type     rec.: {self.nsrctype}")
        print(f"      Number of source by type  rec.: {self.nsrcbytype}")

    def get_coordinates(self):
        nx = self.grid["nx"]; ny = self.grid["ny"]
        dx = self.grid["dx"]; dy = self.grid["dy"]
        x0 = self.grid["x0"]; y0 = self.grid["y0"]

        x = x0 + np.arange(nx) * dx
        y = y0 + np.arange(ny) * dy

        # create 2D meshgrid of coordinates
        X, Y = np.meshgrid(x, y) 

        return X, Y

    def get_gridded_data(self, pollut):
         """
         Returns CALPUFF gridded receptors data.
         """
         if ( not self.has_gridded_receptors ):
             raise ValueError(f"No gridded receptors found in file {self.filepath}")
         if ( pollut in self.species ):
             raise ValueError(f"No pollutant {pollut} found in file {self.filepath}")

         NT=self.run_length
         NX=self.grid['nx']
         NY=self.grid['ny']
         out = np.zeros([NT, NY, NX])

         with open(self.filepath, "rb") as f:

            #Skip header
            _skip_n_lines(f,7+self.NCOM)
            if ( self.nrec > 0):
                _skip_n_lines(f, 2)          # Record #NCOM+6 : discrete receptors
                                             # Record #NCOM+6a: receptors group names
            if (self.nctrec > 0):
                _skip_n_lines(f,1)           # Record #NCOM+7 : complex terrain receptors    
            for i in range(self.nsrctype):
                n=self.nsrcbytype[i]
                if ( n > 0 ):
                    _skip_n_lines(f,1)       # Record #NCOM+8 : source names

            #Main data:
            for t in range(NT): #temporal loop
                #write(io8)nyrab,njulab,nhrab,nsecab,nyre,njule,nhre,nsece
                f.read(4)[0]
                yr1,dy1,hr1,sec1, yr2,dy2,hr2,sec2 = struct.unpack("8i", f.read(32))
                f.read(4) #EOL
            
                #write(io8)ktype,ksource,csrcnam,xmapkm,ymapkm
                f.read(4)[0]
                ktype,ksource,csrcnam,xmapkm,ymapkm = struct.unpack("2i 16s 2f", f.read(32))
                f.read(4)
            
                #Loop on output species
                for sp in range(self.nsp): 

                    # ---    Gridded receptor concentrations
                    if ( pollut == self.species[sp] ):
    
                        if ( self.is_compressed ):
                            f.read(4)
                            ii = struct.unpack("i", f.read(4))[0] # "number of words" in compressed record.
                            f.read(4)
            
                            f.read(4)
                            specie = f.read(15).decode("utf-8").split() 
                            raw = list(struct.unpack(f'{ii}f', f.read(4*ii) ))
                            xdat = _decompress( raw )
                            out[t,:,:] = np.array(xdat).reshape((NY,NX)) #, order='F')
                            f.read(4) #EOL

                        else:
                            f.read(4)
                            specie = f.read(15).decode("utf-8").split() 
                            out[t,:,:] = np.array([[ struct.unpack("f",f.read(4)) for i in range(nyg)] for j in range(nxg)])
                            f.read(4) #EOL
 
                    else:
                        if (self.is_compressed):
                            _skip_n_lines(f, 2)
                        else:
                            _skip_n_lines(f, 1)

                    ## ---    Discrete receptor concentrations
                    if ( self.nrec > 0):
                        _skip_n_lines(f,1) #self.nrec+self.nctrec)

                    ## ---    Discrete CTSG receptor concentrations
                    if ( self.nctrec > 0):
                        _skip_n_lines(f,1) #self.nrec+self.nctrec)
        
         return(out)

    def get_discrete_data(self, pollut):
         """
         Returns CALPUFF discrete receptors data.
         """
        
         if ( self.nrec < 1):
             raise ValueError(f"No discrete receptors found in file {self.filepath}")
         if ( not pollut in self.species ):
             raise ValueError(f"No pollutant {pollut} found in file {self.filepath}")

         NT = self.run_length
         out = np.zeros([NT, self.nrec ])
                                                                                             
         with open( self.filepath, "rb") as f:

            #Skip header
            _skip_n_lines(f,7+self.NCOM)
            if ( self.nrec > 0):
                _skip_n_lines(f, 2)          # Record #NCOM+6 : discrete receptors
                                             # Record #NCOM+6a: receptors group names
            if (self.nctrec > 0):
                _skip_n_lines(f,1)           # Record #NCOM+7 : complex terrain receptors    

            for i in range(self.nsrctype):
                n=self.nsrcbytype[i]
                if ( n > 0 ):
                    _skip_n_lines(f,1)       # Record #NCOM+8 : source names


            #Main data:
            for t in range(NT): #temporal loop
                #write(io8)nyrab,njulab,nhrab,nsecab,nyre,njule,nhre,nsece
                f.read(4)[0]
                yr1,dy1,hr1,sec1, yr2,dy2,hr2,sec2 = struct.unpack("8i", f.read(32))
                f.read(4) #EOL
                
                #write(io8)ktype,ksource,csrcnam,xmapkm,ymapkm
                f.read(4)[0]
                ktype,ksource,csrcnam,xmapkm,ymapkm = struct.unpack("2i 16s 2f", f.read(32))
                f.read(4)

                #Loop on output species                                                                                                    
                for sp in range(self.nsp): 
                                                                                                                                 
                    # ---    Gridded receptor concentrations
                    if ( self.has_gridded_receptors ):
                        if ( pollut == self.species[sp] ):
                                                                                                                                 
                            if ( self.is_compressed ):
                                _skip_n_lines(f,2)

                            else:
                                _skip_n_lines(f,1)
                                                                                                                             
                    ## ---    Discrete receptor concentrations
                    if ( self.species[sp] == pollut ):

                        if ( self.is_compressed ):
                            f.read(4)
                            ii = struct.unpack("i", f.read(4))[0] #
                            f.read(4)
                                                                                                            
                            f.read(4)
                            specie = f.read(15).decode("utf-8").split() 
                            raw = list(struct.unpack(f'{ii}f', f.read(4*ii) ))
                            xdat = _decompress( raw )
                            out[t,:] = np.array(xdat) 
                            f.read(4) #EOL
                        else:
                            f.read(4)
                            specie   = struct.unpack("15s",f.read(15))
                            for i in range(self.nrec):
                                tmp= struct.unpack("f",f.read(4))[0]
                                print(tmp)
                                out[t,i] = tmp
                            f.read(4)
                    else:
                        _skip_n_lines(f,1)

                    ## ---    Discrete CTSG receptor concentrations
                    if ( self.nctrec > 0):
                        _skip_n_lines(f,1) #self.nrec+self.nctrec)

         return(out)



    def time_avg_max(self, pollut, interval):
        """
        Computes, for each grid cell, the maximum of the time-averaged fields.

        Parameters:
            data (np.ndarray): Input array with shape [nt, nx, ny].
            interval (int): Number of time steps to average over.

        Returns:
            np.ndarray: 2D array [nx, ny] with the max of time-averaged values.
        """

        data=self.get_data(pollut)

        nt, nx, ny = data.shape
        if nt % interval != 0:
            data = data[:nt - (nt % interval)]
            print(f"Warning: data trimmed to {data.shape[0]} time steps to match interval")

        # Reshape and average over time chunks
        data_reshaped = data.reshape(-1, interval, nx, ny)
        averaged = data_reshaped.mean(axis=1)  # shape: [n_chunks, nx, ny]

        # Take max over the time-averaged fields
        max_field = averaged.max(axis=0)  # shape: [nx, ny]

        return max_field

def read_file(filepath: str) -> 'CalpuffOutput':
    """
    Opens and reads a CALPUFF output file, parsing all header information
    and initializing a CalpuffOutput object.

    Args:
        filepath: Path to the CALPUFF output file (e.g: CONC.DAT)

    Returns:
        Initialized CalpuffOutput object with all header information

    Raises:
        ValueError: If file format is invalid
        IOError: If file cannot be read
    """
    # Initialize empty CalpuffOutput object
    output = CalpuffOutput()
    output.filepath = filepath

    # Open the binary file
    with open(filepath, "rb") as f:
        print(f"Reading file {filepath}...")
    
        #FORTRAN's UNFORMATTED RECORDS HAS THIS STRUCTURE :
        #    [4-byte length] [binary data] [4-byte length]

        # RECORD #1 --  
        #write(io8) conset,dataver,datamod
        _skip_n_lines(f,1)
        #f.seek(f.read(4)[0]+4, 1)    # skip this record.   f.read(4+96+4)
    
        # RECORD #2 --  NCOM     
        #write(io8) ncom
        f.read(4) #record marker
        NCOM = struct.unpack("i", f.read(4))[0]
        f.read(4) #end of record

        # RECORD #3 to NCOM+2 --  Comments.
        for _ in range(NCOM):
            _skip_n_lines(f,1) #f.seek(f.read(4)[0]+4, 1)#skip record 

        # RECORD #NCOM+3 -- General run parameters
        #write(io8)mmodel,ver,level,ibyr,ibjul,ibhr,ibsec,
        #1 abtz,irlg,iavg,nsecdt,nx,ny,dxkm,dykm,ione,xorigkm,yorigkm,
        #2 nssta(1),ibcomp,iecomp,jbcomp,jecomp,ibsamp,jbsamp,iesamp,
        #3 jesamp,meshdn,nsrctype,msource,nrec,nrgrp,
        #4 nctrec,lsamp,nspout,lcomprs,i2dmet,
        #5 iutmzn,feast,fnorth,rnlat0,relon0,xlat1,xlat2,
        #6 pmap,utmhem,datum,daten,clat0,clon0,clat1,clat2
        f.read(4) 
        MODEL, VERSION, LEVEL          = struct.unpack("12s 12s 12s", f.read(36))
        IBYR, IBJUL, IBHR, IBSEC, ABTZ = struct.unpack("4i 8s",       f.read(24))
        IRLG, IAVG, NSECDT = struct.unpack("3i", f.read(12))
        NX,NY, DXKM,DYKM, IONE, XORIGKM,YORIGKM  = struct.unpack("2i 2f i 2f",f.read(28))

        NSSTA = struct.unpack("i", f.read(4))                            #n ssurface stations
        IBCOMP, IECOMP, JBCOMP, JECOMP = struct.unpack("4i", f.read(16)) #computational grid
        IBSAMP, JBSAMP, IESAMP, JESAMP = struct.unpack("4i", f.read(16)) #sampling grid
        MESHDN = struct.unpack("i", f.read(4))[0]                        #sampling grid spacing factor
    
        nsrctype, msource, nrec, nrgrp, nctrec  = struct.unpack("5i", f.read(20)) #
        lsamp, nspout, lcomprs, i2dmet          = struct.unpack("4i", f.read(16)) # sample grid? # of species. is data compresed?, rel humidity source data?
        iutmzn, feast,fnorth,rnlat0,relon0,xlat1,xlat2 = struct.unpack("i 6f", f.read(28))       # proj parameters.
        pmap,utmhem,datum,daten,clat0,clon0,clat1,clat2 = struct.unpack("8s 4s 8s 12s 16s 16s 16s 16s", f.read(96)) # proj parameters.
        f.read(4) #close record

        #Definitions of Sample Grid (which is the one reported)
        nxg = (IESAMP - IBSAMP )*MESHDN +1
        nyg = (JESAMP - JBSAMP )*MESHDN +1
        dxg = (DXKM/float(MESHDN))*1000.
        dyg = (DYKM/float(MESHDN))*1000.
        x0g = (XORIGKM+IBSAMP*DXKM)*1000
        y0g = (YORIGKM+JBSAMP*DYKM)*1000

        # NON-REPORTED RECORD: how many sources from each type
        #write(io8) (nsrcbytype(n),n=1,nsrctype)
        f.read(4) #record marker
        nsrcbytype = struct.unpack(str(nsrctype)+"i", f.read(nsrctype*4)) #[nsrcbytype(n),n=1,nsrctype]
        f.read(4) #close record
    
        # RECORD #NCOM+4 -- TITLE
        f.read(4) #record marker (record size in bytes)
        TITLE = f.read(80*3).decode("utf-8").strip()
        f.read(4) #EOL
    
        # RECORD #NCOM+5 -- List of species-groups output
        #write(io8)(csout(n),n=1,nspout)
        #for grp in range(ngrp):
        f.read(4)
        csout=[""]*nspout
        for i in range(nspout):
           csout[i] = f.read(15).decode("utf-8").strip()
        f.read(4) #EOL
    
        # RECORD #NCOM+5a-- List of species-groups units 
        #write(io8)(acunit(n),n=1,nspout)
        f.read(4) # record marker
        acunit=[""]*nspout
        for i in range(nspout):
          acunit[i] = f.read(16).decode("utf-8").strip()
        f.read(4) #EOL
    
        # RECORD #NCOM+6 -- Discrete (non-gridded) receptor data
        #write(io8)(tmp3(n),n=1,nrec),
        #1     (tmp4(n2),n2=1,nrec),(elevng(n3),n3=1,nrec),
        #2     (zng(n4),n4=1,nrec), (irgrp(n5),n5=1,nrec)
        if (nrec > 0):
            x=np.zeros(nrec)
            y=np.zeros(nrec) 
            z=np.zeros(nrec) 
            h=np.zeros(nrec)
            igrp=[0]*nrec
            f.read(4) #record marker
            for i in range(nrec):
                x[i], y[i], z[i], h[i], igrp[i] = struct.unpack("4f i", f.read(20))
            f.read(4) #EOL
            
            # RECORD #NCOM+6a -- Receptor-group names
            rgrpnam=[""]*nrgrp
            f.read(4) #record marker
            for i in range(nrgrp):
                rgrpnam[i] = f.read(80).decode("utf-8").strip()
            f.read(4) #EOL

            output.x_r = x*1e3
            output.y_r = y*1e3
            output.z_r = z
            output.h_r = h
            output.igrp = igrp
            output.rgroups = rgrpnam
    
        # RECORD #NCOM+7 -- Complex terrain receptor data
        if (nctrec > 0):
            f.read(4) #record marker
            for i in range(nctrec):
                x[i],y[i],z[i],ihill[i] = struct.unpack("5f", f.read(20))
            f.read(4) #EOL
            
        # RECORD #NCOM+8 -- Source names
        #do itype=1,nsrctype
        #if(nsrcbytype(itype).GT.0) then
        #write(io8) itype,[cnamsrc(n,itype), n=1,nsrcbytype(itype)]
        #print(f"Number of sorces types: {nsrctype}")
        cnamsrc=[""]*nsrctype
        for itype in range(nsrctype):
            n=nsrcbytype[itype] 
            if ( n > 0 ):
                f.read(4) #open record line
                f.read(4) #itype      
                for j in range(n):
                    cnamsrc = f.read(16).decode("utf-8").strip()
                f.read(4) #EOL
    
        # RECORD #NCOM+8  -- Nearest Surface Station for VISIBILITY ONLY
        # RECORD #NCOM+9  -- X coord (UTM) of stations for VISIBILITY ONLY
        # RECORD #NCOM+10 -- Y coord (UTM) of stations for VISIBILITY ONLY

        # Temporal params:
        output.model_version=VERSION
        output.start_date=datetime.strptime(f"{IBYR} {IBJUL} {IBHR} {IBSEC}","%Y %j %H %S")
        output.run_length=IRLG
        output.ave_time =timedelta(hours=IAVG)
        output.time_step=timedelta(seconds=NSECDT)
        output.end_date=output.start_date+output.time_step*output.run_length
                                                                                            
        # Grid params:
        output.grid['nx'] = nxg          ;  output.grid['ny'] = nyg
        output.grid['x0'] = x0g          ;  output.grid['y0'] = y0g         
        output.grid['dx'] = dxg          ;  output.grid['dy'] = dyg      

        # Proj params:
        output.proj['crs']   = pmap.decode("utf-8").strip()     #
        output.proj['datum'] = datum.decode("utf-8").strip()    #
        output.proj['zone']  = iutmzn                            
        output.proj['hemis'] = utmhem.decode("utf-8").strip()   # north/south-ern hemisphere                                                  

        # Species params:
        output.nsp     = nspout
        output.nrgrp   = nrgrp
        output.species = csout
        output.units   = acunit
                                                                                            
        # Receptors params:
        output.has_gridded_receptors = lsamp    # gridded         receptors
        output.nrec    = nrec                   # discrete        receptors
        output.nctrec  = nctrec                 # complex-terrain receptors

        # Sources params:
        output.nsrctype=nsrctype
        output.nsrcbytype=nsrcbytype

        # Internal:
        output.is_compressed = lcomprs
        output.NCOM=NCOM

        return(output)


        
#def get_ctsg_receptors(self):
               ## ---    Discrete CTSG receptor concentrations
               #if ( nctrec > 0):
               #    rec_size=f.read(4)[0]
               #    #for i in range(self.nctrec):
               #    #    CSPECD = struct.unpack("15s",f.read(15))
               #    #    CONCCT(sp,t,i)=
               #    f.read(4+rec_size) #EOL



#  READ(iunit)nyrb,njulb,nhrb,nsecb,nyre,njule,nhre,nsece
#  READ(iunit) istype,isnum,sname,sxkm,sykm
#  ┌─LOOP OVER OUTPUT SPECIES
#  │  
#  │  GRIDDED RECEPTOR CONCENTRATIONS
#  │  IF(LSGRID)READ(iunit)CSPECG,CONCG
#  │  
#  │  DISCRETE RECEPTOR CONCENTRATIONS
#  │  IF(NDREC.GT.0)READ(iunit)CSPECD,CONCD
#  │  
#  │  COMPLEX TERRAIN RECEPTOR CONCENTRATIONS
#  │  IF(NCTREC.GT.0)READ(iunit)CSPECCT,CONCCT
#  │  
#  └─END LOOP OVER OUTPUT SPECIES
#  
#  where the following declarations apply:
#     character*15 CSPECG,CPSECD,CSPECCT
#     character*16 SNAME
#     real CONCG(nxg,nyg),CONCD(NDREC),CONCCT(NCTREC)
#  and
#     nxg = IESAMP - IBSAMP+1
#     nyg = JESAMP - JBSAMP+1


#-----------------------------------------------------------------
#Skip header data:
#Explanaton:  First 4bytes of a fortran's unformated file is the "record-mark" which 
#           is an int telling how long (in bytes) is the current record without counting
#           the "\n" (new line) symbol which lenths an extra 4 bytes.
#f.seek(f.read(4)[0] + 4)        # Record #1
#f.seek(f.read(4)[0] + 4)        # Record #2: ncom
#for _ in range(self.NCOM):
#    f.seek(f.read(4)[0] + 4, 1) # comments
#f.seek(f.read(4)[0] + 4, 1)     # Record #NCOM+3 : general parameters
#f.seek(f.read(4)[0] + 4, 1)     # Record #NCOM+3a: sources counter
#f.seek(f.read(4)[0] + 4, 1)     # Record #NCOM+4 : title
#f.seek(f.read(4)[0] + 4, 1)     # Record #NCOM+5 : species-groups names list
#f.seek(f.read(4)[0] + 4, 1)     # Record #NCOM+5a: species-groups units list
#f.seek(f.read(4)[0] + 4, 1)     # Record #NCOM+6 : discrete receptors
#f.seek(f.read(4)[0] + 4, 1)     # Record #NCOM+6a: receptors group names
#f.seek(f.read(4)[0] + 4, 1)     # Record #NCOM+7 : complex terrain receptors
#f.seek(f.read(4)[0] + 4, 1)     # Record #NCOM+8 : source names
#-----------------------------------------------------------------
