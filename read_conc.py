#!/bin/python3

import struct
import numpy as np
from matplotlib import pyplot as plt
file = "conc.dat"

def decompress(xwork):
    xdat = []
    for value in xwork:
        if value > 0.0:
            xdat.append(value)
        else:
            xdat.extend([0.0] * int(-value))
    return xdat

# Open the binary file
with open(file, "rb") as f:

    # RECORD #1 --  
    #write(io8) conset,dataver,datamod
    f.read(4) # record marker (record size in bytes)
    f.read(96)#FIRST_REC = f.read(96).decode("utf-8").strip()
    f.read(4) #end of line (EOL) record

    # RECORD #2 --  NCOM     
    #write(io8) ncom
    f.read(4) #record marker
    NCOM = struct.unpack("i", f.read(4))[0]
    f.read(4) #end of record
    # RECORD #3 to NCOM+2 --  Comments.
    for i in range(NCOM):
        rec_size=f.read(4)[0] #record marker
        f.read(rec_size) #skip commented lines.
        f.read(4)        #EOL

    # RECORD #NCOM+3 -- General run parameters
    #write(io8)mmodel,ver,level,ibyr,ibjul,ibhr,ibsec,
    #1 abtz,irlg,iavg,nsecdt,nx,ny,dxkm,dykm,ione,xorigkm,yorigkm,
    #2 nssta(1),ibcomp,iecomp,jbcomp,jecomp,ibsamp,jbsamp,iesamp,
    #3 jesamp,meshdn,nsrctype,msource,nrec,nrgrp,
    #4 nctrec,lsamp,nspout,lcomprs,i2dmet,
    #5 iutmzn,feast,fnorth,rnlat0,relon0,xlat1,xlat2,
    #6 pmap,utmhem,datum,daten,clat0,clon0,clat1,clat2
    f.read(4) 
    MODEL, VER, LEVEL        = struct.unpack("12s 12s 12s", f.read(36))
    IBYR, IBJUL, IBHR, IBSEC, ABTZ = struct.unpack("4i 8s", f.read(24))
    IRLG, IAVG, NSECDT = struct.unpack("3i", f.read(12))
    NX,NY, DXKM,DYKM, IONE, XORIGKM,YORIGKM  = struct.unpack("2i 2f i 2f",f.read(28))

    NSSTA = struct.unpack("i", f.read(4))                            #n ssurface stations
    IBCOMP, IECOMP, JBCOMP, JECOMP = struct.unpack("4i", f.read(16)) #computational grid
    IBSAMP, JBSAMP, IESAMP, JESAMP = struct.unpack("4i", f.read(16)) #sampling grid
    MESHDN = struct.unpack("i", f.read(4))                           #sampling grid spacing factor

    nsrctype, msource, nrec, nrgrp, nctrec  = struct.unpack("5i", f.read(20)) #
    lsamp, nspout, lcomprs, i2dmet          = struct.unpack("4i", f.read(16)) # sample grid?. num of species. data compresed?, rel humidity source data?
    iutmzn, feast,fnorth,rnlat0,relon0,xlat1,xlat2 = struct.unpack("i 6f", f.read(28))       # proj parameters.
    proj= f.read(96).decode("utf-8").strip()#pmap,utmhem,datum,daten,clat0,clon0,clat1,clat2 # proj parameters string form.

    f.read(4) #close record

    print(f"Model: {MODEL}, Version: {VER}, Level: {LEVEL}")
    print(f"   Start Date. Year: {IBYR}, Julian Day: {IBJUL}, Hour: {IBHR} ({ABTZ})")
    print(f"   Run Lenght: {IRLG}. Avg. Time: {IAVG}. TimeStep: {NSECDT} secs.")
    print(f"   Grid size: {NX}x{NY}. Cell size: {DXKM}km x  {DYKM}km. Xmin, Ymin: ({XORIGKM},{YORIGKM})")
    print( "   Proj. parameters:",iutmzn,feast,fnorth,rnlat0,relon0,xlat1,xlat2 )
    print(f"   {proj}")
    print(f"   Number of species-groups: {nspout}")
    nxg = IESAMP - IBSAMP+1
    nyg = JESAMP - JBSAMP+1 
    print(f"   Number of gridded         receptors: {nxg} x {nyg}")
    print(f"   Number of discrete        receptors: {nspout}")
    print(f"   Number of complex terrain receptors: {nctrec}")

    # NON-REPORTED RECORD: how many sources from each type
    #write(io8) (nsrcbytype(n),n=1,nsrctype)
    f.read(4) #close record
    nsrcbytype = struct.unpack(str(nsrctype)+"i", f.read(nsrctype*4)) #[nsrcbytype(n),n=1,nsrctype]
    f.read(4) #close record
    print(nsrcbytype)

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
    print(f"   Species groups: {csout}")

    # RECORD #NCOM+5a-- List of species-groups units 
    #write(io8)(acunit(n),n=1,nspout)
    f.read(4) 
    acunit=[""]*nspout
    for i in range(nspout):
      acunit[i] = f.read(16).decode("utf-8").strip()
    f.read(4) #EOL
    print(f"   Species units:  {acunit}")


    # RECORD #NCOM+6 -- Discrete (non-gridded) receptor data
    if (nrec > 0):
        print(f"   Number of discrete receptors: {nrec}")
        f.read(4) #close record
        for i in range(nrec):
            x[i],y[i],z[i],zng[i],irgrp[i] = struct.unpack("5f", f.read(20))
        f.read(4) #EOL
        
        # RECORD #NCOM+6a -- Receptor-group names
        f.read(4) #close record
        for i in range(nrec):
            rgrpnam[i] = f.read(80).decode("utf-8").strip()
        f.read(4) #EOL

    # RECORD #NCOM+7 -- Complex terrain receptor data
    if (nctrec > 0):
        print(f"   Number of complex terrain receptors: {nctrec}")
        f.read(4) #close record
        for i in range(nctrec):
            x[i],y[i],z[i],ihill[i] = struct.unpack("5f", f.read(20))
        f.read(4) #EOL
        print(x,y,z,ihill)

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
                print(cnamsrc)
            f.read(4) #EOL
    print(f"   Source names: {cnamsrc})")

    # RECORD #NCOM+8 -- Nearest Surface Station for VISIBILITY ONLY
    # RECORD #NCOM+9 -- X coord (UTM) of stations for VISIBILITY ONLY
    # RECORD #NCOM+10 -- Y coord (UTM) of stations for VISIBILITY ONLY

    #-----------------------------------------------------------------
    ##MAIN DATA:
    for t in range(IRLG): #temporal loop
        #write(io8)nyrab,njulab,nhrab,nsecab,nyre,njule,nhre,nsece
        f.read(4)[0]
        yr1,dy1,hr1,sec1, yr2,dy2,hr2,sec2 = struct.unpack("8i", f.read(32))
        f.read(4) #EOL

        print(f"Record date: {yr1} - {dy1} - {hr1} hour.")
        print(f"Record date: {yr2} - {dy2} - {hr2} hour.")

        #write(io8)ktype,ksource,csrcnam,xmapkm,ymapkm
        f.read(4)[0]
        ktype,ksource,csrcnam,xmapkm,ymapkm = struct.unpack("2i 16s 2f", f.read(32))
        f.read(4)
        print(f"2nd line: {ktype},{ksource},{csrcnam},{xmapkm},{ymapkm})")       

        #Loop on output species
        for sp in range(nspout): 

            # ---    Gridded receptor concentrations
            if ( lsamp ):
                #if compressed data
                if ( lcomprs ):
                    rec_size = f.read(4)[0]
                    ii = struct.unpack("i", f.read(4))[0]      #header of compressed file, indicating record length "number of words".
                    f.read(4)
                    print(f" (( REC_SIZE {rec_size} )). II = {ii}")

                    f.read(4)
                    CSPECG = f.read(15).decode("utf-8").split() #struct.unpack("15s", f.read(15)) #
                    raw = f.read(4 * ii)
                    xwork = list(struct.unpack(f'{ii}f', raw))
                    xdat = decompress(xwork)
                    CONCG = np.array(xdat, dtype=np.float32).reshape((nyg, nxg), order='F')
                    f.read(4) #EOL

                else:
                    rec_size=f.read(4)[0]
                    CSPECG = f.read(15).decode("utf-8").split() #struct.unpack("15s", f.read(15)) #
                    print(f" (( REC_SIZE {rec_size} )). SPECIE= {CSPECG}")
                    CONCG = np.array([[ struct.unpack("f",f.read(4)) for j in range(nyg)] for i in range(nxg)])
                    f.read(4) #EOL

                print("CONCG =")
                print(CONCG)
                plt.imshow(CONCG)
                plt.show()


            # ---    Discrete receptor concentrations
            if ( nrec > 0 ):
                rec_size=f.read(4)[0]
                #for i in range(NREC):
                #    CSPECD = struc.unpack("15s",f.read(15))
                #    CONCCD(sp,t,i)=
                f.read(4+rec_size) #EOL

            # ---    Discrete CTSG receptor concentrations
            if ( nctrec > 0):
                rec_size=f.read(4)[0]
                #for i in range(NCTREC):
                #    CSPECD = struc.unpack("15s",f.read(15))
                #    CONCCT(sp,t,i)=
                f.read(4+rec_size) #EOL


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

