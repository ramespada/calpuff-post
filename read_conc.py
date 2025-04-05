#!/bin/python3

import struct

file = "conc.dat"

# Open the binary file
with open(file, "rb") as f:

    # RECORD #1 --  DATASET(16),DATAVER(16),DATAMOD(64)
    f.read(4) #record marker (record size in bytes)
    FIRST_REC = f.read(96).decode("utf-8").strip()
    print(f"{FIRST_REC}")
    f.read(4) #end of record

    # RECORD #2 --  NCOM
    f.read(4) #record marker
    NCOM = struct.unpack("i", f.read(4))[0]
    print(f"NCOM: {NCOM}")
    f.read(4) #end of record

    # RECORD #3 to NCOM+2 --  Comments.
    for i in range(NCOM):
        rec_size=f.read(4)[0] #record marker
        f.read(rec_size) #skip commented lines.
        f.read(4)        #end of record

    # RECORD #NCOM+3 -- General run parameters

    f.read(4) #rec_size = struct.unpack("i", f.read(4))[0]  
    #print(rec_size)
    MODEL = f.read(12).decode("utf-8").strip()
    VER   = f.read(12).decode("utf-8").strip()
    LEVEL = f.read(12).decode("utf-8").strip()

    IBYR, IBJUL, IBHR, IBSEC = struct.unpack("4i", f.read(16))
    ABTZ = f.read(8).decode("utf-8").strip()
    IRLG, IAVG, NSECDT = struct.unpack("3i", f.read(12))
    NX,NY = struct.unpack("2i",f.read(8))
    DXKM,DYKM,IONE,XORIGKM,YORIGKM = struct.unpack("5f", f.read(20))

    NSSTA = struct.unpack("i", f.read(4))                            # n ssurface stations
    IBCOMP, IECOMP, JBCOMP, JECOMP = struct.unpack("4i", f.read(16)) #computational grid
    IBSAMP, JBSAMP, IESAMP, JESAMP = struct.unpack("4i", f.read(16)) #sampling grid
    MESHDN = struct.unpack("i", f.read(4))                           #sampling grid spacing factor

    nsrctype, msource, nrec,nrgrp, nctrec = struct.unpack("5i", f.read(20))
    lsamp= struct.unpack("i", f.read(4))
    nspout, lcomprs, i2dmet,iutmzn = struct.unpack("4i", f.read(16))
    feast,fnorth,rnlat0,relon0,xlat1,xlat2 = struct.unpack("6f", f.read(24))
    proj= f.read(96).decode("utf-8").strip()#pmap,utmhem,datum,daten,clat0,clon0,clat1,clat2
    f.read(4) #close record

    print(f"Model: {MODEL}, Version: {VER}, Level: {LEVEL}")
    print(f"  Start Year: {IBYR}, Julian Day: {IBJUL}, Hour: {IBHR} ({ABTZ})")
    print(f"  Length of run: {IRLG}, Avg. Time: {IAVG}, TimeStep: {NSECDT} sec.")
    print(f"  Grid size: {NX}x{NY}. Cell size: {DXKM}km x  {DYKM}km. Xmin, Ymin: ({XORIGKM},{YORIGKM})")

    print("Proj:",iutmzn,feast,fnorth,rnlat0,relon0,xlat1,xlat2 )
    print(proj)

    f.read(4) #close record
    nsrcbytype = struct.unpack(str(nsrctype)+"i", f.read(nsrctype*4)) #[nsrcbytype(n),n=1,nsrctype]
    f.read(4) #close record
    print(nsrcbytype)

    # RECORD #NCOM+4 -- TITLE
    f.read(4) #record marker (record size in bytes)
    TITLE = f.read(80).decode("utf-8").strip()
    print(f"TITLE: {TITLE}")
    f.read(4) #end of record

    # RECORD #NCOM+5 -- List of species-groups output
    #for grp in range(ngrp):
    f.read(4) #close record
    csout=[""]*nspout
    for i in range(nspout):
       csout[i] = f.read(15).decode("utf-8").strip()
    f.read(4) #close record
    print(csout)

    # RECORD #NCOM+5a-- List of species-groups units 
    f.read(4) #close record
    acunit=[""]*nspout
    for i in range(nspout):
      acunit[i] = f.read(16).decode("utf-8").strip()
    f.read(4) #close record
    print(acunit)

    if (nrec > 0):
        # RECORD #NCOM+6 -- Discrete (non-gridded) receptor data
        f.read(4) #close record
        for i in range(nrec):
            x[i],y[i],z[i],zng[i],irgrp[i] = struct.unpack("5f", f.read(20))
        f.read(4) #close record
        
        # RECORD #NCOM+6a-- Receptor-group names
        f.read(4) #close record
        for i in range(nrec):
            rgrpnam[i] = f.read(80).decode("utf-8").strip()
        f.read(4) #close record

    print(nctrec)
    # RECORD #NCOM+7 -- Complex terrain receptor data
    if (nctrec > 0):
        f.read(4) #close record
        for i in range(nctrec):
            x[i],y[i],z[i],ihill[i] = struct.unpack("5f", f.read(20))
        f.read(4) #close record

        print(x,y,z,ihill)

    # RECORD #NCOM+8( -- Source names
    f.read(4) #close record
    for itype in range(nsrctype):
        if ( nsrcbytype(itype) > 0):
            cnamsrc[itype] = f.read(16).decode("utf-8").strip()
    f.read(4) #close record

    # RECORD #NCOM+9 -- Source names



    ##MAIN DATA:
    #write(io8)nyrab,njulab,nhrab,nsecab,nyre,njule,nhre,nsece
    #write(io8)ktype,ksource,csrcnam,xmapkm,ymapkm

    for t in range(IRLG): #temporal loop
        for src in range(nspec):
            for grp in range(ngrup)

                for i in range(
    yr1,dy1,hr1,sec1 = struct.unpack("4i", f.read(16))
    yr2,dy2,hr2,sec2 = struct.unpack("4i", f.read(16))
    f.read(4)

    istype,isnum     = struct.unpack("2i", f.read(8))
    sname            = f.read(16).decode("utf-8").split()
    sxkm,sykm        = struct.unpack("2i", f.read(8))
    f.read(4)

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


#c ---------------------------------------------------------------
#c --- WRITE CONCENTRATIONS TO DISK   (g/m**3,odour_units,Bq/m**3)
#c ---------------------------------------------------------------
#c
#c --- Output date/hour times use 0-23 convention so that hour 24
#c --- of day 12 starts at 23 0000 on day 12 and ends at 00 0000
#c --- on day 13.
#      if(icon.ne.1)go to 492
#c
#c --- Write date/time and source data records
#      write(io8)nyrab,njulab,nhrab,nsecab,nyre,njule,nhre,nsece
#      write(io8)ktype,ksource,csrcnam,xmapkm,ymapkm
#c
#      do 400 ig=1,ngrup
#c --- Identify array storage location for this group
#      i=istore(ig)
#c
#c --- Only species-groups specified are stored on disk
#      if(ioutop(2,ig).eq.1)then
#         cname=cgrup(ig)
#         cname(13:15)='  1'
#c
#c ---    Gridded receptor concentrations
#         if(lsamp) then
#            if(ifull.eq.1)then
#               if(lcomprs)then
#c ---             Write compressed data records
#                  call comprs(chisam(1,1,i),mxnxyg,tmp8,mxnxyg,
#     1              cname,io8)
#               else
#c ---             Write uncompressed data record
#                  call wrdat(io8,cname,chisam(1,1,i),nxsam,nysam)
#               endif
#            else
#               call xtract(chisam(1,1,i),mxnxg,mxnyg,nxsam,nysam,tmp7)
#               if(lcomprs)then
#c ---             Write compressed data records
#                  nwords=nxsam*nysam
#                  call comprs(tmp7,nwords,tmp8,mxnxyg,cname,io8)
#               else
#c ---             Write uncompressed data record
#                  call wrdat(io8,cname,tmp7,nxsam,nysam)
#               endif
#            endif
#         endif
#c
#c ---    Discrete receptor concentrations
#         if(nrec .GT. 0) then
#            if(lcomprs)then
#c ---          Write compressed data records
#               call comprs(chirec(1,i),nrec,tmp3,mxrec,cname,io8)
#            else
#c ---          Write uncompressed data record
#               call wrdat(io8,cname,chirec(1,i),nrec,1)
#            endif
#         endif
#c
#c ---    Discrete CTSG receptor concentrations
#         if(nctrec .GT. 0) then
#            if(lcomprs)then
#c ---          Write compressed data records
#               call comprs(chict(1,i),nctrec,tmp5,mxrect,cname,io8)
#            else
#c ---          Write uncompressed data record
#               call wrdat(io8,cname,chict(1,i),nctrec,1)
#            endif
#         endif
#      endif
#400   continue
#492   continue

