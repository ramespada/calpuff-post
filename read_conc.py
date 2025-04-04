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

    nsrctype, msource, nrec,nrgrp,nctrec = struct.unpack("5i", f.read(20))
    lsamp= struct.unpack("i", f.read(4))
    nspout, lcomprs, i2dmet = struct.unpack("3i", f.read(12))
    iutmzn,feast,fnorth,rnlat0,relon0,xlat1,xlat2 = struct.unpack("7f", f.read(28))
    proj= f.read(96).decode("utf-8").strip()#pmap,utmhem,datum,daten,clat0,clon0,clat1,clat2
    f.read(4) #close record

    print(f"Model: {MODEL}, Version: {VER}, Level: {LEVEL}")
    print(f"  Start Year: {IBYR}, Julian Day: {IBJUL}, Hour: {IBHR}. TZ ({ABTZ})")
    print(f"  Length of run: {IRLG}, Avg. Time: {IAVG}, TimeStep: {NSECDT} sec.")
    print(f"  Grid size: {NX}x{NY}. Cell size: {DXKM}km x  {DYKM}km. Xmin, Ymin: ({XORIGKM},{YORIGKM})")

    f.read(4) #close record
    nsrcbytype = struct.unpack(str(nsrctype)+"i", f.read(nsrctype*4)) #[nsrcbytype(n),n=1,nsrctype]
    f.read(4) #close record
    #print(nsrcbytype)

    # RECORD #NCOM+4 -- TITLE
    f.read(4) #record marker (record size in bytes)
    TITLE = f.read(80).decode("utf-8").strip()
    print(f"TITLE: {TITLE}")
    f.read(4) #end of record


    # RECORD #NCOM+5 -- List of species-groups output
    #for grp in range(ngrp):

    # RECORD #NCOM+5a-- List of species-groups units 
    # RECORD #NCOM+6 -- Discrete (non-gridded) receptor data
    # RECORD #NCOM+6a-- Receptor-group names
    # RECORD #NCOM+7 -- Complex terrain receptor data
    # RECORD #NCOM+( -- Source names

