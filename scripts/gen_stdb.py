#!/usr/bin/env python
# encoding: utf-8
''' 
    Program: gen_stdb.py

    Description:
    Create Station Database dictionary based on an input file in one of two
    formats:
    1) chS csv:          -------Start---------   --------END----------
     NET,STA,locID,CHAN,YYYY-MM-DD,HH:MM:SS.SSS,YYYY-MM-DD,HH:MM:SS.SSS,lat,lon
    
    2) IPO SPC                 --Start--- ---END----
     NET STA CHAN lat lon elev YYYY-MM-DD YYYY-MM-DD

    The station dictionary contains keys which are named NET.STA.CHAN, where CHAN
    is a two character representation of the desired channel (ex, BH, HH, LH).
    Within each KEY is the set of data used in later programs to define the 
    station information. The data is stored in the pickle file as a dictionary,
    with each key having a list:
        stdb[stkey] = [NET, STA, lat, lon, CHAN, UTCDateTime(Start Year), UTCDateTime(End Year)
    
    When the station database pickle is read in by subsequent programs, the 
    contents of each key are converted into a struct:
        stdb[stkey]:
                .network
                .station
                .stla
                .stlo
                .cha
                .dstart
                .dend

'''


try:
    import pickle as pickle
except:
    import pickle
import sys
import os.path as osp
from obspy.core import UTCDateTime


if __name__=='__main__':

    #-- Check if called with inputs
    if  len(sys.argv)!=2:
        print("Error: Must provide a station list as input.\n")
        print("    Description:")
        print("    Create Station Database dictionary based on an input file in one of two")
        print("    formats:")
        print("    1) chS csv:          -------Start---------   --------END----------")
        print("     NET,STA,locID,CHAN,YYYY-MM-DD,HH:MM:SS.SSS,YYYY-MM-DD,HH:MM:SS.SSS,lat,lon")
        print("    ")
        print("    2) IPO SPC                 --Start--- ---END----")
        print("     NET STA CHAN lat lon elev YYYY-MM-DD YYYY-MM-DD")
        print("")
        print("    The station dictionary contains keys which are named NET.STA.CHAN, where CHAN")
        print("    is a two character representation of the desired channel (ex, BH, HH, LH).")
        print("    Within each KEY is the set of data used in later programs to define the ")
        print("    station information. The data is stored in the pickle file as a dictionary,")
        print("    with each key having a list:")
        print("      stdb[stkey] = [NET, STA, lat, lon, CHAN, UTCDateTime(Start Year), UTCDateTime(End Year)")
        print("    ")
        print("    When the station database pickle is read in by subsequent programs, the ")
        print("    contents of each key are converted into a struct:")
        print("      stdb[stkey]:")
        print("          .network")
        print("          .station")
        print("          .stla")
        print("          .stlo")
        print("          .cha")
        print("          .dstart")
        print("          .dend")
        sys.exit()

    else:
        if not osp.exists(sys.argv[1]):
            print("Input File Does not exist...")
            sys.exit()
        else:
            ext=osp.splitext(sys.argv[1])[1]

            if ext!=".pkl":
                #-- Station List...Pickle it.
                print("Parse Station List "+sys.argv[1])
                pklfile=sys.argv[1]+".pkl"
                
                fin=open(sys.argv[1],'r')
                stations={}
                for line in fin:
                    line=line.strip()
                    if len(line.split(','))>6:
                        line=line.split(','); net=line[0]; stn=line[1]; loc=line[2];chn=line[3][0:2]
                        stdt=line[4]; sttm=line[5]; eddt=line[6]; edtm=line[7]
                        lat=float(line[8]); lon=float(line[9])
                    elif len(line.split())>6:
                        line=line.split(); net=line[0]; stn=line[1]; chn=line[2][0:2]
                        stdt=line[6]; eddt=line[7]
                        lat=float(line[3]); lon=float(line[4])

                    #-- Now Add lines to station Dictionary
                    key="{0:s}.{1:s}.{2:2s}".format(net.strip(),stn.strip(),chn.strip())
                    if key not in stations:
                        stations[key]=[net,stn,lat,lon,chn,UTCDateTime(stdt),UTCDateTime(eddt)]
                        print("Adding key: "+key)
                    else:
                        print("Warning: Key "+key+" already exists...Skip")

                #-- Save and Pickle the station database
                print("  Pickling to "+pklfile)
                pickle.dump(stations, open(pklfile, 'wb'))
            else:
                print("Error: Must supply a station list, not a Pickle File")
                sys.exit()
