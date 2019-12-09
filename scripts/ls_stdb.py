#!/usr/bin/env python
# encoding: utf-8
''' 
    Program: ls_stdb.py

    Description:
    List Station Database Dictionary contained in pickle file.

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
        print("Error: Must provide a station database pickle file as input.")
        sys.exit()

    else:
        if not osp.exists(sys.argv[1]):
            print("Input File Does not exist...")
            sys.exit()
        else:
            ext=osp.splitext(sys.argv[1])[1]
            if ext==".pkl":
                #-- Pickle Already Created...
                print("Listing Station Pickle: ")
                db = pickle.load(open(sys.argv[1], 'rb'))
                keys=list(db.keys())
                keys.sort()
                ikey=0
                for key in keys:
                    ikey+=1
                    print("--------------------------------------------------------------------------")
                    print("{0:.0f}) {1:s}".format(ikey,key))
                    print("   Station:   ",db[key][0]+"."+db[key][1])
                    print("   Channel:   ",db[key][4])
                    print("   Lon, Lat:  ",db[key][3],db[key][2])
                    print("   TimeRange: ",db[key][5]," - ",db[key][6])
                    print("")
            else:
                print("Error: Must Enter a .pkl station database pickle file")
                sys.exit()
