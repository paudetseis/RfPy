#!/usr/bin/env python

# Copyright 2019 Pascal Audet
#
# This file is part of RfPy.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


# Import modules and functions
import numpy as np
import os.path
import pickle
import glob
import stdb
from rfpy import options
from rfpy import RFData


def main():

    print()
    print("########################################################")
    print("#                                                      #")
    print("#        __                                    _       #")
    print("#  _ __ / _|_ __  _   _     _ __ ___  ___ __ _| | ___  #")
    print("# | '__| |_| '_ \| | | |   | '__/ _ \/ __/ _` | |/ __| #")
    print("# | |  |  _| |_) | |_| |   | | |  __/ (_| (_| | | (__  #")
    print("# |_|  |_| | .__/ \__, |___|_|  \___|\___\__,_|_|\___| #")
    print("#          |_|    |___/_____|                          #")
    print("#                                                      #")
    print("########################################################")
    print()

    # Run Input Parser
    (opts, indb) = options.get_recalc_options()

    # Load Database
    db = stdb.io.load_db(fname=indb)

    # Construct station key loop
    allkeys = db.keys()
    sorted(allkeys)

    # Extract key subset
    if len(opts.stkeys) > 0:
        stkeys = []
        for skey in opts.stkeys:
            stkeys.extend([s for s in allkeys if skey in s])
    else:
        stkeys = db.keys()

    # Loop over station keys
    for stkey in list(stkeys):

        # Extract station information from dictionary
        sta = db[stkey]

        # Define path to see if it exists
        if opts.phase in ['P', 'PP', 'allP']:
            datapath = 'P_DATA/' + stkey
        elif opts.phase in ['S', 'SKS', 'allS']:
            datapath = 'S_DATA/' + stkey
        if not os.path.isdir(datapath):
            print('Path to ' + datapath + ' doesn`t exist - continuing')
            continue

        # Temporary print locations
        tlocs = sta.location
        if len(tlocs) == 0:
            tlocs = ['']
        for il in range(0, len(tlocs)):
            if len(tlocs[il]) == 0:
                tlocs[il] = "--"
        sta.location = tlocs

        # Update Display
        print(" ")
        print(" ")
        print("|===============================================|")
        print("|===============================================|")
        print("|                   {0:>8s}                    |".format(
            sta.station))
        print("|===============================================|")
        print("|===============================================|")
        print("|  Station: {0:>2s}.{1:5s}                            |".format(
            sta.network, sta.station))
        print("|      Channel: {0:2s}; Locations: {1:15s}  |".format(
            sta.channel, ",".join(tlocs)))
        print("|      Lon: {0:7.2f}; Lat: {1:6.2f}                |".format(
            sta.longitude, sta.latitude))
        print("|-----------------------------------------------|")

        for folder in os.listdir(datapath):

            # Re-initialize RFData object
            rfdata = RFData(sta)

            # Load meta data
            if not os.path.isfile(datapath+"/"+folder+"/Meta_Data.pkl"):
                continue
            rfdata.meta = pickle.load(open(
                datapath+"/"+folder+"/Meta_Data.pkl",'rb'))

            # Skip data not in list of phases
            if rfdata.meta.phase not in opts.listphase:
                continue

            if opts.verb:
                print("* Station: {0}; folder: {1}".format(stkey,folder))

            # if not hasattr(rfdata.meta, 'phase'):
            #     rfdata.meta.phase = 'P'

            # Load ZNE data
            rfdata.data = pickle.load(open(
                datapath+"/"+folder+"/ZNE_Data.pkl",'rb'))

            # Remove rotated flag and snr flag
            rfdata.meta.rotated = False
            rfdata.meta.snr = None

            # Rotate from ZNE to 'align' ('ZRT', 'LQT', or 'PVH')
            rfdata.rotate(vp=opts.vp, vs=opts.vs, align=opts.align)

            # Calculate SNR
            rfdata.calc_snr(dt=opts.dt_snr, fmin=opts.fmin, fmax=opts.fmax)
    
            if opts.verb:
                print("* SNR: {}".format(rfdata.meta.snr))

            # Deconvolve data
            rfdata.deconvolve(
                vp=opts.vp, vs=opts.vs, 
                align=opts.align, method=opts.method,
                gfilt=opts.gfilt, wlevel=opts.wlevel,
                pre_filt=opts.pre_filt)

            # Get cross-correlation QC
            rfdata.calc_cc()

            if opts.verb:
                print("* CC: {}".format(rfdata.meta.cc))

            # Convert to Stream
            rfstream = rfdata.to_stream()

            # Save RF Traces
            pickle.dump(rfstream, open(
                datapath+"/"+folder+"/RF_Data.pkl", "wb"))

            # Update
            if opts.verb:
                print("* Output files written")
                print("**************************************************")


if __name__ == "__main__":

    # Run main program
    main()
