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
from obspy.clients.fdsn import Client
from obspy.core import Stream, UTCDateTime
from rfpy import options, binning, plotting
from rfpy import HkStack

# Main function


def main():

    # Run Input Parser
    (opts, indb) = options.get_hk_options()

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
        sorted(stkeys)

    # Loop over station keys
    for stkey in list(stkeys):

        # Extract station information from dictionary
        sta = db[stkey]

        # Define path to see if it exists
        datapath = 'DATA/' + stkey
        if not os.path.isdir(datapath):
            raise(Exception('Path to '+datapath+' doesn`t exist - aborting'))

        # Get search start time
        if opts.startT is None:
            tstart = sta.startdate
        else:
            tstart = opts.startT

        # Get search end time
        if opts.endT is None:
            tend = sta.enddate
        else:
            tend = opts.endT

        if tstart > sta.enddate or tend < sta.startdate:
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
        print("|      Start time: {0:19s}          |".format(
            sta.startdate.strftime("%Y-%m-%d %H:%M:%S")))
        print("|      End time:   {0:19s}          |".format(
            sta.enddate.strftime("%Y-%m-%d %H:%M:%S")))
        print("|-----------------------------------------------|")

        rfRstream = Stream()

        for folder in os.listdir(datapath):

            date = folder.split('_')[0]
            year = date[0:4]
            month = date[4:6]
            day = date[6:8]
            dateUTC = UTCDateTime(year+'-'+month+'-'+day)

            if dateUTC > tstart and dateUTC < tend:

                file = open(datapath+"/"+folder+"/RF_Data.pkl", "rb")
                rfdata = pickle.load(file)
                rfRstream.append(rfdata)

            else:
                continue

        # plotting.wiggle(rfRstream, sort='baz')

        # Try binning if specified
        if opts.nbin is not None:
            rf_tmp = binning.bin(rfRstream, typ='slow', nbin=opts.nbin)
            rfRstream = rf_tmp

        # Get a copy of the radial component and filter
        if opts.copy:
            rfRstream_copy = rfRstream.copy()
            rfRstream_copy.filter('bandpass', freqmin=opts.freqs_copy[0],
                                  freqmax=opts.freqs_copy[1], corners=2,
                                  zerophase=True)

        # Filter original stream
        rfRstream.filter('bandpass', freqmin=opts.freqs[0],
                         freqmax=opts.freqs[1], corners=2,
                         zerophase=True)

        # Initialize the HkStack object
        try:
            hkstack = HkStack(rfRstream, rfV2=rfRstream_copy,
                              strike=opts.strike, dip=opts.dip, vp=opts.vp)
        except:
            hkstack = HkStack(rfRstream,
                              strike=opts.strike, dip=opts.dip, vp=opts.vp)

        # Update attributes
        hkstack.hbound = opts.hbound
        hkstack.kbound = opts.kbound
        hkstack.dh = opts.dh
        hkstack.dk = opts.dk
        hkstack.weights = opts.weights

        # Stack with or without dip 
        if opts.calc_dip:
            hkstack.stack_dip()
        else:
            hkstack.stack()

        # Average stacks
        hkstack.average(typ=opts.typ)

        hkstack.plot()

        # Save the hkstack object to file.
        # Add check at beginning to see if file is present.
        # If it is (and overwrite is specified), load it
        # and add the option to simply try another stacking method, weights,
        # and/or plotting

if __name__ == "__main__":

    # Run main program
    main()
