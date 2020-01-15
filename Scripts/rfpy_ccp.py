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
from rfpy import CCPimage


def main():

    # Run Input Parser
    (opts, indb) = options.get_ccp_options()

    # Load Database
    db = stdb.io.load_db(fname=indb)

    # Construct station key loop
    allkeys = db.keys()

    # Extract key subset
    if len(opts.stkeys) > 0:
        stkeys = []
        for skey in opts.stkeys:
            stkeys.extend([s for s in allkeys if skey in s])
    else:
        stkeys = db.keys()

    # Initialize CCPimage object
    ccpimage = CCPimage(coord_start=opts.coord_start,
                        coord_end=opts.coord_end)

    # Loop over station keys
    for stkey in list(stkeys):

        # Extract station information from dictionary
        sta = db[stkey]

        # Define path to see if it exists
        datapath = 'DATA/' + stkey
        if not os.path.isdir(datapath):
            raise(Exception('Path to '+datapath+' doesn`t exist - aborting'))

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

        rfRstream = Stream()

        for folder in os.listdir(datapath):

            filename = datapath+"/"+folder+"/RF_Data.pkl"
            if os.path.isfile(filename):
                file = open(filename, "rb")
                rfdata = pickle.load(file)
                if rfdata[1].stats.snr > opts.snr:
                    rfRstream.append(rfdata[1])
                file.close()

        ccpimage.add_rfstream(rfRstream)

    ccpimage.prep_data(f1=opts.f1, f2ps=opts.f2ps,
                       f2pps=opts.f2pps, f2pss=opts.f2pss,
                       n_depth=opts.n_depth, nbaz=opts.nbaz,
                       nslow=opts.nslow)

    ccpimage.prestack(cell_length=opts.cell_length)

    if opts.run_ccp:
        ccpimage.ccp()
        ccpimage.stack_ccp()
        ccpimage.plot_ccp()

    if opts.run_gccp:
        ccpimage.gccp()
        ccpimage.pws_gccp()
        ccpimage.plot_gccp()


if __name__ == "__main__":

    # Run main program
    main()
