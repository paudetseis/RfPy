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
from rfpy import options
from rfpy import RFData

# Main function


def main():

    # Run Input Parser
    (opts, indb) = options.get_calc_options()

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
            print('Path to '+datapath+' doesn`t exist - creating it')
            os.makedirs(datapath)

        # Establish client for data
        if len(opts.UserAuth) == 0:
            data_client = Client(opts.Server)
        else:
            data_client = Client(opts.Server, user=opts.UserAuth[0],
                                 password=opts.UserAuth[1])

        # Establish client for events
        event_client = Client()

        # Get catalogue search start time
        if opts.startT is None:
            tstart = sta.startdate
        else:
            tstart = opts.startT

        # Get catalogue search end time
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
        print("| Searching Possible events:                    |")
        print("|   Start: {0:19s}                  |".format(
            tstart.strftime("%Y-%m-%d %H:%M:%S")))
        print("|   End:   {0:19s}                  |".format(
            tend.strftime("%Y-%m-%d %H:%M:%S")))
        if opts.maxmag is None:
            print("|   Mag:   >{0:3.1f}", format(opts.minmag) +
                  "                                 |")
        else:
            msg = "|   Mag:   {0:3.1f}".format(opts.minmag) + \
                " - {0:3.1f}".format(opts.maxmag) + \
                "                            |"
            print(msg)

        print("| ...                                           |")

        # Get catalogue using deployment start and end
        cat = event_client.get_events(
            starttime=tstart, endtime=tend,
            minmagnitude=opts.minmag, maxmagnitude=opts.maxmag)

        # Total number of events in Catalogue
        nevK = 0
        nevtT = len(cat)
        print(
            "|  Found {0:5d}".format(nevtT) +
            " possible events                  |")
        ievs = range(0, nevtT)

        # Get Local Data Availabilty
        if len(opts.localdata) > 0:
            print("|-----------------------------------------------|")
            print("| Cataloging Local Data...                      |")
            if opts.useNet:
                stalcllist = options.list_local_data_stn(
                    lcldrs=opts.localdata, sta=sta.station,
                    net=sta.network, altnet=sta.altnet)
                print("|   {0:>2s}.{1:5s}: {2:6d}".format(
                    sta.network, sta.station, len(stalcllist)) +
                    " files                      |")
                print(stalcllist[0:10])
            else:
                stalcllist = options.list_local_data_stn(
                    lcldrs=opts.localdata, sta=sta.station)
                print("|   {0:5s}: {1:6d} files                " +
                      "        |".format(
                          sta.station, len(stalcllist)))
        else:
            stalcllist = []
        print("|===============================================|")

        # Select order of processing
        if opts.reverse:
            ievs = range(0, nevtT)
        else:
            ievs = range(nevtT-1, -1, -1)

        # Read through catalogue
        for iev in ievs:

            # Extract event
            ev = cat[iev]

            # Initialize RF object with station info
            rfdata = RFData(sta)

            # Add event to Split object
            accept = rfdata.add_event(
                ev, gacmin=opts.mindist, gacmax=opts.maxdist, returned=True)

            # Define time stamp
            yr = str(rfdata.meta.time.year).zfill(4)
            jd = str(rfdata.meta.time.julday).zfill(3)
            hr = str(rfdata.meta.time.hour).zfill(2)

            # If distance between mindist and maxdist deg:
            if accept:

                # Display Event Info
                nevK = nevK + 1
                if opts.reverse:
                    inum = iev + 1
                else:
                    inum = nevtT - iev + 1
                print(" ")
                print("**************************************************")
                print("* #{0:d} ({1:d}/{2:d}):  {3:13s}".format(
                    nevK, inum, nevtT, rfdata.meta.time.strftime(
                        "%Y%m%d_%H%M%S")))
                print("*   Origin Time: " +
                      rfdata.meta.time.strftime("%Y-%m-%d %H:%M:%S"))
                print(
                    "*   Lat: {0:6.2f};        Lon: {1:7.2f}".format(
                        rfdata.meta.lat, rfdata.meta.lon))
                print(
                    "*   Dep: {0:6.2f} km;     Mag: {1:3.1f}".format(
                        rfdata.meta.dep, rfdata.meta.mag))
                print("*   Dist: {0:7.2f} km;".format(rfdata.meta.epi_dist) +
                      "   Epi dist: {0:6.2f} deg\n".format(rfdata.meta.gac) +
                      "*   Baz:  {0:6.2f} deg;".format(rfdata.meta.baz) +
                      "   Az: {0:6.2f} deg".format(rfdata.meta.az))

                # Event Folder
                timekey = rfdata.meta.time.strftime("%Y%m%d_%H%M%S")
                evtdir = datapath + "/" + timekey

                # Check if RF data already exist and overwrite has been set
                if os.path.isdir(evtdir):
                    if os.path.isfile(evtdir+"/RF_Data.pkl"):
                        if not opts.ovr:
                            continue

                # Get data
                has_data = rfdata.download_data(
                    client=data_client, dts=opts.dts, stdata=stalcllist,
                    ndval=opts.ndval, new_sr=opts.new_sampling_rate,
                    returned=True)

                if not has_data:
                    continue

                # Create Folder if it doesn't exist
                if not os.path.isdir(evtdir):
                    os.makedirs(evtdir)

                # Save ZNE Traces
                pickle.dump(rfdata.data, open(
                    evtdir + "/ZNE_Data.pkl", "wb"))

                # Rotate from ZNE to 'align' ('ZRT', 'LQT', or 'PVH')
                rfdata.rotate(vp=opts.vp, vs=opts.vs, align=opts.align)

                # Calculate snr over dt_snr seconds
                rfdata.calc_snr(
                    dt=opts.dt_snr, fmin=opts.fmin, fmax=opts.fmax)
                print("* SNR: {}".format(rfdata.meta.snr))

                # Make sure no processing happens for NaNs
                if np.isnan(rfdata.meta.snr):
                    print("* SNR NaN...Skipping")
                    print("**************************************************")
                    continue

                # Deconvolve data
                rfdata.deconvolve(
                    twin=opts.twin, vp=opts.vp, vs=opts.vs,
                    align=opts.align, method=opts.method)

                # Get cross-correlation QC
                rfdata.get_QC()
                print("* CC: {}".format(rfdata.meta.cc))

                # Convert to Stream
                rfstream = rfdata.to_stream()

                # Save event meta data
                pickle.dump(rfdata.meta, open(
                    evtdir + "/Meta_Data.pkl", "wb"))

                # Save Station Data
                pickle.dump(rfdata.sta, open(
                    evtdir + "/Station_Data.pkl", "wb"))

                # Save RF Traces
                pickle.dump(rfstream, open(
                    evtdir + "/RF_Data.pkl", "wb"))

                # Update
                print("* Wrote Output Files to: ")
                print("*     "+evtdir)
                print("**************************************************")


if __name__ == "__main__":

    # Run main program
    main()
