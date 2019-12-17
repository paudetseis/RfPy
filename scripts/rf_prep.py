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


"""
"""

# Import modules and functions
import numpy as np
import os.path
import pickle
import glob
import stdb
from obspy.taup import TauPyModel
from obspy.clients.fdsn import Client
from rfpy import options
from rfpy import RFData

# Main function


def main():

    # Run Input Parser
    (opts, indb) = options.get_prep_options()

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

    # Initialize Taup Model
    tpmodel = TauPyModel(model='iasp91')

    # Loop over station keys
    for stkey in list(stkeys):

        # Extract station information from dictionary
        sta = db[stkey]

        # Initialize RF object with station info
        rfdata = RFData(sta, opts.vp, opts.vs, opts.align)

        # Define path to see if it exists
        datapath = 'DATA/' + stkey
        if not os.path.isdir(datapath):
            print('Path to '+datapath+' doesn`t exist - creating it')
            os.makedirs(datapath)

        # Establish client
        if len(opts.UserAuth) == 0:
            client = Client(opts.Server)
        else:
            client = Client(
                opts.Server,
                user=opts.UserAuth[0],
                password=opts.UserAuth[1])

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
        print("|  Station: {0:>2s}.{1:5s}                      " +
              "      |".format(
                  sta.network, sta.station))
        print("|      Channel: {0:2s}; Locations: {1:15s}  |".format(
            sta.channel, ",".join(tlocs)))
        print("|      Lon: {0:7.2f}; Lat: {1:6.2f}              " +
              "  |".format(
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
            print("|   Mag:   >{0:3.1f}                         " +
                  "        |".format(
                      opts.minmag))
        else:
            print(
                "|   Mag:   {0:3.1f} - {1:3.1f}                 " +
                "           |".format(opts.minmag, opts.maxmag))
        print("| ...                                           |")

        # Get catalogue using deployment start and end
        cat = client.get_events(
            starttime=tstart, endtime=tend,
            minmagnitude=opts.minmag, maxmagnitude=opts.maxmag)

        # Total number of events in Catalogue
        nevK = 0
        nevtT = len(cat)
        print(
            "|  Found {0:5d} possible events                  |".format(nevtT))
        ievs = range(0, nevtT)

        # Get Local Data Availabilty
        if len(opts.localdata) > 0:
            print("|-----------------------------------------------|")
            print("| Cataloging Local Data...                      |")
            if opts.useNet:
                stalcllist = rfdata.io.list_local_data_stn(
                    lcldrs=opts.localdata, sta=sta.station,
                    net=sta.network, altnet=sta.altnet)
                print("|   {0:>2s}.{1:5s}: {2:6d} files         " +
                      "             |".format(
                          sta.network, sta.station, len(stalcllist)))
            else:
                stalcllist = rfdata.io.list_local_data_stn(
                    lcldrs=opts.localdata, sta=sta.station)
                print("|   {0:5s}: {1:6d} files                 " +
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

            # Add event to Split object
            rfdata.add_event(ev)

            # Define time stamp
            yr = str(rfdata.meta.time.year).zfill(4)
            jd = str(rfdata.meta.time.julday).zfill(3)
            hr = str(rfdata.meta.time.hour).zfill(2)

            # If distance between mindist and maxdist deg:
            if (rfdata.meta.gac > opts.mindist and
                    rfdata.meta.gac < opts.maxdist):

                # Display Event Info
                nevK = nevK + 1
                if opts.reverse:
                    inum = iev + 1
                else:
                    inum = nevtT - iev + 1
                print(" ")
                print("****************************************************")
                print("* #{0:d} ({1:d}/{2:d}):  {3:13s}".format(
                    nevK, inum, nevtT,
                    rfdata.meta.time.strftime("%Y%m%d_%H%M%S")))
                print("*   Origin Time: " +
                      rfdata.meta.time.strftime("%Y-%m-%d %H:%M:%S"))
                print(
                    "*   Lat: {0:6.2f}; Lon: {1:7.2f}".format(
                        rfdata.meta.lat, rfdata.meta.lon))
                print(
                    "*   Dep: {0:6.2f}; Mag: {1:3.1f}".format(
                        rfdata.meta.dep/1000., rfdata.meta.mag))
                print("*     {0:5s} -> Ev: {1:7.2f} km; {2:7.2f} deg; " +
                      "{3:6.2f}; {4:6.2f}".format(
                          rfdata.sta.station, rfdata.meta.epi_dist,
                          rfdata.meta.gac, rfdata.meta.baz, rfdata.meta.az))

                # Get Travel times (Careful: here dep is in meters)
                tt = tpmodel.get_travel_times(
                    distance_in_degree=rfdata.meta.gac,
                    source_depth_in_km=rfdata.meta.dep/1000.)

                # Loop over all times in tt
                for t in tt:

                    # Extract P arrival
                    if t.name == 'P':

                        # Add P phase to Split
                        rfdata.add_phase(t)

                        # Break out of loop
                        break

                # Get data
                rfdata.get_data_NEZ(
                    client=client, dts=opts.dts, stdata=stalcllist,
                    ndval=opts.ndval, new_sr=opts.new_sampling_rate)

                if rfdata.err:
                    continue

                # Rotate from ZEN to LQT (Longitudinal, Radial, Transverse)
                rfdata.rotate()

                # Calculate snr over dt_snr seconds
                rfdata.calc_snrz(dt=opts.dt_snr)

                # Make sure no processing happens for NaNs
                if np.isnan(rfdata.snrz):
                    print(
                        "* SNR NaN...Skipping")
                    print(
                        "****************************************************")
                    continue

                # Create Event Folder
                evtdir = datapath + "/" + \
                    rfdata.meta.time.strftime("%Y%m%d_%H%M%S")

                # Create Folder
                if not os.path.isdir(evtdir):
                    os.makedirs(evtdir)

                # Event Data
                pickle.dump(rfdata.meta, open(evtdir + "/Meta_Data.pkl", "wb"))

                # Station Data
                pickle.dump(rfdata.sta, open(
                    evtdir + "/Station_Data.pkl", "wb"))

                # Trace Filenames
                trpref = rfdata.meta.time.strftime(
                    "%Y%m%d_%H%M%S") + "_" + sta.network + "." + sta.station

                # Raw Trace files
                pickle.dump(
                    [rfdata.data.trN, rfdata.data.trE, rfdata.data.trZ],
                    open(evtdir + "/NEZ_Data.pkl", "wb"))
                rfdata.data.trN.write(os.path.join(
                    evtdir, trpref + ".N.mseed"), format='MSEED')
                rfdata.data.trE.write(os.path.join(
                    evtdir, trpref + ".E.mseed"), format='MSEED')
                rfdata.data.trZ.write(os.path.join(
                    evtdir, trpref + ".Z.mseed"), format='MSEED')

                # LQT Traces
                pickle.dump(
                    [rfdata.data.trL, rfdata.data.trQ, rfdata.data.trT],
                    open(evtdir + "/LQT_Data.pkl", "wb"))
                rfdata.data.trL.write(os.path.join(
                    evtdir, trpref + ".L.mseed"), format='MSEED')
                rfdata.data.trQ.write(os.path.join(
                    evtdir, trpref + ".Q.mseed"), format='MSEED')
                rfdata.data.trT.write(os.path.join(
                    evtdir, trpref + ".T.mseed"), format='MSEED')

                # Update
                print("* Wrote Output Files to: ")
                print("*     "+evtdir)
                print("****************************************************")


if __name__ == "__main__":

    # Run main program
    main()
