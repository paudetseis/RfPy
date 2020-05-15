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
import pickle
import stdb
from obspy.clients.fdsn import Client
from rfpy import arguments, utils
from rfpy import RFData
from pathlib import Path


def main():

    print()
    print("##############################################")
    print("#        __                          _       #")
    print("#  _ __ / _|_ __  _   _     ___ __ _| | ___  #")
    print("# | '__| |_| '_ \| | | |   / __/ _` | |/ __| #")
    print("# | |  |  _| |_) | |_| |  | (_| (_| | | (__  #")
    print("# |_|  |_| | .__/ \__, |___\___\__,_|_|\___| #")
    print("#          |_|    |___/_____|                #")
    print("#                                            #")
    print("##############################################")
    print()

    # Run Input Parser
    args = arguments.get_calc_arguments()

    # Load Database
    db = stdb.io.load_db(fname=args.indb)

    # Construct station key loop
    allkeys = db.keys()
    sorted(allkeys)

    # Extract key subset
    if len(args.stkeys) > 0:
        stkeys = []
        for skey in args.stkeys:
            stkeys.extend([s for s in allkeys if skey in s])
    else:
        stkeys = db.keys()
        sorted(stkeys)

    # Loop over station keys
    for stkey in list(stkeys):

        # Extract station information from dictionary
        sta = db[stkey]

        # Define path to see if it exists
        if args.phase in ['P', 'PP']:
            datapath = Path('P_DATA') / stkey
        elif args.phase in ['S', 'SKS']:
            datapath = Path('S_DATA') / stkey
        if not datapath.is_dir():
            print('Path to '+str(datapath)+' doesn`t exist - creating it')
            datapath.mkdir()

        # Establish client for data
        if len(args.UserAuth) == 0:
            data_client = Client(args.Server)
        else:
            data_client = Client(args.Server, user=args.UserAuth[0],
                                 password=args.UserAuth[1])

        # Establish client for events
        event_client = Client()

        # Get catalogue search start time
        if args.startT is None:
            tstart = sta.startdate
        else:
            tstart = args.startT

        # Get catalogue search end time
        if args.endT is None:
            tend = sta.enddate
        else:
            tend = args.endT
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
        if args.maxmag is None:
            print("|   Mag:   >{0:3.1f}", format(args.minmag) +
                  "                                 |")
        else:
            msg = "|   Mag:   {0:3.1f}".format(args.minmag) + \
                " - {0:3.1f}".format(args.maxmag) + \
                "                            |"
            print(msg)

        print("| ...                                           |")

        # Get catalogue using deployment start and end
        cat = event_client.get_events(
            starttime=tstart, endtime=tend,
            minmagnitude=args.minmag, maxmagnitude=args.maxmag)

        # Total number of events in Catalogue
        nevK = 0
        nevtT = len(cat)
        print(
            "|  Found {0:5d}".format(nevtT) +
            " possible events                  |")
        ievs = range(0, nevtT)

        # Get Local Data Availabilty
        if len(args.localdata) > 0:
            print("|-----------------------------------------------|")
            print("| Cataloging Local Data...                      |")
            if args.useNet:
                stalcllist = utils.list_local_data_stn(
                    lcldrs=args.localdata, sta=sta.station,
                    net=sta.network, altnet=sta.altnet)
                print("|   {0:>2s}.{1:5s}: {2:6d}".format(
                    sta.network, sta.station, len(stalcllist)) +
                    " files                      |")
                print(stalcllist[0:10])
            else:
                stalcllist = utils.list_local_data_stn(
                    lcldrs=args.localdata, sta=sta.station)
                print("|   {0:5s}: {1:6d} files                " +
                      "        |".format(
                          sta.station, len(stalcllist)))
        else:
            stalcllist = []
        print("|===============================================|")

        # Select order of processing
        if args.reverse:
            ievs = range(0, nevtT)
        else:
            ievs = range(nevtT-1, -1, -1)

        # Read through catalogue
        for iev in ievs:

            # Extract event
            ev = cat[iev]

            # Initialize RF object with station info
            rfdata = RFData(sta)

            # Add event to rfdata object
            accept = rfdata.add_event(
                ev, gacmin=args.mindist, gacmax=args.maxdist,
                phase=args.phase, returned=True)

            # Define time stamp
            yr = str(rfdata.meta.time.year).zfill(4)
            jd = str(rfdata.meta.time.julday).zfill(3)
            hr = str(rfdata.meta.time.hour).zfill(2)

            # If event is accepted (data exists)
            if accept:

                # Display Event Info
                nevK = nevK + 1
                if args.reverse:
                    inum = iev + 1
                else:
                    inum = nevtT - iev + 1
                print(" ")
                print("**************************************************")
                print("* #{0:d} ({1:d}/{2:d}):  {3:13s} {4}".format(
                    nevK, inum, nevtT, rfdata.meta.time.strftime(
                        "%Y%m%d_%H%M%S"), stkey))
                if args.verb:
                    print("*   Phase: {}".format(args.phase))
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
                evtdir = datapath / timekey
                RFfile = evtdir / 'RF_Data.pkl'
                ZNEfile = evtdir / 'ZNE_Data.pkl'
                metafile = evtdir / 'Meta_Data.pkl'
                stafile = evtdir / 'Station_Data.pkl'

                # Check if RF data already exist and overwrite has been set
                if evtdir.exists():
                    if RFfile.exists():
                        if not args.ovr:
                            continue

                # Get data
                has_data = rfdata.download_data(
                    client=data_client, dts=args.dts, stdata=stalcllist,
                    ndval=args.ndval, new_sr=args.new_sampling_rate,
                    returned=True, verbose=args.verb)

                if not has_data:
                    continue

                # Create Folder if it doesn't exist
                if evtdir.exists():
                    evtdir.mkdir()

                # Save ZNE Traces
                pickle.dump(rfdata.data, open(ZNEfile, "wb"))

                # Rotate from ZNE to 'align' ('ZRT', 'LQT', or 'PVH')
                rfdata.rotate(vp=args.vp, vs=args.vs, align=args.align)

                # Calculate snr over dt_snr seconds
                rfdata.calc_snr(
                    dt=args.dt_snr, fmin=args.fmin, fmax=args.fmax)
                if args.verb:
                    print("* SNR: {}".format(rfdata.meta.snr))

                # Make sure no processing happens for NaNs
                if np.isnan(rfdata.meta.snr):
                    if args.verb:
                        print("* SNR NaN...Skipping")
                    print("**************************************************")
                    continue

                # Deconvolve data
                rfdata.deconvolve(
                    vp=args.vp, vs=args.vs,
                    align=args.align, method=args.method,
                    gfilt=args.gfilt, wlevel=args.wlevel,
                    pre_filt=args.pre_filt)

                # Get cross-correlation QC
                rfdata.calc_cc()
                if args.verb:
                    print("* CC: {}".format(rfdata.meta.cc))

                # Convert to Stream
                rfstream = rfdata.to_stream()

                # Save event meta data
                pickle.dump(rfdata.meta, open(metafile, "wb"))

                # Save Station Data
                pickle.dump(rfdata.sta, open(stafile, "wb"))

                # Save RF Traces
                pickle.dump(rfstream, open(RFfile, "wb"))

                # Update
                if args.verb:
                    print("* Wrote Output Files to: ")
                    print("*     "+str(evtdir))
                print("**************************************************")


if __name__ == "__main__":

    # Run main program
    main()
