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
from obspy.core import Stream, UTCDateTime
from rfpy import arguments, binning, plotting
from rfpy import Harmonics
from pathlib import Path


def main():

    print()
    print("################################################################################")
    print("#        __                 _                                      _           #")
    print("#  _ __ / _|_ __  _   _    | |__   __ _ _ __ _ __ ___   ___  _ __ (_) ___ ___  #")
    print("# | '__| |_| '_ \| | | |   | '_ \ / _` | '__| '_ ` _ \ / _ \| '_ \| |/ __/ __| #")
    print(
        "# | |  |  _| |_) | |_| |   | | | | (_| | |  | | | | | | (_) | | | | | (__\__ \ #")
    print("# |_|  |_| | .__/ \__, |___|_| |_|\__,_|_|  |_| |_| |_|\___/|_| |_|_|\___|___/ #")
    print("#          |_|    |___/_____|                                                  #")
    print("#                                                                              #")
    print("################################################################################")
    print()

    # Run Input Parser
    args = arguments.get_harmonics_arguments()

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
        if args.phase in ['P', 'PP', 'allP']:
            datapath = Path('P_DATA') / stkey
        elif args.phase in ['S', 'SKS', 'allS']:
            datapath = Path('S_DATA') / stkey
        if not datapath.is_dir():
            print('Path to ' + str(datapath) + ' doesn`t exist - continuing')
            continue

        # Get search start time
        if args.startT is None:
            tstart = sta.startdate
        else:
            tstart = args.startT

        # Get search end time
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

        rfRstream = Stream()
        rfTstream = Stream()

        datafiles = [x for x in datapath.iterdir() if x.is_dir()]
        for folder in datafiles:

            # Skip hidden folders
            if folder.name.startswith('.'):
                continue

            date = folder.name.split('_')[0]
            year = date[0:4]
            month = date[4:6]
            day = date[6:8]
            dateUTC = UTCDateTime(year+'-'+month+'-'+day)

            if dateUTC > tstart and dateUTC < tend:

                filename = folder / "RF_Data.pkl"
                if filename.is_file():
                    file = open(filename, "rb")
                    rfdata = pickle.load(file)
                    if rfdata[0].stats.snrh > args.snrh and rfdata[0].stats.snr and \
                            rfdata[0].stats.cc > args.cc:

                        rfRstream.append(rfdata[1])
                        rfTstream.append(rfdata[2])

                    file.close()

            else:
                continue

        if args.no_outl:
            # Remove outliers wrt variance
            varR = np.array([np.var(tr.data) for tr in rfRstream])

            # Calculate outliers
            medvarR = np.median(varR)
            madvarR = 1.4826*np.median(np.abs(varR-medvarR))
            robustR = np.abs((varR-medvarR)/madvarR)
            outliersR = np.arange(len(rfRstream))[robustR > 2.]
            for i in outliersR[::-1]:
                rfRstream.remove(rfRstream[i])
                rfTstream.remove(rfTstream[i])

            # Do the same for transverse
            varT = np.array([np.var(tr.data) for tr in rfTstream])
            medvarT = np.median(varT)
            madvarT = 1.4826*np.median(np.abs(varT-medvarT))
            robustT = np.abs((varT-medvarT)/madvarT)
            outliersT = np.arange(len(rfTstream))[robustT > 2.]
            for i in outliersT[::-1]:
                rfRstream.remove(rfRstream[i])
                rfTstream.remove(rfTstream[i])

        # Try binning if specified
        if args.nbin is not None:
            rf_tmp = binning.bin(rfRstream, rfTstream,
                                 typ='baz', nbin=args.nbin+1)
            rfRstream = rf_tmp[0]
            rfTstream = rf_tmp[1]

        # Filter original streams
        rfRstream.filter('bandpass', freqmin=args.bp[0],
                         freqmax=args.bp[1], corners=2,
                         zerophase=True)
        rfTstream.filter('bandpass', freqmin=args.bp[0],
                         freqmax=args.bp[1], corners=2,
                         zerophase=True)

        # Initialize the HkStack object
        harmonics = Harmonics(rfRstream, rfTstream)

        # Stack with or without dip
        if args.find_azim:
            harmonics.dcomp_find_azim(xmin=args.trange[0], xmax=args.trange[1])
            print("Optimal azimuth for trange between " +
                  str(args.trange[0])+" and "+str(args.trange[1]) +
                  "is: "+str(harmonics.azim))
        else:
            harmonics.dcomp_fix_azim(azim=args.azim)

        if args.plot:
            harmonics.plot(args.ymax, args.scale,
                           args.save_plot, args.title, args.form)

        if args.save:
            filename = datapath / (hkstack.hstream[0].stats.station +
                                   ".harmonics.pkl")
            harmonics.save()


if __name__ == "__main__":

    # Run main program
    main()
