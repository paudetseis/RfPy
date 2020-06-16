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
from rfpy import HkStack
from pathlib import Path


def main():

    print()
    print("#########################################")
    print("#        __                 _     _     #")
    print("#  _ __ / _|_ __  _   _    | |__ | | __ #")
    print("# | '__| |_| '_ \| | | |   | '_ \| |/ / #")
    print("# | |  |  _| |_) | |_| |   | | | |   <  #")
    print("# |_|  |_| | .__/ \__, |___|_| |_|_|\_\ #")
    print("#          |_|    |___/_____|           #")
    print("#                                       #")
    print("#########################################")
    print()

    # Run Input Parser
    args = arguments.get_hk_arguments()

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

        # Define save path
        if args.save:
            savepath = Path('HK_DATA') / stkey
            if not savepath.is_dir():
                print('Path to '+str(savepath)+' doesn`t exist - creating it')
                savepath.mkdir(parents=True)

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

                # Load meta data
                metafile = folder / "Meta_Data.pkl"
                if not metafile.is_file():
                    continue
                meta = pickle.load(open(metafile, 'rb'))

                # Skip data not in list of phases
                if meta.phase not in args.listphase:
                    continue

                # QC Thresholding
                if meta.snrh < args.snrh:
                    continue
                if meta.snr < args.snr:
                    continue
                if meta.cc < args.cc:
                    continue

                # # Check bounds on data
                # if meta.slow < args.slowbound[0] and meta.slow > args.slowbound[1]:
                #     continue
                # if meta.baz < args.bazbound[0] and meta.baz > args.bazbound[1]:
                #     continue

                # If everything passed, load the RF data
                filename = folder / "RF_Data.pkl"
                if filename.is_file():
                    file = open(filename, "rb")
                    rfdata = pickle.load(file)
                    rfRstream.append(rfdata[1])
                    file.close()
                if rfdata[0].stats.npts != 1451:
                    print(folder)

        if len(rfRstream) == 0:
            continue

        if args.no_outl:
            t1 = 0.
            t2 = 30.

            varR = []
            for i in range(len(rfRstream)):
                taxis = rfRstream[i].stats.taxis
                tselect = (taxis > t1) & (taxis < t2)
                varR.append(np.var(rfRstream[i].data[tselect]))
            varR = np.array(varR)

            # Remove outliers wrt variance within time range
            medvarR = np.median(varR)
            madvarR = 1.4826*np.median(np.abs(varR-medvarR))
            robustR = np.abs((varR-medvarR)/madvarR)
            outliersR = np.arange(len(rfRstream))[robustR > 2.5]
            for i in outliersR[::-1]:
                rfRstream.remove(rfRstream[i])

        print('')
        print("Number of radial RF data: " + str(len(rfRstream)))
        print('')

        # Try binning if specified
        if args.calc_dip:
            rf_tmp = binning.bin_baz_slow(rfRstream,
                                          nbaz=args.nbaz+1,
                                          nslow=args.nslow+1, 
                                          pws=args.pws)
            rfRstream = rf_tmp[0]
        else:
            rf_tmp = binning.bin(rfRstream,
                                 typ='slow',
                                 nbin=args.nslow+1, 
                                 pws=args.pws)
            rfRstream = rf_tmp[0]

        # Get a copy of the radial component and filter
        if args.copy:
            rfRstream_copy = rfRstream.copy()
            rfRstream_copy.filter('bandpass', freqmin=args.bp_copy[0],
                                  freqmax=args.bp_copy[1], corners=2,
                                  zerophase=True)

        # Check bin counts:
        for tr in rfRstream:
            if (tr.stats.nbin < args.binlim):
                rfRstream.remove(tr)

        # Continue if stream is too short
        if len(rfRstream) < 5:
            continue

        if args.save_plot and not Path('HK_PLOTS').is_dir():
            Path('HK_PLOTS').mkdir(parents=True)

        print('')
        print("Number of radial RF bins: " + str(len(rfRstream)))
        print('')

        # Filter original stream
        rfRstream.filter('bandpass', freqmin=args.bp[0],
                         freqmax=args.bp[1], corners=2,
                         zerophase=True)

        # Initialize the HkStack object
        try:
            hkstack = HkStack(rfRstream, rfV2=rfRstream_copy,
                              strike=args.strike, dip=args.dip, vp=args.vp)
        except:
            hkstack = HkStack(rfRstream,
                              strike=args.strike, dip=args.dip, vp=args.vp)

        # Update attributes
        hkstack.hbound = args.hbound
        hkstack.kbound = args.kbound
        hkstack.dh = args.dh
        hkstack.dk = args.dk
        hkstack.weights = args.weights

        # Stack with or without dip
        if args.calc_dip:
            hkstack.stack_dip()
        else:
            hkstack.stack()

        # Average stacks
        hkstack.average(typ=args.typ)

        if args.plot:
            hkstack.plot(args.save_plot, args.title, args.form)

        if args.save:
            filename = savepath / (hkstack.rfV1[0].stats.station + \
                ".hkstack."+args.typ+".pkl")

            hkstack.save(file=filename)


if __name__ == "__main__":

    # Run main program
    main()
