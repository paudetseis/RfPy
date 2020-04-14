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
        if opts.phase in ['P', 'PP', 'allP']:
            datapath = 'P_DATA/' + stkey
        elif opts.phase in ['S', 'SKS', 'allS']:
            datapath = 'S_DATA/' + stkey
        if not os.path.isdir(datapath):
            print('Path to ' + datapath + ' doesn`t exist - continuing')
            continue

        # Define save path
        if opts.save:
            savepath = 'HK_DATA/' + stkey
            if not os.path.isdir(savepath):
                print('Path to '+savepath+' doesn`t exist - creating it')
                os.makedirs(savepath)

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

            if folder.startswith('.'):
                continue
                
            date = folder.split('_')[0]
            year = date[0:4]
            month = date[4:6]
            day = date[6:8]
            dateUTC = UTCDateTime(year+'-'+month+'-'+day)

            if dateUTC > tstart and dateUTC < tend:

                # Load meta data
                if not os.path.isfile(datapath+"/"+folder+"/Meta_Data.pkl"):
                    continue
                meta = pickle.load(open(
                    datapath + "/" + folder + "/Meta_Data.pkl",'rb'))

                # Skip data not in list of phases
                if meta.phase not in opts.listphase:
                    continue

                # QC Thresholding
                if meta.snrh < opts.snrh:
                    continue
                if meta.snr < opts.snr:
                    continue
                if meta.cc < opts.cc:
                    continue

                # # Check bounds on data
                # if meta.slow < opts.slowbound[0] and meta.slow > opts.slowbound[1]:
                #     continue
                # if meta.baz < opts.bazbound[0] and meta.baz > opts.bazbound[1]:
                #     continue

                # If everything passed, load the RF data
                filename = datapath + "/" + folder + "/RF_Data.pkl"
                if os.path.isfile(filename):
                    file = open(filename, "rb")
                    rfdata = pickle.load(file)
                    rfRstream.append(rfdata[1])
                    file.close()
                if rfdata[0].stats.npts != 1451:
                    print(folder)

        if len(rfRstream) == 0:
            continue

        if opts.no_outl:
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
        if opts.calc_dip:
            rf_tmp = binning.bin_baz_slow(
                rfRstream, nbaz=opts.nbaz+1, nslow=opts.nslow+1)
            rfRstream = rf_tmp[0]
        else:
            rf_tmp = binning.bin(
                rfRstream, typ='slow', nbin=opts.nslow+1)
            rfRstream = rf_tmp[0]

        # Get a copy of the radial component and filter
        if opts.copy:
            rfRstream_copy = rfRstream.copy()
            rfRstream_copy.filter('bandpass', freqmin=opts.bp_copy[0],
                                  freqmax=opts.bp_copy[1], corners=2,
                                  zerophase=True)

        # Check bin counts:
        for tr in rfRstream:
            if (tr.stats.nbin < opts.binlim):
                rfRstream.remove(tr)

        if opts.save_plot and not os.path.isdir('HK_PLOTS'):
            os.makedirs('HK_PLOTS')

        print('')
        print("Number of radial RF bins: " + str(len(rfRstream)))
        print('')

        # Filter original stream
        rfRstream.filter('bandpass', freqmin=opts.bp[0],
                         freqmax=opts.bp[1], corners=2,
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

        if opts.plot:
            hkstack.plot(opts.save_plot, opts.title, opts.form)

        if opts.save:
            filename = savepath + "/" + hkstack.rfV1[0].stats.station + \
                ".hkstack."+opts.typ+".pkl"

            hkstack.save(file=filename)


if __name__ == "__main__":

    # Run main program
    main()
