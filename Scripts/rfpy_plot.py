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
from obspy import Stream
from rfpy import options, binning, plotting


def main():

    print()
    print("#################################################")
    print("#        __                        _       _    #")
    print("#  _ __ / _|_ __  _   _      _ __ | | ___ | |_  #")
    print("# | '__| |_| '_ \| | | |    | '_ \| |/ _ \| __| #")
    print("# | |  |  _| |_) | |_| |    | |_) | | (_) | |_  #")
    print("# |_|  |_| | .__/ \__, |____| .__/|_|\___/ \__| #")
    print("#          |_|    |___/_____|_|                 #")
    print("#                                               #")
    print("#################################################")
    print()

    # Run Input Parser
    (opts, indb) = options.get_plot_options()

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
        print("|===============================================|")
        print("|                   {0:>8s}                    |".format(
            sta.station))
        print("|===============================================|")
        print("|  Station: {0:>2s}.{1:5s}                            |".format(
            sta.network, sta.station))
        print("|      Channel: {0:2s}; Locations: {1:15s}  |".format(
            sta.channel, ",".join(tlocs)))
        print("|      Lon: {0:7.2f}; Lat: {1:6.2f}                |".format(
            sta.longitude, sta.latitude))
        print("|-----------------------------------------------|")

        rfRstream = Stream()
        rfTstream = Stream()

        for folder in os.listdir(datapath):

            # Skip hidden folders
            if folder.startswith('.'):
                continue

            # Load meta data
            filename = datapath+"/"+folder+"/Meta_Data.pkl"
            if not os.path.isfile(filename):
                continue
            meta = pickle.load(open(filename, 'rb'))

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

            # Check bounds on data
            if meta.slow < opts.slowbound[0] or meta.slow > opts.slowbound[1]:
                continue
            if meta.baz < opts.bazbound[0] or meta.baz > opts.bazbound[1]:
                continue

            # If everything passed, load the RF data
            filename = datapath + "/" + folder + "/RF_Data.pkl"
            if os.path.isfile(filename):
                file = open(filename, "rb")
                rfdata = pickle.load(file)
                if opts.phase in ['P', 'PP', 'allP']:
                    Rcmp = 1
                    Tcmp = 2
                elif opts.phase in ['S', 'SKS', 'allS']:
                    Rcmp = 1
                    Tcmp = 2
                rfRstream.append(rfdata[Rcmp])
                rfTstream.append(rfdata[Tcmp])
                file.close()
            # print(folder)
            # print(rfdata[1].stats.npts)

        if len(rfRstream) == 0:
            continue

        if opts.no_outl:

            varR = []
            for i in range(len(rfRstream)):
                taxis = rfRstream[i].stats.taxis
                tselect = (taxis > opts.trange[0]) & (taxis < opts.trange[1])
                varR.append(np.var(rfRstream[i].data[tselect]))
            varR = np.array(varR)

            # Remove outliers wrt variance within time range
            medvarR = np.median(varR)
            madvarR = 1.4826*np.median(np.abs(varR-medvarR))
            robustR = np.abs((varR-medvarR)/madvarR)
            outliersR = np.arange(len(rfRstream))[robustR > 2.5]
            for i in outliersR[::-1]:
                rfRstream.remove(rfRstream[i])
                rfTstream.remove(rfTstream[i])

            # Do the same for transverse
            varT = []
            for i in range(len(rfRstream)):
                taxis = rfRstream[i].stats.taxis
                tselect = (taxis > opts.trange[0]) & (taxis < opts.trange[1])
                varT.append(np.var(rfTstream[i].data[tselect]))
            varT = np.array(varT)

            medvarT = np.median(varT)
            madvarT = 1.4826*np.median(np.abs(varT-medvarT))
            robustT = np.abs((varT-medvarT)/madvarT)
            outliersT = np.arange(len(rfTstream))[robustT > 2.5]
            for i in outliersT[::-1]:
                rfRstream.remove(rfRstream[i])
                rfTstream.remove(rfTstream[i])

        # Filter
        if opts.bp:
            rfRstream.filter('bandpass', freqmin=opts.bp[0],
                             freqmax=opts.bp[1], corners=2,
                             zerophase=True)
            rfTstream.filter('bandpass', freqmin=opts.bp[0],
                             freqmax=opts.bp[1], corners=2,
                             zerophase=True)

        if opts.saveplot and not os.path.isdir('RF_PLOTS'):
            os.makedirs('RF_PLOTS')

        print('')
        print("Number of radial RF data: " + str(len(rfRstream)))
        print("Number of transverse RF data: " + str(len(rfTstream)))
        print('')

        if opts.nbaz:
            # Bin according to BAZ
            rf_tmp = binning.bin(rfRstream, rfTstream,
                                 typ='baz', nbin=opts.nbaz+1,
                                 pws=opts.pws)

        elif opts.nslow:
            # Bin according to slowness
            rf_tmp = binning.bin(rfRstream, rfTstream,
                                 typ='slow', nbin=opts.nslow+1,
                                 pws=opts.pws)

        # Check bin counts:
        for tr in rf_tmp[0]:
            if (tr.stats.nbin < opts.binlim):
                rf_tmp[0].remove(tr)
        for tr in rf_tmp[1]:
            if (tr.stats.nbin < opts.binlim):
                rf_tmp[1].remove(tr)

        # Show a stacked trace on top OR normalize option specified
        if opts.stacked or opts.norm:
            st_tmp = binning.bin_all(rf_tmp[0], rf_tmp[1], pws=opts.pws)
            tr1 = st_tmp[0]
            tr2 = st_tmp[1]
            # Find normalization constant
            normR = np.amax(np.abs(
                tr1.data[(taxis > opts.trange[0]) & (taxis < opts.trange[1])]))
            normT = np.amax(np.abs(
                tr2.data[(taxis > opts.trange[0]) & (taxis < opts.trange[1])]))
            norm = np.max([normR, normT])
        else:
            norm = None
            tr1 = None
            tr2 = None

        # Now plot
        if opts.nbaz:
            plotting.wiggle_bins(rf_tmp[0], rf_tmp[1], tr1=tr1, tr2=tr2,
                                 btyp='baz', scale=opts.scale,
                                 tmin=opts.trange[0], tmax=opts.trange[1],
                                 norm=norm, save=opts.saveplot,
                                 title=opts.titleplot, form=opts.form)
        elif opts.nslow:
            plotting.wiggle_bins(rf_tmp[0], rf_tmp[1], tr1=tr1, tr2=tr2,
                                 btyp='slow', scale=opts.scale,
                                 tmin=opts.trange[0], tmax=opts.trange[1],
                                 norm=norm, save=opts.saveplot,
                                 title=opts.titleplot, form=opts.form)

        # Event distribution
        if opts.plot_event_dist:
            plotting.event_dist(rfRstream, phase=opts.phase, save=opts.saveplot,
                                title=opts.titleplot, form=opts.form)

if __name__ == "__main__":

    # Run main program
    main()
