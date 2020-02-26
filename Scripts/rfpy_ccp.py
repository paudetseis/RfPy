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

    print()
    print("############################################")
    print("#        __                                #")
    print("#  _ __ / _|_ __  _   _     ___ ___ _ __   #")
    print("# | '__| |_| '_ \| | | |   / __/ __| '_ \  #")
    print("# | |  |  _| |_) | |_| |  | (_| (__| |_) | #")
    print("# |_|  |_| | .__/ \__, |___\___\___| .__/  #")
    print("#          |_|    |___/_____|      |_|     #")
    print("#                                          #")
    print("############################################")
    print()

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

    if opts.load:

        # Check if CCPimage object exists and whether overwrite has been set
        load_file = 'CCP_load.pkl'
        if os.path.isfile(load_file) and not opts.ovr:
            ccpimage = pickle.load(open(load_file,"rb"))

        else:

            print()
            print("|-----------------------------------------------|")
            print("|  Loading data                                 |")
            print("|-----------------------------------------------|")
            print("| Gridding: ")
            print("|     start    = {0:5.1f},{1:6.1f}".format(
                opts.coord_start[0],opts.coord_start[1]))
            print("|     end      = {0:5.1f},{1:6.1f}".format(
                opts.coord_end[0],opts.coord_end[1]))
            print("|     dz    = {0} (km)".format(str(opts.dz)))
            print("|     dx    = {0} (km)".format(str(opts.dx)))
            print()

            # Initialize CCPimage object
            ccpimage = CCPimage(coord_start=opts.coord_start,
                                coord_end=opts.coord_end,
                                dz=opts.dz, dx=opts.dx)

            # Loop over station keys
            for stkey in list(stkeys):

                # Extract station information from dictionary
                sta = db[stkey]

                # Define path to see if it exists
                datapath = 'DATA/' + stkey
                if not os.path.isdir(datapath):
                    print('Path to '+datapath+' doesn`t exist - continuing')
                    continue

                # Temporary print locations
                tlocs = sta.location
                if len(tlocs) == 0:
                    tlocs = ['']
                for il in range(0, len(tlocs)):
                    if len(tlocs[il]) == 0:
                        tlocs[il] = "--"
                sta.location = tlocs

                rfRstream = Stream()
                for folder in os.listdir(datapath):

                    filename = datapath+"/"+folder+"/RF_Data.pkl"
                    if os.path.isfile(filename):
                        file = open(filename, "rb")
                        rfdata = pickle.load(file)
                        if rfdata[0].stats.snrh > opts.snrh and rfdata[0].stats.snr and \
                                rfdata[0].stats.cc > opts.cc:
                            rfRstream.append(rfdata[1])
                        file.close()

                if opts.no_outl:
                    # Remove outliers wrt variance
                    var = np.array([np.var(tr.data) for tr in rfRstream])
                    medvar = np.median(var)
                    madvar = 1.4826*np.median(np.abs(var-medvar))
                    robust = np.abs((var-medvar)/madvar)
                    outliers = np.arange(len(rfRstream))[robust>2.]
                    for i in outliers[::-1]:
                        rfRstream.remove(rfRstream[i])      

                print("Station: {0:>2s}.{1:5s} -  {2} traces loaded".format(
                    sta.network, sta.station, len(rfRstream)))
                if len(rfRstream)==0:
                    continue

                ccpimage.add_rfstream(rfRstream)

            if len(ccpimage.radialRF) > 0:
                ccpimage.save("CCP_load.pkl")
                ccpimage.is_ready_for_prep = True
                print()
                print("CCPimage saved to 'CCP_load.pkl'")
            else:
                ccpimage.is_ready_for_prep = False
    else:
        pass

    if opts.prep:

        prep_file = "CCP_prep.pkl"
        if os.path.isfile(prep_file) and not opts.ovr:
            ccpimage = pickle.load(open(prep_file, "rb"))
        else:
            load_file = 'CCP_load.pkl'
            if not os.path.isfile(load_file):
                raise(Exception("No CCP_load.pkl file available - aborting"))
            else:
                print()
                print("|-----------------------------------------------|")
                print("|  Preparing data before stacking               |")
                print("|-----------------------------------------------|")
                print("| Frequencies: ")
                print("|     f1    = {0:4.2f} (Hz)".format(opts.f1))
                print("|     f2ps  = {0:4.2f} (Hz)".format(opts.f2ps))
                print("|     f2pps = {0:4.2f} (Hz)".format(opts.f2pps))
                print("|     f2pss = {0:4.2f} (Hz)".format(opts.f2pss))
                print("| Binning: ")
                print("|     nbaz  = {0}".format(str(opts.nbaz)))
                print("|     nslow = {0}".format(str(opts.nslow)))
                print()

                ccpimage = pickle.load(open(load_file,"rb"))
                ccpimage.prep_data(f1=opts.f1, f2ps=opts.f2ps,
                                   f2pps=opts.f2pps, f2pss=opts.f2pss,
                                   nbaz=opts.nbaz, nslow=opts.nslow)
                ccpimage.is_ready_for_prestack = True
                ccpimage.save(prep_file)
                print()
                print("CCPimage saved to {0}".format(prep_file))

    else:
        pass


    if opts.prestack:

        prestack_file = "CCP_prestack.pkl"
        if os.path.isfile(prestack_file) and not opts.ovr:
            ccpimage = pickle.load(open(prestack_file, "rb"))
        else:
            prep_file = "CCP_prep.pkl"
            if not os.path.isfile(prep_file):
                raise(Exception("No CCP_prep.pkl file available - aborting"))
            else:
                print()
                print("|-----------------------------------------------|")
                print("|  CCP pre-stacking each phase                  |")
                print("|-----------------------------------------------|")
                print()

                ccpimage = pickle.load(open(prep_file, "rb"))
                ccpimage.prestack()
                ccpimage.save(prestack_file)
                print()
                print("CCPimage saved to {0}".format(prestack_file))

    else:
        pass

    if opts.ccp:

        ccp_file = "CCP_stack.pkl"
        if os.path.isfile(ccp_file) and not opts.ovr:
            ccpimage = pickle.load(open(ccp_file, "rb"))
        else:
            prestack_file = "CCP_prestack.pkl"
            if not os.path.isfile(prestack_file):
                raise(Exception("No CCP_prestack.pkl file available - aborting"))
            else:
                if opts.linear:
                    print()
                    print("|-----------------------------------------------|")
                    print("|  Linear CCP stack - all phases                |")
                    print("|-----------------------------------------------|")
                    print()
                elif opts.phase:
                    print()
                    print("|-----------------------------------------------|")
                    print("|  Phase-weighted CCP stack - all phases        |")
                    print("|-----------------------------------------------|")
                    print()

                ccpimage = pickle.load(open(prestack_file, "rb"))
                ccpimage.ccp()
                if opts.linear:
                    ccpimage.linear_stack(typ='ccp')
                elif opts.phase:
                    ccpimage.phase_weighted_stack(typ='ccp')
                ccpimage.save(ccp_file)
                print()
                print("CCPimage saved to {0}".format(ccp_file))

        if opts.ccp_figure:
            ccpimage.plot_ccp(save=opts.save_figure, fmt=opts.fmt,
                vmin=-1.*opts.cbound, vmax=opts.cbound, title=opts.title)

    else:
        pass

    if opts.gccp:

        gccp_file = "GCCP_stack.pkl"
        if os.path.isfile(gccp_file) and not opts.ovr:
            ccpimage = pickle.load(open(gccp_file, "rb"))
        else:
            prestack_file = "CCP_prestack.pkl"
            if not os.path.isfile(prestack_file):
                raise(Exception("No CCP_prestack.pkl file available - aborting"))
            else:
                if opts.linear:
                    print()
                    print("|-----------------------------------------------|")
                    print("|  Linear GCCP stack - all phases               |")
                    print("|-----------------------------------------------|")
                    print()
                elif opts.phase:
                    print()
                    print("|-----------------------------------------------|")
                    print("|  Phase-weighted GCCP stack - all phases       |")
                    print("|-----------------------------------------------|")
                    print()

                ccpimage = pickle.load(open(prestack_file, "rb"))
                ccpimage.weights = [1., 3., -3.]
                ccpimage.gccp(wlen=opts.wlen)
                if opts.linear:
                    ccpimage.linear_stack(typ='gccp')
                elif opts.phase:
                    ccpimage.phase_weighted_stack(typ='gccp')
                ccpimage.save(gccp_file)
                print()
                print("CCPimage saved to {0}".format(gccp_file))

        if opts.ccp_figure:
            ccpimage.plot_gccp(save=opts.save_figure, fmt=opts.fmt,
                vmin=-1.*opts.cbound, vmax=opts.cbound, title=opts.title)

    else:
        pass


if __name__ == "__main__":

    # Run main program
    main()
