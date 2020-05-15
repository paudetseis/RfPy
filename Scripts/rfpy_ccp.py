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
from rfpy import CCPimage
from pathlib import Path


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
    args = arguments.get_ccp_arguments()

    # Load Database
    db = stdb.io.load_db(fname=args.indb)

    # Construct station key loop
    allkeys = db.keys()

    # Extract key subset
    if len(args.stkeys) > 0:
        stkeys = []
        for skey in args.stkeys:
            stkeys.extend([s for s in allkeys if skey in s])
    else:
        stkeys = db.keys()

    if args.load:

        # Check if CCPimage object exists and whether overwrite has been set
        load_file = Path('CCP_load.pkl')
        if load_file.is_file() and not args.ovr:
            ccpfile = open(load_file, "rb")
            ccpimage = pickle.load(ccpfile)
            ccpfile.close()

        else:

            print()
            print("|-----------------------------------------------|")
            print("|  Loading data                                 |")
            print("|-----------------------------------------------|")
            print("| Gridding: ")
            print("|     start    = {0:5.1f},{1:6.1f}".format(
                args.coord_start[0],args.coord_start[1]))
            print("|     end      = {0:5.1f},{1:6.1f}".format(
                args.coord_end[0],args.coord_end[1]))
            print("|     dz    = {0} (km)".format(str(args.dz)))
            print("|     dx    = {0} (km)".format(str(args.dx)))
            print()

            # Initialize CCPimage object
            ccpimage = CCPimage(coord_start=args.coord_start,
                                coord_end=args.coord_end,
                                dz=args.dz, dx=args.dx)

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

                # Temporary print locations
                tlocs = sta.location
                if len(tlocs) == 0:
                    tlocs = ['']
                for il in range(0, len(tlocs)):
                    if len(tlocs[il]) == 0:
                        tlocs[il] = "--"
                sta.location = tlocs

                rfRstream = Stream()

                datafiles = [x for x in datapath.iterdir() if x.is_dir()]
                for folder in datafiles:

                    # Skip hidden folders
                    if folder.name.startswith('.'):
                        continue

                    # Load meta data
                    filename = folder / "Meta_Data.pkl"
                    if not filename.is_file():
                        continue
                    metafile = open(filename, 'rb')
                    meta = pickle.load(metafile)
                    metafile.close()

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

                    # If everything passed, load the RF data
                    filename = folder / "RF_Data.pkl"
                    if filename.is_file():
                        file = open(filename, "rb")
                        rfdata = pickle.load(file)
                        rfRstream.append(rfdata[1])
                        file.close()

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

    if args.prep:

        prep_file = Path("CCP_prep.pkl")
        if prep_file.is_file() and not args.ovr:
            ccpfile = open(prep_file, 'rb')
            ccpimage = pickle.load(ccpfile)
            ccpfile.close()
        else:
            load_file = Path('CCP_load.pkl')
            if not load_file.is_file():
                raise(Exception("No CCP_load.pkl file available - aborting"))
            else:
                print()
                print("|-----------------------------------------------|")
                print("|  Preparing data before stacking               |")
                print("|-----------------------------------------------|")
                print("| Frequencies: ")
                print("|     f1    = {0:4.2f} (Hz)".format(args.f1))
                print("|     f2ps  = {0:4.2f} (Hz)".format(args.f2ps))
                print("|     f2pps = {0:4.2f} (Hz)".format(args.f2pps))
                print("|     f2pss = {0:4.2f} (Hz)".format(args.f2pss))
                print("| Binning: ")
                print("|     nbaz  = {0}".format(str(args.nbaz)))
                print("|     nslow = {0}".format(str(args.nslow)))
                print()

                ccpfile = open(load_file,"rb")
                ccpimage = pickle.load(ccpfile)
                ccpfile.close()
                ccpimage.prep_data(f1=args.f1, f2ps=args.f2ps,
                                   f2pps=args.f2pps, f2pss=args.f2pss,
                                   nbaz=args.nbaz, nslow=args.nslow)
                ccpimage.is_ready_for_prestack = True
                ccpimage.save(prep_file)
                print()
                print("CCPimage saved to {0}".format(str(prep_file)))

    else:
        pass


    if args.prestack:

        prestack_file = Path("CCP_prestack.pkl")
        if prestack_file.is_file() and not args.ovr:
            ccpfile = open(prestack_file, 'rb')
            ccpimage = pickle.load(ccpfile)
            ccpfile.close()
        else:
            prep_file = Path("CCP_prep.pkl")
            if not prep_file.is_file():
                raise(Exception("No CCP_prep.pkl file available - aborting"))
            else:
                print()
                print("|-----------------------------------------------|")
                print("|  CCP pre-stacking each phase                  |")
                print("|-----------------------------------------------|")
                print()

                ccpfile = open(prep_file, 'rb')
                ccpimage = pickle.load(ccpfile)
                ccpfile.close()
                ccpimage.prestack()
                ccpimage.save(prestack_file)
                print()
                print("CCPimage saved to {0}".format(str(prestack_file)))

    else:
        pass

    if args.ccp:

        ccp_file = Path("CCP_stack.pkl")
        if ccp_file.is_file() and not args.ovr:
            ccpfile = open(ccp_file, 'rb')
            ccpimage = pickle.load(ccpfile)
            ccpfile.close()
        else:
            prestack_file = Path("CCP_prestack.pkl")
            if not prestack_file.is_file():
                raise(Exception("No CCP_prestack.pkl file available - aborting"))
            else:
                if args.linear:
                    print()
                    print("|-----------------------------------------------|")
                    print("|  Linear CCP stack - all phases                |")
                    print("|-----------------------------------------------|")
                    print()
                elif args.pws:
                    print()
                    print("|-----------------------------------------------|")
                    print("|  Phase-weighted CCP stack - all phases        |")
                    print("|-----------------------------------------------|")
                    print()

                ccpfile = open(prestack_file, 'rb')
                ccpimage = pickle.load(ccpfile)
                ccpfile.close()
                ccpimage.ccp()
                if args.linear:
                    if args.weights:
                        ccpimage.weights = args.weights
                    ccpimage.linear_stack(typ='ccp')
                elif args.pws:
                    if args.weights:
                        ccpimage.weights = args.weights
                    ccpimage.phase_weighted_stack(typ='ccp')
                ccpimage.save(ccp_file)
                print()
                print("CCPimage saved to {0}".format(str(ccp_file)))

        if args.ccp_figure:
            ccpimage.plot_ccp(save=args.save_figure, fmt=args.fmt,
                vmin=-1.*args.cbound, vmax=args.cbound, title=args.title)

    else:
        pass

    if args.gccp:

        gccp_file = Path("GCCP_stack.pkl")
        if gccp_file.is_file() and not args.ovr:
            ccpfile = open(gccp_file, 'rb')
            ccpimage = pickle.load(ccpfile)
            ccpfile.close()
        else:
            prestack_file = Path("CCP_prestack.pkl")
            if not prestack_file.is_file():
                raise(Exception("No CCP_prestack.pkl file available - aborting"))
            else:
                if args.linear:
                    print()
                    print("|-----------------------------------------------|")
                    print("|  Linear GCCP stack - all phases               |")
                    print("|-----------------------------------------------|")
                    print()
                elif args.pws:
                    print()
                    print("|-----------------------------------------------|")
                    print("|  Phase-weighted GCCP stack - all phases       |")
                    print("|-----------------------------------------------|")
                    print()

                ccpfile = open(prestack_file, 'rb')
                ccpimage = pickle.load(ccpfile)
                ccpfile.close()
                ccpimage.gccp(wlen=args.wlen)
                if args.linear:
                    if args.weights:
                        ccpimage.weights = args.weights
                    ccpimage.linear_stack(typ='gccp')
                elif args.pws:
                    if args.weights:
                        ccpimage.weights = args.weights
                    ccpimage.phase_weighted_stack(typ='gccp')
                ccpimage.save(gccp_file)
                print()
                print("CCPimage saved to {0}".format(str(gccp_file)))

        if args.ccp_figure:
            ccpimage.plot_gccp(save=args.save_figure, fmt=args.fmt,
                vmin=-1.*args.cbound, vmax=args.cbound, title=args.title)

    else:
        pass


if __name__ == "__main__":

    # Run main program
    main()
