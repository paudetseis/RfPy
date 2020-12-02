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
from rfpy import binning, plotting, CCPimage
from pathlib import Path
from argparse import ArgumentParser
from os.path import exists as exist
from numpy import nan


def get_ccp_arguments(argv=None):
    """
    Get Options from :class:`~optparse.OptionParser` objects.

    This function is used for data processing on-the-fly (requires web connection)

    """

    parser = ArgumentParser(
        usage="%(prog)s [arguments] <station database>",
        description="Script used to process receiver function data " +
        "for common-conversion-point (CCP) imaging.")

    # General Settings
    parser.add_argument(
        "indb",
        help="Station Database to process from.",
        type=str)
    parser.add_argument(
        "--keys",
        action="store",
        type=str,
        dest="stkeys",
        default="",
        help="Specify a comma separated list of station keys for " +
        "which to perform the analysis. These must be " +
        "contained within the station database. Partial keys will " +
        "be used to match against those in the dictionary. For " +
        "instance, providing IU will match with all stations in " +
        "the IU network [Default processes all stations in the database]")
    parser.add_argument(
        "-v", "-V", "--verbose",
        action="store_true",
        dest="verb",
        default=False,
        help="Specify to increase verbosity.")
    parser.add_argument(
        "-O", "--overwrite",
        action="store_true",
        dest="ovr",
        default=False,
        help="Force the overwriting of pre-existing data. " +
        "[Default False]")
    parser.add_argument(
        "-L", "--long-name",
        action="store_true",
        dest="lkey",
        default=False,
        help="Force folder names to use long-key form (NET.STN.CHN). " +
        "Default behaviour uses short key form (NET.STN) for the folder " +
        "names, regardless of the key type of the database."
    )

    LineGroup = parser.add_argument_group(
        title='Line Geometry Settings',
        description="Options for defining the line along which to " +
        "produce the CCP image")
    LineGroup.add_argument(
        "--start",
        action="store",
        type=str,
        dest="coord_start",
        default=None,
        help="Specify a list of two floats with the latitude and longitude " +
        "of the start point, in this respective order. [Exception raised " +
        "if not specified]")
    LineGroup.add_argument(
        "--end",
        action="store",
        dest="coord_end",
        type=str,
        default=None,
        help="Specify a list of two floats with the latitude and longitude" +
        "of the end point, in this respective order. [Exception raised " +
        "if not specified]")
    LineGroup.add_argument(
        "--dz",
        action="store",
        dest="dz",
        type=int,
        default=1.,
        help="Specify vertical cell size in km. " +
        "[Default 1.]")
    LineGroup.add_argument(
        "--dx",
        action="store",
        dest="dx",
        type=float,
        default=2.5,
        help="Specify horizontal cell size in km. " +
        "[Default 2.5]")

    PreGroup = parser.add_argument_group(
        title='Pre-processing Settings',
        description="Options for pre-processing of receiver function " +
        "data for CCP stacking")
    PreGroup.add_argument(
        "--snr",
        action="store",
        type=float,
        dest="snr",
        default=-9999.,
        help="Specify the SNR threshold for extracting receiver functions. " +
        "[Default None]")
    PreGroup.add_argument(
        "--snrh",
        action="store",
        type=float,
        dest="snrh",
        default=-9999,
        help="Specify the horizontal component SNR threshold for " +
        "extracting receiver functions. [Default None]")
    PreGroup.add_argument(
        "--cc",
        action="store",
        type=float,
        dest="cc",
        default=-1.,
        help="Specify the CC threshold for extracting receiver functions. " +
        "[Default None]")
    PreGroup.add_argument(
        "--no-outlier",
        action="store_true",
        dest="no_outl",
        default=False,
        help="Set this option to delete outliers based on the MAD " +
        "on the variance. [Default False]")
    PreGroup.add_argument(
        "--binlim",
        action="store",
        type=float,
        dest="binlim",
        default=1,
        help="Specify the minimum number of RFs in each bin. [Default 3]")
    PreGroup.add_argument(
        "--f1",
        action="store",
        type=float,
        dest="f1",
        default=0.05,
        help="Specify the low frequency corner for the bandpass filter " +
        "for all phases (Hz). [Default [0.05]]")
    PreGroup.add_argument(
        "--f2ps",
        action="store",
        type=float,
        dest="f2ps",
        default=0.75,
        help="Specify the high frequency corner for the bandpass filter " +
        "for the Ps phase (Hz). [Default [0.75]]")
    PreGroup.add_argument(
        "--f2pps",
        action="store",
        type=float,
        dest="f2pps",
        default=0.36,
        help="Specify the high frequency corner for the bandpass filter " +
        "for the Pps phase (Hz). [Default [0.36]]")
    PreGroup.add_argument(
        "--f2pss",
        action="store",
        type=float,
        dest="f2pss",
        default=0.3,
        help="Specify the high frequency corner for the bandpass filter " +
        "for the Pss phase (Hz). [Default [0.3]]")
    PreGroup.add_argument(
        "--nbaz",
        action="store",
        dest="nbaz",
        type=int,
        default=36,
        help="Specify integer number of back-azimuth bins to consider. " +
        "[Default 36]")
    PreGroup.add_argument(
        "--nslow",
        action="store",
        dest="nslow",
        type=int,
        default=40,
        help="Specify integer number of slowness bins to consider. " +
        "[Default 40]")
    PreGroup.add_argument(
        "--wlen",
        action="store",
        dest="wlen",
        type=float,
        default=35.,
        help="Specify wavelength of P-wave as sensitivity (km). " +
        "[Default 35.]")
    PreGroup.add_argument(
        "--phase",
        action="store",
        type=str,
        dest="phase",
        default='allP',
        help="Specify the phase name to plot.  " +
        "Options are 'P', 'PP', 'allP', 'S', 'SKS' or 'allS'. " +
        "[Default 'allP']")

    CCPGroup = parser.add_argument_group(
        title='CCP Settings',
        description="Options for specifying the type of CCP stacking " +
        "to perform")
    CCPGroup.add_argument(
        "--load",
        action="store_true",
        dest="load",
        default=False,
        help="Step 1. Set this option to load rfstreams into CCPimage " +
        "object. [Default False]")
    CCPGroup.add_argument(
        "--prep",
        action="store_true",
        dest="prep",
        default=False,
        help="Step 2. Set this option to prepare CCPimage before " +
        "pre-stacking. [Default False]")
    CCPGroup.add_argument(
        "--prestack",
        action="store_true",
        dest="prestack",
        default=False,
        help="Step 3. Set this option to prestack all phases before " +
        "CCP averaging. [Default False]")
    CCPGroup.add_argument(
        "--ccp",
        action="store_true",
        dest="ccp",
        default=False,
        help="Step 4a. Set this option for standard CCP stacking with " +
        "multiples. [Default False]")
    CCPGroup.add_argument(
        "--gccp",
        action="store_true",
        dest="gccp",
        default=False,
        help="Step 4b. Set this option for Gaussian-weighted " +
        "CCP stacking with multiples. [Default False]")
    CCPGroup.add_argument(
        "--linear",
        action="store_true",
        dest="linear",
        default=False,
        help="Step 5a. Set this option to produce a linear, weighted " +
        "stack for the final [G]CCP image. [Default True unless " +
        "--phase is set]")
    CCPGroup.add_argument(
        "--pws",
        action="store_true",
        dest="pws",
        default=False,
        help="Step 5b. Set this option to produce a phase weighted stack " +
        "for the final [G]CCP image. [Default False]")
    CCPGroup.add_argument(
        "--weights",
        action="store",
        dest="weights",
        default=None,
        help="Option to define weights for each of the three phases: " +
        "Ps, Pps and Pss, by specifying three comma-separated floats. " +
        "[Default 1., 3., -3.]")

    FigGroup = parser.add_argument_group(
        title='Figure Settings',
        description="Options for specifying the settings for the final figure")
    FigGroup.add_argument(
        "--figure",
        action="store_true",
        dest="ccp_figure",
        default=False,
        help="Set this option to plot the final [G]CCP figure. " +
        "[Default False]")
    FigGroup.add_argument(
        "--cbound",
        action="store",
        dest="cbound",
        type=float,
        default=None,
        help="Set the maximum value for the color palette. " +
        "[Default 0.05 for --ccp or 0.015 for --gccp]")
    FigGroup.add_argument(
        "--save-fig",
        action="store_true",
        dest="save_figure",
        default=False,
        help="Set this option to save the final [G]CCP figure. " +
        "This option can only be set if --figure is also set." +
        "[Default False]")
    FigGroup.add_argument(
        "--title",
        action="store",
        dest="title",
        type=str,
        default="",
        help="Set Figure title. [Default None]")
    FigGroup.add_argument(
        "--format",
        action="store",
        dest="fmt",
        type=str,
        default='png',
        help="Set format of figure. You can choose among " +
        "'png', 'jpg', 'eps', 'pdf'. [Default 'png']")

    args = parser.parse_args(argv)

    # Check inputs
    if not exist(args.indb):
        parser.error("Input file " + args.indb + " does not exist")

    # create station key list
    if len(args.stkeys) > 0:
        args.stkeys = args.stkeys.split(',')

    if args.phase not in ['P', 'PP', 'allP', 'S', 'SKS', 'allS']:
        parser.error(
            "Error: choose between 'P', 'PP', 'allP', 'S', 'SKS' and 'allS'.")
    if args.phase == 'allP':
        args.listphase = ['P', 'PP']
    elif args.phase == 'allS':
        args.listphase = ['S', 'SKS']
    else:
        args.listphase = [args.phase]

    if args.weights is None:
        args.weights = [1., 3., -3.]
    else:
        args.weights = [float(val) for val in args.weights.split(',')]
        if (len(args.weights)) != 3:
            parser.error(
                "Error: --weights should contain 3 " +
                "comma-separated floats")

    if args.load and args.coord_start is None:
        parser.error("--start=lon,lat is required")
    elif args.load and args.coord_start is not None:
        args.coord_start = [float(val) for val in args.coord_start.split(',')]
        if (len(args.coord_start)) != 2:
            parser.error(
                "Error: --start should contain 2 " +
                "comma-separated floats using '--load'")

    if args.load and args.coord_end is None:
        parser.error("--end=lon,lat is required")
    elif args.load and args.coord_end is not None:
        args.coord_end = [float(val) for val in args.coord_end.split(',')]
        if (len(args.coord_end)) != 2:
            parser.error(
                "Error: --end should contain 2 " +
                "comma-separated floats using '--load'")

    if not (args.load or args.prep or args.prestack or args.ccp
            or args.gccp):
        parser.error(
            "Error: needs at least one CCP Setting (--load, --prep, " +
            "--prestack, --ccp or --gccp")

    if args.linear and args.pws:
        parser.error(
            "Error: cannot use --linear and --pws at the same time")

    if args.ccp and not args.linear and not args.pws:
        args.linear = True
    if args.gccp and not args.linear and not args.pws:
        args.pws = True

    if args.ccp or args.gccp:
        if (args.save_figure or args.cbound or args.fmt) and not args.ccp_figure:
            print("Warning: Figure will not be produced since --figure " +
                  "has not been set.")

    if args.ccp_figure and not (args.ccp or args.gccp):
        parser.error(
            "Error: Cannot produce Figure without specifying the " +
            "type of stacking [--ccp or --gccp].")

    if not args.cbound and args.gccp:
        args.cbound = 0.015
    elif not args.cbound and args.ccp:
        args.cbound = 0.05

    return args


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
    args = get_ccp_arguments()

    # Load Database
    db, stkeys = stdb.io.load_db(fname=args.indb, keys=args.stkeys)

    # Track processed folders
    procfold = []

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
                args.coord_start[0], args.coord_start[1]))
            print("|     end      = {0:5.1f},{1:6.1f}".format(
                args.coord_end[0], args.coord_end[1]))
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

                # Construct Folder Name
                stfld = stkey
                if not args.lkey:
                    stfld = stkey.split('.')[0]+"."+stkey.split('.')[1]

                # Check for folder already processed
                if stfld in procfold:
                    print('{0} already processed...skipping   '.format(stfld))
                    continue

                # Define path to see if it exists
                if args.phase in ['P', 'PP', 'allP']:
                    datapath = Path('P_DATA') / stfld
                elif args.phase in ['S', 'SKS', 'allS']:
                    datapath = Path('S_DATA') / stfld
                if not datapath.is_dir():
                    print('Path to ' + str(datapath) +
                          ' doesn`t exist - continuing')
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
                if len(rfRstream) == 0:
                    continue

                ccpimage.add_rfstream(rfRstream)

                # Update processed folders
                procfold.append(stfld)

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

                ccpfile = open(load_file, "rb")
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
