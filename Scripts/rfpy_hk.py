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
from rfpy import binning, plotting, HkStack
from pathlib import Path
from argparse import ArgumentParser
from os.path import exists as exist
from numpy import nan


def get_hk_arguments(argv=None):
    """
    Get Options from :class:`~optparse.OptionParser` objects.

    This function is used for data processing on-the-fly (requires web connection)

    """

    parser = ArgumentParser(
        usage="%(prog)s [arguments] <station database>",
        description="Script used to process receiver function data " +
        "for H-k stacking.")

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

    # Event Selection Criteria
    TimeGroup = parser.add_argument_group(
        title="Time Settings",
        description="Settings associated with refining " +
        "the times to include in searching for receiver function data")
    TimeGroup.add_argument(
        "--start",
        action="store",
        type=str,
        dest="startT",
        default="",
        help="Specify a UTCDateTime compatible string representing " +
        "the start time for the search. This will override any " +
        "station start times. [Default start date of station]")
    TimeGroup.add_argument(
        "--end",
        action="store",
        type=str,
        dest="endT",
        default="",
        help="Specify a UTCDateTime compatible string representing " +
        "the end time for the search. This will override any " +
        "station end times [Default end date of station]")

    PreGroup = parser.add_argument_group(
        title='Pre-processing Settings',
        description="Options for pre-processing of receiver function " +
        "data prior to H-k stacking")
    PreGroup.add_argument(
        "--binlim",
        action="store",
        type=float,
        dest="binlim",
        default=1,
        help="Specify the minimum number of RFs in each bin. [Default 3]")
    PreGroup.add_argument(
        "--bp",
        action="store",
        type=str,
        dest="bp",
        default=None,
        help="Specify the corner frequencies for the bandpass filter. " +
        "[Default 0.05,0.5]")
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
        "--slowbound",
        action="store",
        dest="slowbound",
        type=str,
        default=None,
        help="Specify a list of two floats with minimum and maximum" +
        "bounds on slowness (s/km). [Default [0.04, 0.08]]")
    PreGroup.add_argument(
        "--bazbound",
        action="store",
        dest="bazbound",
        type=str,
        default=None,
        help="Specify a list of two floats with minimum and maximum" +
        "bounds on back azimuth (degrees). [Default [0, 360]]")
    PreGroup.add_argument(
        "--pws",
        action="store_true",
        dest="pws",
        default=False,
        help="Set this option to use phase-weighted stacking during binning " +
        " [Default False]")
    PreGroup.add_argument(
        "--phase",
        action="store",
        type=str,
        dest="phase",
        default='allP',
        help="Specify the phase name to plot.  " +
        "Options are 'P', 'PP', 'allP', 'S', 'SKS' or 'allS'. " +
        "[Default 'allP']")
    PreGroup.add_argument(
        "--copy",
        action="store_true",
        dest="copy",
        default=False,
        help="Set this option to use a copy of the radial component " +
        "filtered at different corners for the Pps and Pss phases. " +
        "[Default False]")
    PreGroup.add_argument(
        "--bp-copy",
        action="store",
        dest="bp_copy",
        type=str,
        default=None,
        help="Specify a list of two floats with minimum and maximum" +
        "frequency for the copies stream (Hz). [Default [0.05, 0.35]]")

    HKGroup = parser.add_argument_group(
        title='Settings for H-k Stacking',
        description="Specify parameters of H-k search, including" +
        "bounds on search, weights, type of stacking, etc.")
    HKGroup.add_argument(
        "--hbound",
        action="store",
        type=str,
        dest="hbound",
        default=None,
        help="Specify a list of two floats with minimum and maximum" +
        "bounds on Moho depth (H, in km). [Default [20., 50.]]")
    HKGroup.add_argument(
        "--dh",
        action="store",
        type=float,
        dest="dh",
        default=0.5,
        help="Specify search interval for H (km). [Default 0.5]")
    HKGroup.add_argument(
        "--kbound",
        action="store",
        type=str,
        dest="kbound",
        default=None,
        help="Specify a list of two floats with minimum and maximum" +
        "bounds on Vp/Vs (k). [Default [1.56, 2.1]]")
    HKGroup.add_argument(
        "--dk",
        action="store",
        type=float,
        dest="dk",
        default=0.02,
        help="Specify search interval for k. [Default 0.02]")
    HKGroup.add_argument(
        "--weights",
        action="store",
        type=str,
        dest="weights",
        default=None,
        help="Specify a list of three floats with for Ps, Pps and Pass " +
        "weights in final stack. [Default [0.5, 2., -1.]]")
    HKGroup.add_argument(
        "--type",
        action="store",
        type=str,
        dest="typ",
        default="sum",
        help="Specify type of final stacking. Options are: 'sum' for " +
        "a weighted average (using weights), or 'product' for the product " +
        "of positive values in stacks. [Default 'sum']")
    HKGroup.add_argument(
        "--save",
        action="store_true",
        dest="save",
        default=False,
        help="Set this option to save the HkStack object to file. " +
        "[Default doesn't save]")

    # Constants Settings
    ModelGroup = parser.add_argument_group(
        title='Model Settings',
        description="Miscellaneous default values and settings")
    ModelGroup.add_argument(
        "--vp",
        action="store",
        type=float,
        dest="vp",
        default=6.0,
        help="Specify mean crustal Vp (km/s). [Default 6.0]")
    ModelGroup.add_argument(
        "--strike",
        action="store",
        type=float,
        dest="strike",
        default=None,
        help="Specify the strike of dipping Moho. [Default None]")
    ModelGroup.add_argument(
        "--dip",
        action="store",
        type=float,
        dest="dip",
        default=None,
        help="Specify the dip of dipping Moho. [Default None]")

    PlotGroup = parser.add_argument_group(
        title='Settings for plotting results',
        description="Specify parameters for plotting the H-k stacks.")
    PlotGroup.add_argument(
        "--plot",
        action="store_true",
        dest="plot",
        default=False,
        help="Set this option to produce a plot of the stacks [Default " +
        "does not produce plot]")
    PlotGroup.add_argument(
        "--save-plot",
        action="store_true",
        dest="save_plot",
        default=False,
        help="Set this option to save the plot [Default doesn't save]")
    PlotGroup.add_argument(
        "--title",
        action="store",
        type=str,
        dest="title",
        default="",
        help="Specify plot title [Default has no title]")
    PlotGroup.add_argument(
        "--format",
        action="store",
        type=str,
        dest="form",
        default="png",
        help="Specify format of figure. Can be any one of the valid" +
        "matplotlib formats: 'png', 'jpg', 'eps', 'pdf'. [Default 'png']")

    args = parser.parse_args(argv)

    # Check inputs
    if not exist(args.indb):
        parser.error("Input file " + args.indb + " does not exist")

    # create station key list
    if len(args.stkeys) > 0:
        args.stkeys = args.stkeys.split(',')

    # construct start time
    if len(args.startT) > 0:
        try:
            args.startT = UTCDateTime(args.startT)
        except:
            parser.error(
                "Cannot construct UTCDateTime from start time: " +
                args.startT)
    else:
        args.startT = None

    # construct end time
    if len(args.endT) > 0:
        try:
            args.endT = UTCDateTime(args.endT)
        except:
            parser.error(
                "Cannot construct UTCDateTime from end time: " +
                args.endT)
    else:
        args.endT = None

    if args.strike is None and args.dip is None:
        args.calc_dip = False
        args.nbaz = None
    elif args.strike is None or args.dip is None:
        parser.error("Specify both strike and dip for this type " +
                     "of analysis")
    else:
        args.calc_dip = True

    if args.bp is None:
        args.bp = [0.05, 0.5]
    else:
        args.bp = [float(val) for val in args.bp.split(',')]
        args.bp = sorted(args.bp)
        if (len(args.bp)) != 2:
            parser.error(
                "Error: --bp should contain 2 " +
                "comma-separated floats")

## JMG ##
    if args.slowbound is None:
        args.slowbound = [0.04, 0.08]
    else:
        args.slowbound = [float(val) for val in args.slowbound.split(',')]
        args.slowbound = sorted(args.slowbound)
        if (len(args.slowbound)) != 2:
            parser.error(
                "Error: --slowbound should contain 2 " +
                "comma-separated floats")

    if args.bazbound is None:
        args.bazbound = [0.0, 360.0]
    else:
        args.bazbound = [float(val) for val in args.bazbound.split(',')]
        args.bazbound = sorted(args.bazbound)
        if (len(args.bazbound)) != 2:
            parser.error(
                "Error: --bazbound should contain 2 " +
                "comma-separated floats")
## JMG ##

    if args.phase not in ['P', 'PP', 'allP', 'S', 'SKS', 'allS']:
        parser.error(
            "Error: choose between 'P', 'PP', 'allP', 'S', 'SKS' and 'allS'.")
    if args.phase == 'allP':
        args.listphase = ['P', 'PP']
    elif args.phase == 'allS':
        args.listphase = ['S', 'SKS']
    else:
        args.listphase = [args.phase]

    if args.typ not in ['sum', 'product']:
        parser.error(
            "Error: choose between 'sum' and 'product'")

    if args.copy:
        if args.bp_copy is None:
            args.bp_copy = [0.05, 0.35]
        else:
            args.bp_copy = [float(val)
                            for val in args.bp_copy.split(',')]
            args.bp_copy = sorted(args.bp_copy)
            if (len(args.bp_copy)) != 2:
                parser.error(
                    "Error: --bp_copy should contain 2 " +
                    "comma-separated floats")

    if args.hbound is None:
        args.hbound = [20., 50.]
    else:
        args.hbound = [float(val) for val in args.hbound.split(',')]
        args.hbound = sorted(args.hbound)
        if (len(args.hbound)) != 2:
            parser.error(
                "Error: --hbound should contain 2 " +
                "comma-separated floats")

    if args.kbound is None:
        args.kbound = [1.56, 2.1]
    else:
        args.kbound = [float(val) for val in args.kbound.split(',')]
        args.kbound = sorted(args.kbound)
        if (len(args.kbound)) != 2:
            parser.error(
                "Error: --kbound should contain 2 " +
                "comma-separated floats")

    if args.weights is None:
        args.weights = [0.5, 2.0, -1.0]
    else:
        args.weights = [float(val) for val in args.weights.split(',')]
        if (len(args.weights)) != 3:
            parser.error(
                "Error: --weights should contain 3 " +
                "comma-separated floats")

    return args


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
    args = get_hk_arguments()

    # Load Database
    db, stkeys = stdb.io.load_db(fname=args.indb, keys=args.stkeys)

    # Track processed folders
    procfold = []

    # Loop over station keys
    for stkey in list(stkeys):

        # Extract station information from dictionary
        sta = db[stkey]

        # Construct Folder Name
        stfld = stkey
        if not args.lkey:
            stfld = stkey.split('.')[0]+"."+stkey.split('.')[1]

        # Define path to see if it exists
        if args.phase in ['P', 'PP', 'allP']:
            datapath = Path('P_DATA') / stfld
        elif args.phase in ['S', 'SKS', 'allS']:
            datapath = Path('S_DATA') / stfld
        if not datapath.is_dir():
            print('Path to ' + str(datapath) + ' doesn`t exist - continuing')
            continue

        # Define save path
        if args.save:
            savepath = Path('HK_DATA') / stfld
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

        # Check for folder already processed
        if stfld in procfold:
            print('  {0} already processed...skipping   '.format(stfld))
            continue

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

                ''' # Check bounds on data
                # if meta.slow < args.slowbound[0] and meta.slow > args.slowbound[1]:
                #     continue
                # if meta.baz < args.bazbound[0] and meta.baz > args.bazbound[1]:
                #     continue
                '''

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
            filename = savepath / (hkstack.rfV1[0].stats.station +
                                   ".hkstack."+args.typ+".pkl")

            hkstack.save(file=filename)

        # Update processed folders
        procfold.append(stfld)


if __name__ == "__main__":

    # Run main program
    main()
