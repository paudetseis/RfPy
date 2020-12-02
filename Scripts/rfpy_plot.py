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
from obspy import Stream, UTCDateTime
from rfpy import binning, plotting
from pathlib import Path
from argparse import ArgumentParser
from os.path import exists as exist
from numpy import nan


def get_plot_arguments(argv=None):
    """
    Get Options from :class:`~optparse.OptionParser` objects.

    This function is used for data processing on-the-fly (requires web connection)

    """

    parser = ArgumentParser(
        usage="%(prog)s [arguments] <station database>",
        description="Script used to plot receiver function data ")

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
        help="Force the overwriting of pre-existing figures. " +
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

    PreGroup = parser.add_argument_group(
        title='Pre-processing Settings',
        description="Options for pre-processing of receiver function " +
        "data before plotting")
    PreGroup.add_argument(
        "--snr",
        action="store",
        type=float,
        dest="snr",
        default=-9999.,
        help="Specify the vertical component SNR threshold for " +
        "extracting receiver functions. [Default 5.]")
    PreGroup.add_argument(
        "--snrh",
        action="store",
        type=float,
        dest="snrh",
        default=-9999.,
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
        "--bp",
        action="store",
        type=str,
        dest="bp",
        default=None,
        help="Specify the corner frequencies for the bandpass filter. " +
        "[Default no filtering]")
    PreGroup.add_argument(
        "--pws",
        action="store_true",
        dest="pws",
        default=False,
        help="Set this option to use phase-weighted stacking during binning " +
        " [Default False]")
    PreGroup.add_argument(
        "--nbaz",
        action="store",
        dest="nbaz",
        type=int,
        default=None,
        help="Specify integer number of back-azimuth bins to consider " +
        "(typically 36 or 72). If not None, the plot will show receiver " +
        "functions sorted by back-azimuth values. [Default None]")
    PreGroup.add_argument(
        "--nslow",
        action="store",
        dest="nslow",
        type=int,
        default=None,
        help="Specify integer number of slowness bins to consider " +
        "(typically 20 or 40). If not None, the plot will show receiver " +
        "functions sorted by slowness values. [Default None]")
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
        "--phase",
        action="store",
        type=str,
        dest="phase",
        default='allP',
        help="Specify the phase name to plot.  " +
        "Options are 'P', 'PP', 'allP', 'S', 'SKS' or 'allS'. " +
        "[Default 'allP']")

    PlotGroup = parser.add_argument_group(
        title='Plot Settings',
        description="Options for plot format")
    PlotGroup.add_argument(
        "--scale",
        action="store",
        dest="scale",
        default=None,
        type=float,
        help="Specify the scaling factor for the amplitude of the " +
        "receiver functions in the wiggle plots. [Default 100. for " +
        "a back-azimuth plot, 0.02 for a slowness plot]")
    PlotGroup.add_argument(
        "--normalize",
        action="store_true",
        dest="norm",
        default=False,
        help="Set this option to produce receiver functions normalized " +
        "by the max amplitude of stacked RFs. [Default False]")
    PlotGroup.add_argument(
        "--trange",
        action="store",
        default=None,
        type=str,
        dest="trange",
        help="Specify the time range for the x-axis (sec). Negative times " +
        "are allowed [Default 0., 30.]")
    PlotGroup.add_argument(
        "--stacked",
        action="store_true",
        dest="stacked",
        default=False,
        help="Set this option to plot a stack of all traces in top panel. " +
        "[Default does not plot stacked traces]")
    PlotGroup.add_argument(
        "--save",
        action="store_true",
        dest="saveplot",
        default=False,
        help="Set this option if you wish to save the figure. [Default " +
        "does not save figure]")
    PlotGroup.add_argument(
        "--title",
        action="store",
        dest="titleplot",
        type=str,
        default='',
        help="Specify title of figure. [Default None]")
    PlotGroup.add_argument(
        "--format",
        action="store",
        type=str,
        dest="form",
        default="png",
        help="Specify format of figure. Can be any one of the valid" +
        "matplotlib formats: 'png', 'jpg', 'eps', 'pdf'. [Default 'png']")
    PlotGroup.add_argument(
        "--plot-event-dist",
        action="store_true",
        dest="plot_event_dist",
        default=False,
        help="Plot distribution of events on map. Other Plotting Options " +
        "will be applied to this figure (title, save, etc.). " +
        "[Default no plot]")

    args = parser.parse_args(argv)

    # Check inputs
    if not exist(args.indb):
        parser.error("Input file " + args.indb + " does not exist")

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

    if args.phase not in ['P', 'PP', 'allP', 'S', 'SKS', 'allS']:
        parser.error(
            "Error: choose between 'P', 'PP', 'allP', 'S', 'SKS' and 'allS'.")
    if args.phase == 'allP':
        args.listphase = ['P', 'PP']
    elif args.phase == 'allS':
        args.listphase = ['S', 'SKS']
    else:
        args.listphase = [args.phase]

    if args.bp is not None:
        args.bp = [float(val) for val in args.bp.split(',')]
        args.bp = sorted(args.bp)
        if (len(args.bp)) != 2:
            parser.error(
                "Error: --bp should contain 2 " +
                "comma-separated floats")

    if args.nbaz is None and args.nslow is None:
        args.nbaz = 72
        print("'nbaz' or 'nslow' not specified - plotting using " +
              "'nbaz=72'")
    elif args.nbaz is not None and args.nslow is not None:
        parser.error(
            "Error: Cannot specify both 'nbaz' and 'nslow'")

    if args.trange is None:
        args.trange = [0., 30.]
    else:
        args.trange = [float(val) for val in args.trange.split(',')]
        args.trange = sorted(args.trange)
        if (len(args.trange)) != 2:
            parser.error(
                "Error: --trange should contain 2 " +
                "comma-separated floats")

    # create station key list
    if len(args.stkeys) > 0:
        args.stkeys = args.stkeys.split(',')

    if args.nbaz is None and args.nslow is None:
        parser.error("Specify at least one of --nbaz or --nslow")
    elif args.nbaz is not None and args.nslow is not None:
        parser.error("Specify only one of --nbaz or --nslow")

    return args


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
    args = get_plot_arguments()

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

        # Check for folder already processed
        if stfld in procfold:
            print('  {0} already processed...skipping   '.format(stfld))
            continue

        rfRstream = Stream()
        rfTstream = Stream()

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

            # Check bounds on data
            if meta.slow < args.slowbound[0] or meta.slow > args.slowbound[1]:
                continue
            if meta.baz < args.bazbound[0] or meta.baz > args.bazbound[1]:
                continue

            # If everything passed, load the RF data
            filename = folder / "RF_Data.pkl"
            if filename.is_file():
                file = open(filename, "rb")
                rfdata = pickle.load(file)
                if args.phase in ['P', 'PP', 'allP']:
                    Rcmp = 1
                    Tcmp = 2
                elif args.phase in ['S', 'SKS', 'allS']:
                    Rcmp = 1
                    Tcmp = 2
                rfRstream.append(rfdata[Rcmp])
                rfTstream.append(rfdata[Tcmp])
                file.close()

        if len(rfRstream) == 0:
            continue

        if args.no_outl:

            varR = []
            for i in range(len(rfRstream)):
                taxis = rfRstream[i].stats.taxis
                tselect = (taxis > args.trange[0]) & (taxis < args.trange[1])
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
                tselect = (taxis > args.trange[0]) & (taxis < args.trange[1])
                varT.append(np.var(rfTstream[i].data[tselect]))
            varT = np.array(varT)

            medvarT = np.median(varT)
            madvarT = 1.4826*np.median(np.abs(varT-medvarT))
            robustT = np.abs((varT-medvarT)/madvarT)
            outliersT = np.arange(len(rfTstream))[robustT > 2.5]
            for i in outliersT[::-1]:
                rfRstream.remove(rfRstream[i])
                rfTstream.remove(rfTstream[i])

        else:
            taxis = rfRstream[0].stats.taxis

        # Filter
        if args.bp:
            rfRstream.filter('bandpass', freqmin=args.bp[0],
                             freqmax=args.bp[1], corners=2,
                             zerophase=True)
            rfTstream.filter('bandpass', freqmin=args.bp[0],
                             freqmax=args.bp[1], corners=2,
                             zerophase=True)

        if args.saveplot and not Path('RF_PLOTS').is_dir():
            Path('RF_PLOTS').mkdir(parents=True)

        print('')
        print(datapath)
        print("Number of radial RF data: " + str(len(rfRstream)))
        print("Number of transverse RF data: " + str(len(rfTstream)))
        print('')

        if args.nbaz:
            # Bin according to BAZ
            rf_tmp = binning.bin(rfRstream, rfTstream,
                                 typ='baz', nbin=args.nbaz+1,
                                 pws=args.pws)

        elif args.nslow:
            # Bin according to slowness
            rf_tmp = binning.bin(rfRstream, rfTstream,
                                 typ='slow', nbin=args.nslow+1,
                                 pws=args.pws)

        # Check bin counts:
        for tr in rf_tmp[0]:
            if (tr.stats.nbin < args.binlim):
                rf_tmp[0].remove(tr)
        for tr in rf_tmp[1]:
            if (tr.stats.nbin < args.binlim):
                rf_tmp[1].remove(tr)

        # Show a stacked trace on top OR normalize option specified
        if args.stacked or args.norm:
            st_tmp = binning.bin_all(rf_tmp[0], rf_tmp[1], pws=args.pws)
            tr1 = st_tmp[0]
            tr2 = st_tmp[1]
            if args.norm:
                # Find normalization constant
                # st_tmp = binning.bin_all(rf_tmp[0], rf_tmp[1], pws=args.pws)
                # tr1 = st_tmp[0]
                # tr2 = st_tmp[1]
                # tmp1 = tr1.data[(taxis > args.trange[0]) & (
                #     taxis < args.trange[1])]
                # tmp2 = tr2.data[(taxis > args.trange[0]) & (
                #     taxis < args.trange[1])]
                # normR = np.amax(np.abs(tmp1))
                # normT = np.amax(np.abs(tmp2))
                # norm = np.max([normR, normT])
                tmp1 = np.array([tr.data[(taxis > args.trange[0]) & (
                    taxis < args.trange[1])] for tr in rf_tmp[0]])
                tmp2 = np.array([tr.data[(taxis > args.trange[0]) & (
                    taxis < args.trange[1])] for tr in rf_tmp[1]])
                normR = np.amax(np.abs(tmp1))
                normT = np.amax(np.abs(tmp2))
                norm = np.max([normR, normT])
            else:
                norm = None
        else:
            norm = None
            tr1 = None
            tr2 = None

        # Now plot
        if args.nbaz:
            plotting.wiggle_bins(rf_tmp[0], rf_tmp[1], tr1=tr1, tr2=tr2,
                                 btyp='baz', scale=args.scale,
                                 tmin=args.trange[0], tmax=args.trange[1],
                                 norm=norm, save=args.saveplot,
                                 title=args.titleplot, form=args.form)
        elif args.nslow:
            plotting.wiggle_bins(rf_tmp[0], rf_tmp[1], tr1=tr1, tr2=tr2,
                                 btyp='slow', scale=args.scale,
                                 tmin=args.trange[0], tmax=args.trange[1],
                                 norm=norm, save=args.saveplot,
                                 title=args.titleplot, form=args.form)

        # Update processed folders
        procfold.append(stfld)

        # Event distribution
        if args.plot_event_dist:
            plotting.event_dist(
                rfRstream, phase=args.phase, save=args.saveplot,
                title=args.titleplot, form=args.form)


if __name__ == "__main__":

    # Run main program
    main()
