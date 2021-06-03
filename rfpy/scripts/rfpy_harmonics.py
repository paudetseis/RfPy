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
from rfpy import binning, plotting, Harmonics
from pathlib import Path
from argparse import ArgumentParser
from os.path import exists as exist
from numpy import nan


def get_harmonics_arguments(argv=None):
    """
    Get Options from :class:`~optparse.OptionParser` objects.

    This function is used for data processing on-the-fly (requires web connection)

    """

    parser = ArgumentParser(
        usage="%(prog)s [arguments] <station database>",
        description="Script used to process receiver function data " +
        "for harmonic decomposition.")

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
        "data prior to harmonic decomposition")
    PreGroup.add_argument(
        "--bp",
        action="store",
        type=str,
        dest="bp",
        default=None,
        help="Specify the corner frequencies for the bandpass filter. " +
        "[Default 0.05,0.5]")
    PreGroup.add_argument(
        "--bin",
        action="store",
        dest="nbin",
        type=int,
        default=None,
        help="Specify integer number of back-azimuth bins to consider " +
        "(typically 36 or 72). [Default does not bin data]")
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
        "--phase",
        action="store",
        type=str,
        dest="phase",
        default='allP',
        help="Specify the phase name to plot.  " +
        "Options are 'P', 'PP', 'allP', 'S', 'SKS' or 'allS'. " +
        "[Default 'allP']")

    HarmonicGroup = parser.add_argument_group(
        title='Settings for harmonic decomposition',
        description="Specify parameters for the decomposition, e.g. " +
        "a fixed azimuth, depth range for finding the optimal azimuth, etc.")
    HarmonicGroup.add_argument(
        "--azim",
        action="store",
        type=float,
        dest="azim",
        default=None,
        help="Specify the azimuth angle along with to perform the " +
        "decomposition. [Default 0.]")
    HarmonicGroup.add_argument(
        "--find-azim",
        action="store_true",
        dest="find_azim",
        default=False,
        help="Set this option to calculate the optimal azimuth. [Default " +
        "uses the '--azim' value]")
    HarmonicGroup.add_argument(
        "--trange",
        action="store",
        type=str,
        dest="trange",
        default=None,
        help="Specify a list of two floats with minimum and maximum" +
        "bounds on time range for finding the optimal azimuth (sec). " +
        "[Default [0., 10.] when '--find-azim' is set]")
    HarmonicGroup.add_argument(
        "--save",
        action="store_true",
        dest="save",
        default=False,
        help="Set this option to save the Harmonics object " +
        "to a pickled file. [Default does not save object]")

    PlotGroup = parser.add_argument_group(
        title='Settings for plotting results',
        description="Specify parameters for plotting the back-azimuth " +
        "harmonics.")
    PlotGroup.add_argument(
        "--plot",
        action="store_true",
        dest="plot",
        default=False,
        help="Set this option to produce a plot of the back-azimuth harmonics")
    PlotGroup.add_argument(
        "--ymax",
        action="store",
        type=float,
        dest="ymax",
        default=30.,
        help="Specify the maximum y axis value for the plot in units of the" +
        "dependent variable (e.g., sec). [Default 30.]")
    PlotGroup.add_argument(
        "--scale",
        action="store",
        type=float,
        dest="scale",
        default=30.,
        help="Specify the scaling value that multiplies the amplitude " +
        "of the harmonic components. [Default 10.]")
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

    if args.phase not in ['P', 'PP', 'allP', 'S', 'SKS', 'allS']:
        parser.error(
            "Error: choose between 'P', 'PP', 'allP', 'S', 'SKS' and 'allS'.")
    if args.phase == 'allP':
        args.listphase = ['P', 'PP']
    elif args.phase == 'allS':
        args.listphase = ['S', 'SKS']
    else:
        args.listphase = [args.phase]

    if args.bp is None:
        args.bp = [0.05, 0.5]
    else:
        args.bp = [float(val) for val in args.bp.split(',')]
        args.bp = sorted(args.bp)
        if (len(args.bp)) != 2:
            parser.error(
                "Error: --bp should contain 2 " +
                "comma-separated floats")

    if args.azim is not None and args.find_azim:
        print("Warning: Setting both '--azim' and '--find-azim' is " +
              "conflictual. Ignoring '--find-azim'")
        args.find_azim = False
    elif args.azim is None and not args.find_azim:
        args.azim = 0.
    if args.find_azim:
        if args.trange is None:
            args.trange = [0., 10.]
        else:
            args.trange = [float(val) for val in args.trange.split(',')]
            args.trange = sorted(args.trange)
            if (len(args.trange)) != 2:
                parser.error(
                    "Error: --trange should contain 2 " +
                    "comma-separated floats")

    return args


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
    args = get_harmonics_arguments()

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
            print('Path to ' + str(datapath) +
                  ' doesn`t exist - continuing')
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
                    if rfdata[0].stats.snrh > args.snrh and \
                            rfdata[0].stats.snr and \
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

        # Initialize the Harmonics object
        harmonics = Harmonics(rfRstream, rfTstream)

        # Stack with or without dip
        if args.find_azim:
            harmonics.dcomp_find_azim(xmin=args.trange[0], xmax=args.trange[1])
            print("Optimal azimuth for trange between " +
                  str(args.trange[0])+" and "+str(args.trange[1]) +
                  " seconds is: "+str(harmonics.azim))
        else:
            harmonics.dcomp_fix_azim(azim=args.azim)

        if args.save_plot and not Path('FIGURES').is_dir():
             Path('FIGURES').mkdir(parents=True)

        if args.plot:
            harmonics.plot(args.ymax, args.scale,
                           args.save_plot, args.title, args.form)

        if args.save:
            filename = datapath / (harmonics.hstream[0].stats.station +
                                   ".harmonics.pkl")
            harmonics.save(filename)

        # Update processed folders
        procfold.append(stfld)


if __name__ == "__main__":

    # Run main program
    main()
