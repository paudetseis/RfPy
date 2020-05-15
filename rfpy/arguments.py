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
"""

Module containing the main utility functions used in the `RfPy` scripts
that accompany this package.

"""

# -*- coding: utf-8 -*-
from argparse import ArgumentParser
from os.path import exists as exist
from obspy import UTCDateTime
from numpy import nan


def get_calc_arguments(argv=None):
    """
    Get Options from :class:`~optparse.OptionParser` objects.

    This function is used for data processing on-the-fly (requires web connection)

    """

    parser = ArgumentParser(
        usage="Usage: %prog [options] <station database>",
        description="Script used to download and pre-process " +
        "three-component (Z, N, and E), seismograms for individual " +
        "events and calculate teleseismic P-wave receiver functions" +
        "This version requests data on the fly for a given date " +
        "range. Data are requested from the internet using the " +
        "client services framework. The stations are processed one " +
        "by one and the data are stored to disk.")

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

    # Server Settings
    ServerGroup = parser.add_argument_group(
        title="Server Settings",
        description="Settings associated with which "
        "datacenter to log into.")
    ServerGroup.add_argument(
        "-S", "--Server",
        action="store",
        type=str,
        dest="Server",
        default="IRIS",
        help="Specify the server to connect to. Options include: " +
        "BGR, ETH, GEONET, GFZ, INGV, IPGP, IRIS, KOERI, " +
        "LMU, NCEDC, NEIP, NERIES, ODC, ORFEUS, RESIF, SCEDC, USGS, USP. " +
        "[Default IRIS]")
    ServerGroup.add_argument(
        "-U", "--User-Auth",
        action="store",
        type=str,
        dest="UserAuth",
        default="",
        help="Enter your IRIS Authentification Username and Password " +
        "(--User-Auth='username:authpassword') to " +
        "access and download restricted data. " +
        "[Default no user and password]")

    # Database Settings
    DataGroup = parser.add_argument_group(
        title="Local Data Settings",
        description="Settings associated with defining " +
        "and using a local data base of pre-downloaded " +
        "day-long SAC files.")
    DataGroup.add_argument(
        "--local-data",
        action="store",
        type=str,
        dest="localdata",
        default=None,
        help="Specify a comma separated list of paths containing " +
        "day-long sac files of data already downloaded. " +
        "If data exists for a seismogram is already present on disk, " +
        "it is selected preferentially over downloading " +
        "the data using the Client interface")
    DataGroup.add_argument(
        "--no-data-zero",
        action="store_true",
        dest="ndval",
        default=False,
        help="Specify to force missing data to be set as zero, rather " +
        "than default behaviour which sets to nan.")
    DataGroup.add_argument(
        "--no-local-net",
        action="store_false",
        dest="useNet",
        default=True,
        help="Specify to prevent using the Network code in the " +
        "search for local data (sometimes for CN stations " +
        "the dictionary name for a station may disagree with that " +
        "in the filename. [Default Network used]")

    # Event Selection Criteria
    EventGroup = parser.add_argument_group(
        title="Event Settings",
        description="Settings associated with refining " +
        "the events to include in matching event-station pairs")
    EventGroup.add_argument(
        "--start",
        action="store",
        type=str,
        dest="startT",
        default="",
        help="Specify a UTCDateTime compatible string representing " +
        "the start time for the event search. This will override any " +
        "station start times. [Default start date of station]")
    EventGroup.add_argument(
        "--end",
        action="store",
        type=str,
        dest="endT",
        default="",
        help="Specify a UTCDateTime compatible string representing " +
        "the end time for the event search. This will override any " +
        "station end times [Default end date of station]")
    EventGroup.add_argument(
        "--reverse", "-R",
        action="store_true",
        dest="reverse",
        default=False,
        help="Reverse order of events. Default behaviour starts at " +
        "oldest event and works towards most recent. Specify reverse " +
        "order and instead the program will start with the most recent " +
        "events and work towards older")
    EventGroup.add_argument(
        "--minmag",
        action="store",
        type=float,
        dest="minmag",
        default=6.0,
        help="Specify the minimum magnitude of event for which to search. " +
        "[Default 6.0]")
    EventGroup.add_argument(
        "--maxmag",
        action="store",
        type=float,
        dest="maxmag",
        default=9.0,
        help="Specify the maximum magnitude of event for which to search. " +
        "[Default None, i.e. no limit]")

    # Geometry Settings
    PhaseGroup = parser.add_argument_group(
        title="Geometry Settings",
        description="Settings associatd with the "
        "event-station geometries for the specified phase")
    PhaseGroup.add_argument(
        "--phase",
        action="store",
        type=str,
        dest="phase",
        default='P',
        help="Specify the phase name to use. Be careful with the distance. " +
        "setting. Options are 'P' or 'PP'. [Default 'P']")
    PhaseGroup.add_argument(
        "--mindist",
        action="store",
        type=float,
        dest="mindist",
        default=None,
        help="Specify the minimum great circle distance (degrees) between " +
        "the station and event. [Default depends on phase]")
    PhaseGroup.add_argument(
        "--maxdist",
        action="store",
        type=float,
        dest="maxdist",
        default=None,
        help="Specify the maximum great circle distance (degrees) between " +
        "the station and event. [Default depends on phase]")

    # Constants Settings
    ConstGroup = parser.add_argument_group(
        title='Parameter Settings',
        description="Miscellaneous default values and settings")
    ConstGroup.add_argument(
        "--sampling-rate",
        action="store",
        type=float,
        dest="new_sampling_rate",
        default=10.,
        help="Specify new sampling rate in Hz. [Default 10.]")
    ConstGroup.add_argument(
        "--dts",
        action="store",
        type=float,
        dest="dts",
        default=150.,
        help="Specify the window length in sec (symmetric about arrival " +
        "time). [Default 150.]")
    ConstGroup.add_argument(
        "--align",
        action="store",
        type=str,
        dest="align",
        default=None,
        help="Specify component alignment key. Can be either " +
        "ZRT, LQT, or PVH. [Default ZRT]")
    ConstGroup.add_argument(
        "--vp",
        action="store",
        type=float,
        dest="vp",
        default=6.0,
        help="Specify near-surface Vp to use with --align=PVH (km/s). " +
        "[Default 6.0]")
    ConstGroup.add_argument(
        "--vs",
        action="store",
        type=float,
        dest="vs",
        default=3.5,
        help="Specify near-surface Vs to use with --align=PVH (km/s). " +
        "[Default 3.5]")
    ConstGroup.add_argument(
        "--dt-snr",
        action="store",
        type=float,
        dest="dt_snr",
        default=30.,
        help="Specify the window length over which to calculate " +
        "the SNR in sec. [Default 30.]")
    ConstGroup.add_argument(
        "--pre-filt",
        action="store",
        type=str,
        dest="pre_filt",
        default=None,
        help="Specify two floats with low and high frequency corners for " +
        "pre-filter (before deconvolution). [Default None]")
    ConstGroup.add_argument(
        "--fmin",
        action="store",
        type=float,
        dest="fmin",
        default=0.05,
        help="Specify the minimum frequency corner for SNR and CC " +
        "filter (Hz). [Default 0.05]")
    ConstGroup.add_argument(
        "--fmax",
        action="store",
        type=float,
        dest="fmax",
        default=1.0,
        help="Specify the maximum frequency corner for SNR and CC " +
        "filter (Hz). [Default 1.0]")

    # Constants Settings
    DeconGroup = parser.add_argument_group(
        title='Deconvolution Settings',
        description="Parameters for deconvolution")
    DeconGroup.add_argument(
        "--method",
        action="store",
        dest="method",
        type=str,
        default="wiener",
        help="Specify the deconvolution method. Available methods " +
        "include 'wiener', 'water' and 'multitaper'. [Default 'wiener']")
    DeconGroup.add_argument(
        "--gfilt",
        action="store",
        dest="gfilt",
        type=float,
        default=None,
        help="Specify the Gaussian filter width in Hz. " +
        "[Default None]")
    DeconGroup.add_argument(
        "--wlevel",
        action="store",
        dest="wlevel",
        type=float,
        default=0.01,
        help="Specify the water level, used in the 'water' method. " +
        "[Default 0.01]")

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

    # Parse User Authentification
    if not len(args.UserAuth) == 0:
        tt = args.UserAuth.split(':')
        if not len(tt) == 2:
            parser.error(
                "Error: Incorrect Username and Password " +
                "Strings for User Authentification")
        else:
            args.UserAuth = tt
    else:
        args.UserAuth = []

    # Parse Local Data directories
    if args.localdata is not None:
        args.localdata = args.localdata.split(',')
    else:
        args.localdata = []

    # Check NoData Value
    if args.ndval:
        args.ndval = 0.0
    else:
        args.ndval = nan

    # Check distances for selected phase
    if args.phase not in ['P', 'PP', 'S', 'SKS']:
        parser.error(
            "Error: choose between 'P', 'PP', 'S' and 'SKS'.")
    if args.phase == 'P':
        if not args.mindist:
            args.mindist = 30.
        if not args.maxdist:
            args.maxdist = 100.
        if args.mindist < 30. or args.maxdist > 100.:
            parser.error(
                "Distances should be between 30 and 100 deg. for " +
                "teleseismic 'P' waves.")
    elif args.phase == 'PP':
        if not args.mindist:
            args.mindist = 100.
        if not args.maxdist:
            args.maxdist = 180.
        if args.mindist < 100. or args.maxdist > 180.:
            parser.error(
                "Distances should be between 100 and 180 deg. for " +
                "teleseismic 'PP' waves.")
    elif args.phase == 'S':
        if not args.mindist:
            args.mindist = 55.
        if not args.maxdist:
            args.maxdist = 85.
        if args.mindist < 55. or args.maxdist > 85.:
            parser.error(
                "Distances should be between 55 and 85 deg. for " +
                "teleseismic 'S' waves.")
    elif args.phase == 'SKS':
        if not args.mindist:
            args.mindist = 85.
        if not args.maxdist:
            args.maxdist = 115.
        if args.mindist < 85. or args.maxdist > 115.:
            parser.error(
                "Distances should be between 85 and 115 deg. for " +
                "teleseismic 'SKS' waves.")

    if args.pre_filt is not None:
        args.pre_filt = [float(val) for val in args.pre_filt.split(',')]
        args.pre_filt = sorted(args.pre_filt)
        if (len(args.pre_filt)) != 2:
            parser.error(
                "Error: --pre-filt should contain 2 " +
                "comma-separated floats")

    # Check alignment options
    if args.align is None:
        args.align = 'ZRT'
    elif args.align not in ['ZRT', 'LQT', 'PVH']:
        parser.error(
            "Error: Incorrect alignment specifier. Should be " +
            "either 'ZRT', 'LQT', or 'PVH'.")

    if args.dt_snr > args.dts:
        args.dt_snr = args.dts - 10.
        print("SNR window > data window. Defaulting to data " +
              "window minus 10 sec.")

    if args.method not in ['wiener', 'water', 'multitaper']:
        parser.error(
            "Error: 'method' should be either 'wiener', 'water' or " +
            "'multitaper'")

    return args


def get_recalc_arguments(argv=None):
    """
    Get Options from :class:`~optparse.OptionParser` objects.

    This function is used for data processing on-the-fly (requires web connection)

    """

    parser = ArgumentParser(
        usage="Usage: %prog [options] <station database>",
        description="Script used to download and pre-process " +
        "three-component (Z, N, and E), seismograms for individual " +
        "events and calculate teleseismic P-wave receiver functions" +
        "This version requests data on the fly for a given date " +
        "range. Data are requested from the internet using the " +
        "client services framework. The stations are processed one " +
        "by one and the data are stored to disk.")

    # General Settings
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

    # Constants Settings
    ConstGroup = parser.add_argument_group(
        title='Parameter Settings',
        description="Miscellaneous default values and settings")
    ConstGroup.add_argument(
        "--phase",
        action="store",
        type=str,
        dest="phase",
        default='allP',
        help="Specify the phase name to use. Be careful with the distance. " +
        "setting. Options are 'P', 'PP', 'allP', 'S', 'SKS' or 'allS'. " +
        "[Default 'allP']")
    ConstGroup.add_argument(
        "--align",
        action="store",
        type=str,
        dest="align",
        default=None,
        help="Specify component alignment key. Can be either " +
        "ZRT, LQT, or PVH. [Default ZRT]")
    ConstGroup.add_argument(
        "--vp",
        action="store",
        type=float,
        dest="vp",
        default=6.0,
        help="Specify near-surface Vp to use with --align=PVH (km/s). " +
        "[Default 6.0]")
    ConstGroup.add_argument(
        "--vs",
        action="store",
        type=float,
        dest="vs",
        default=3.5,
        help="Specify near-surface Vs to use with --align=PVH (km/s). " +
        "[Default 3.5]")
    ConstGroup.add_argument(
        "--dt-snr",
        action="store",
        type=float,
        dest="dt_snr",
        default=30.,
        help="Specify the window length over which to calculate " +
        "the SNR in sec. [Default 30.]")
    ConstGroup.add_argument(
        "--pre-filt",
        action="store",
        type=str,
        dest="pre_filt",
        default=None,
        help="Specify two floats with low and high frequency corners for " +
        "pre-filter (before deconvolution). [Default None]")
    ConstGroup.add_argument(
        "--fmin",
        action="store",
        type=float,
        dest="fmin",
        default=0.05,
        help="Specify the minimum frequency corner for SNR " +
        "filter (Hz). [Default 0.05]")
    ConstGroup.add_argument(
        "--fmax",
        action="store",
        type=float,
        dest="fmax",
        default=1.0,
        help="Specify the maximum frequency corner for SNR " +
        "filter (Hz). [Default 1.0]")

    # Constants Settings
    DeconGroup = parser.add_argument_group(
        title='Deconvolution Settings',
        description="Parameters for deconvolution")
    DeconGroup.add_argument(
        "--method",
        action="store",
        dest="method",
        type=str,
        default="wiener",
        help="Specify the deconvolution method. Available methods " +
        "include 'wiener', 'water' and 'multitaper'. [Default 'wiener']")
    DeconGroup.add_argument(
        "--gfilt",
        action="store",
        dest="gfilt",
        type=float,
        default=None,
        help="Specify the Gaussian filter width in Hz. " +
        "[Default None]")
    DeconGroup.add_argument(
        "--wlevel",
        action="store",
        dest="wlevel",
        type=float,
        default=0.01,
        help="Specify the water level, used in the 'water' method. " +
        "[Default 0.01]")

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

    if args.align is None:
        args.align = 'ZRT'
    elif args.align not in ['ZRT', 'LQT', 'PVH']:
        parser.error(
            "Error: Incorrect alignment specifier. Should be " +
            "either 'ZRT', 'LQT', or 'PVH'.")

    if args.method not in ['wiener', 'water', 'multitaper']:
        parser.error(
            "Error: 'method' should be either 'wiener', 'water' or " +
            "'multitaper'")

    if args.pre_filt is not None:
        args.pre_filt = [float(val) for val in args.pre_filt.split(',')]
        args.pre_filt = sorted(args.pre_filt)
        if (len(args.pre_filt)) != 2:
            parser.error(
                "Error: --pre-filt should contain 2 " +
                "comma-separated floats")

    return args


def get_hk_arguments(argv=None):
    """
    Get Options from :class:`~optparse.OptionParser` objects.

    This function is used for data processing on-the-fly (requires web connection)

    """

    parser = ArgumentParser(
        usage="Usage: %prog [options] <station database>",
        description="Script used to process receiver function data " +
        "for H-k stacking.")

    # General Settings
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
        default=3,
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


def get_harmonics_arguments(argv=None):
    """
    Get Options from :class:`~optparse.OptionParser` objects.

    This function is used for data processing on-the-fly (requires web connection)

    """

    parser = ArgumentParser(
        usage="Usage: %prog [options] <station database>",
        description="Script used to process receiver function data " +
        "for harmonic decomposition.")

    # General Settings
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
            print(args.trange.split(','))
            args.trange = [float(val) for val in args.trange.split(',')]
            args.trange = sorted(args.trange)
            if (len(args.trange)) != 2:
                parser.error(
                    "Error: --trange should contain 2 " +
                    "comma-separated floats")

    return args


def get_ccp_arguments(argv=None):
    """
    Get Options from :class:`~optparse.OptionParser` objects.

    This function is used for data processing on-the-fly (requires web connection)

    """

    parser = ArgumentParser(
        usage="Usage: %prog [options] <station database>",
        description="Script used to process receiver function data " +
        "for common-conversion-point (CCP) imaging.")

    # General Settings
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
        default=3,
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
        help="Option to define weights for each of the three phases: "+
        "Ps, Pps and Pss, by specifying three comma-separated floats. "+
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
                "comma-separated floats")

    if args.load and args.coord_end is None:
        parser.error("--end=lon,lat is required")
    elif args.load and args.coord_end is not None:
        args.coord_end = [float(val) for val in args.coord_end.split(',')]
        if (len(args.coord_end)) != 2:
            parser.error(
                "Error: --end should contain 2 " +
                "comma-separated floats")

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


def get_plot_arguments(argv=None):
    """
    Get Options from :class:`~optparse.OptionParser` objects.

    This function is used for data processing on-the-fly (requires web connection)

    """

    parser = ArgumentParser(
        usage="Usage: %prog [options] <station database>",
        description="Script used to plot receiver function data ")

    # General Settings
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
        help="Specify the vertical component SNR threshold for extracting receiver functions. " +
        "[Default 5.]")
    PreGroup.add_argument(
        "--snrh",
        action="store",
        type=float,
        dest="snrh",
        default=-9999.,
        help="Specify the horizontal component SNR threshold for extracting receiver functions. " +
        "[Default None]")
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
        default=3,
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
        help="Plot distribution of events on map. Other Plotting Options "+
        "will be applied to this figure (title, save, etc.). "+
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
