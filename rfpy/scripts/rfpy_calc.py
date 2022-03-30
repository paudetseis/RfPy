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
from obspy import Catalog, UTCDateTime
from http.client import IncompleteRead
from rfpy import utils, RFData
from pathlib import Path
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
        usage="%(prog)s [arguments] <station database>",
        description="Script used to download and pre-process " +
        "three-component ('Z', 'N', and 'E'), seismograms for individual " +
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
    parser.add_argument(
        "-L", "--long-name",
        action="store_true",
        dest="lkey",
        default=False,
        help="Force folder names to use long-key form (NET.STN.CHN). " +
        "Default behaviour uses short key form (NET.STN) for the folder " +
        "names, regardless of the key type of the database."
    )

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
        "day-long SAC or MSEED files.")
    DataGroup.add_argument(
        "--local-data",
        action="store",
        type=str,
        dest="localdata",
        default=None,
        help="Specify a comma separated list of paths containing " +
        "day-long sac or mseed files of data already downloaded. " +
        "If data exists for a seismogram is already present on disk, " +
        "it is selected preferentially over downloading " +
        "the data using the Client interface")
    DataGroup.add_argument(
    	"--dtype",
    	action="store",
    	type=str,
    	dest="dtype",
    	default='SAC',
    	help="Specify the data archive file type, either SAC " +
    	" or MSEED. Note the default behaviour is to search for " +
    	"SAC files. Local archive files must have extensions of '.SAC' "+
    	" or '.MSEED. These are case dependent, so specify the correct case"+
        "here.")
    DataGroup.add_argument(
        "--no-data-zero",
        action="store_true",
        dest="ndval",
        default=False,
        help="Specify to force missing data to be set as zero, rather " +
        "than default behaviour which sets to nan. Note this is applied " +
        "only to the SAC data")
    DataGroup.add_argument(
        "--no-local-net",
        action="store_false",
        dest="useNet",
        default=True,
        help="Specify to prevent using the Network code in the " +
        "search for local data (sometimes for CN stations " +
        "the dictionary name for a station may disagree with that " +
        "in the filename. [Default Network used]")
    DataGroup.add_argument(
        "--save-Z12",
        action="store_true",
        dest="saveZ12",
        default=False,
        help="Specify to save Z12 (un-rotated) components. [Default " +
        "False]")

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

    # Check Datatype specification
    if (not args.dtype.upper() == 'MSEED') and  (not args.dtype.upper() == 'SAC'):
    	parser.error(
    		"Error: Local Data Archive must be of types 'SAC'" +
    		"or MSEED. These must match the file extensions for " +
    		" the archived data.")

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

    # Check alignment arguments
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


def main():

    print()
    print("##############################################")
    print("#        __                          _       #")
    print("#  _ __ / _|_ __  _   _     ___ __ _| | ___  #")
    print("# | '__| |_| '_ \| | | |   / __/ _` | |/ __| #")
    print("# | |  |  _| |_) | |_| |  | (_| (_| | | (__  #")
    print("# |_|  |_| | .__/ \__, |___\___\__,_|_|\___| #")
    print("#          |_|    |___/_____|                #")
    print("#                                            #")
    print("##############################################")
    print()

    # Run Input Parser
    args = get_calc_arguments()

    # Load Database
    db, stkeys = stdb.io.load_db(fname=args.indb, keys=args.stkeys)

    # Loop over station keys
    for stkey in list(stkeys):

        # Extract station information from dictionary
        sta = db[stkey]

        # Construct Folder Name
        stfld = stkey
        if not args.lkey:
            stfld = stkey.split('.')[0]+"."+stkey.split('.')[1]

        # Define path to see if it exists
        if args.phase in ['P', 'PP']:
            datapath = Path('P_DATA') / stfld
        elif args.phase in ['S', 'SKS']:
            datapath = Path('S_DATA') / stfld
        if not datapath.exists():
            print('Path to '+str(datapath)+' doesn`t exist - creating it')
            datapath.mkdir(parents=True)

        # Establish client for data
        if len(args.UserAuth) == 0:
            data_client = Client(args.Server)
        else:
            data_client = Client(args.Server, user=args.UserAuth[0],
                                 password=args.UserAuth[1])

        # Establish client for events
        event_client = Client()

        # Get catalogue search start time
        if args.startT is None:
            tstart = sta.startdate
        else:
            tstart = args.startT

        # Get catalogue search end time
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
        print("| Searching Possible events:                    |")
        print("|   Start: {0:19s}                  |".format(
            tstart.strftime("%Y-%m-%d %H:%M:%S")))
        print("|   End:   {0:19s}                  |".format(
            tend.strftime("%Y-%m-%d %H:%M:%S")))
        if args.maxmag is None:
            print("|   Mag:   >{0:3.1f}", format(args.minmag) +
                  "                                 |")
        else:
            msg = "|   Mag:   {0:3.1f}".format(args.minmag) + \
                " - {0:3.1f}".format(args.maxmag) + \
                "                            |"
            print(msg)

        print("| ...                                           |")

        # Get catalogue using deployment start and end
        try:
            cat = event_client.get_events(
                starttime=tstart, endtime=tend,
                minmagnitude=args.minmag, maxmagnitude=args.maxmag)
        except IncompleteRead:
            # See http.client.IncompleteRead
            # https://github.com/obspy/obspy/issues/3028#issue-1179808237
            
            # Read yearly chunks
            print('| Warning: Unable to get entire catalog at once |')
            print('| Trying to get one year at a time              |')
            print('|                                               |')
            
            chunk = 365 * 86400  # Query length in seconds
            cat = Catalog()
            tend = min(tend, UTCDateTime.now())
            while tstart < tend:
                print("| Start:   {0:19s}                  |".format(
                    tstart.strftime("%Y-%m-%d %H:%M:%S")))
                cat += event_client.get_events(
                    starttime=tstart, endtime=tstart + chunk,
                    minmagnitude=args.minmag, maxmagnitude=args.maxmag)

                # Make sure that we go all the way to tend
                if tstart + chunk > tend:
                    chunk = tend - tstart

                    # But do not get caught up in an infinite loop due to
                    # rounding errors
                    if chunk <= 1:
                        break

                tstart += chunk

        # Total number of events in Catalogue
        nevK = 0
        nevtT = len(cat)
        print(
            "|  Found {0:5d}".format(nevtT) +
            " possible events                  |")
        ievs = range(0, nevtT)

        # Get Local Data Availabilty
        if len(args.localdata) > 0:
            print("|-----------------------------------------------|")
            print("| Cataloging Local Data...                      |")
            if args.useNet:
                stalcllist = utils.list_local_data_stn(
                    lcldrs=args.localdata, sta=sta.station,
                    net=sta.network, dtype=args.dtype, altnet=sta.altnet)
                print("|   {0:>2s}.{1:5s}: {2:6d}".format(
                    sta.network, sta.station, len(stalcllist)) +
                    " files                      |")
                #print(stalcllist[0:10])
            else:
                stalcllist = utils.list_local_data_stn(
                    lcldrs=args.localdata, sta=sta.station, dtype=args.dtype)
                print("|   {0:5s}: {1:6d} files                " +
                      "        |".format(
                          sta.station, len(stalcllist)))
        else:
            stalcllist = []
        print("|===============================================|")

        # Select order of processing
        if args.reverse:
            ievs = range(0, nevtT)
        else:
            ievs = range(nevtT-1, -1, -1)

        # Read through catalogue
        for iev in ievs:

            # Extract event
            ev = cat[iev]

            # Initialize RF object with station info
            rfdata = RFData(sta)

            # Add event to rfdata object
            accept = rfdata.add_event(
                ev, gacmin=args.mindist, gacmax=args.maxdist,
                phase=args.phase, returned=True)

            # Define time stamp
            yr = str(rfdata.meta.time.year).zfill(4)
            jd = str(rfdata.meta.time.julday).zfill(3)
            hr = str(rfdata.meta.time.hour).zfill(2)

            # If event is accepted (data exists)
            if accept:

                # Display Event Info
                nevK = nevK + 1
                if args.reverse:
                    inum = iev + 1
                else:
                    inum = nevtT - iev + 1
                print(" ")
                print("*"*50)
                print("* #{0:d} ({1:d}/{2:d}):  {3:13s} {4}".format(
                    nevK, inum, nevtT, rfdata.meta.time.strftime(
                        "%Y%m%d_%H%M%S"), stkey))
                if args.verb:
                    print("*   Phase: {}".format(args.phase))
                    print("*   Origin Time: " +
                          rfdata.meta.time.strftime("%Y-%m-%d %H:%M:%S"))
                    print(
                        "*   Lat: {0:6.2f};        Lon: {1:7.2f}".format(
                            rfdata.meta.lat, rfdata.meta.lon))
                    print(
                        "*   Dep: {0:6.2f} km;     Mag: {1:3.1f}".format(
                            rfdata.meta.dep, rfdata.meta.mag))
                    print(
                        "*   Dist: {0:7.2f} km;".format(rfdata.meta.epi_dist) +
                        "   Epi dist: {0:6.2f} deg\n".format(rfdata.meta.gac) +
                        "*   Baz:  {0:6.2f} deg;".format(rfdata.meta.baz) +
                        "   Az: {0:6.2f} deg".format(rfdata.meta.az))

                # Event Folder
                timekey = rfdata.meta.time.strftime("%Y%m%d_%H%M%S")
                evtdir = datapath / timekey
                RFfile = evtdir / 'RF_Data.pkl'
                ZNEfile = evtdir / 'ZNE_Data.pkl'
                metafile = evtdir / 'Meta_Data.pkl'
                stafile = evtdir / 'Station_Data.pkl'

                # Check if RF data already exist and overwrite has been set
                if evtdir.exists():
                    if RFfile.exists():
                        if not args.ovr:
                            continue

                # Get data
                has_data = rfdata.download_data(
                    client=data_client, dts=args.dts, stdata=stalcllist,
                    ndval=args.ndval, dtype=args.dtype, new_sr=args.new_sampling_rate,
                    returned=True, verbose=args.verb)

                if not has_data:
                    continue

                # Create Folder if it doesn't exist
                if not evtdir.exists():
                    evtdir.mkdir(parents=True)

                # Save ZNE Traces
                pickle.dump(rfdata.data, open(ZNEfile, "wb"))

                # Save Z12 if components exist
                if hasattr(rfdata, "dataZ12"):
                    Z12file = evtdir / 'Z12_Data.pkl'
                    pickle.dump(rfdata.dataZ12, open(Z12file, "wb"))

                # Rotate from ZNE to 'align' ('ZRT', 'LQT', or 'PVH')
                rfdata.rotate(vp=args.vp, vs=args.vs, align=args.align)

                # Calculate snr over dt_snr seconds
                rfdata.calc_snr(
                    dt=args.dt_snr, fmin=args.fmin, fmax=args.fmax)
                if args.verb:
                    print("* SNR: {}".format(rfdata.meta.snr))

                # Make sure no processing happens for NaNs
                if np.isnan(rfdata.meta.snr):
                    if args.verb:
                        print("* SNR NaN...Skipping")
                    print("*"*50)
                    continue

                # Deconvolve data
                rfdata.deconvolve(
                    vp=args.vp, vs=args.vs,
                    align=args.align, method=args.method,
                    gfilt=args.gfilt, wlevel=args.wlevel,
                    pre_filt=args.pre_filt)

                # Get cross-correlation QC
                rfdata.calc_cc()
                if args.verb:
                    print("* CC: {}".format(rfdata.meta.cc))

                # Convert to Stream
                rfstream = rfdata.to_stream()

                # Save event meta data
                pickle.dump(rfdata.meta, open(metafile, "wb"))

                # Save Station Data
                pickle.dump(rfdata.sta, open(stafile, "wb"))

                # Save RF Traces
                pickle.dump(rfstream, open(RFfile, "wb"))

                # Update
                if args.verb:
                    print("* Wrote Output Files to: ")
                    print("*     "+str(evtdir))
                print("*"*50)


if __name__ == "__main__":

    # Run main program
    main()
