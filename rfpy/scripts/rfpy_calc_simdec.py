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
from argparse import ArgumentParser
from os.path import exists as exist
import numpy as np
import pickle
import stdb
from obspy import Stream
from rfpy import RFData, binning, plotting
from pathlib import Path


def get_simdec_arguments(argv=None):
    """
    Get Options from :class:`~optparse.OptionParser` objects.

    This function is used for data processing on-the-fly (requires web connection)

    """

    parser = ArgumentParser(
        usage="%(prog)s [arguments] <station database>",
        description="Script used to calculate receiver functions " +
        "using simultaneous deconvolution. The stations are processed one " +
        "by one and the data are stored to disk. ")

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
        "-L", "--long-name",
        action="store_true",
        dest="lkey",
        default=False,
        help="Force folder names to use long-key form (NET.STN.CHN). " +
        "Default behaviour uses short key form (NET.STN) for the folder " +
        "names, regardless of the key type of the database."
    )

    # Constants Settings
    ConstGroup = parser.add_argument_group(
        title='Parameter Settings',
        description="Miscellaneous default values and settings")
    ConstGroup.add_argument(
        "--Z12",
        action="store_true",
        dest="Z12",
        default=False,
        help="Use Z12 data if available. [Default uses ZNE data]")
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
        "--resample",
        action="store",
        type=float,
        dest="resample",
        default=None,
        help="Specify the new sampling-rate for the receiver functions. " +
        "Note the sampling rate of the original data (ZNE or Z12) stored " +
        "on disk is unchanged. [Default None]")
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

    # Deconvolution Settings
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

    PreGroup = parser.add_argument_group(
        title='Pre-processing Settings',
        description="Options for pre-processing " +
        "data before stacking")
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
        "--nbaz",
        action="store",
        dest="nbaz",
        type=int,
        default=36,
        help="Specify integer number of back-azimuth bins to consider " +
        "(typically 36 or 72). If not None, the plot will show receiver " +
        "functions sorted by back-azimuth values. [Default 36]")
    PreGroup.add_argument(
        "--nslow",
        action="store",
        dest="nslow",
        type=int,
        default=20,
        help="Specify integer number of slowness bins to consider " +
        "(range of slowness values is 0.04 to 0.08). " + 
        "If not None, the plot will show receiver " +
        "functions sorted by slowness values. [Default 20]")
    PreGroup.add_argument(
        "--bp",
        action="store",
        type=str,
        dest="bp",
        default=None,
        help="Specify the corner frequencies for the bandpass filter. " +
        "[Default no filtering]")

    PlotGroup = parser.add_argument_group(
        title='Plot Settings',
        description="Options for plot format")
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

    if args.pre_filt is not None:
        args.pre_filt = [float(val) for val in args.pre_filt.split(',')]
        args.pre_filt = sorted(args.pre_filt)
        if (len(args.pre_filt)) != 2:
            parser.error(
                "Error: --pre-filt should contain 2 " +
                "comma-separated floats")

    return args


def main():

    print()
    print("##################################################################################")
    print("#                                                                                #")
    print("#        __                          _              _               _            #")
    print("#  _ __ / _|_ __  _   _     ___ __ _| | ___     ___(_)_ __ ___   __| | ___  ___  #")
    print("# | '__| |_| '_ \| | | |   / __/ _` | |/ __|   / __| | '_ ` _ \ / _` |/ _ \/ __| #")
    print("# | |  |  _| |_) | |_| |  | (_| (_| | | (__    \__ \ | | | | | | (_| |  __/ (__  #")
    print("# |_|  |_| | .__/ \__, |___\___\__,_|_|\___|___|___/_|_| |_| |_|\__,_|\___|\___| #")
    print("#          |_|    |___/_____|             |_____|                                #")
    print("#                                                                                #")
    print("##################################################################################")
    print()

    # Run Input Parser
    args = get_simdec_arguments()

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

        # # Check for folder already processed
        # if stfld in procfold:
        #     print('  {0} already processed...skipping   '.format(stfld))
        #     continue

        LLstreams = Stream()
        QLstreams = Stream()
        TLstreams = Stream()
        NNstreams = Stream()

        datafiles = [x for x in datapath.iterdir() if x.is_dir()]
        for folder in datafiles:

            # Skip hidden folders
            if folder.name.startswith('.'):
                continue

            # Re-initialize RFData object
            rfdata = RFData(sta)

            # Load meta data
            metafile = folder / "Meta_Data.pkl"
            if not metafile.is_file():
                continue
            rfdata.meta = pickle.load(open(metafile, 'rb'))

            # Skip data not in list of phases
            if rfdata.meta.phase not in args.listphase:
                continue

            if args.verb:
                print("* Station: {0}; folder: {1}".format(stkey, folder))

            if args.Z12:
                try:
                    # Try loading Z12_Data
                    Z12file = folder / "Z12_Data.pkl"
                    rfdata.data = pickle.load(open(Z12file, 'rb'))
                    # Remove rotated flag and snr flag
                    rfdata.meta.rotated = False
                    rfdata.rotate(align='ZNE')
                except:

                    print("Z12_Data.pkl not available - using ZNE_Data.pkl")
                    # Load ZNE data
                    ZNEfile = folder / "ZNE_Data.pkl"
                    rfdata.data = pickle.load(open(ZNEfile, 'rb'))
            else:
                # Load ZNE data
                ZNEfile = folder / "ZNE_Data.pkl"
                rfdata.data = pickle.load(open(ZNEfile, 'rb'))

            # Resample if requested
            if args.resample:
                rfdata.data.resample(
                    args.resample, no_filter=False)

            # Remove rotated flag and snr flag
            rfdata.meta.rotated = False
            rfdata.meta.snr = None

            # Rotate from ZNE to 'align' ('ZRT', 'LQT', or 'PVH')
            rfdata.rotate(vp=args.vp, vs=args.vs, align=args.align)

            # Calculate SNR
            rfdata.calc_snr(dt=args.dt_snr, fmin=args.fmin, fmax=args.fmax)

            # QC Thresholding
            if rfdata.meta.snrh < args.snrh:
                continue
            if rfdata.meta.snr < args.snr:
                continue

            # Calculate spectra
            rfdata.calc_spectra(
                    vp=args.vp, vs=args.vs,
                    align=args.align, method=args.method,
                    wavelet='complete', envelope_threshold=0.05, 
                    time=5, norm=True)

            # Convert to Stream
            rfstream = rfdata.to_stream(store='specs')

            # Add to Stream objects
            LLstreams.append(rfstream[0])
            QLstreams.append(rfstream[1])
            TLstreams.append(rfstream[2])
            NNstreams.append(rfstream[3])

        # Bin into baz and slowness bins
        LL_tmp, = binning.bin_baz_slow(LLstreams, nbaz=args.nbaz+1, nslow=args.nslow+1,
            pws=False, include_empty=False)
        QL_tmp, = binning.bin_baz_slow(QLstreams, nbaz=args.nbaz+1, nslow=args.nslow+1,
            pws=False, include_empty=False)
        TL_tmp, = binning.bin_baz_slow(TLstreams, nbaz=args.nbaz+1, nslow=args.nslow+1,
            pws=False, include_empty=False)
        NN_tmp, = binning.bin_baz_slow(NNstreams, nbaz=args.nbaz+1, nslow=args.nslow+1,
            pws=False, include_empty=False)

        # Initialize new empty streams
        RFL = Stream()
        RFQ = Stream()
        RFT = Stream()

        # Cycle through binned spectra
        for i in range(len(LL_tmp)):

            # Initialize new RFData object with empty Meta header
            rfdata = RFData(sta)

            # Extract the specs from the binned streams
            rfdata.specs = Stream(traces=[LL_tmp[i], QL_tmp[i], TL_tmp[i], NN_tmp[i]])

            # Deconvolve
            rfdata.deconvolve(
                align=args.align, method=args.method,
                gfilt=args.gfilt, wlevel=args.wlevel)

            # Append RFs to streams
            RFL.append(rfdata.rf[0])
            RFQ.append(rfdata.rf[1])
            RFT.append(rfdata.rf[2])

        # Filter
        if args.bp:
            RFQ.filter('bandpass', freqmin=args.bp[0],
                             freqmax=args.bp[1], corners=2,
                             zerophase=True)
            RFT.filter('bandpass', freqmin=args.bp[0],
                             freqmax=args.bp[1], corners=2,
                             zerophase=True)

        plotting.panel(RFQ, RFT, tmin=args.trange[0], tmax=args.trange[1],
            normalize=True, save=args.saveplot,
                title=args.titleplot, form=args.form)
        # # Update processed folders
        # procfold.append(stfld)


if __name__ == "__main__":

    # Run main program
    main()
