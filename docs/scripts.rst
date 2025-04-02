Scripts
=======

There are several Python scripts that accompany :mod:`~rfpy`, which can be used
in bash scripts to automate data processing. These include scripts to download 
three-component seismogram data and calculate receiver functions, and perform 
post-processing for `H-k` stacking and harmonic decomposition. All of them use 
a station database provided as a :class:`~stdb.StDb` dictionary. 


``rfpy_calc``
++++++++++++++++

Description
-----------

Downloads three-component ('Z', 'N' and 'E') seismograms based
on available times of earthquakes and performs `P`-wave receiver function
calculation. Station selection is specified by a network and 
station code. The database is provided as a :class:`~stdb.StDb` dictionary.

Usage
-----

.. code-block::

    $ rfpy_calc -h

    ##############################################
    #        __                          _       #
    #  _ __ / _|_ __  _   _     ___ __ _| | ___  #
    # | '__| |_| '_ \| | | |   / __/ _` | |/ __| #
    # | |  |  _| |_) | |_| |  | (_| (_| | | (__  #
    # |_|  |_| | .__/ \__, |___\___\__,_|_|\___| #
    #          |_|    |___/_____|                #
    #                                            #
    ##############################################

    usage: rfpy_calc [arguments] <station database>

    Script used to download and pre-process three-component ('Z', 'N', and 'E'),
    seismograms for individual events and calculate teleseismic P-wave receiver
    functionsThis version requests data on the fly for a given date range. Data
    are requested from the internet using the client services framework. The
    stations are processed one by one and the data are stored to disk.

    options:
      -h, --help            show this help message and exit
      --keys STKEYS         Specify a comma separated list of station keys for
                            which to perform the analysis. These must be contained
                            within the station database. Partial keys will be used
                            to match against those in the dictionary. For
                            instance, providing IU will match with all stations in
                            the IU network [Default processes all stations in the
                            database]
      -v, -V, --verbose     Specify to increase verbosity.
      -O, --overwrite       Force the overwriting of pre-existing data. [Default
                            False]
      --zcomp ZCOMP         Specify the Vertical Component Channel Identifier.
                            [Default Z].
      -L, --long-name       Force folder names to use long-key form (NET.STN.CHN).
                            Default behaviour uses short key form (NET.STN) for
                            the folder names, regardless of the key type of the
                            database.

    Server Settings:
      Settings associated with which datacenter to log into.

  --server SERVER       Base URL of FDSN web service compatible server (e.g.
                        “http://service.iris.edu”) or key string for
                        recognized server (one of 'AUSPASS', 'BGR',
                        'EARTHSCOPE', 'EIDA', 'EMSC', 'ETH', 'GEOFON',
                        'GEONET', 'GFZ', 'ICGC', 'IESDMC', 'INGV', 'IPGP',
                        'IRIS', 'IRISPH5', 'ISC', 'KNMI', 'KOERI', 'LMU',
                        'NCEDC', 'NIEP', 'NOA', 'NRCAN', 'ODC', 'ORFEUS',
                        'RASPISHAKE', 'RESIF', 'RESIFPH5', 'SCEDC', 'TEXNET',
                        'UIB-NORSAR', 'USGS', 'USP'). [Default 'IRIS']
  --user-auth USERAUTH  Authentification Username and Password for the
                        waveform server (--user-auth='username:authpassword')
                        to access and download restricted data. [Default no
                        user and password]
  --eida-token TOKENFILE
                        Token for EIDA authentication mechanism, see
                        http://geofon.gfz-
                        potsdam.de/waveform/archive/auth/index.php. If a token
                        is provided, argument --user-auth will be ignored.
                        This mechanism is only available on select EIDA nodes.
                        The token can be provided in form of the PGP message
                        as a string, or the filename of a local file with the
                        PGP message in it. [Default None]

    Local Data Settings:
      Settings associated with defining and using a local data base of pre-
      downloaded day-long SAC files.

      --local-data LOCALDATA
                            Specify a comma separated list of paths containing
                            day-long sac files of data already downloaded. If data
                            exists for a seismogram is already present on disk, it
                            is selected preferentially over downloading the data
                            using the Client interface
      --dtype DTYPE         Specify the data archive file type, either SAC or
                            MSEED. Note the default behaviour is to search for SAC
                            files. Local archive files must have extensions of
                            '.SAC' or '.MSEED. These are case dependent, so
                            specify the correct casehere.
      --no-data-zero        Specify to force missing data to be set as zero,
                            rather than default behaviour which sets to nan.
      --no-local-net        Specify to prevent using the Network code in the
                            search for local data (sometimes for CN stations the
                            dictionary name for a station may disagree with that
                            in the filename. [Default Network used]
      --save-Z12            Specify to save Z12 (un-rotated) components. [Default
                            False]

    Event Settings:
      Settings associated with refining the events to include in matching event-
      station pairs

      --start STARTT        Specify a UTCDateTime compatible string representing
                            the start time for the event search. This will
                            override any station start times. [Default start date
                            of station]
      --end ENDT            Specify a UTCDateTime compatible string representing
                            the end time for the event search. This will override
                            any station end times [Default end date of station]
      --reverse, -R         Reverse order of events. Default behaviour starts at
                            oldest event and works towards most recent. Specify
                            reverse order and instead the program will start with
                            the most recent events and work towards older
      --minmag MINMAG       Specify the minimum magnitude of event for which to
                            search. [Default 6.0]
      --maxmag MAXMAG       Specify the maximum magnitude of event for which to
                            search. [Default None, i.e. no limit]

    Geometry Settings:
      Settings associatd with the event-station geometries for the specified
      phase

      --phase PHASE         Specify the phase name to use. Be careful with the
                            distance. setting. Options are 'P' or 'PP'. [Default
                            'P']
      --mindist MINDIST     Specify the minimum great circle distance (degrees)
                            between the station and event. [Default depends on
                            phase]
      --maxdist MAXDIST     Specify the maximum great circle distance (degrees)
                            between the station and event. [Default depends on
                            phase]

    Parameter Settings:
      Miscellaneous default values and settings

      --sampling-rate NEW_SAMPLING_RATE
                            Specify new sampling rate in Hz. [Default 10.]
      --dts DTS             Specify the window length in sec (symmetric about
                            arrival time). [Default 150.]
      --align ALIGN         Specify component alignment key. Can be either ZRT,
                            LQT, or PVH. [Default ZRT]
      --vp VP               Specify near-surface Vp to use with --align=PVH
                            (km/s). [Default 6.0]
      --vs VS               Specify near-surface Vs to use with --align=PVH
                            (km/s). [Default 3.5]
      --dt-snr DT_SNR       Specify the window length over which to calculate the
                            SNR in sec. [Default 30.]
      --pre-filt PRE_FILT   Specify two floats with low and high frequency corners
                            for pre-filter (before deconvolution). [Default None]
      --fmin FMIN           Specify the minimum frequency corner for SNR and CC
                            filter (Hz). [Default 0.05]
      --fmax FMAX           Specify the maximum frequency corner for SNR and CC
                            filter (Hz). [Default 1.0]

    Deconvolution Settings:
      Parameters for deconvolution

      --method METHOD       Specify the deconvolution method. Available methods
                            include 'wiener', 'water' and 'multitaper'. [Default
                            'wiener']
      --gfilt GFILT         Specify the Gaussian filter width in Hz. [Default
                            None]
      --wlevel WLEVEL       Specify the water level, used in the 'water' method.
                            [Default 0.01]


``rfpy_recalc``
++++++++++++++++

Description
-----------

Looks for available receiver functions on disk and re-calculates them
using different processing options. Station selection is specified by 
a network and station code. The database is provided as a :class:`~stdb.StDb` 
dictionary.

Usage
-----

.. code-block::

    $ rfpy_recalc -h

    ########################################################
    #                                                      #
    #        __                                    _       #
    #  _ __ / _|_ __  _   _     _ __ ___  ___ __ _| | ___  #
    # | '__| |_| '_ \| | | |   | '__/ _ \/ __/ _` | |/ __| #
    # | |  |  _| |_) | |_| |   | | |  __/ (_| (_| | | (__  #
    # |_|  |_| | .__/ \__, |___|_|  \___|\___\__,_|_|\___| #
    #          |_|    |___/_____|                          #
    #                                                      #
    ########################################################

    usage: rfpy_recalc [arguments] <station database>

    Script used to re-calculate receiver functions that already exist on disk, but
    using different processing options. The stations are processed one by one and
    the data are over-written to disk. 

    positional arguments:
      indb                 Station Database to process from.

    optional arguments:
      -h, --help           show this help message and exit
      --keys STKEYS        Specify a comma separated list of station keys for
                           which to perform the analysis. These must be contained
                           within the station database. Partial keys will be used
                           to match against those in the dictionary. For instance,
                           providing IU will match with all stations in the IU
                           network [Default processes all stations in the
                           database]
      -v, -V, --verbose    Specify to increase verbosity.
      -L, --long-name      Force folder names to use long-key form (NET.STN.CHN).
                           Default behaviour uses short key form (NET.STN) for the
                           folder names, regardless of the key type of the
                           database.


    Parameter Settings:
      Miscellaneous default values and settings

      --Z12                Use Z12 data if available. [Default uses ZNE data]
      --phase PHASE        Specify the phase name to use. Be careful with the
                           distance. setting. Options are 'P', 'PP', 'allP', 'S',
                           'SKS' or 'allS'. [Default 'allP']
      --resample RESAMPLE  Specify the new sampling-rate for the receiver
                           functions. Note the sampling rate of the original data
                           (ZNE or Z12) stored on disk is unchanged. [Default
                           None]      
      --align ALIGN        Specify component alignment key. Can be either ZRT,
                           LQT, or PVH. [Default ZRT]
      --vp VP              Specify near-surface Vp to use with --align=PVH (km/s).
                           [Default 6.0]
      --vs VS              Specify near-surface Vs to use with --align=PVH (km/s).
                           [Default 3.5]
      --dt-snr DT_SNR      Specify the window length over which to calculate the
                           SNR in sec. [Default 30.]
      --pre-filt PRE_FILT  Specify two floats with low and high frequency corners
                           for pre-filter (before deconvolution). [Default None]
      --fmin FMIN          Specify the minimum frequency corner for SNR filter
                           (Hz). [Default 0.05]
      --fmax FMAX          Specify the maximum frequency corner for SNR filter
                           (Hz). [Default 1.0]

    Deconvolution Settings:
      Parameters for deconvolution

      --method METHOD      Specify the deconvolution method. Available methods
                           include 'wiener', 'water' and 'multitaper'. [Default
                           'wiener']
      --gfilt GFILT        Specify the Gaussian filter width in Hz. [Default None]
      --wlevel WLEVEL      Specify the water level, used in the 'water' method.
                           [Default 0.01]


``rfpy_plot``
++++++++++++++++

Description
-----------

Script used to make plots of receiver function panels sorted by
back-azimuth (averaging all slowness information) or by slowness
(averaging all back-azimuth information).

Usage
-----

.. code-block::

    $ rfpy_plot -h

    #################################################
    #        __                        _       _    #
    #  _ __ / _|_ __  _   _      _ __ | | ___ | |_  #
    # | '__| |_| '_ \| | | |    | '_ \| |/ _ \| __| #
    # | |  |  _| |_) | |_| |    | |_) | | (_) | |_  #
    # |_|  |_| | .__/ \__, |____| .__/|_|\___/ \__| #
    #          |_|    |___/_____|_|                 #
    #                                               #
    #################################################

    usage: rfpy_plot [arguments] <station database>

    Script used to plot receiver function data

    positional arguments:
      indb                  Station Database to process from.

    optional arguments:
      -h, --help            show this help message and exit
      --keys STKEYS         Specify a comma separated list of station keys for
                            which to perform the analysis. These must be contained
                            within the station database. Partial keys will be used
                            to match against those in the dictionary. For
                            instance, providing IU will match with all stations in
                            the IU network [Default processes all stations in the
                            database]
      -v, -V, --verbose     Specify to increase verbosity.
      -O, --overwrite       Force the overwriting of pre-existing figures.
                            [Default False]
      -L, --long-name       Force folder names to use long-key form (NET.STN.CHN).
                            Default behaviour uses short key form (NET.STN) for
                            the folder names, regardless of the key type of the
                            database.

    Pre-processing Settings:
      Options for pre-processing of receiver function data before plotting

      --snr SNR             Specify the vertical component SNR threshold for
                            extracting receiver functions. [Default 5.]
      --snrh SNRH           Specify the horizontal component SNR threshold for
                            extracting receiver functions. [Default None]
      --cc CC               Specify the CC threshold for extracting receiver
                            functions. [Default None]
      --no-outlier          Set this option to delete outliers based on the MAD on
                            the variance. [Default False]
      --binlim BINLIM       Specify the minimum number of RFs in each bin.
                            [Default 1]
      --bp BP               Specify the corner frequencies for the bandpass
                            filter. [Default no filtering]
      --pws                 Set this option to use phase-weighted stacking during
                            binning [Default False]
      --nbaz NBAZ           Specify integer number of back-azimuth bins to
                            consider (typically 36 or 72). If not None, the plot
                            will show receiver functions sorted by back-azimuth
                            values. [Default None]
      --nslow NSLOW         Specify integer number of slowness bins to consider
                            (typically 20 or 40). If not None, the plot will show
                            receiver functions sorted by slowness values. [Default
                            None]
      --slowbound SLOWBOUND
                            Specify a list of two floats with minimum and
                            maximumbounds on slowness (s/km). [Default [0.04,
                            0.08]]
      --bazbound BAZBOUND   Specify a list of two floats with minimum and
                            maximumbounds on back azimuth (degrees). [Default [0,
                            360]]
      --phase PHASE         Specify the phase name to plot. Options are 'P', 'PP',
                            'allP', 'S', 'SKS' or 'allS'. [Default 'allP']

    Plot Settings:
      Options for plot format

      --scale SCALE         Specify the scaling factor for the amplitude of the
                            receiver functions in the wiggle plots. [Default 100.
                            for a back-azimuth plot, 0.02 for a slowness plot]
      --normalize           Set this option to produce receiver functions
                            normalized by the max amplitude of stacked RFs.
                            [Default False]
      --trange TRANGE       Specify the time range for the x-axis (sec). Negative
                            times are allowed [Default 0., 30.]
      --stacked             Set this option to plot a stack of all traces in top
                            panel. [Default does not plot stacked traces]
      --save                Set this option if you wish to save the figure.
                            [Default does not save figure]
      --title TITLEPLOT     Specify title of figure. [Default None]
      --format FORM         Specify format of figure. Can be any one of the
                            validmatplotlib formats: 'png', 'jpg', 'eps', 'pdf'.
                            [Default 'png']
      --plot-event-dist     Plot distribution of events on map. Other Plotting
                            Options will be applied to this figure (title, save,
                            etc.). [Default no plot]


``rfpy_hk``
++++++++++++++

Description
-----------

Loads radial-component receiver function data available on disk
and calculates Moho depth ('H') and Vp/Vs ('k') of the assumed 1D
crustal structure. Station selection is specified by a network and 
station code. The database is provided as a :class:`~stdb.StDb` dictionary.

Usage
-----

.. code-block::

    $ rfpy_hk -h

    #########################################
    #        __                 _     _     #
    #  _ __ / _|_ __  _   _    | |__ | | __ #
    # | '__| |_| '_ \| | | |   | '_ \| |/ / #
    # | |  |  _| |_) | |_| |   | | | |   <  #
    # |_|  |_| | .__/ \__, |___|_| |_|_|\_\ #
    #          |_|    |___/_____|           #
    #                                       #
    #########################################

    usage: rfpy_hk [arguments] <station database>

    Script used to process receiver function data for H-k stacking.

    positional arguments:
      indb                  Station Database to process from.

    optional arguments:
      -h, --help            show this help message and exit
      --keys STKEYS         Specify a comma separated list of station keys for
                            which to perform the analysis. These must be contained
                            within the station database. Partial keys will be used
                            to match against those in the dictionary. For
                            instance, providing IU will match with all stations in
                            the IU network [Default processes all stations in the
                            database]
      -v, -V, --verbose     Specify to increase verbosity.
      -O, --overwrite       Force the overwriting of pre-existing data. [Default
                            False]
      -L, --long-name       Force folder names to use long-key form (NET.STN.CHN).
                            Default behaviour uses short key form (NET.STN) for
                            the folder names, regardless of the key type of the
                            database.

    Time Settings:
      Settings associated with refining the times to include in searching for
      receiver function data

      --start STARTT        Specify a UTCDateTime compatible string representing
                            the start time for the search. This will override any
                            station start times. [Default start date of station]
      --end ENDT            Specify a UTCDateTime compatible string representing
                            the end time for the search. This will override any
                            station end times [Default end date of station]

    Pre-processing Settings:
      Options for pre-processing of receiver function data prior to H-k stacking

      --binlim BINLIM       Specify the minimum number of RFs in each bin.
                            [Default 3]
      --bp BP               Specify the corner frequencies for the bandpass
                            filter. [Default 0.05,0.5]
      --nbaz NBAZ           Specify integer number of back-azimuth bins to
                            consider. [Default 36]
      --nslow NSLOW         Specify integer number of slowness bins to consider.
                            [Default 40]
      --snr SNR             Specify the SNR threshold for extracting receiver
                            functions. [Default None]
      --snrh SNRH           Specify the horizontal component SNR threshold for
                            extracting receiver functions. [Default None]
      --cc CC               Specify the CC threshold for extracting receiver
                            functions. [Default None]
      --no-outlier          Set this option to delete outliers based on the MAD on
                            the variance. [Default False]
      --slowbound SLOWBOUND
                            Specify a list of two floats with minimum and
                            maximumbounds on slowness (s/km). [Default [0.04,
                            0.08]]
      --bazbound BAZBOUND   Specify a list of two floats with minimum and
                            maximumbounds on back azimuth (degrees). [Default [0,
                            360]]
      --pws                 Set this option to use phase-weighted stacking during
                            binning [Default False]
      --phase PHASE         Specify the phase name to plot. Options are 'P', 'PP',
                            'allP', 'S', 'SKS' or 'allS'. [Default 'allP']
      --copy                Set this option to use a copy of the radial component
                            filtered at different corners for the Pps and Pss
                            phases. [Default False]
      --bp-copy BP_COPY     Specify a list of two floats with minimum and
                            maximumfrequency for the copied stream (Hz). [Default
                            [0.05, 0.35]]

    Settings for H-k Stacking:
      Specify parameters of H-k search, includingbounds on search, weights, type
      of stacking, etc.

      --hbound HBOUND       Specify a list of two floats with minimum and
                            maximumbounds on Moho depth (H, in km). [Default [20.,
                            50.]]
      --dh DH               Specify search interval for H (km). [Default 0.5]
      --kbound KBOUND       Specify a list of two floats with minimum and
                            maximumbounds on Vp/Vs (k). [Default [1.56, 2.1]]
      --dk DK               Specify search interval for k. [Default 0.02]
      --weights WEIGHTS     Specify a list of three floats with for Ps, Pps and
                            Pass weights in final stack. [Default [0.5, 2., -1.]]
      --type TYP            Specify type of final stacking. Options are: 'sum' for
                            a weighted average (using weights), or 'product' for
                            the product of positive values in stacks. [Default
                            'sum']
      --save                Set this option to save the HkStack object to file.
                            [Default doesn't save]

    Model Settings:
      Miscellaneous default values and settings

      --vp VP               Specify mean crustal Vp (km/s). [Default 6.0]
      --strike STRIKE       Specify the strike of dipping Moho. [Default None]
      --dip DIP             Specify the dip of dipping Moho. [Default None]

    Settings for plotting results:
      Specify parameters for plotting the H-k stacks.

      --plot                Set this option to produce a plot of the stacks
                            [Default does not produce plot]
      --save-plot           Set this option to save the plot [Default doesn't
                            save]
      --title TITLE         Specify plot title [Default has no title]
      --format FORM         Specify format of figure. Can be any one of the
                            validmatplotlib formats: 'png', 'jpg', 'eps', 'pdf'.
                            [Default 'png']


``rfpy_harmonics``
+++++++++++++++++++++

Description
-----------

Loads radial and transverse component receiver function data available on disk
and decomposes them into back-azimuth harmonics. Station selection is specified 
by a network and station code. The database is provided as a :class:`~stdb.StDb` 
dictionary.

Usage
-----

.. code-block::

    $ rfpy_harmonics -h

    ################################################################################
    #        __                 _                                      _           #
    #  _ __ / _|_ __  _   _    | |__   __ _ _ __ _ __ ___   ___  _ __ (_) ___ ___  #
    # | '__| |_| '_ \| | | |   | '_ \ / _` | '__| '_ ` _ \ / _ \| '_ \| |/ __/ __| #
    # | |  |  _| |_) | |_| |   | | | | (_| | |  | | | | | | (_) | | | | | (__\__ \ #
    # |_|  |_| | .__/ \__, |___|_| |_|\__,_|_|  |_| |_| |_|\___/|_| |_|_|\___|___/ #
    #          |_|    |___/_____|                                                  #
    #                                                                              #
    ################################################################################

    usage: rfpy_harmonics [arguments] <station database>

    Script used to process receiver function data for harmonic decomposition.

    positional arguments:
      indb               Station Database to process from.

    optional arguments:
      -h, --help         show this help message and exit
      --keys STKEYS      Specify a comma separated list of station keys for which
                         to perform the analysis. These must be contained within
                         the station database. Partial keys will be used to match
                         against those in the dictionary. For instance, providing
                         IU will match with all stations in the IU network
                         [Default processes all stations in the database]
      -v, -V, --verbose  Specify to increase verbosity.
      -O, --overwrite    Force the overwriting of pre-existing data. [Default
                         False]
      -L, --long-name    Force folder names to use long-key form (NET.STN.CHN).
                         Default behaviour uses short key form (NET.STN) for the
                         folder names, regardless of the key type of the database.

    Time Settings:
      Settings associated with refining the times to include in searching for
      receiver function data

      --start STARTT     Specify a UTCDateTime compatible string representing the
                         start time for the search. This will override any station
                         start times. [Default start date of station]
      --end ENDT         Specify a UTCDateTime compatible string representing the
                         end time for the search. This will override any station
                         end times [Default end date of station]

    Pre-processing Settings:
      Options for pre-processing of receiver function data prior to harmonic
      decomposition

      --bp BP            Specify the corner frequencies for the bandpass filter.
                         [Default 0.05,0.5]
      --bin NBIN         Specify integer number of back-azimuth bins to consider
                         (typically 36 or 72). [Default does not bin data]
      --snr SNR          Specify the SNR threshold for extracting receiver
                         functions. [Default None]
      --snrh SNRH        Specify the horizontal component SNR threshold for
                         extracting receiver functions. [Default None]
      --cc CC            Specify the CC threshold for extracting receiver
                         functions. [Default None]
      --no-outlier       Set this option to delete outliers based on the MAD on
                         the variance. [Default False]
      --phase PHASE      Specify the phase name to plot. Options are 'P', 'PP',
                         'allP', 'S', 'SKS' or 'allS'. [Default 'allP']

    Settings for harmonic decomposition:
      Specify parameters for the decomposition, e.g. a fixed azimuth, depth
      range for finding the optimal azimuth, etc.

      --azim AZIM        Specify the azimuth angle along with to perform the
                         decomposition. [Default 0.]
      --find-azim        Set this option to calculate the optimal azimuth.
                         [Default uses the '--azim' value]
      --trange TRANGE    Specify a list of two floats with minimum and
                         maximumbounds on time range for finding the optimal
                         azimuth (sec). [Default [0., 10.] when '--find-azim' is
                         set]
      --save             Set this option to save the Harmonics object to a pickled
                         file. [Default does not save object]

    Settings for plotting results:
      Specify parameters for plotting the back-azimuth harmonics.

      --plot             Set this option to produce a plot of the back-azimuth
                         harmonics
      --ymax YMAX        Specify the maximum y axis value for the plot in units of
                         thedependent variable (e.g., sec). [Default 30.]
      --scale SCALE      Specify the scaling value that multiplies the amplitude
                         of the harmonic components. [Default 10.]
      --save-plot        Set this option to save the plot [Default doesn't save]
      --title TITLE      Specify plot title [Default has no title]
      --format FORM      Specify format of figure. Can be any one of the
                         validmatplotlib formats: 'png', 'jpg', 'eps', 'pdf'.
                         [Default 'png']


``rfpy_ccp``
+++++++++++++++++++++

Description
-----------

Loads radial component receiver function data available on disk
and processes them for Common Conversion Point stacking along a linear
profile. The three CCP phase stacks (Ps, Pps and Pss) are averaged
using a weighted sum, or using phase-weighted stacking to downweight
incoherent signal across all stacks. The phase stacks can be further 
smoothed using a Gaussian kernel that simulates P-wave sensitivity.
Station selection is specified by a network and station code. 
The database is provided as a :class:`~stdb.StDb` dictionary.

.. note::

    The start and end coordinates (latitude, longitude) of the profile 
    must be supplied as `--start=` and `--end=` parameters. The CCP
    stacks will be projected along the line, regardless of station distance
    normal to the line. 

Usage
-----

.. code-block::

    $ rfpy_ccp -h

    ############################################
    #        __                                #
    #  _ __ / _|_ __  _   _     ___ ___ _ __   #
    # | '__| |_| '_ \| | | |   / __/ __| '_ \  #
    # | |  |  _| |_) | |_| |  | (_| (__| |_) | #
    # |_|  |_| | .__/ \__, |___\___\___| .__/  #
    #          |_|    |___/_____|      |_|     #
    #                                          #
    ############################################

    usage: rfpy_ccp [arguments] <station database>

    Script used to process receiver function data for common-conversion-point
    (CCP) imaging.

    positional arguments:
      indb                 Station Database to process from.

    optional arguments:
      -h, --help           show this help message and exit
      --keys STKEYS        Specify a comma separated list of station keys for
                           which to perform the analysis. These must be contained
                           within the station database. Partial keys will be used
                           to match against those in the dictionary. For instance,
                           providing IU will match with all stations in the IU
                           network [Default processes all stations in the
                           database]
      -v, -V, --verbose    Specify to increase verbosity.
      -O, --overwrite      Force the overwriting of pre-existing data. [Default
                           False]
      -L, --long-name      Force folder names to use long-key form (NET.STN.CHN).
                           Default behaviour uses short key form (NET.STN) for the
                           folder names, regardless of the key type of the
                           database.
                       
    Line Geometry Settings:
      Options for defining the line along which to produce the CCP image

      --start COORD_START  Specify a list of two floats with the latitude and
                           longitude of the start point, in this respective order.
                           [Exception raised if not specified]
      --end COORD_END      Specify a list of two floats with the latitude and
                           longitudeof the end point, in this respective order.
                           [Exception raised if not specified]
      --dz DZ              Specify vertical cell size in km. [Default 1.]
      --dx DX              Specify horizontal cell size in km. [Default 2.5]

    Pre-processing Settings:
      Options for pre-processing of receiver function data for CCP stacking

      --snr SNR            Specify the SNR threshold for extracting receiver
                           functions. [Default None]
      --snrh SNRH          Specify the horizontal component SNR threshold for
                           extracting receiver functions. [Default None]
      --cc CC              Specify the CC threshold for extracting receiver
                           functions. [Default None]
      --no-outlier         Set this option to delete outliers based on the MAD on
                           the variance. [Default False]
      --binlim BINLIM      Specify the minimum number of RFs in each bin. [Default
                           3]
      --f1 F1              Specify the low frequency corner for the bandpass
                           filter for all phases (Hz). [Default [0.05]]
      --f2ps F2PS          Specify the high frequency corner for the bandpass
                           filter for the Ps phase (Hz). [Default [0.75]]
      --f2pps F2PPS        Specify the high frequency corner for the bandpass
                           filter for the Pps phase (Hz). [Default [0.36]]
      --f2pss F2PSS        Specify the high frequency corner for the bandpass
                           filter for the Pss phase (Hz). [Default [0.3]]
      --nbaz NBAZ          Specify integer number of back-azimuth bins to
                           consider. [Default 36]
      --nslow NSLOW        Specify integer number of slowness bins to consider.
                           [Default 40]
      --wlen WLEN          Specify wavelength of P-wave as sensitivity (km).
                           [Default 35.]
      --phase PHASE        Specify the phase name to plot. Options are 'P', 'PP',
                           'allP', 'S', 'SKS' or 'allS'. [Default 'allP']

    CCP Settings:
      Options for specifying the type of CCP stacking to perform

      --load               Step 1. Set this option to load rfstreams into CCPimage
                           object. [Default False]
      --prep               Step 2. Set this option to prepare CCPimage before pre-
                           stacking. [Default False]
      --prestack           Step 3. Set this option to prestack all phases before
                           CCP averaging. [Default False]
      --ccp                Step 4a. Set this option for standard CCP stacking with
                           multiples. [Default False]
      --gccp               Step 4b. Set this option for Gaussian-weighted CCP
                           stacking with multiples. [Default False]
      --linear             Step 5a. Set this option to produce a linear, weighted
                           stack for the final [G]CCP image. [Default True unless
                           --phase is set]
      --pws                Step 5b. Set this option to produce a phase weighted
                           stack for the final [G]CCP image. [Default False]
      --weights WEIGHTS    Option to define weights for each of the three phases:
                           Ps, Pps and Pss, by specifying three comma-separated
                           floats. [Default 1., 3., -3.]

    Figure Settings:
      Options for specifying the settings for the final figure

      --figure             Set this option to plot the final [G]CCP figure.
                           [Default False]
      --cbound CBOUND      Set the maximum value for the color palette. [Default
                           0.05 for --ccp or 0.015 for --gccp]
      --save-fig           Set this option to save the final [G]CCP figure. This
                           option can only be set if --figure is also set.[Default
                           False]
      --title TITLE        Set Figure title. [Default None]
      --format FMT         Set format of figure. You can choose among 'png',
                           'jpg', 'eps', 'pdf'. [Default 'png']