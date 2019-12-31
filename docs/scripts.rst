Scripts
=======

There are several Python scripts that accompany :mod:`~rfpy`. These can be used
in bash scripts to automate data processing. These include scripts to download 
three-component seismogram data and calculate receiver functions, and perform 
post-processing for `H-k` stacking and harmonic decomposition. All of them use 
a station database provided as a :class:`~stdb.StDb` dictionary. 


``rfpy_calc.py``
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

    $ rfpy_calc.py -h
    Usage: rfpy_calc.py [options] <station database>

    Script used to download and pre-process three-component (Z, N, and E),
    seismograms for individual events and calculate teleseismic P-wave receiver
    functionsThis version requests data on the fly for a given date range. Data
    are requested from the internet using the client services framework. The
    stations are processed one by one and the data are stored to disk.

    Options:
      -h, --help            show this help message and exit
      --keys=STKEYS         Specify a comma separated list of station keys for
                            which to perform the analysis. These must be contained
                            within the station database. Partial keys will be used
                            to match against those in the dictionary. For
                            instance, providing IU will match with all stations in
                            the IU network [Default processes all stations in the
                            database]
      -v, -V, --verbose     Specify to increase verbosity.
      -O, --overwrite       Force the overwriting of pre-existing data. [Default
                            False]

      Server Settings:
        Settings associated with which datacenter to log into.

        -S SERVER, --Server=SERVER
                            Specify the server to connect to. Options include:
                            BGR, ETH, GEONET, GFZ, INGV, IPGP, IRIS, KOERI, LMU,
                            NCEDC, NEIP, NERIES, ODC, ORFEUS, RESIF, SCEDC, USGS,
                            USP. [Default IRIS]
        -U USERAUTH, --User-Auth=USERAUTH
                            Enter your IRIS Authentification Username and Password
                            (--User-Auth='username:authpassword') to access and
                            download restricted data. [Default no user and
                            password]

      Local Data Settings:
        Settings associated with defining and using a local data base of pre-
        downloaded day-long SAC files.

        --local-data=LOCALDATA
                            Specify a comma separated list of paths containing
                            day-long sac files of data already downloaded. If data
                            exists for a seismogram is already present on disk, it
                            is selected preferentially over downloading the data
                            using the Client interface
        --no-data-zero      Specify to force missing data to be set as zero,
                            rather than default behaviour which sets to nan.

      Event Settings:
        Settings associated with refining the events to include in matching
        station pairs

        --start=STARTT      Specify a UTCDateTime compatible string representing
                            the start time for the event search. This will
                            override any station start times. [Default start date
                            of station]
        --end=ENDT          Specify a UTCDateTime compatible string representing
                            the end time for the event search. This will override
                            any station end times [Default end date of station]
        -R, --reverse       Reverse order of events. Default behaviour starts at
                            oldest event and works towards most recent. Specify
                            reverse order and instead the program will start with
                            the most recent events and work towards older
        --minmag=MINMAG     Specify the minimum magnitude of event for which to
                            search. [Default 6.0]
        --maxmag=MAXMAG     Specify the maximum magnitude of event for which to
                            search. [Default None, i.e. no limit]
        --dts=DTS           Specify the window length in sec (symmetric about
                            arrival time). [Default 120.]

      Geometry Settings:
        Settings associatd with the event-station geometries

        --mindist=MINDIST   Specify the minimum great circle distance (degrees)
                            between the station and event. [Default 30.]
        --maxdist=MAXDIST   Specify the maximum great circle distance (degrees)
                            between the station and event. [Default 120.]

      Parameter Settings:
        Miscellaneous default values and settings

        --sampling-rate=NEW_SAMPLING_RATE
                            Specify new sampling rate in Hz. [Default 5.]
        --align=ALIGN       Specify component alignment key. Can be either ZRT,
                            LQT, or PVH. [Default ZRT]
        --vp=VP             Specify near-surface Vp (km/s). [Default 6.0]
        --vs=VS             Specify near-surface Vs (km/s). [Default 3.6]
        --dt_snr=DT_SNR     Specify the window length over which to calculate the
                            SNR in sec. [Default 30.]
        --fmin=FMIN         Specify the minimum frequency corner for SNR filter
                            (Hz). [Default 0.1]
        --fmax=FMAX         Specify the maximum frequency corner for SNR filter
                            (Hz). [Default 1.0]
        --twin=TWIN         Specify the source time duration for deconvolution
                            (sec). [Default 30.]


``rfpy_hk.py``
++++++++++++++

Description
-----------

Loads radial-component receiver function data available on disk
and calculates Moho depth (H) and Vp/Vs (k) of the assumed 1D
crustal structure. Station selection is specified by a network and 
station code. The database is provided as a :class:`~stdb.StDb` dictionary.

Usage
-----

.. code-block::

    $ rfpy_hk.py -h
    Usage: rfpy_hk.py [options] <station database>

    Script used to process receiver function data for H-k stacking.

    Options:
      -h, --help            show this help message and exit
      --keys=STKEYS         Specify a comma separated list of station keys for
                            which to perform the analysis. These must be contained
                            within the station database. Partial keys will be used
                            to match against those in the dictionary. For
                            instance, providing IU will match with all stations in
                            the IU network [Default processes all stations in the
                            database]
      -v, -V, --verbose     Specify to increase verbosity.
      -O, --overwrite       Force the overwriting of pre-existing data. [Default
                            False]

      Time Settings:
        Settings associated with refining the times to include in searching
        for receiver function data

        --start=STARTT      Specify a UTCDateTime compatible string representing
                            the start time for the search. This will override any
                            station start times. [Default start date of station]
        --end=ENDT          Specify a UTCDateTime compatible string representing
                            the end time for the search. This will override any
                            station end times [Default end date of station]

      Pre-processing Settings:
        Options for pre-processing of receiver function data prior to H-k
        stacking

        --freqs=FREQS       Specify a list of two floats with the minimum and
                            maximum frequency corner for the bandpass filter (Hz).
                            [Default [0.05, 0.5]]
        --bin=NBIN          Specify integer number of slowness bins to consider.
                            Use realistic bin number around 20 to start. [Default
                            does not bin data]
        --copy              Set this option to use a copy of the radial component
                            filtered at different corners for the Pps and Pss
                            phases. [Default False]
        --freqs_copy=FREQS_COPY
                            Specify a list of two floats with minimum and
                            maximumfrequency for the copies stream (Hz). [Default
                            [0.05, 0.35]]

      Settings for H-k Stacking:
        Specify parameters of H-k search, includingbounds on search, weights,
        type of stacking, etc.

        --hbound=HBOUND     Specify a list of two floats with minimum and
                            maximumbounds on Moho depth (H; km). [Default [20.,
                            50.]]
        --dh=DH             Specify interval in H for search (km). [Default 0.5]
        --kbound=KBOUND     Specify a list of two floats with minimum and
                            maximumbounds on Moho depth (H; km). [Default [1.56,
                            2.1]]
        --dk=DK             Specify interval in k for search. [Default 0.02]
        --weights=WEIGHTS   Specify a list of three floats with for Ps, Pps and
                            Pass weights in final stack. [Default [0.5, 2., -1.]]
        --type=TYP          Specify type of final stacking. Options are: 'sum' for
                            a weighted average (using weights), or 'prod' for the
                            product of positive values in stacks. [Default 'sum']

      Model Settings:
        Miscellaneous default values and settings

        --vp=VP             Specify mean crustal Vp (km/s). [Default 6.0]
        --strike=STRIKE     Specify the strike of dipping Moho. [Default None]
        --dip=DIP           Specify the dip of dipping Moho. [Default None]


``rfpy_harmonics.py``
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

    $ rfpy_harmonics.py -h
    Usage: rfpy_harmonics.py [options] <station database>

    Script used to process receiver function data for harmonic decomposition.

    Options:
      -h, --help         show this help message and exit
      --keys=STKEYS      Specify a comma separated list of station keys for which
                         to perform the analysis. These must be contained within
                         the station database. Partial keys will be used to match
                         against those in the dictionary. For instance, providing
                         IU will match with all stations in the IU network
                         [Default processes all stations in the database]
      -v, -V, --verbose  Specify to increase verbosity.
      -O, --overwrite    Force the overwriting of pre-existing data. [Default
                         False]

      Time Settings:
        Settings associated with refining the times to include in searching
        for receiver function data

        --start=STARTT   Specify a UTCDateTime compatible string representing the
                         start time for the search. This will override any station
                         start times. [Default start date of station]
        --end=ENDT       Specify a UTCDateTime compatible string representing the
                         end time for the search. This will override any station
                         end times [Default end date of station]

      Pre-processing Settings:
        Options for pre-processing of receiver function data prior to harmonic
        decomposition

        --freqs=FREQS    Specify a list of two floats with the minimum and maximum
                         frequency corner for the bandpass filter (Hz). [Default
                         [0.05, 0.5]]
        --bin=NBIN       Specify integer number of back-azimuth bins to consider
                         (typically 36 or 72). [Default does not bin data]

      Settings for harmonic decomposition:
        Specify parameters for the decomposition, e.g. a fixed azimuth, depth
        range for finding the optimal azimuth, etc.

        --azim=AZIM      Specify the azimuth angle along with to perform the
                         decomposition. [Default 0.]
        --find-azim      Set this option to calculate the optimal azimuth.
                         [Default uses the '--azim' value]
        --trange=TRANGE  Specify a list of two floats with minimum and
                         maximumbounds on time range for finding the optimal
                         azimuth (sec). [Default [0., 10.] when '--find-azim' is
                         set]
