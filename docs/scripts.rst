Scripts
=======

There are several Python scripts that accompany :mod:`~rfpy`. These can be used
in bash scripts to automate data processing. These include scripts to download 
three-component seismogram data and calculate receiver functions, and perform 
post-processing for `H-k` stacking and harmonic decomposition. All of them use 
a station database provided as a :class:`~stdb.StDb` dictionary. These scripts are:

- rfpy_calc.py


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

