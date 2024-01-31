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

# -*- coding: utf-8 -*-
from math import ceil
import pickle
import numpy as np
from obspy import Trace, Stream, UTCDateTime
from rfpy import utils


class Meta(object):
    """
    A Meta object contains attributes associated with the metadata
    for a single receiver function analysis.

    Parameters
    ----------
    time : :class:`~obspy.core.UTCDateTime`
        Origin time of earthquake
    dep : float
        Depth of hypocenter (km)
    lon : float
        Longitude coordinate of epicenter
    lat : float
        Latitude coordinate of epicenter
    mag : float
        Magnitude of earthquake
    gac : float
        Great arc circle between station and epicenter (degrees)
    epi_dist : float
        Epicentral distance between station and epicenter (km)
    baz : float
        Back-azimuth - pointing to earthquake from station (degrees)
    az : float
        Azimuth - pointing to station from earthquake (degrees)
    ttime : float
        Predicted arrival time (sec) 
    ph : str
        Phase name 
    slow : float
        Horizontal slowness of phase 
    inc : float
        Incidence angle of phase at surface 
    vp : float
        P-wave velocity at surface (km/s)
    vs : float
        S-wave velocity at surface (km/s)
    align : str
        Alignment of coordinate system for rotation
        ('ZRT', 'LQT', or 'PVH')
    rotated : bool
        Whether or not data have been rotated to ``align``
        coordinate system

    """

    def __init__(self, sta, event, gacmin=30., gacmax=90., phase='P'):

        from obspy.geodetics.base import gps2dist_azimuth as epi
        from obspy.geodetics import kilometer2degrees as k2d
        from obspy.taup import TauPyModel

        # Extract event 4D parameters
        self.time = event.origins[0].time
        self.lon = event.origins[0].longitude
        self.lat = event.origins[0].latitude
        self.dep = event.origins[0].depth

        # Check if depth is valid type
        if self.dep is not None:
            if self.dep > 1000.:
                self.dep = self.dep/1000.
        else:
            self.dep = 10.

        # Magnitude
        self.mag = event.magnitudes[0].mag
        if self.mag is None:
            self.mag = -9.

        # Calculate epicentral distance
        self.epi_dist, self.az, self.baz = epi(
            self.lat, self.lon, sta.latitude, sta.longitude)
        self.epi_dist /= 1000
        self.gac = k2d(self.epi_dist)

        if self.gac > gacmin and self.gac < gacmax:

            # Get travel time info
            tpmodel = TauPyModel(model='iasp91')

            # Get Travel times (Careful: here dep is in meters)
            arrivals = tpmodel.get_travel_times(
                distance_in_degree=self.gac,
                source_depth_in_km=self.dep,
                phase_list=[phase])
            if len(arrivals) > 1:
                print("arrival has many entries: ", len(arrivals))
            elif len(arrivals) == 0:
                print("no arrival found")
                self.accept = False
                return

            arrival = arrivals[0]

            # Attributes from parameters
            self.ttime = arrival.time
            self.slow = arrival.ray_param_sec_degree/111.
            self.inc = arrival.incident_angle
            self.phase = phase
            self.accept = True
        else:
            self.ttime = np.nan
            self.slow = np.nan
            self.inc = np.nan
            self.phase = None
            self.accept = False

        # Defaults for non - station-event geometry attributes
        self.vp = 6.0
        self.vs = 3.5
        self.align = 'ZRT'

        # Attributes that get updated as analysis progresses
        self.rotated = False
        self.snr = np.nan
        self.snrh = np.nan
        self.cc = np.nan


class RFData(object):
    """
    A RFData object contains Class attributes that associate
    station information with a single event (i.e., earthquake) 
    metadata, corresponding raw and rotated seismograms, spectra and 
    receiver functions.

    Note
    ----
    The object is initialized with the ``sta`` field only, and
    other attributes are added to the object as the analysis proceeds.

    Parameters
    ----------
    sta : object
        Object containing station information - from :mod:`~stdb` database.
        When RFData is initialized without sta, functions requiring station info
        coordinate will not be working.
    meta : :class:`~rfpy.rfdata.Meta`
        Object of metadata information for single event (initially set to None)
    data : :class:`~obspy.core.Stream`
        Stream object containing the three-component seismograms (either
        un-rotated or rotated by the method :func:`~rfpy.rfdata.rotate`)
    rf : :class:`~obspy.core.Stream`
        Stream object containing the receiver functions
        (created when executing :func:`~rfpy.rfdata.deconvlove`)
    specs : :class:`~obspy.core.Stream`
        Stream object containing the spectra of the data
        (created when executing :func:`~rfpy.rfdata.calc_spectra`)
    """

    def __init__(self, sta=None):

        # Load example data if initializing empty object
        if sta == 'demo' or sta == 'Demo':
            print("Uploading demo data - station NY.MMPY")
            import os
            sta = pickle.load(
                open(os.path.join(
                    os.path.dirname(__file__),
                    "examples/data", "MMPY.pkl"), 'rb'))['NY.MMPY']

        # Attributes from parameters
        if sta:
            self.sta = sta
        else:
            self.sta = None

        # Initialize meta and data objects as None
        self.meta = None
        self.data = None

    def add_event(self, event, gacmin=30., gacmax=90., phase='P',
                  returned=False):
        """
        Adds event metadata to RFData object, including travel time info 
        of P wave. 

        Parameters
        ----------
        event : :class:`~obspy.core.event`
            Event XML object
        returned : bool
            Whether or not to return the ``accept`` attribute

        Returns
        -------
        accept : bool
            Whether or not the object is accepted for further analysis

        """
        from obspy.geodetics.base import gps2dist_azimuth as epi
        from obspy.geodetics import kilometer2degrees as k2d
        from obspy.taup import TauPyModel
        from obspy.core.event.event import Event

        if event == 'demo' or event == 'Demo':
            from obspy.clients.fdsn import Client
            from obspy.core import UTCDateTime
            client = Client()
            # Get catalogue using deployment start and end
            event = client.get_events(
                starttime=UTCDateTime('2014-06-30T19:00:00'),
                endtime=UTCDateTime('2014-06-30T21:00:00'),
                minmagnitude=6.0,
                maxmagnitude=6.5)[0]
            print(event.short_str())

        if not isinstance(event, Event):
            raise(Exception("Event has incorrect type"))

        # Store as object attributes
        self.meta = Meta(sta=self.sta, event=event,
                         gacmin=gacmin, gacmax=gacmax,
                         phase=phase)

        if returned:
            return self.meta.accept

    def add_data(self, stream, returned=False, new_sr=5.):
        """
        Adds stream as object attribute

        Parameters
        ----------
        stream : :class:`~obspy.core.Stream`
            Stream container for NEZ seismograms
        returned : bool
            Whether or not to return the ``accept`` attribute

        Attributes
        ----------
        zne_data : :class:`~obspy.core.Stream`
            Stream container for NEZ seismograms

        Returns
        -------
        accept : bool
            Whether or not the object is accepted for further analysis

        """

        if not self.meta:
            raise(Exception("No meta data available - aborting"))

        if not self.meta.accept:
            return

        # Load demo data
        if stream == 'demo' or stream == 'Demo':
            import os
            file = open(os.path.join(
                os.path.dirname(__file__),
                "examples/data", "ZNE_Data.pkl"), "rb")
            stream = pickle.load(file)
            print(stream)

        if not isinstance(stream, Stream):
            raise(Exception("Event has incorrect type"))

        try:
            self.data = stream

            if not np.allclose(
                    [tr.stats.npts for tr in stream[1:]], stream[0].stats.npts):
                self.meta.accept = False

            # Filter Traces
            if not stream[0].stats.sampling_rate == new_sr:
                self.data.filter('lowpass', freq=0.5*new_sr,
                                 corners=2, zerophase=True)
                self.data.resample(new_sr, no_filter=True)

        except:
            print("Error: Not all channels are available")
            self.meta.accept = False

        if returned:
            return self.meta.accept

    def download_data(self, client, stdata=[], dtype='SAC', ndval=np.nan,
                       new_sr=5.,dts=120., remove_response=False,
                       local_response_dir='',
                       returned=False, verbose=False):
        """
        Downloads seismograms based on event origin time and
        P phase arrival.

        Parameters
        ----------
        client : :class:`~obspy.client.fdsn.Client`
            Client object
        ndval : float
            Fill in value for missing data
        new_sr : float
            New sampling rate (Hz)
        dts : float
            Time duration (sec)
        stdata : List
            Station list
        remove_response : bool
            Remove instrument response from seismogram and resitute to true ground
            velocity (m/s) using obspy.core.trace.Trace.remove_response()
        returned : bool
            Whether or not to return the ``accept`` attribute
        verbose : bool
            Output diagnostics to screen

        Returns
        -------
        accept : bool
            Whether or not the object is accepted for further analysis

        Attributes
        ----------
        data : :class:`~obspy.core.Stream`
            Stream containing :class:`~obspy.core.Trace` objects

        """

        def _resample(tr, new_sr):

             if tr.stats.sampling_rate > new_sr:
                 # only filter when downsampling
                 tr.filter('lowpass', freq=0.5*new_sr, corners=2, zerophase=True)

             tr.resample(new_sr, no_filter=True)
             return tr


        if self.meta is None:
            raise(Exception("Requires event data as attribute - aborting"))

        if not self.meta.accept:
            return

        # Define start time for request
        tstart = self.meta.time + self.meta.ttime - dts
        tend = self.meta.time + self.meta.ttime + dts

        # Get waveforms
        print("* Requesting Waveforms: ")
        print("*    Startime: " + tstart.strftime("%Y-%m-%d %H:%M:%S"))
        print("*    Endtime:  " + tend.strftime("%Y-%m-%d %H:%M:%S"))

        # Download data
        err, stream = utils.download_data(
            client=client, sta=self.sta, start=tstart, end=tend,
            stdata=stdata, dtype=dtype, ndval=ndval, new_sr=new_sr,
            remove_response=remove_response,
            local_response_dir=local_response_dir,
            verbose=verbose)

        # Store as attributes with traces in dictionary
        try:
            trE = stream.select(component='E')[0]
            trN = stream.select(component='N')[0]
            trZ = stream.select(component='Z')[0]
            self.data = Stream(traces=[trZ, trN, trE])

            # Filter Traces when downsampling
            if self.data[0].stats.sampling_rate > new_sr:
                self.data.filter('lowpass', freq=0.5*new_sr,
                                 corners=2, zerophase=True)
            # Resample
            self.data.resample(new_sr, no_filter=True)

        # If there is no ZNE, perhaps there is Z12?
        except Exception as e:
            if verbose:
                print('* Met exception: ')
                print('* ' + str(e))
                print('* Trying to access Z12 data.')

            try:
                tr1 = stream.select(component='1')[0]
                tr2 = stream.select(component='2')[0]
                trZ = stream.select(component='Z')[0]
                self.data = Stream(traces=[trZ, tr1, tr2])

                # Filter Traces when downsampling
                if self.data[0].stats.sampling_rate > new_sr:
                    self.data.filter('lowpass', freq=0.5*new_sr,
                                     corners=2, zerophase=True)
                # Resample
                self.data.resample(new_sr, no_filter=True)

                # Save Z12 components in case it's necessary for later
                self.dataZ12 = self.data.copy()

                # Rotate from Z12 to ZNE using StDb azcorr attribute
                self.rotate(align='ZNE')

            except Exception as e:
                if verbose:
                    print('* Met exception: ')
                    print('* ' + str(e))
                    print('* Cannot process')
                self.meta.accept = False

        if returned:
            return self.meta.accept

    def rotate(self, vp=None, vs=None, align=None):
        """
        Rotates 3-component seismograms from vertical (Z),
        east (E) and north (N) to longitudinal (L), 
        radial (Q) and tangential (T) components of motion.
        Note that the method 'rotate' from ``obspy.core.stream.Stream``
        is used for the rotation ``'ZNE->ZRT'`` and ``'ZNE->LQT'``.
        Rotation ``'ZNE->PVH'`` is implemented separately here 
        due to different conventions.

        Parameters
        ----------
        vp : float
            P-wave velocity at surface (km/s)
        vs : float
            S-wave velocity at surface (km/s)
        align : str
            Alignment of coordinate system for rotation
            ('ZRT', 'LQT', or 'PVH')

        Returns
        -------
        rotated : bool
            Whether or not the object has been rotated

        Note
        ----
        'PVH' following Equation 18 in:
        Kennett (1991). The removal of free surface interactions from three-component
        seismograms. Geophysical Journal International
        """

        if not self.meta.accept:
            return

        if self.meta.rotated:
            print("Data have been rotated already - continuing")
            return

        # Use default values from meta data if arguments are not specified
        if not align:
            align = self.meta.align

        if align == 'ZNE':
            from obspy.signal.rotate import rotate2zne

            # Copy traces
            trZ = self.data.select(component='Z')[0].copy()
            trN = self.data.select(component='1')[0].copy()
            trE = self.data.select(component='2')[0].copy()

            azim = self.sta.azcorr

            # Try with left handed system
            Z, N, E = rotate2zne(trZ.data, 0., -90., trN.data,
                                 azim, 0., trE.data, azim+90., 0.)

            # Z, N, E = rotate2zne(trZ.data, 0., -90., trN.data,
            #                      azim, 0., trE.data, azim+90., 0.)
            trN.data = N
            trE.data = E

            # Update stats of streams
            trN.stats.channel = trN.stats.channel[:-1] + 'N'
            trE.stats.channel = trE.stats.channel[:-1] + 'E'

            self.data = Stream(traces=[trZ, trN, trE])

        elif align == 'ZRT':
            self.data.rotate('NE->RT',
                             back_azimuth=self.meta.baz)
            self.meta.align = align
            self.meta.rotated = True

        elif align == 'LQT':
            self.data.rotate('ZNE->LQT',
                             back_azimuth=self.meta.baz,
                             inclination=self.meta.inc)
            for tr in self.data:
                if tr.stats.channel.endswith('Q'):
                    tr.data = -tr.data
            self.meta.align = align
            self.meta.rotated = True

        elif align == 'PVH':

            # First rotate to ZRT
            self.data.rotate('NE->RT',
                             back_azimuth=self.meta.baz)

            # Copy traces
            trP = self.data.select(component='Z')[0].copy()
            trV = self.data.select(component='R')[0].copy()
            trH = self.data.select(component='T')[0].copy()

            slow = self.meta.slow
            if not vp:
                vp = self.meta.vp
            if not vs:
                vs = self.meta.vs

            # Vertical slownesses
            # P vertical slowness
            qp = np.sqrt(1./vp/vp-slow*slow)
            # S vertical slowness
            qs = np.sqrt(1./vs/vs-slow*slow)

            # Elements of rotation matrix
            m11 = slow*vs*vs/vp
            m12 = -(1.-2.*vs*vs*slow*slow)/(2.*vp*qp)
            m21 = (1.-2.*vs*vs*slow*slow)/(2.*vs*qs)
            m22 = slow*vs

            # Rotation matrix
            rot = np.array([[-m11, m12], [-m21, m22]])

            # Vector of Radial and Vertical
            r_z = np.array([trV.data, trP.data])

            # Rotation
            vec = np.dot(rot, r_z)

            # Extract P and SV, SH components
            trP.data = vec[0, :]
            trV.data = vec[1, :]
            trH.data = -trH.data/2.

            # Update stats of streams
            trP.stats.channel = trP.stats.channel[:-1] + 'P'
            trV.stats.channel = trV.stats.channel[:-1] + 'V'
            trH.stats.channel = trH.stats.channel[:-1] + 'H'

            # Over-write data attribute
            self.data = Stream(traces=[trP, trV, trH])
            self.meta.align = align
            self.meta.rotated = True

        else:
            raise(Exception("incorrect 'align' argument"))


    def calc_snr(self, dt=30., fmin=0.05, fmax=1.):
        """
        Calculates signal-to-noise ratio on either Z, L or P component

        Parameters
        ----------
        dt : float
            Duration (sec)
        fmin : float
            Minimum frequency corner for SNR filter (Hz)
        fmax : float
            Maximum frequency corner for SNR filter (Hz)

        Attributes
        ----------
        snr : float
            Signal-to-noise ratio on vertical component (dB)
        snrh : float
            Signal-to-noise ratio on radial component (dB)

        """

        if not self.meta.accept:
            return

        if self.meta.snr:  # False if None (for backward compatibility)
            if np.isfinite(self.meta.snr):  # False if nan (None throws error)
                print("SNR already calculated - continuing")
                return

        t1 = self.meta.time + self.meta.ttime

        # SNR for dominant component ('Z', 'L' or 'P')
        comp = self.meta.align[0]

        # Copy trace to signal and noise traces
        trSig = self.data.select(component=comp)[0].copy()
        trNze = self.data.select(component=comp)[0].copy()

        trSig.detrend().taper(max_percentage=0.05)
        trNze.detrend().taper(max_percentage=0.05)

        # Filter between 0.1 and 1.0 (dominant P wave frequencies)
        trSig.filter('bandpass', freqmin=fmin, freqmax=fmax,
                     corners=2, zerophase=True)
        trNze.filter('bandpass', freqmin=fmin, freqmax=fmax,
                     corners=2, zerophase=True)

        # Trim around P-wave arrival
        trSig.trim(t1, t1 + dt)
        trNze.trim(t1 - dt, t1)

        # Calculate root mean square (RMS)
        srms = np.sqrt(np.mean(np.square(trSig.data)))
        nrms = np.sqrt(np.mean(np.square(trNze.data)))

        # Calculate signal/noise ratio in dB
        self.meta.snr = 10*np.log10(srms*srms/nrms/nrms)

        # SNR for radial component ('R', 'Q' or 'V')
        comp = self.meta.align[1]

        # Copy trace to signal and noise traces
        trSig = self.data.select(component=comp)[0].copy()
        trNze = self.data.select(component=comp)[0].copy()

        trSig.detrend().taper(max_percentage=0.05)
        trNze.detrend().taper(max_percentage=0.05)

        # Filter between 0.1 and 1.0 (dominant P wave frequencies)
        trSig.filter('bandpass', freqmin=fmin, freqmax=fmax,
                     corners=2, zerophase=True)
        trNze.filter('bandpass', freqmin=fmin, freqmax=fmax,
                     corners=2, zerophase=True)

        # Trim around P-wave arrival
        trSig.trim(t1, t1 + dt)
        trNze.trim(t1 - dt, t1)

        # Calculate root mean square (RMS)
        srms = np.sqrt(np.mean(np.square(trSig.data)))
        nrms = np.sqrt(np.mean(np.square(trNze.data)))

        # Calculate signal/noise ratio in dB
        self.meta.snrh = 10*np.log10(srms*srms/nrms/nrms)

    def calc_spectra(self, phase='P', vp=None, vs=None,
                   align=None, method='wiener', wavelet='complete',
                   envelope_threshold=0.05, time=5., pre_filt=None,
                   writeto=False, norm=False, length=None):
        """
        Calculate the numerators and denominator of the spectral division for receiver
        function calculation using one component as the source wavelet.
        The source component is always taken as the dominant compressional 
        component, which can be either 'Z', 'L', or 'P'. 

        Parameters
        ----------
        vp : float
            P-wave velocity at surface (km/s)
        vs : float
            S-wave velocity at surface (km/s)
        align : str
            Alignment of coordinate system for rotation
            ('ZRT', 'LQT', or 'PVH')
        method : str
            Method for deconvolution. Options are 'wiener', 'water', 'water2' or 
            'multitaper'
        wavelet : str
            Type of wavelet for deconvolution. Options are 'complete', 'time' or 
            'envelope'
        envelope_threshold : float
            Threshold [0-1] used in ``wavelet='envelope'``.
        time : float
            Window length used in ``wavelet='time'``.
            Minimum window length for ``wavelet='envelope'``.
        pre_filt : list of 2 floats
            Low and High frequency corners of bandpass filter applied
            before deconvolution
        gfilt : float
            Center frequency of Gaussian filter (Hz). 
        wlevel : float
            Water level used in ``method='water'`` and ``method='water2'`` 
        writeto : str or None
            Write wavelets for deconvolution to file.
        length : float or None
            Time length (s) of the numerator. If None, use entire traces.

        Attributes
        ----------
        specs : :class:`~obspy.core.Stream`
            Stream containing the cross spectral quantities

        """

        if not self.meta.accept:
            return

        if not align:
            try:
                align = self.meta.align
            except AttributeError:
                msg = "No alignment found. "
                msg += "Please supply via align keyword."
                raise ValueError(msg)

        def _npow2(x):
            return 1 if x == 0 else 2**(x-1).bit_length()

        def _pad(array, n):
            tmp = np.zeros(n)
            tmp[:array.shape[0]] = array
            return tmp

        def _gauss_filt(dt, nft, f0):
            df = 1./(nft*dt)
            nft21 = int(0.5*nft + 1)
            f = df*np.arange(nft21)
            w = 2.*np.pi*f
            gauss = np.zeros(nft)
            gauss[:nft21] = np.exp(-0.25*(w/f0)**2.)/dt
            gauss[nft21:] = np.flip(gauss[1:nft21-1])
            return gauss

        def _Pwavelet(parent, method='complete', overhang=5.,
                envelope_threshold=0.05, time=5.):

            """
            Select wavelet from the parent function for deconvolution using method.
            parent: obspy.Trace
                wavefrom to extract the wavelet from
            method: str
                'complete' use complete parent signal after P arrival  (current
                    implementation)
                'envelope' use only the part of the parent signal after the
                    P arrival where
                    envelope > envelope_threshold*max(envelope)
                    fall back to 'complete' if condition not reached
                'time' use only this many seconds after P arrival
                    fall back to 'complete' if longer than parent
                'noise' use the time before the P arrival
            overhang: float
                seconds before start and after end of wavelet to be used for
                tapering
            envelope_threshold: float
                fraction of the envelope that defines wavelet (for
                method='envelope')
            time: float
                window (seconds) that defines the wavelet (for method='time')
                minimum time (seconds) of the wavelet (for method='envelope')

            Return:
            left, right: (obspy.UTCDateTime) start and end time of wavelet
            """

            import obspy.signal.filter as osf

            if method not in ['complete', 'envelope', 'time', 'noise']:
                msg = 'Unknow method for wavelet extraction: ' + method
                raise NotImplementedError(msg)

            errflag = False

            if method == 'envelope':
                split = int((self.meta.time + self.meta.ttime +
                    time - parent.stats.starttime ) *
                    parent.stats.sampling_rate)
                env = osf.envelope(parent.data)
                env /= max(env)  # normalize
                env = env[split:]  # only look after P + time
                try:
                    i = np.nonzero(
                            np.diff(
                                np.array(
                                    env > envelope_threshold, dtype=int))==-1)[0][0]
                except IndexError:
                    i = len(parent.data)-1
                dts = i * parent.stats.delta + time
                left = self.meta.time + self.meta.ttime - overhang
                right = self.meta.time + self.meta.ttime + dts + overhang
                if left < parent.stats.starttime or right > parent.stats.endtime:
                    print('Envelope wavelet longer than trace.')
                    print('Falling back to "complete" wavelet')
                    errflag = True

            if method == 'time':
                left = self.meta.time + self.meta.ttime - overhang
                right = self.meta.time + self.meta.ttime + time + overhang
                if left < parent.stats.starttime or right > parent.stats.endtime:
                    print('Time wavelet longer than trace.')
                    print('Falling back to "complete" wavelet')
                    errflag = True

            if method == 'complete' or errflag:
                dts = len(parent.data)*parent.stats.delta/2.
                left = self.meta.time + self.meta.ttime - overhang
                right = self.meta.time + self.meta.ttime + dts - 2*overhang
            
            if method == 'noise':
                dts = len(parent.data)*parent.stats.delta/2.
                left = self.meta.time + self.meta.ttime - dts
                right = self.meta.time + self.meta.ttime - overhang

            return left, right
            
        def _calc_specs(parent, daughter1, daughter2, noise, nn, method):

            # Sampling interval
            dt = parent.stats.delta

            # Wiener or Water level deconvolution
            if method == 'wiener' or method == 'water' or method == 'water2':

                # npad = _npow2(nn*2)
                npad = nn
                freqs = np.fft.fftfreq(npad, d=dt)

                # Fourier transform
                Fp = np.fft.fft(parent.data, n=npad)
                Fd1 = np.fft.fft(daughter1.data, n=npad)
                Fd2 = np.fft.fft(daughter2.data, n=npad)
                Fn = np.fft.fft(noise.data, n=npad)

                # Auto and cross spectra
                Spp = np.real(Fp*np.conjugate(Fp))
                Sd1p = Fd1*np.conjugate(Fp)
                Sd2p = Fd2*np.conjugate(Fp)
                Snn = np.real(Fn*np.conjugate(Fn))

            # Multitaper deconvolution
            elif method == 'multitaper':

                from spectrum import dpss

                npad = nn
                # Re-check length and pad with zeros if necessary
                if not np.allclose(
                        [tr.stats.npts for tr in [parent,
                                                  daughter1,
                                                  daughter2,
                                                  noise]], npad):
                    parent.data = _pad(parent.data, npad)
                    daughter1.data = _pad(daughter1.data, npad)
                    daughter2.data = _pad(daughter2.data, npad)
                    noise.data = _pad(noise.data, npad)

                freqs = np.fft.fftfreq(npad, d=dt)

                NW = 2.5
                Kmax = int(NW*2-2)
                [tapers, eigenvalues] = dpss(npad, NW, Kmax)

                # Get multitaper spectrum of data
                Fp = np.fft.fft(np.multiply(tapers.transpose(),
                                            parent.data))
                Fd1 = np.fft.fft(np.multiply(tapers.transpose(),
                                             daughter1.data))
                Fd2 = np.fft.fft(np.multiply(tapers.transpose(),
                                             daughter2.data))
                Fn = np.fft.fft(np.multiply(tapers.transpose(),
                                            noise.data))

                # Auto and cross spectra
                Spp = np.sum(np.real(Fp*np.conjugate(Fp)), axis=0)
                Sd1p = np.sum(Fd1*np.conjugate(Fp), axis=0)
                Sd2p = np.sum(Fd2*np.conjugate(Fp), axis=0)
                Snn = np.sum(np.real(Fn*np.conjugate(Fn)), axis=0)

            else:
                print("Method not implemented")
                pass

            return Spp, Sd1p, Sd2p, Snn


        if not self.meta.rotated:
            print("Warning: Data have not been rotated yet - rotating now")
            self.rotate(vp=vp, vs=vs, align=align)

        #  v--True if None      v--True if nan, error if None
        if not self.meta.snr or not np.isfinite(self.meta.snr):
            print("Warning: SNR has not been calculated - " +
                  "calculating now using default")
            self.calc_snr()

        if hasattr(self, 'rf'):
            print("Warning: Data have been deconvolved already - passing")
            return

        # Get the name of components (order is critical here)
        cL = self.meta.align[0]
        cQ = self.meta.align[1]
        cT = self.meta.align[2]

        # Define signal and noise
        trL = self.data.select(component=cL)[0].copy()
        trQ = self.data.select(component=cQ)[0].copy()
        trT = self.data.select(component=cT)[0].copy()
        trNl = self.data.select(component=cL)[0].copy()
        trNq = self.data.select(component=cQ)[0].copy()

        trL.stats.channel = 'WV' + self.meta.align[0]
        trQ.stats.channel = 'WV' + self.meta.align[1]
        trT.stats.channel = 'WV' + self.meta.align[2]
        trNl.stats.channel = 'WV' + self.meta.align[0]
        trNq.stats.channel = 'WV' + self.meta.align[1]

        if phase == 'P' or 'PP':

            # Get signal length (i.e., seismogram to deconvolve) from trace length
            over = 5
            dtsqt = len(trL.data) * trL.stats.delta / 2.

            if length and length > dtsqt:
                msg = f"Warning: Requested length {length} longer than data {dtsqt} - passing"
                raise ValueError(msg)

            if length is not None:
                dtsqt = length

            # Traces will be zero-paded to this length (samples)
            nn = int(round((dtsqt + over) * trL.stats.sampling_rate)) + 1

            sig_left, sig_right  = _Pwavelet(trL, method=wavelet,
                envelope_threshold=envelope_threshold, time=time, overhang=over)

            # Trim wavelet
            trL.trim(sig_left, sig_right, nearest_sample=False, pad=True,
                fill_value=0.)

            # Signal window (-5. to dtsqt-10 sec)
            if length is None:
                sig_left, sig_right  = _Pwavelet(trQ, method='complete', overhang=over)
            else:
                sig_left, sig_right  = _Pwavelet(trQ, method='time', time=length, overhang=over)

            # Trim signal traces
            [tr.trim(sig_left, sig_right, nearest_sample=False, 
                pad=True, fill_value=0.) for tr in [trQ, trT]]

            # Noise window (-dtsqt to -5. sec)
            noise_left, noise_right  = _Pwavelet(trQ, method='noise', overhang=over)

            # Trim noise traces
            [tr.trim(noise_left, noise_right, nearest_sample=False, 
                pad=True, fill_value=0.) for tr in [trNl, trNq]]

        # Can't deal with S phases right now ####
        # elif phase == 'S' or 'SKS':

        #     # Get signal length (i.e., seismogram to deconvolve) from trace length
        #     dts = len(trL.data)*trL.stats.delta/2.

        #     # Trim signal traces (-5. to dts-10 sec)
        #     trL.trim(self.meta.time+self.meta.ttime+25.-dts/2.,
        #              self.meta.time+self.meta.ttime+25.)
        #     trQ.trim(self.meta.time+self.meta.ttime+25.-dts/2.,
        #              self.meta.time+self.meta.ttime+25.)
        #     trT.trim(self.meta.time+self.meta.ttime+25.-dts/2.,
        #              self.meta.time+self.meta.ttime+25.)
            
        #     # Trim noise traces (-dts to -5 sec)
        #     trNl.trim(self.meta.time+self.meta.ttime-dts,
        #               self.meta.time+self.meta.ttime-dts/2.)
        #     trNq.trim(self.meta.time+self.meta.ttime-dts,
        #               self.meta.time+self.meta.ttime-dts/2.)

        # Taper traces - only necessary processing after trimming
        # TODO: What does this to the multitaper method
        [tr.taper(max_percentage=0.05, max_length=2.) 
         for tr in [trL, trQ, trT, trNl, trNq]]

        # Pre-filter waveforms before calculating spectra
        if pre_filt:
            [tr.filter('bandpass', freqmin=pre_filt[0], freqmax=pre_filt[1],
                       corners=2, zerophase=True) 
             for tr in [trL, trQ, trT, trNl, trNq]]

        if norm:
            # Standardize all trace data - useful in simdec
            meanN = np.mean(trNl.data)
            stdN = np.std(trNl.data)
            [(tr.data - meanN)/stdN for tr in [trL, trQ, trT, trNl, trNq]]

        if writeto:
            with open(writeto, 'wb') as f:
                pickle.dump(Stream(traces=[trL, trQ, trT]), f)

        # Calculate spectra
        if phase == 'P' or 'PP':
            spp, sd1p, sd2p, snn = _calc_specs(trL, trQ, trT, trNl, nn, method)

        # Can't deal with S phases right now ####
        # elif phase == 'S' or 'SKS':
        #     rfQ, rfL, rfT = _calc_specs(trQ, trL, trT, trNq, nn, method)

        # Copy traces
        sLL = trL.copy()
        sQL = trQ.copy()
        sTL = trT.copy()
        sNN = trNl.copy()
        sLL.data = spp
        sQL.data = sd1p
        sTL.data = sd2p
        sNN.data = snn

        # Update stats of streams
        sLL.stats.channel = 'S' + self.meta.align[0] + self.meta.align[0]
        sQL.stats.channel = 'S' + self.meta.align[1] + self.meta.align[0]
        sTL.stats.channel = 'S' + self.meta.align[2] + self.meta.align[0]
        sNN.stats.channel = 'SNN'

        self.specs = Stream(traces=[sLL, sQL, sTL, sNN])


    def deconvolve(self, align=None, method='wiener',
                   gfilt=None, wlevel=0.01):
        """
        Deconvolves three-component data using one component as the source wavelet.
        The source component is always taken as the dominant compressional 
        component, which can be either 'Z', 'L', or 'P'. 

        Parameters
        ----------
        align : str
            Alignment of coordinate system for rotation
            ('ZRT', 'LQT', or 'PVH')
        method : str
            Method for deconvolution. Options are 'wiener', 'water' (aka 'wlevel'),
            'water2' (aka 'wlevel2'), or 'multitaper'
        gfilt : float
            Center frequency of Gaussian filter (Hz). 
        wlevel : float
            Water level used in ``method='water'``.

        Attributes
        ----------
        rf : :class:`~obspy.core.Stream`
            Stream containing the receiver function traces

        """

        if not hasattr(self, 'specs'):
            msg = "Spectra have not been calculated."
            msg += "Call RFData.calc_spectra() first."
            raise Exception(msg)

        if not align:
            try:
                align = self.meta.align
            except AttributeError:
                msg = "No alignment found. "
                msg += "Please supply via align keyword."
                raise ValueError(msg)

        # Make everything explicit
        SLL = self.specs[0].copy()
        SQL = self.specs[1].copy()
        STL = self.specs[2].copy()
        SNN = self.specs[3].copy()

        dt = SLL.stats.delta
        npad = SLL.stats.npts
        freqs = np.fft.fftfreq(npad, d=dt)

        # Wiener or Water level deconvolution
        if method in ['wiener', 'multitaper']:
            # Denominator (Spp + Snn)
            Sdenom = SLL.data + SNN.data

        elif method == 'water' or method == 'wlevel':
            phi = np.amax(SLL.data)*wlevel
            Sdenom = SLL.data
            Sdenom[Sdenom < phi] = phi

        elif method == 'water2' or method == 'wlevel2':
            # David Gubbins
            # Time Series Analysis and Inverse Therory for Geophysicists
            # "Wiener Filter", Eq. 10.21
            beta = np.amax(SLL.data)*wlevel
            Sdenom = SLL.data + beta

        else:
            raise ValueError('Unknown method: ' + method)

        # Apply Gaussian filter?
        if gfilt:
            gauss = _gauss_filt(dt, npad, gfilt)
            gnorm = np.sum(gauss)*(freqs[1]-freqs[0])*dt
        else:
            gauss = np.ones(npad)
            gnorm = 1.

        # Copy traces
        rfL = SLL.copy()
        rfQ = SQL.copy()
        rfT = STL.copy()

        # Spectral division and inverse transform
        rfL.data = np.fft.ifftshift(np.real(np.fft.ifft(
            gauss*SLL.data/Sdenom))/gnorm)
        rfQ.data = np.fft.ifftshift(np.real(np.fft.ifft(
            gauss*SQL.data/Sdenom))/gnorm)
        rfT.data = np.fft.ifftshift(np.real(np.fft.ifft(
            gauss*STL.data/Sdenom))/gnorm)

        # Update stats of streams
        rfL.stats.channel = 'RF' + align[0]
        rfQ.stats.channel = 'RF' + align[1]
        rfT.stats.channel = 'RF' + align[2]

        self.rf = Stream(traces=[rfL, rfQ, rfT])


    def calc_cc(self):

        if not self.meta.accept:
            return

        if not hasattr(self, 'rf'):
            raise(Exception("Warning: Receiver functions are not available"))

        obs_L = self.data[0].copy()
        obs_Q = self.data[1].copy()
        obs_rfQ = self.rf[1].copy()
        sr = obs_L.stats.sampling_rate

        # Filter using SNR bandpass
        obs_L.detrend().taper(max_percentage=0.05, max_length=2.)
        obs_Q.detrend().taper(max_percentage=0.05, max_length=2.)
        obs_L.filter('bandpass', freqmin=0.05, freqmax=1., corners=2, zerophase=True)
        obs_Q.filter('bandpass', freqmin=0.05, freqmax=1., corners=2, zerophase=True)
        obs_rfQ.filter('bandpass', freqmin=0.05, freqmax=1., corners=2, zerophase=True)

        # Convolve L with rfQ to obtain predicted Q
        pred_Q = obs_Q.copy()
        pred_Q.stats.channel = 'PRR'
        ind1 = int((len(obs_rfQ.data)/2.))
        ind2 = ind1+len(obs_L.data)
        pred_Q.data = np.convolve(
            obs_L.data, obs_rfQ.data, mode='full')[ind1:ind2]

        # trim all traces from 0 to 20. sec following P-wave (fftshift first)
        obs_L.data = np.fft.fftshift(obs_L.data)[0:int(sr*10.)]
        obs_Q.data = np.fft.fftshift(obs_Q.data)[0:int(sr*10.)]
        pred_Q.data = np.fft.fftshift(pred_Q.data)[0:int(sr*10.)]

        # Get cross correlation coefficient between observed and predicted Q
        self.meta.cc = np.corrcoef(obs_Q.data, pred_Q.data)[0][1]


    def to_stream(self, store=None):
        """
        Method to switch from RFData object to Stream object.
        This allows easier manipulation of the receiver functions
        for post-processing.

        """

        if not self.meta.accept:
            return

        def _add_stats(trace, store):
            trace.stats.snr = self.meta.snr
            trace.stats.snrh = self.meta.snrh
            trace.stats.cc = self.meta.cc
            trace.stats.slow = self.meta.slow
            trace.stats.baz = self.meta.baz
            trace.stats.gac = self.meta.gac
            trace.stats.stlo = self.sta.longitude
            trace.stats.stla = self.sta.latitude
            trace.stats.evlo = self.meta.lon
            trace.stats.evla = self.meta.lat
            trace.stats.vp = self.meta.vp
            trace.stats.vs = self.meta.vs
            trace.stats.phase = self.meta.phase
            trace.stats.align = self.meta.align
            if store == 'rf':
                trace.stats.is_rf = True
                nn = self.rf[0].stats.npts
                sr = self.rf[0].stats.sampling_rate
            elif store == 'specs':
                trace.stats.is_specs = True
                nn = self.specs[0].stats.npts
                sr = self.specs[0].stats.sampling_rate
            else:
                nn = self.data[0].stats.npts
                sr = self.data[0].stats.sampling_rate
            trace.stats.taxis = np.fft.fftshift(np.fft.fftfreq(nn, sr)*nn)

            return trace

        if store == None or store == 'data':
            stream = self.data
            for tr in stream:
                tr = _add_stats(tr, store)

        if store == 'rf':
            if not hasattr(self, 'rf'):
                raise(Exception("Warning: Receiver functions are not available"))
            else:
                stream = self.rf
                for tr in stream:
                    tr = _add_stats(tr, store)

        elif store == 'specs':
            if not hasattr(self, 'specs'):
                raise(Exception("Warning: Spectra are not available"))
            else:
                stream = self.specs
                for tr in stream:
                    tr = _add_stats(tr, store)

        return stream

    def save(self, file):
        """
        Saves RFData object to file

        Parameters
        ----------
        file : str
            File name for RFData object

        """

        output = open(file, 'wb')
        pickle.dump(self, output)
        output.close()

