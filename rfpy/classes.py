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
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
"""

:mod:`~rfpy` defines the following base classes:

- :class:`~rfpy.classes.RFData`

The class :class:`~rfpy.classes.RFData` contains attributes
and methods for the analysis of teleseismic receiver functions 
from three-component seismograms. 

"""

# -*- coding: utf-8 -*-
import numpy as np
import rfpy
from obspy.core import Trace
from obspy.geodetics.base import gps2dist_azimuth as epi
from obspy.geodetics import kilometer2degrees as k2d
from obspy.signal.rotate import rotate_ne_rt
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gspec

class Meta(object):
    """
    A Meta object contains attributes associated with the metadata
    for a single receiver function analysis. 
    
    Attributes
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
    align : str
        Alignment of coordinate system ('ZRT', 'LQT', or 'PVH')

    """

    def __init__(self, time, dep, lon, lat, mag, gac, epi_dist, baz, az):

        self.time = time
        self.dep = dep
        self.lon = lon
        self.lat = lat
        self.mag = mag
        self.gac = gac
        self.epi_dist = epi_dist
        self.baz = baz
        self.az = az
        self.ttime = None
        self.ph = None
        self.slow = None
        self.inc = None
        self.rotate = None

class Data(object):
    """
    A Data object contains three-component raw (ZNE) and rotated (ZRT, LQT, PVH) 
    waveforms centered on the arrival time of interest.
    
    Attributes
    ----------

    trN : :class:`~obspy.core.Trace`
        Trace of North component of motion
    trE : :class:`~obspy.core.Trace`
        Trace of East component of motion
    trZ : :class:`~obspy.core.Trace` 
        Trace of Vertical component of motion
    trL : :class:`~obspy.core.Trace`
        Trace of longitudinal/vertical component of motion
    trQ : :class:`~obspy.core.Trace`
        Trace of radial/SV component of motion
    trT : :class:`~obspy.core.Trace`
        Trace of tangential/transverse/SH component of motion
    rfL : :class:`~obspy.core.Trace`
        Trace of longitudinal/vertical component of receiver functions
    rfQ : :class:`~obspy.core.Trace`
        Trace of radial/SV component of receiver functions
    rfT : :class:`~obspy.core.Trace`
        Trace of tangential/transverse/SH component of receiver functions

    """

    def __init__(self, trE, trN, trZ):

        self.trE = trE
        self.trN = trN
        self.trZ = trZ
        self.trL = Trace()
        self.trQ = Trace()
        self.trT = Trace()
        self.rfL = Trace()
        self.rfQ = Trace()
        self.rfT = Trace()

class RFData(object):
    """
    A RFData object contains Class attributes that associate
    station information with a single event (i.e., earthquake) 
    metadata, corresponding raw and rotated seismograms and 
    receiver functions.

    Note
    ----
    The object is initialized with the ``sta`` field only, and
    other attributes are added to the object as the analysis proceeds.

    Attributes
    ----------
    sta : object
        Object containing station information - from :mod:`~stdb` database.
    meta : :class:`~splitpy.classes.Meta`
        Object of metadata information for single event.
    data : :class:`~splitpy.classes.Data`
        Object containing trace data in :class:`~obspy.core.Trace` format
    err : bool
        Whether or not the `get_data_NEZ` successfully retrieved waveforms
    snrq : float
        Signal-to-noise ratio for radial (Q) component
    snrt : float
        Signal-to-noise ratio for tangential (T) component
    maxdt : float
        Max delay time between slow and fast axes in search
    ddt : float
        Sampling distance (in sec) for delay time search
    dphi : float
        Sampling distance (in degrees) for azimuth search 

    """

    def __init__(self, sta, vp=6.0, vs=3.6, align='ZRT'):

        self.sta = sta
        self.vp = vp
        self.vs = vs
        self.align = align

    def add_event(self, event):
        """
        Adds event metadata to RFData object. 

        Parameters
        ----------
        event : :class:`~obspy.core.event`
            Event metadata

        Attributes
        ----------
        meta :
            Object containing metadata information

        """

        time = event.origins[0].time
        dep = event.origins[0].depth
        lon = event.origins[0].longitude
        lat = event.origins[0].latitude

        # Problem with mag
        mag = event.magnitudes[0].mag
        if mag is None: mag = -9.

        # Calculate epicentral distance
        epi_dist, az, baz = epi(lat, lon, self.sta.latitude, self.sta.longitude)
        epi_dist /= 1000
        gac = k2d(epi_dist)

        # Store as object attributes
        self.meta = Meta(time, dep, lon, lat, mag, gac, epi_dist, baz, az)

    def add_phase(self, t):
        """
        Adds phase information for P arrival from taup model

        Parameters
        ----------
        t : :class:`~obspy.taup.TaupPyModel`
            Travel-time table metadata

        Attributes
        ----------
        meta.ttime : float
            Travel time between earthquake and station (sec)
        meta.ph : str
            Phase name ('SKS')
        meta.slow : float
            Horizontal slowness of phase
        meta.inc : float
            Incidence angle of phase at surface

        """

        # Store as attributes
        self.meta.ttime = t.time
        self.meta.ph = t.name
        self.meta.slow = t.ray_param_sec_degree/111.
        self.meta.inc = np.arcsin(self.vp*self.meta.slow)*180./np.pi

    def add_NEZ(self, stream):
        """
        Adds seismograms from available stream

        Parameters
        ----------
        stream : :class:`~obspy.core.Stream`
            Stream container for NEZ seismograms

        Attributes
        ----------
        data : :class:`~rfpy.classes.Data`
            Object containing :class:`obspy.core.Trace` objects

        """

        self.data = Data(stream[0], stream[1], stream[2])

    def add_LQT(self, stream):
        """
        Adds seismograms from available stream

        Parameters
        ----------
        stream : :class:`~obspy.core.Stream`
            Stream container for NEZ seismograms

        Attributes
        ----------
        data : :class:`~rfpy.classes.Data`
            Object containing :class:`~obspy.core.Trace` objects

        """

        self.data.trL = stream[0]
        self.data.trQ = stream[1]
        self.data.trT = stream[2]

    def get_data_NEZ(self, client, dts, stdata, ndval, new_sr):
        """
        Downloads seismograms based on event origin time and
        P phase arrival.

        Parameters
        ----------
        client : :class:`~obspy.client.fdsn.Client`
            Client object
        dts : float
            Time duration (?)
        stdata : :class:`stdb.classes.StDbElement`
            Station metadata
        ndval : float
            Fill in value for missing data

        Attributes
        ----------

        data : :class:`~rfpy.classes.Data`
            Object containing :class:`~obspy.core.Trace` objects

        """

        # Define start and end times for requests
        tstart = self.meta.time + self.meta.ttime - dts
        tend = self.meta.time + self.meta.ttime + dts

        # Get waveforms
        print ("* Requesting Waveforms: ")
        print ("*    Startime: " + tstart.strftime("%Y-%m-%d %H:%M:%S"))
        print ("*    Endtime:  " + tend.strftime("%Y-%m-%d %H:%M:%S"))

        err, trN, trE, trZ = rfpy.utils.get_data_NEZ(client=client, \
            sta=self.sta, start=tstart, \
            end=tend, stdata=stdata, ndval=ndval, new_sr=new_sr)

        # Store as attributes with traces in dictionay
        self.err = err
        self.data = Data(trN, trE, trZ)

    def rotate(self):
        """
        Rotates 3-component seismograms from vertical (Z),
        east (E) and north (N) to longitudinal (L), 
        radial (Q) and tangential (T) components of motion

        Attributes
        ----------

        data : :class:`~splitpy.classes.Data`
            Object containing :class:`~obspy.core.Trace` objects

        """

        if self.align=='ZRT':
            self.data.trL = self.data.trZ.copy()
            self.data.trQ.data, self.data.trT.data = rotate_ne_rt(\
                self.data.trN.data, self.data.trE.data, self.meta.baz)

        elif self.align=='LQT':
            inc = self.meta.inc*np.pi/180.
            baz = self.meta.baz*np.pi/180.

            M = np.zeros((3,3))
            M[0,0] = np.cos(inc)
            M[0,1] = -np.sin(inc) * np.sin(baz)
            M[0,2] = -np.sin(inc) * np.cos(baz)
            M[1,0] = np.sin(inc)
            M[1,1] = np.cos(inc) * np.sin(baz)
            M[1,2] = np.cos(inc) * np.cos(baz)
            M[2,0] = 0.
            M[2,1] = -np.cos(baz)
            M[2,2] = np.sin(baz)

            # Perform 3-D rotation
            LQT = np.dot(np.array(M), np.array(\
                [self.data.trZ.data, self.data.trE.data, self.data.trN.data]))

            # Store into traces and add as new items in attribute dictionary
            self.data.trL = Trace(data=LQT[0], header=self.data.trZ.stats)
            self.data.trQ = Trace(data=LQT[1], header=self.data.trN.stats)
            self.data.trT = Trace(data=LQT[2], header=self.data.trE.stats)

        elif self.align=='PVH':
            # First rotate to ZRT
            self.data.trL = self.data.trZ.copy()
            self.data.trQ.data, self.data.trT.data = rotate_ne_rt(\
                self.data.trN.data, self.data.trE.data, self.meta.baz)

            # Copy traces
            trP = self.data.trL.copy()
            trV = self.data.trQ.copy()
            trH = self.data.trT.copy()

            # Vertical slownesses
            qp = np.sqrt(1/vp/vp-self.meta.slow*self.meta.slow)          # P vertical slowness
            qs = np.sqrt(1/vs/vs-self.meta.slow*self.meta.slow)          # S vertical slowness

            # Elements of rotation matrix
            m11 = self.meta.slow*vs*vs/vp
            m12 = -(1-2*vs*vs*self.meta.slow*self.meta.slow)/(2*vp*qp)
            m21 = (1-2*vs*vs*self.meta.slow*self.meta.slow)/(2*vs*qs)
            m22 = self.meta.slow*vs

            # Rotation matrix
            rot = np.array([[-m11, m12], [-m21, m22]])
            
            # Vector of Radial and Vertical
            r_z = np.array([trR.data,trZ.data])
                        
            # Rotation
            vec = np.dot(rot, r_z)
                        
            # Extract P and SV, SH components - store as attributes
            trP.data = vec[0,:]
            trV.data = vec[1,:]
            trH.data = -trT.data/2.
            self.data.trL = trP
            self.data.trQ = trV
            self.data.trT = trH

        self.meta.align = self.align


    def calc_snrz(self, t1=None, dt=30.):
        """
        Calculates signal-to-noise ration on vertical (L) component

        Parameters
        ----------
        t1 : :class:`~obspy.core.utcdatetime.UTCDateTime`
            Predicted arrival time of phase
        dt : float
            Duration (sec)

        Attributes
        ----------
        snr : float
            Signal-to-noise ratio  (dB)

        """

        if t1 is None:
            t1 = self.meta.time + self.meta.ttime - 5.

        # Copy Z trace to signal and noise traces
        trSig = self.data.trL.copy()
        trNze = self.data.trL.copy()

        # Filter between 0.1 and 1.0 (dominant P wave frequencies)
        trSig.filter('bandpass',freqmin=0.1,freqmax=1.,corners=2,zerophase=True)
        trNze.filter('bandpass',freqmin=0.1,freqmax=1.,corners=2,zerophase=True)

        # Trim twin seconds around P-wave arrival
        trSig.trim(t1, t1 + dt)
        trNze.trim(t1 - dt, t1)

        # Calculate root mean square (RMS)
        srms = np.sqrt(np.mean(np.square(trSig.data)))
        nrms = np.sqrt(np.mean(np.square(trNze.data)))

        # Calculate signal/noise ratio in dB
        self.snrz = 10*np.log10(srms*srms/nrms/nrms)

    def deconvolve(self, twin=30.):
        """
        Calculates the shear-wave splitting parameters based 
        on two alternative method: the Rotation-Correlation (RC)
        method and the Silver-Chan (SC) method. Each set of results
        is stored in a Dictionary as attributes of the split object.

        Parameters
        ----------
        t1 : :class:`~obspy.core.utcdatetime.UTCDateTime`
            Start time of picking window
        t2 : :class:`~obspy.core.utcdatetime.UTCDateTime`
            End time of picking window

        Attributes
        ----------
        RC_res : :class:`~splitpy.classes.Result`
            Object containing results of Rotation-Correlation method
        SC_res : :class:`~splitpy.classes.Result`
            Object containing results of Silver-Chan method

        """

        def _taper(nt, ns):

            tap = np.ones(nt)
            win = np.hanning(2*ns)
            tap[0:ns] = win[0:ns]
            tap[nt-ns:nt] = win[ns:2*ns]

        # Define source and noise
        trL = self.data.trL.copy()
        trQ = self.data.trQ.copy()
        trT = self.data.trT.copy()
        trS = self.data.trL.copy()
        trNl = self.data.trL.copy() # Noise on L
        trNq = self.data.trQ.copy() # Noise on Q

        # trim traces 115 sec in each direction
        trL.trim(self.meta.time+self.meta.ttime-5., self.meta.time+self.meta.ttime+110.)
        trQ.trim(self.meta.time+self.meta.ttime-5., self.meta.time+self.meta.ttime+110.)
        trT.trim(self.meta.time+self.meta.ttime-5., self.meta.time+self.meta.ttime+110.)
        trS.trim(self.meta.time+self.meta.ttime-5., self.meta.time+self.meta.ttime+110.)
        trNl.trim(self.meta.time+self.meta.ttime-120., self.meta.time+self.meta.ttime-5.)
        trNq.trim(self.meta.time+self.meta.ttime-120., self.meta.time+self.meta.ttime-5.)

        # Taper trS 
        window = np.zeros(len(trS.data))
        tap = _taper(int(twin/trS.stats.delta), int(2./trS.stats.delta))
        window[0:int(twin/trS.stats.delta)] = tap
        trS.data *= window

        # Taper other traces
        window = np.zeros(len(trL.data))
        tap = _taper(len(trL.data), int(2./trL.stats.delta))
        window[0:len(trL.data)] = tap
        
        # Some checks
        lwin = len(window)
        if not (lwin==len(trL.data) and lwin==len(trQ.data) and lwin==len(trT.data)\
                and lwin==len(trNl.data) and lwin==len(trNq.data)):
            print('problem with lwin')
            self.data.rfL = Trace()
            self.data.rfQ = Trace()
            self.data.rfT = Trace()

        # Apply taper
        trL.data *=window
        trQ.data *=window
        trT.data *=window
        trNl.data *=window
        trNq.data *=window

        # Fourier transform
        Fl = np.fft.fft(trL.data)
        Fq = np.fft.fft(trQ.data)
        Ft = np.fft.fft(trT.data)
        Fs = np.fft.fft(trS.data)
        Fnl = np.fft.fft(trNl.data)
        Fnq = np.fft.fft(trNq.data)
        
        # Auto and cross spectra
        Sl = Fl*np.conjugate(Fs)
        Sq = Fq*np.conjugate(Fs)
        St = Ft*np.conjugate(Fs)
        Ss = Fs*np.conjugate(Fs)
        Snl = Fnl*np.conjugate(Fnl)
        Snq = Fnq*np.conjugate(Fnq)
        Snlq = Fnq*np.conjugate(Fnl)

        # Denominator
        Sdenom = 0.25*(Snl+Snq)+0.5*np.abs(Snlq)
       
        # Copy traces
        rfL = trL.copy()
        rfQ = trQ.copy()
        rfT = trT.copy()
        
        # Spectral division and inverse transform
        rfL.data = np.real(np.fft.ifft(Sl/(Ss+Sdenom)))
        rfQ.data = np.real(np.fft.ifft(Sq/(Ss+Sdenom))/np.amax(rfL.data))
        rfT.data = np.real(np.fft.ifft(St/(Ss+Sdenom))/np.amax(rfL.data))

        self.data.rfL = rfL
        self.data.rfQ = rfQ
        self.data.rfT = rfT

        
    def save(self, file):
        """
        Saves Split object to file

        Parameters
        ----------
        file : str
            File name for split object

        """
        
        import pickle
        output = open(file, 'wb')
        pickle.dump(self, output)
        output.close()

