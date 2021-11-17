def RFStreams(object):

    def _init__(self, rfdata=None):

        self.rfdata = []

        if isinstance(rfdata, RFData):
            rfdata = [rfdata]
        if rfdata:
            self.rfdata.extend(rfdata)

    def __add__(self, other):
        """
        Add two `:class:`~plateflex.classes.Grid` objects or a 
        :class:`~plateflex.classes.Project` object with a single grid.

        """
        if isinstance(other, RFData):
            other = RFStreams([other])
        if not isinstance(other, RFStreams):
            raise TypeError
        rfdata = self.rfdata + other.rfdata
        return self.__class__(rfdata=rfdata)

    def __iter__(self):
        """
        Return a robust iterator for :class:`~plateflex.classes.Grid` 
        objects

        """
        return list(self.rfdata).__iter__()

    def append(self, rfdata):
        """
        Append a single :class:`~plateflex.classes.Grid` object to the 
        current `:class:`~plateflex.classes.Project` object.

        :type grid: :class:`~plateflex.classes.Grid`
        :param grid: object to append to project

        .. rubric:: Example

        >>> import numpy as np
        >>> from plateflex import Grid, Project
        >>> nn = 200; dd = 10.
        >>> grid = Grid(np.random.randn(nn, nn), dd, dd)
        >>> project = Project()
        >>> project.append(grid)
        """

        if isinstance(rfdata, RFData):
            self.rfdata.append(rfdata)
        else:
            msg = 'Append only supports a single RFData object as an argument.'
            raise TypeError(msg)

        return self

    def extend(self, rfdata_list):
        """
        Extend the current Project object with a list of Grid objects.

        :param trace_list: list of :class:`~plateflex.classes.Grid` objects or
            :class:`~plateflex.classes.Project`.

        .. rubric:: Example

        >>> import numpy as np
        >>> from plateflex import Grid, Project
        >>> nn = 200; dd = 10.
        >>> grid1 = Grid(np.random.randn(nn, nn), dd, dd)
        >>> grid2 = Grid(np.random.randn(nn, nn), dd, dd)
        >>> project = Project()
        >>> project.extend(grids=[grid1, grid2])

        """
        if isinstance(rfdata_list, list):
            for _i in rfdata_list:
                # Make sure each item in the list is a Grid object.
                if not isinstance(_i, RFData):
                    msg = 'Extend only accepts a list of RFData objects.'
                    raise TypeError(msg)
            self.rfdata.extend(rfdata_list)
        elif isinstance(rfdata_list, RFStreams):
            self.rfdata.extend(rfdata_list.rfdata)
        else:
            msg = 'Extend only supports a list of RFData objects as argument.'
            raise TypeError(msg)
        return self


    def calc_spectra(self, wavelet='complete',):

        def _npow2(x):
            return 1 if x == 0 else 2**(x-1).bit_length()

        def _pad(array, n):
            tmp = np.zeros(n)
            tmp[:array.shape[0]] = array
            return tmp

        def _Pwavelet(parent, method='complete', overhang=5,
                envelope_threshold=0.05, time=5):

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
                'time' use only this many seonds after P arrival`
                    fall back to 'complete' if longer than parent
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

        def _specs(trL, trQ, trT, trN, nn, method):

            # Get length, zero padding parameters and frequencies
            dt = parent.stats.delta

            # npad = _npow2(nn*2)
            npad = nn
            freqs = np.fft.fftfreq(npad, d=dt)

            # Fourier transform
            Fl = np.fft.fft(trL.data, n=npad)
            Fq = np.fft.fft(trQ.data, n=npad)
            Ft = np.fft.fft(trT.data, n=npad)
            Fn = np.fft.fft(trN.data, n=npad)

            # Copy traces
            Sll = trL.copy()
            Sql = trQ.copy()
            Stl = trT.copy()
            Snn = trN.copy()

            # Auto and cross spectra
            Sll.data = np.real(Fl*np.conjugate(Fl))
            Sql.data = Fq*np.conjugate(Fl)
            Stl.data = Ft*np.conjugate(Fl)
            Snn.data = np.real(Fn*np.conjugate(Fn))

            return Sll, Sql, Stl, Snn


        for rfdata in self.rfdata:

            # Get the name of components (order is critical here)
            cL = stream.meta.align[0]
            cQ = stream.meta.align[1]
            cT = stream.meta.align[2]

            # Define signal and noise
            trL = stream.data.select(component=cL)[0].copy()
            trQ = stream.data.select(component=cQ)[0].copy()
            trT = stream.data.select(component=cT)[0].copy()
            trN = stream.data.select(component=cL)[0].copy()

            # Standardize all trace data
            meanN = np.mean(trN.data)
            stdN = np.std(trN.data)
            [tr.data = (tr.data - meanN)/stdN for tr in [trL, trQ, trT, trN]]

            trL.stats.channel = 'WV' + stream.meta.align[0]
            trQ.stats.channel = 'WV' + stream.meta.align[1]
            trT.stats.channel = 'WV' + stream.meta.align[2]
            trN.stats.channel = 'WV' + stream.meta.align[0]

            if phase == 'P' or 'PP':

                # Get signal length (i.e., seismogram to deconvolve) from trace length
                over = 5
                dtsqt = len(trL.data)*trL.stats.delta/2.

                # Traces will be zero-paded to this length (samples)
                nn = int(round((dtsqt+over)*trL.stats.sampling_rate)) + 1

                sig_left, sig_right  = _Pwavelet(trL, method=wavelet,
                    envelope_threshold=envelope_threshold, time=time, overhang=over)

                # Trim wavelet
                trL.trim(sig_left, sig_right, nearest_sample=False, pad=True,
                    fill_value=0.)

                # Signal window (-5. to dtsqt-10 sec)
                sig_left, sig_right  = _Pwavelet(trQ, method='complete', overhang=over)

                # Trim signal traces
                [tr.trim(sig_left, sig_right, nearest_sample=False, 
                    pad=True, fill_value=0.) for tr in [trQ, trT]]

                # Noise window (-dtsqt to -5. sec)
                noise_left, noise_right  = _Pwavelet(trQ, method='noise', overhang=over)

                # Trim noise trace
                trN.trim(noise_left, noise_right, nearest_sample=False, 
                    pad=True, fill_value=0.)

            elif phase == 'S' or 'SKS':

                # Get signal length (i.e., seismogram to deconvolve) from trace length
                dts = len(trL.data)*trL.stats.delta/2.

                # Trim signal traces (-5. to dts-10 sec)
                trL.trim(self.meta.time+self.meta.ttime+25.-dts/2.,
                         self.meta.time+self.meta.ttime+25.)
                trQ.trim(self.meta.time+self.meta.ttime+25.-dts/2.,
                         self.meta.time+self.meta.ttime+25.)
                trT.trim(self.meta.time+self.meta.ttime+25.-dts/2.,
                         self.meta.time+self.meta.ttime+25.)
                
                # Trim noise traces (-dts to -5 sec)
                trN.trim(self.meta.time+self.meta.ttime-dts,
                          self.meta.time+self.meta.ttime-dts/2.)

            # Taper traces - only necessary processing after trimming
            # TODO: What does this to the multitaper method
            [tr.taper(max_percentage=0.05, max_length=2.) 
             for tr in [trL, trQ, trT, trN]]

            # Spectra
            if phase == 'P' or 'PP':
                sll, sql, stl, snn = _specs(trL, trQ, trT, trN, nn)

            rfdata.specs = Stream(traces=[sll, sql, stl, snn])


        # elif phase == 'S' or 'SKS':
        #     rfQ, rfL, rfT = _decon(trQ, trL, trT, trNq, nn, method)


    def bin_baz_slow_specs(self, nbaz=72, nslow=1):

            """ 
            Function to stack auto- and cross-component spectra into back-azimuth 
            and slowness bins.

            Parameters
            ----------
            nbaz : int
                Number of back-azimuth samples in bins
            nslow : int
                Number of slowness samples in bins
            include_empty : bool
                Return empty bins as null traces (default omits them)

            Returns
            -------
            stack : :class:`~obspy.core.Stream`
                Stream containing one or two stacked traces,
                depending on the number of input streams

            """

            # Define back-azimuth and slowness bins
            baz_bins = np.linspace(0, 360, nbaz)
            slow_bins = np.linspace(0.04, 0.08, nslow)

            # Extract baz and slowness
            baz = [stream[i].stats.baz for i in range(len(stream))]
            slow = [stream[i].stats.slow for i in range(len(stream))]

            # Digitize baz and slowness
            ibaz = np.digitize(baz, baz_bins)
            islow = np.digitize(slow, slow_bins)

            final_stream = []

            for stream in self.streams:
                try:
                    # Define empty streams
                    binned_stream = Stream()

                    # Loop through baz_bins
                    for i in range(nbaz):
                        for j in range(nslow):

                            nbin = 0
                            array = np.zeros(len(stream[0].specs.data))

                            # Loop all traces
                            for k, tr in enumerate(stream):

                                # If index of baz_bins is equal to ibaz
                                if i == ibaz[k] and j == islow[k]:

                                    nbin += 1
                                    array += tr.data
                                    
                                    continue

                            if nbin > 0 or include_empty:

                                # Average and update stats
                                array /= nbin
                                weight = np.real(abs(weight/nbin))

                                trace = Trace(header=stream[0].stats)
                                trace.stats.baz = baz_bins[i]
                                trace.stats.slow = slow_bins[j]
                                trace.stats.nbin = nbin

                                if not pws:
                                    weight = np.ones(len(stream[0].data))
                                trace.data = weight*array
                                binned_stream.append(trace)

                    final_stream.append(binned_stream)

                except:
                    continue

            return final_stream


    def deconvolve(self, phase='P', vp=None, vs=None,
                   align=None, method='wiener', 
                   envelope_threshold=0.05, time=5, pre_filt=None,
                   gfilt=None, wlevel=0.01, writeto=None):
        """
        Deconvolves three-component data using one component as the source wavelet.
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
            Method for deconvolution. Options are 'wiener', 'water' or 
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
            Water level used in ``method='water'``.
        writeto : str or None
            Write wavelets for deconvolution to file.

        Attributes
        ----------
        rf : :class:`~obspy.core.Stream`
            Stream containing the receiver function traces

        """

        if not self.meta.accept:
            return


        def _gauss_filt(dt, nft, f0):
            df = 1./(nft*dt)
            nft21 = int(0.5*nft + 1)
            f = df*np.arange(nft21)
            w = 2.*np.pi*f
            gauss = np.zeros(nft)
            gauss[:nft21] = np.exp(-0.25*(w/f0)**2.)/dt
            gauss[nft21:] = np.flip(gauss[1:nft21-1])
            return gauss

            
        def _decon(parent, daughter1, daughter2, noise, nn, method):


                # Final processing depends on method
                if method == 'wiener':
                    Sdenom = Spp + Snn
                elif method == 'water':
                    phi = np.amax(Spp)*wlevel
                    Sdenom = Spp
                    Sdenom[Sdenom < phi] = phi

            # Apply Gaussian filter?
            if gfilt:
                gauss = _gauss_filt(dt, npad, gfilt)
                gnorm = np.sum(gauss)*(freqs[1]-freqs[0])*dt
            else:
                gauss = np.ones(npad)
                gnorm = 1.

            # Is this correct?
            #parent de-noised = np.fft.ifftshift(np.real(np.fft.ifft(gauss*Sdenom))/gnorm)
            #daughter1 de-noised = np.fft.ifftshift(np.real(np.fft.ifft( gauss*Sd1p))/gnorm)
            #daughter2 de-noised = np.fft.ifftshift(np.real(np.fft.ifft( gauss*Sd2p))/gnorm)

            # Copy traces
            rfp = parent.copy()
            rfd1 = daughter1.copy()
            rfd2 = daughter2.copy()

            # Spectral division and inverse transform
            rfp.data = np.fft.ifftshift(np.real(np.fft.ifft(
                gauss*Spp/Sdenom))/gnorm)
            rfd1.data = np.fft.ifftshift(np.real(np.fft.ifft(
                gauss*Sd1p/Sdenom))/gnorm)
            rfd2.data = np.fft.ifftshift(np.real(np.fft.ifft(
                gauss*Sd2p/Sdenom))/gnorm)

            return rfp, rfd1, rfd2

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


        # Pre-filter waveforms before deconvolution
        if pre_filt:
            [tr.filter('bandpass', freqmin=pre_filt[0], freqmax=pre_filt[1],
                       corners=2, zerophase=True) 
             for tr in [trL, trQ, trT, trNl, trNq]]

        if writeto:
            with open(writeto, 'wb') as f:
                pickle.dump(Stream(traces=[trL, trQ, trT]), f)

        # Deconvolve
        if phase == 'P' or 'PP':
            rfL, rfQ, rfT = _decon(trL, trQ, trT, trNl, nn, method)

        elif phase == 'S' or 'SKS':
            rfQ, rfL, rfT = _decon(trQ, trL, trT, trNq, nn, method)

        # Update stats of streams
        rfL.stats.channel = 'RF' + self.meta.align[0]
        rfQ.stats.channel = 'RF' + self.meta.align[1]
        rfT.stats.channel = 'RF' + self.meta.align[2]

        self.rf = Stream(traces=[rfL, rfQ, rfT])


    def to_stream(self):
        """
        Method to switch from RFData object to Stream object.
        This allows easier manipulation of the receiver functions
        for post-processing.

        """

        if not self.meta.accept:
            return

        def _add_specstats(trace):
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
            trace.stats.is_rf = True
            nn = self.rf[0].stats.npts
            sr = self.rf[0].stats.sampling_rate
            trace.stats.taxis = np.fft.fftshift(np.fft.fftfreq(nn, sr)*nn)
            return trace

        if not hasattr(self, 'rf'):
            raise(Exception("Warning: Receiver functions are not available"))

        stream = self.rf
        for tr in stream:
            tr = _add_rfstats(tr)

        return stream
