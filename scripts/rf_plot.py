#!/usr/bin/env python 
# encoding: utf-8
'''
PROGRAM rf_plot.py

Plots single station P-receiver functions.

Options to plot receiver functions by back-azimuth vs time/depth,
or slowness vs time. 

See module rf_bins.py for details.

'''

# Import modules and functions
import os.path
import fnmatch
import matplotlib.pyplot as plt
from sys import exit
from rf import rf_bins as rfb
from rf import rf_wiggle as rfw
from rf.io_utils import load_db,get_RFs



################################################################################
# Main plotting function
################################################################################
def plot_rfs(db, sta_key, snr=5, f1=0.05, f2=0.5, typ='baztime', rfpath='RF_DATA', save=False, figpath=None, figname=None):
    # Extract station information from dictionary
    sta = db[sta_key]

    if (not typ=="baztime") and (not typ=="slowtime") and (not typ=="bazdepth"):
        print("Error: Incorrect Type")
        exit()

    #--
    if not os.path.isdir('RF_PLOTS'):
        os.makedirs('RF_PLOTS')

    #-- Construct Figure name
    if figname==None:
        savename="."+sta.network+"_SNR"+str(snr)+"_"+str(f1)+"-"+str(f2)+"_"+typ
    else:
        savename="."+sta.network+"_"+figname

    # Get data from database of RFs
    print("  Loading RFs")
    rfV, rfH = get_RFs(os.path.join(rfpath,sta.network+"."+sta.station), sta.network, sta.station, snr)

    # Filter streams
    rfV.filter('bandpass',freqmin=f1, freqmax=f2, corners=4, zerophase=True)
    rfH.filter('bandpass',freqmin=f1, freqmax=f2, corners=4, zerophase=True)

    print("  ---------------------------")
    print("    Number RFs: ",len(rfV))
    print("    Type: "+typ)
    print("    Threshold: "+str(snr))
    print("    Frequency: "+str(f1)+" - "+str(f2)+" Hz")
    print("    Save: ",save)
    if save:
        print("      Name: "+sta.station+"."+savename)
    print("  ---------------------------")

    if len(rfV)==0:
        print("  Skipping...")
        print(" ")
        return

    if typ=='baztime':
            
        #----------------------
        # Plotting baz bins as a function of time
        
        # Get back-azimuth bins
        rfVbin, rfHbin = rfb.rf_bins(rfV, rfH, 'baz', 36+1)
        
        # Stack all traces
        rfVstack, rfHstack = rfb.rf_stack_all(rfV, rfH)
        
        # Plot as wiggles
        rfw.rf_wiggle_bins(rfVbin, rfHbin, rfVstack, rfHstack,  sta.station, 'baz', 30, 'time', save=save, title=savename)

    elif typ=='slowtime':
    
        #----------------------
        # Plotting slow bins as a function of time
        
        # Get slowness bins
        rfVbin, rfHbin = rfb.rf_bins(rfV, rfH, 'slow', 20+1)
        
        # Stack all traces
        rfVstack, rfHstack = rfb.rf_stack_all(rfV, rfH)
        
        # Plot as wiggles
        rfw.rf_wiggle_bins(rfVbin, rfHbin, rfVstack, rfHstack, sta.station, 'slow', 30, 'time', save=save, title=savename)

    elif typ=='bazdepth':

        #----------------------
        # Plotting baz bins as a function of depth
        
        # Migrate to depth
        rfVdep, rfHdep = rfb.rf_migrate(rfV, rfH, 100+1)
        
        # Get back-azimuth bins
        rfVbin, rfHbin = rfb.rf_bins(rfVdep, rfHdep, 'baz', 36+1)
        
        # Stack all traces
        rfVstack, rfHstack = rfb.rf_stack_all(rfVbin, rfHbin)
        
        # Plot as wiggles
        rfw.rf_wiggle_bins(rfVbin, rfHbin, rfVstack, rfHstack, sta.station, 'baz', 100, 'depth', save=save, title=savename)
################################################################################



################################################################################
#-- Parse Commandline Options
################################################################################
def get_options():
    from optparse import OptionParser,OptionGroup
    parser = OptionParser(usage="Usage: %prog [options] <Station Pickle File>",description="RF Plotting Wrapper")
    parser.add_option("-q","-Q","--quiet",action="store_false",dest="verb",default=True,help="Disable output (default Verbose)")

    statGroup=OptionGroup(parser,"Stations")
    statGroup.add_option("--Stations",action="store",type="string",dest="sta_list",default="",help="Specify a comma-separated list of stations to process. By default, all stations in database are processed.")

    foldGroup=OptionGroup(parser,"Folder Options")
    foldGroup.add_option("--RFDir",action="store",type="string",dest="rfpath",default="RF_DATA",help="Specify the RF Directory (default RF_DATA).")
    foldGroup.add_option("--PlotDIR",action="store",type="string",dest="pltdir",default="RF_Plots",help="Specify the directory to store plot files (default RF_Plots)")

    threshGroup=OptionGroup(parser,"Threshold Options")
    threshGroup.add_option("-S","--SNR",action="store",type="float",dest="SNR",default=5.,help="Specify minimum threshold for SNR (default 5.0)")
    threshGroup.add_option("--f1","--minFreq",action="store",type="float",dest="minFreq",default=0.05,help="Specify the minimum filtering frequency (default 0.05Hz)")
    threshGroup.add_option("--f2","--maxFreq",action="store",type="float",dest="maxFreq",default=0.5,help="Specify the maximum filtering frequency (default 0.5Hz)")

    plotGroup=OptionGroup(parser,"Plotting Options")
    plotGroup.add_option("--type",action="store",type="string",dest="type",default='baztime',help="Specify RF plot type. Options are baztime, slowtime, bazdepth (default baztime).")
    plotGroup.add_option("--save",action="store_true",dest="save",default=False,help="Save Figures?")
    plotGroup.add_option("--FigName",action="store",type="string",dest="figname",default=None,help="Specify Figure Name suffix. Prefix is always Net.Sta. Leave blank to include type, filter, and threshold in filename (ex: TA.EPYK.baztime_0.05-0.5_SNR5.png)")

    #-- Add Parsing Groups
    parser.add_option_group(statGroup)
    parser.add_option_group(foldGroup)
    parser.add_option_group(threshGroup)
    parser.add_option_group(plotGroup)

    #-- Parse Inputs
    (opts, args) = parser.parse_args()

    #-- Check for Pickle
    if len(args)!=1: parser.error("Need Station Database pickle file")
    if not os.path.exists(args[0]):
        parser.error("Error: Station Pickle File "+args[0]+"doesn't exist")
    dbpkl=args[0]

    #-- Parser Sta list
    if len(opts.sta_list)==0:
        opts.sta_list=list()
    else:
        opts.sta_list=opts.sta_list.split(',')

    return (opts,dbpkl)
################################################################################




if __name__=='__main__':

    (opts,dbfile)=get_options()

    # Load it
    stationdb = load_db(dbfile)

    # Select station to process
    if len(opts.sta_list)==0:
        sta_keys=list(stationdb.keys())
    else:
        sta_keys = opts.sta_list
    sta_keys.sort()

    dbkeys=list(stationdb.keys())
    dbkeys.sort()


    # Loop over station keys and calculate RFs
    ik=1
    for skey in sta_keys:
            pkeys=[s for s in dbkeys if skey in s]
            for key in pkeys:
                print(" ")
                print("==============================")
                print("{0:.0f}) {1:5s}".format(ik,key))
                try:
                    if key not in stationdb:
                        print("  Station Not in Database...Skip")
                        continue
                except:
                    if key not in stationdb:
                        print("  Station Not in Database...Skip")
                        continue

                plot_rfs(stationdb, key, opts.SNR, opts.minFreq, opts.maxFreq, opts.type.lower(), opts.rfpath, opts.save, opts.pltdir, opts.figname)
                ik+=1
                print("==============================")
                print(" ")
###############################
