import obspy
import matplotlib.pyplot as plt
import numpy as np
class STpair(object):
    def __init__(self, st1, st2):
        # Parameters for first iteration
        if isinstance(st1, obspy.Stream):
            self.st1    = st1
        elif isinstance(st1, obspy.Trace):
            self.st1    = obspy.Stream()
            self.st1.append(st1)
        else:
            raise ValueError('st1 should be obspy stream or trace!')
        if isinstance(st2, obspy.Stream):
            self.st2    = st2
        elif isinstance(st2, obspy.Trace):
            self.st2    = obspy.Stream()
            self.st2.append(st2)
        else:
            raise ValueError('st2 should be obspy stream or trace!')
        return
    
    def compare(self, fmin=None, fmax=None, vmin=None, vmax=None):
        chanLst = ['BXE', 'BXN', 'BXZ']
        fig = plt.figure()
        i=1
        for chan in chanLst:
            try:
                tr1 = self.st1.select(channel=chan)[0]
                tr2 = self.st2.select(channel=chan)[0]
            except:
                continue
            if fmin !=None and fmax !=None:
                tr1.filter('bandpass', freqmin=fmin, freqmax=fmax)
                tr2.filter('bandpass', freqmin=fmin, freqmax=fmax)
            elif fmin == None and fmax !=None:
                tr1.filter('highpass', freq=fmax)
                tr2.filter('highpass', freq=fmax)
            elif fmin != None and fmax == None:
                tr1.filter('lowpass', freq=fmin)
                tr2.filter('lowpass', freq=fmin)
            time = np.arange(tr1.stats.npts, dtype=float)*tr1.stats.delta
            ax1 = fig.add_subplot(3, 1, i)
            plt.title(chan)
            ax1.plot(time, tr1.data, 'k-', lw=2)
            ax1.plot(time, tr2.data, 'r--', lw=2)
            ax1.plot(time, (tr1.data - tr2.data)*2., 'g-', lw=2)
            tmin    = -1e10; tmax = 1e10
            if vmin!=None: tmax = tr1.stats.sac.dist/vmin; 
            if vmax!=None: tmin = tr1.stats.sac.dist/vmax;
            print tmin, tmax, (tr1.stats.npts-1.)*tr1.stats.delta
            tmin = max(tmin, 0.)
            tmax = min(tmax, (tr1.stats.npts-1.)*tr1.stats.delta)
            # print tmax
            plt.xlim(tmin, tmax)
            i+=1
        plt.show()
        return
            
            
            