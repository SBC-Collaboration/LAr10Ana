'''Shared t0 / scintillation-pulse helpers used by more than one tab.

The Slow DAQ and Piezo tabs both draw the same set of t0 references and the same
per-CAEN-trigger pulse overlay, so the colour convention and the recon lookups live
here rather than in either tab. These cannot live in ped.py: ped.py imports the tabs,
so a tab importing back from ped would be circular.

All three t0 values exposed through get_t0_ms are in the same frame: milliseconds
relative to the trigger. See ana/ScintT0.py, where
    pressureT0_in_corrected_time = latch_time_corrected + pt0_rel_to_trig
establishes that latch_time_corrected is the acoustic trigger expressed in the CAEN
corrected clock, so subtracting the latch puts CAEN times in that same frame.
'''

import os

import numpy as np

from GetEvent import GetEvent
from sbcbinaryformat import Streamer

# Colour per t0 type, shared so a given t0 looks the same on every tab.
T0_COLORS = {'pressure': 'r', 'bubble': 'g', 'scint': 'b'}

# CAEN scintillation digitizer clock, ticks/s (matches ana/ScintT0.py)
SCINT_SAMPLE_RATE = 125e6

# Approximate photo-electron scale for the channel-summed hit_area.
#
# Treat the resulting numbers as "order of a p.e.", not as calibrated phd:
#  - hit_area is NOT in mV. ana/SiPMPulses.py builds it as hit_trace_V.sum(axis=0)
#    over a baseline-subtracted trace, but that trace is in ADC counts (the recon
#    file stores baseline ~2116 and rms ~9 alongside it), so hit_area is really
#    ADC counts x samples. The 160 was supplied as a rough mV/phd figure, so the
#    division is an empirical rescale rather than a unit conversion.
#  - Summing over channels before dividing assumes the SiPM channel gains are
#    roughly matched, which is only approximately true.
# It is only ever used for the log colour/height scale of the pulse overlay, where
# a constant factor cancels out, so an imprecise value costs nothing visually.
SIPM_AREA_PER_PHD = 160.0


def unwrap_caen_timestamp(ts, max_ts=2**31):
    # The CAEN trigger clock is 31 bits and rolls over mid-event. Mirrors
    # ana/ScintT0.py so pulse times match the t0s computed there; kept local
    # to avoid importing that module's sklearn dependency into the GUI.
    ts = np.asarray(ts, dtype=np.int64)
    rollovers = np.diff(ts, axis=-1, prepend=0) < 0
    return ts + np.cumsum(rollovers, axis=-1) * max_ts


class ScintPulses:
    '''Mixin on Application: t0 lookup plus a per-run cache of CAEN pulse times/areas.

    Mixed in once, so every tab that draws the overlay shares one cache and one set
    of disk reads.
    '''

    def __init__(self):
        self._scint_pulse_run = None     # run whose recon file is cached
        self._scint_pulse_all = None     # recon scintillation.sbc for the whole run
        self._scint_pulse_cache = {}     # ev -> (t_ms, area) or None

    def get_t0_ms(self, key, ev):
        # Trigger-relative t0 [ms] for one of T0_COLORS, or None. May be NaN: the
        # underlying fits fail often enough that callers must check np.isfinite.
        if key == 'pressure':
            return self.pressure_t0.get(ev)
        elif key == 'scint':
            return self.scint_t0.get(ev)
        elif key == 'bubble':
            return self.bubble_t0_ms(ev)

    def load_scint_pulses(self, ev):
        # Per-CAEN-trigger (time, summed pulse area) for one event, or None.
        # Time comes from the raw file's TriggerTimeTag and area from the recon
        # file; recon carries no timestamps and raw carries no pulse areas, so
        # both are needed. Recon row N is raw trigger N.
        if ev in self._scint_pulse_cache:
            return self._scint_pulse_cache[ev]

        result = None
        try:
            if self._scint_pulse_run != self.run or self._scint_pulse_all is None:
                self._scint_pulse_all = None
                self._scint_pulse_cache = {}  # event numbers repeat across runs
                self._scint_pulse_run = self.run
                path = self._find_recon('scintillation.sbc', self.run)
                if path is not None:
                    self._scint_pulse_all = Streamer(path).data

            recon = self._scint_pulse_all
            if recon is not None:
                if 'ev' in recon.dtype.names:
                    recon = recon[recon['ev'] == ev]
                if len(recon):
                    # NaN marks a channel with no hit, so nansum gives 0 for a
                    # trigger with nothing on any channel.
                    area = np.nansum(np.asarray(recon['hit_area']), axis=1)

                    # TriggerTimeTag is read eagerly even in lazy mode, so this
                    # never pulls waveforms off disk.
                    run_path = os.path.join(self.raw_directory, self.run)
                    event = GetEvent(run_path, ev, 'run_control', 'scintillation',
                                     strictMode=False, lazy_load_scintillation=True)
                    scint = event['scintillation']
                    if scint.get('loaded'):
                        tt = unwrap_caen_timestamp(scint['TriggerTimeTag'])
                        t_ms = (tt - tt[0]) / SCINT_SAMPLE_RATE * 1000.0
                        n = min(len(t_ms), len(area))
                        if n:
                            result = (np.asarray(t_ms[:n], dtype=float),
                                      np.asarray(area[:n], dtype=float))
        except Exception as e:
            self.logger.error(
                'failed to load scint pulses for ev {}: {}'.format(ev, e))
            result = None

        self._scint_pulse_cache[ev] = result
        return result

    def scint_pulses_unavailable(self, ev):
        # Reason the overlay can't be drawn, or None if it can. Cheap checks only
        # (dict lookups + a path probe) so this is safe to call on every redraw.
        if ev is None:
            return 'no event selected'
        try:
            # reco_directory is set from the config after the tabs are constructed
            if self._find_recon('scintillation.sbc', self.run) is None:
                return 'no recon scintillation.sbc for this run'
        except AttributeError:
            return 'reco directory not configured yet'
        latch = self.scint_latch.get(ev)
        if latch is None or not np.isfinite(latch):
            return 'no scint t0 for this event (missing scint_t0.sbc, or Failed)'
        return None
