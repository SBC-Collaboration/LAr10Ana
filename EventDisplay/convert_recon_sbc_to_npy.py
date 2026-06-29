# Builds a global, (run, ev)-keyed reco table for the EventDisplay cuts
#
# This does NOT touch raw_events.npy or per run npy files
# EventDisplay filters this table based on cuts and what is locally-available
# preventing rendering of events that are not in the unpacked run directory.
#
# python convert_recon_sbc_to_npy.py <recon-root> <npy-output-dir>
# python convert_recon_sbc_to_npy.py /path/SBC-25-recon /path/to/eventdisplay/npy/SBC-25
from glob import glob
import numpy as np
import os
import re
import sys
import time

try:
    from sbcbinaryformat import Streamer
except ImportError:
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                    '..', '..', 'SBCBinaryFormat', 'python'))
    from sbcbinaryformat import Streamer

# Per-event recon columns to expose for cuts, grouped by the .sbc file
# Each listed file must be one row per event. To surface more variables,
# just add an entry here.
# Columns from a file that is missing for a given run are left as NaN or ''
EVENT_BINARIES = {
    'event.sbc': {'pset_lo': 'f4', 'pset_hi': 'f4', 'trigger_source': 'U100'},
    # 'pressure_t0.sbc': {'t0_fitting': 'f4', 't0_chi_sq': 'f4'},
    # 't_expansion.sbc': {'expansion_time': 'f4'},
    # 'run.sbc':         {'source1_ID': 'U100', 'source1_location': 'U100'},
}

# event.sbc defines which events exist in a run so it must be present.
CANONICAL_FILE = 'event.sbc'

# Fixed output schema, so every run yields the same dtype and can be concatenated.
OUTPUT_DTYPE = np.dtype(
    [('run', 'U12'), ('ev', 'i4')]
    + [(name, dt) for cols in EVENT_BINARIES.values() for name, dt in cols.items()])


def natural_sort(things):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(things, key=alphanum_key)


def _align(canonical_ev, source_ev):
    """Indices into a source file's rows that line them up with canonical_ev (-1 where
    an event is absent from the source)."""
    pos = {int(e): i for i, e in enumerate(source_ev)}
    return np.array([pos.get(int(e), -1) for e in canonical_ev], dtype=int)


def load_run_events(run, run_dir):
    """Read the per-run recon .sbc files in run_dir and return a (run, ev)-keyed
    structured array (OUTPUT_DTYPE). Output is determined by EVENT_BINARIES"""
    canonical = Streamer(os.path.join(run_dir, CANONICAL_FILE)).to_dict()
    ev = np.asarray(canonical['ev'], dtype='i4')
    n = len(ev)

    rows = np.zeros(n, dtype=OUTPUT_DTYPE)
    for name in OUTPUT_DTYPE.names:
        if OUTPUT_DTYPE[name].kind == 'f':
            rows[name] = np.nan # unfilled numeric columns read as NaN, not 0
    rows['run'] = run
    rows['ev'] = ev

    for fname, cols in EVENT_BINARIES.items():
        path = os.path.join(run_dir, fname)
        data = canonical if fname == CANONICAL_FILE else (
            Streamer(path).to_dict() if os.path.isfile(path) else None)
        if data is None:
            print('  {}: no {}, leaving {} empty'.format(run, fname, list(cols)))
            continue
        idx = _align(ev, np.asarray(data['ev'], dtype='i4'))
        present = idx >= 0
        for name in cols:
            src = np.asarray(data[name])
            rows[name][present] = src[idx[present]]

    return rows


def main():
    if len(sys.argv) != 3:
        print('Usage: python convert_recon_sbc_to_npy.py <recon-root> <npy-output-dir>')
        sys.exit(1)

    recon_root = sys.argv[1]
    npy_dir = sys.argv[2]
    dev_output = os.path.join(recon_root, 'dev-output')

    print('Reading per-run recon from: {}'.format(dev_output))
    start = time.time()

    run_dirs = natural_sort(glob(os.path.join(dev_output, '*' + os.sep)))
    all_rows = []
    n_runs = 0
    for run_dir in run_dirs:
        run = os.path.basename(run_dir.rstrip(os.sep))
        if not os.path.isfile(os.path.join(run_dir, CANONICAL_FILE)):
            print('  skip {}: no {}'.format(run, CANONICAL_FILE))
            continue
        try:
            rows = load_run_events(run, run_dir)
        except Exception as e:
            print('  WARNING: failed to read recon for {}: {}'.format(run, e))
            continue
        all_rows.append(rows)
        n_runs += 1
        print('  {}: {} events'.format(run, len(rows)))

    if not all_rows:
        print('No {} files found under {}. Aborting.'.format(CANONICAL_FILE, dev_output))
        sys.exit(1)

    reco_all = np.concatenate(all_rows)

    os.makedirs(npy_dir, exist_ok=True)
    out_path = os.path.join(npy_dir, 'reco_events.npy')
    np.save(out_path, reco_all)

    print('Wrote {} events from {} runs to {}'.format(len(reco_all), n_runs, out_path))
    print('Fields: {}'.format(reco_all.dtype.names))
    print('finished in {:.0f} seconds'.format(time.time() - start))


if __name__ == '__main__':
    main()
