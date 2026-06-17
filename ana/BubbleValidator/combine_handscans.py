import csv
import glob
import os
import sys
from datetime import datetime

# scan_source mappings come from EventDisplay code
SOURCE_LABELS = {0: "bulk", 1: "wall", 2: "dome", 3: "bellows", 4: "other"}
WEEKDAYS = {"Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun"}

COLUMNS = [
    "run", "ev", "scanner",
    "scan_source", "scan_source_label",
    "scan_nbub", "scan_trigger", "scan_crosshairsgood",
    "scan_comment", "source_file", "scan_time",
]


def parse_scan_time(filename):
    """scan_<run>_<scanner>_<Day>_<Mon>_<DD>_<H>_<M>_<S>_<YYYY>.txt to ISO time."""
    tokens = os.path.splitext(os.path.basename(filename))[0].split("_")
    for i, tok in enumerate(tokens):
        if tok in WEEKDAYS and len(tokens[i:]) == 7:
            try:
                return datetime.strptime(" ".join(tokens[i:]), "%a %b %d %H %M %S %Y").isoformat()
            except ValueError:
                return ""
    return ""


def parse_row(line):
    """A scan row has >=7 fields with an integer event number; else it's a header."""
    parts = line.split(None, 7)
    if len(parts) < 7:
        return None
    try:
        ev, src, nbub, trig, cross = (int(parts[i]) for i in (1, 3, 4, 5, 6))
    except ValueError:
        return None

    comment = parts[7].strip() if len(parts) == 8 else ""
    if len(comment) >= 2 and comment[0] == comment[-1] == "'":
        comment = comment[1:-1]

    return {
        "run": parts[0],
        "ev": ev,
        "scanner": parts[2],
        "scan_source": src,
        "scan_source_label": SOURCE_LABELS.get(src, f"unknown_{src}"),
        "scan_nbub": nbub,
        "scan_trigger": trig,
        "scan_crosshairsgood": cross,
        "scan_comment": comment,
    }


def main(indir, outpath):
    rows = []
    for path in sorted(glob.glob(os.path.join(indir, "scan_*.txt"))):
        scan_time = parse_scan_time(path)
        with open(path) as fh:
            for line in fh:
                row = parse_row(line)
                if row:
                    row["source_file"] = os.path.basename(path)
                    row["scan_time"] = scan_time
                    rows.append(row)

    with open(outpath, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS)
        w.writeheader()
        w.writerows(rows)
    print(f"{len(rows)} rows -> {outpath}")


if __name__ == "__main__":
    here = os.path.dirname(os.path.abspath(__file__))
    indir = sys.argv[1] if len(sys.argv) > 1 else "/exp/e961/data/SBC-25-handscan"
    outpath = sys.argv[2] if len(sys.argv) > 2 else os.path.join(here, "all_handscans.csv")
    main(indir, outpath)
