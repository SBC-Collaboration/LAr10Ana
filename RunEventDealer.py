from EventDealer import ProcessSingleRun
import sys

if __name__ == "__main__":
    if len(sys.argv) > 1:
        ProcessSingleRun(
            rundir=sys.argv[1],
            recondir=sys.argv[2],
            process_list = ["event", "exposure", "scintillation", "acoustic"])
    else:
        ProcessSingleRun(
            rundir="/exp/e961/data/SBC-25-daqdata/20250611_1",
            recondir="/exp/e961/data/users/gputnam/test-sbcdaq", # Use your own directory for testing~
            process_list = ["event", "exposure", "scintillation", "acoustic"])
