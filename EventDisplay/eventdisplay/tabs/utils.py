import sys
import os
from json import load
import numpy as np
from collections import OrderedDict
from importlib.metadata import version, PackageNotFoundError
from pathlib import Path


class Streamer:
    """
    This class manages opening a sbc binary file. It reads the header and 
    saves data into a dictionary of numpy arrays.
    """
    def __init__(self, file, max_size=1000000000):
        self.system_endianess = sys.byteorder
        self.file_size = os.path.getsize(file)
        self.file = file
        self.is_all_in_ram = self.file_size < max_size

        # Read header info
        with open(file, "rb") as f:
            # Read and check endianess
            file_endianess_val = np.fromfile(f, dtype=np.uint32, count=1)[0]
            if file_endianess_val == 0x01020304:
                self.file_endianess = "little"
            elif file_endianess_val == 0x04030201:
                self.file_endianess = "big"
            else:
                raise OSError(f"Endianess not supported: {file_endianess_val}")

            # Read header length and header string
            self.header_length = int(np.fromfile(f, dtype=np.uint16, count=1)[0])
            header = f.read(self.header_length).decode('ascii')
            header_items = header.split(';')[:-1]  # last element is empty

            if len(header_items) % 3 != 0:
                raise OSError("Header format error: items not in multiples of 3")

            num_columns = len(header_items) // 3
            self.columns = []
            self.dtypes = []
            self.sizes = []
            for i in range(num_columns):
                name = header_items[i * 3]
                type_str = header_items[i * 3 + 1]
                size_str = header_items[i * 3 + 2]
                self.columns.append(name)
                self.dtypes.append(sbcstring_to_type(type_str, self.file_endianess))
                self.sizes.append(list(map(int, size_str.split(','))))

            # Read the expected number of elements (can be 0 if unknown)
            self.expected_num_elems = np.fromfile(f, dtype=np.int32, count=1)[0]
            self.header_size = f.tell()

        # Create a structured dtype from the header info
        fields = []
        for col, dtype, sizes in zip(self.columns, self.dtypes, self.sizes):
            # If sizes is [1], store as scalar; otherwise as a subarray.
            if sizes == [1]:
                fields.append((col, dtype))
            else:
                fields.append((col, dtype, tuple(sizes)))
        self.row_dtype = np.dtype(fields)

        # Calculate the number of rows by comparing the payload size to the row size.
        self.payload_size = self.file_size - self.header_size
        if self.payload_size % self.row_dtype.itemsize != 0:
            raise OSError("File payload size does not match structured row size")
        self.num_elems = self.payload_size // self.row_dtype.itemsize
        if self.expected_num_elems and self.expected_num_elems != self.num_elems:
            raise OSError("Expected number of elements does not match calculated number")

        # Read the data.
        self.data = np.fromfile(file, dtype=self.row_dtype, offset=self.header_size)

    def __getitem__(self, idx):
        return self.data[idx]

    def to_dict(self):
        # Convert the structured array into a dictionary for easier access.
        return OrderedDict({name: self.data[name] for name in self.columns})

def sbcstring_to_type(type_str, endianess):
    out_type_str = ""
    if endianess == 'little':
        out_type_str += '<'
    elif endianess == 'big':
        out_type_str += '>'

    if type_str.startswith('string'):
        return np.dtype(out_type_str+type_str.replace("string", "U"))

    string_to_type = {'char': 'i1',
                      'int8': 'i1',
                      'int16': 'i2',
                      'int32': 'i4',
                      'int64': 'i8',
                      'uint8': 'u1',
                      'uint16': 'u2',
                      'uint32': 'u4',
                      'uint64': 'u8',
                      'single': 'f',
                      'float32': 'f',
                      'double': 'd',
                      'float64': 'd',
                      'float128': 'f16'}

    return np.dtype(out_type_str+string_to_type[type_str])


def type_to_sbcstring(sbc_type_str):
    if sbc_type_str.startswith("U"):
        return sbc_type_str.replace("U", "string")

    string_to_type = {'i1': 'int8',
                      'i2': 'int16',
                      'i4': 'int32',
                      'i8': 'int64',
                      'u1': 'uint8',
                      'u2': 'uint16',
                      'u4': 'uint32',
                      'u8': 'uint64',
                      'f': 'float32',
                      'd': 'double',
                      'f16': 'float128'}

    return string_to_type[sbc_type_str]

def GetEvent(path, ev):
    event_dir = Path(path) / str(ev)
    if not event_dir.is_dir():
        return None

    matches = list(event_dir.glob("acoustics*.sbc.bin"))
    if not matches:
        return None

    file_path = matches[0]

    # Try to get sample rate from rc.json one level up
    try:
        rc_path = file_path.parent.parent / "rc.json"

        with open(rc_path, "r") as f:
            rc_data = load(f)

        acous_config = rc_data.get("acous", {})
        sample_rate_str = acous_config.get("sample_rate", "").strip().upper()

        if "MS/S" in sample_rate_str:
            sample_rate = float(sample_rate_str.replace("MS/S", "").strip()) * 1_000_000
        elif "KS/S" in sample_rate_str:
            sample_rate = float(sample_rate_str.replace("KS/S", "").strip()) * 1_000
        elif "S/S" in sample_rate_str:
            sample_rate = float(sample_rate_str.replace("S/S", "").strip())
        else:
            raise ValueError(f"Unrecognized sample rate format: '{sample_rate_str}'")
    except Exception as e:
        print(f"Warning: failed to read sample rate from rc.json: {e}")
        sample_rate = None

    # Try to load the binary data
    try:
        data = Streamer(str(file_path)).to_dict()
        data["sample_rate"] = sample_rate  # Add sample rate inside the 'acoustic' dict
        return {"acoustic": data}
    except Exception as e:
        print(f"Error reading binary file: {e}")
        return None