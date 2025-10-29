from sbcbinaryformat import Streamer
import tarfile
import numpy as np
from io import BytesIO
import warnings

class TarStreamer(Streamer):
    """
    This class manages opening a sbc binary file from a tar archive. 
    It reads the header and saves data into a dictionary of numpy arrays. 
    For large files, it uses block-based loading to handle partial reads 
    efficiently.
    
    :param tar_file: Path to the .tar archive containing a SBC binary file
    :type tar_file: str
    :param member: Path to SBC binary file inside tar_file
    :type member: str
    :param max_size: Maximum size in MB to load at once (0 = load all)
    :type max_size: int
    :param block_len: Number of rows per block for streaming mode
    :type block_len: int
    :raises OSError: If the file format is invalid or unsupported endianness
    :raises ValueError: If the header format is incorrect or incompatible
    :raises IndexError: If the requested range is out of bounds
    :raises NotImplementedError: If step slicing is attempted in streaming mode
    :raises ValueError: If the 'end' and 'length' parameters are used together
    :raises ValueError: If the 'start' parameter is negative and exceeds file length
    :raises ValueError: If the 'length' parameter is not positive
    :raises ValueError: If the 'end' parameter is negative and exceeds file length
    """

    def __init__(self, tar_file, member, max_size=0, block_len=100):
        self.tar_file = tar_file
        self.member = member
        with tarfile.open(self.tar_file, "r") as tf:
            self.file_size = tf.getmember(self.member).size

        self.max_size_bytes = max_size * 1024 * 1024
        self.block_len = block_len

        self._parse_header()
        self.single_read = (self.file_size < self.max_size_bytes or max_size == 0)

        if self.single_read:
            with tarfile.open(self.tar_file, "r") as tf:
                with tf.extractfile(self.member) as f:
                    self.data = np.frombuffer(f.read(), dtype=self.row_dtype, offset=self.header_size)
            self.tar_file_handle = None
            self.file_handle = None
        else:
            self.tar_file_handle = tarfile.open(self.tar_file, "r")
            self.file_handle = self.tar_file_handle.extractfile(self.member)
            self.data = None

    def _parse_header(self):
        with tarfile.open(self.tar_file, "r") as tf:
            with tf.extractfile(self.member) as f:
                endian_val = np.frombuffer(f.read(np.dtype(np.uint32).itemsize), dtype=np.uint32, count=1)[0]
                if endian_val == 0x01020304:
                    self.endianness = "little"
                elif endian_val == 0x04030201:
                    self.endianness = "big"
                else:
                    raise OSError(f"Unspported endianness: {endian_val}")

                # Read header
                header_length = np.frombuffer(f.read(np.dtype(np.uint16).itemsize), dtype=np.uint16, count=1)[0]
                header = f.read(header_length).decode('ascii')
                header_items = header.split(';')[:-1]  # Remove empty last element
                
                if len(header_items) % 3 != 0:
                    raise OSError("Invalid header format")

                # Parse columns
                num_columns = len(header_items) // 3
                self.columns = []
                fields = []
                
                for i in range(num_columns):
                    name = header_items[i * 3]
                    type_str = header_items[i * 3 + 1]
                    size_str = header_items[i * 3 + 2]
                    
                    self.columns.append(name)
                    dtype = self._parse_dtype(type_str)
                    sizes = list(map(int, size_str.split(',')))
                    
                    # Build numpy dtype field
                    if sizes == [1]:
                        fields.append((name, dtype))
                    else:
                        fields.append((name, dtype, tuple(sizes)))

                self.row_dtype = np.dtype(fields)
                # Skip expected number of elements (we'll calculate it)
                f.seek(4, 1)
                self.header_size = f.tell()
                    
                # Calculate actual number of elements
                payload_size = self.file_size - self.header_size
                if payload_size % self.row_dtype.itemsize != 0:
                    # Instead of raising an error, warn the user and proceed.
                    warnings.warn(
                        f"File '{self.tar_file}' may be truncated or corrupted. "
                        "The last partial row will be ignored.",
                        UserWarning
                    )

                self.num_elems = payload_size // self.row_dtype.itemsize
                self.block_len = min(self.block_len, self.max_size_bytes // self.row_dtype.itemsize)

    def _get_block_data(self, start, stop):
        """Get data for a single block"""
        # Load block from file
        block_start = start
        block_end = min(block_start + self.block_len, self.num_elems)
        
        byte_offset = self.header_size + block_start * self.row_dtype.itemsize
        self.file_handle.seek(byte_offset)
        
        block_data = np.frombuffer(
            self.file_handle.read((block_end - block_start)*self.row_dtype.itemsize),
            dtype=self.row_dtype, 
            count=block_end - block_start
        )
        
        # Return requested slice
        local_start = start - block_start
        local_stop = local_start + (stop - start)
        return block_data[local_start:local_stop]
    

    def __exit__(self, exc_type, exc_value, traceback):
        if self.file_handle:
            self.file_handle.close()
        if self.self.tar_file_handle:
            self.self.tar_file_handle.close()
