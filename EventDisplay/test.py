from cProfile import run
from tkinter import Image
import zipfile
import linecache
from eventdisplay.PICOcode.DataHandling.GetEvent import Event
from PIL import Image

path_to_zip = '.../data/2l-15-data/'

test_path = '../data/30l-16-data/'

file_name = '20160322_0.zip'
test_file = '20160912_0/0/PLClog.txt'

archive = zipfile.ZipFile(path_to_zip + file_name, 'r')

with archive.open('20160322_0/5/PLClog.txt', 'r') as file:
    entries = file.readlines()[6].split()
    entries = [entry.decode() for entry in entries]
    print(entries)


print(linecache.getline(test_path + test_file, 7).split())
