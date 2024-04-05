from numpy import *
import math
import os.path
import numpy as np
import pandas as pd
from itertools import islice
from scipy.fft import fft, ifft,fftfreq

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.interpolate import make_interp_spline,BSpline
import matplotlib.patches as m_patches


ifile="./pubchem_organo_fluorine.csv"

# Read the first line of the CSV file to get the column labels
with open(ifile, 'r') as f:
    first_line = f.readline().strip()

# Split the first line into a list of column labels
column_labels = first_line.split(',')
print(column_labels)
print(len(column_labels))
# Read the CSV file into a DataFrame, skipping the first row and setting the column labels
df = pd.read_csv(ifile, skiprows=1, names=column_labels)

#  Select the columns  and write them to a new CSV file
df[['cid', 'isosmiles', 'charge','exactmass']].to_csv('output_file.csv', index=False)
