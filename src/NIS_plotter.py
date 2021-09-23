import re
import numpy as np
import matplotlib.pyplot as plt


string1 = 'NIS Radar'
string2 = 'NIS LiDAR'

nis_radar = np.array
nis_lidar = np.array

# opening a text file to read
file1 = open("../build/UKF_log.txt", "r")

# opening text files to wite
# f = open("myfile.txt", "w")

# setting flag and index to
index = 0

# Loop through the file line by line
for line in file1:
    # checking string is present in line or not
    if string1 in line:
        x = line[11:]
        print type(int(x))
        np.append = (nis_radar, x)
        print nis_radar

# print np.size(nis_radar)
# plt.plot(nis_radar)
# plt.ylabel('some numbers')
# plt.show()

# closing text file
file1.close()
