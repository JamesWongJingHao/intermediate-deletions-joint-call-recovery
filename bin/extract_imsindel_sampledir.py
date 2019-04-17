#/home/james/IMSindel_RK272_20171226
#['', 'home', 'james', 'IMSindel_RK064_20171206', '']

import sys
import re

dir_name1 = re.split("\/", sys.argv[1])
#print(dir_name1)
dir_name2 = re.split("_", dir_name1[3])
samplename = dir_name2[0] + "_" + dir_name2[1]
print(samplename)

#samplename = dir_name[0] + "_" + dir_name[1]
#print(samplename)
