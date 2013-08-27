#from pyraf import iraf
from numpy import *
import string
import sys
import functions
import os

def write_data(data,outfile_name):
    output = open(outfile_name,"w")
    output.write(str(data))
    output.close()

### Check of number is NaN
def isnan(num):
    return num == num

### Check if string is number
def is_number(s):
    test = False
    try:
        float(s)
        test = True
    except ValueError:
        test = False

    if test:
        test = isnan(s)
    return test

### takes raw ccdlist output and turns it into format of 
### [file_location,object]
def ccdlist_extract(ccdlist_info):
    for i in range(len(ccdlist_info)):
        ccdlist_info[i] = string.split(ccdlist_info[i],":")
        ccdlist_info[i][0] = string.split(ccdlist_info[i][0],"[")
        ### Identify the actual file name, excluding the path
        ccdlist_info_i = string.split(ccdlist_info[i][0][0],"/")
        ccdlist_info_file_name = ccdlist_info_i[len(ccdlist_info_i)-1]
        ccdlist_info_file_name = string.split(ccdlist_info_file_name,'.')
        ccdlist_info_file_name = ccdlist_info_file_name[0] #we don't want the .fits part
        ccdlist_info[i] = [ccdlist_info[i][0][0],ccdlist_info[i][1],ccdlist_info_file_name]
    return ccdlist_info

### takes output from ccdlist_extract and identifies 
### images of objects matching string target
### outputs list location of those images
def ccdlist_identify(ccdlist_info,target):
    match = []
    for i in range(len(ccdlist_info)):
        #if target == ccdlist_info[i][1]:
        if target in ccdlist_info[i][1]:
            match.append(i)
    return match

### takes output from ccdlist_extract and identifies 
### objects starting with HD
def ccdlist_identify_HD(ccdlist_info):
    match = []
    for i in range(len(ccdlist_info)):
        if ccdlist_info[i][1][:2] == "HD":
            match.append(i)
    return match

### Multiline input
def multiline_input(multiline):
    input_line = ""
    
    multiline = string.split(multiline,"\n")
    for i in range(len(multiline)):
        input_line = input_line + multiline[i]

    eval(input_line)

### Read ascii file function
def read_ascii(file_location):
    ascii_file_temp = []
    ascii_file = open(file_location).read()
    ascii_file = string.split(ascii_file,"\n")
    ascii_file = ascii_file[:len(ascii_file)]
    for i in range(len(ascii_file)):
        if not ascii_file[i] == "":
            if not ascii_file[i][0] == "#":
                ascii_file_temp.append(ascii_file[i])
    return ascii_file_temp

### Tables are passed on from read_ascii to read_table
def read_table(input_list):
    for i in range(len(input_list)):
        input_list[i] = string.split(input_list[i])
        input_list_temp = []
        for j in range(len(input_list[i])):
            if not input_list[i][j] == "":
                input_list_temp.append(input_list[i][j])
        input_list[i] = input_list_temp
        for j in range(len(input_list[i])):
            if is_number(input_list[i][j]):
                input_list[i][j] = float(input_list[i][j])
    return input_list

### Write table to output_file
def write_table(table,output_file):
    for i in range(len(table)):
        line_i = ""
        for j in range(len(table[i])):
            line_i = line_i + str(table[i][j]) + " "
        line_i = line_i + "\n"
        output_file.write(line_i)

### Write string to output file
def write_string_to_file(input_string,output_file):
    output_file = open(output_file,"w")
    output_file.write(input_string)
    output_file.close()

### Change list to string
def list_to_string(input_list):
    output_string = ""
    for i in input_list:
        output_string = output_string + str(i) + " "
    output_string = output_string + "\n"
    return output_string

### Read field in config file
def read_config_file(field):
    os.system("grep " + field + " config_file | awk '{print $2}' > temp")
    field_entry = open("temp").read()
    field_entry = string.split(field_entry)[0]
    os.system("rm temp")
    return field_entry

### Sort two lists according to list 1
def sort_lists(list1,list2):
    sorted_list1 = sorted(list1)
    sorted_list2 = []
    for i in range(len(sorted_list1)):
        for j in range(len(list1)):
            if sorted_list1[i] == list1[j]:
                sorted_list2.append(list2[j])
                break
    return sorted_list1,sorted_list2

### Sort array according to a specific column
def sort_array(input_array,column):
    column_to_sort = transpose(input_array)[column]
    sorted_indicies = column_to_sort.argsort()
    temp_array = []
    for index in sorted_indicies:
        temp_array.append(input_array[index])
    return temp_array

# Convert HH:MM:SS.SSS into Degrees :
def convHMS(ra):
   try :
      sep1 = ra.find(':')
      hh=int(ra[0:sep1])
      sep2 = ra[sep1+1:].find(':')
      mm=int(ra[sep1+1:sep1+sep2+1])
      ss=float(ra[sep1+sep2+2:])
   except:
      raise
   else:
      pass
   
   return(hh*15.+mm/4.+ss/240.)

# Convert +DD:MM:SS.SSS into Degrees :
def convDMS(dec):

   Csign=dec[0]
   if Csign=='-':
      sign=-1.
      off = 1
   elif Csign=='+':
      sign= 1.
      off = 1
   else:
      sign= 1.
      off = 0

   try :
      sep1 = dec.find(':')
      deg=int(dec[off:sep1])
      sep2 = dec[sep1+1:].find(':')
      arcmin=int(dec[sep1+1:sep1+sep2+1])
      arcsec=float(dec[sep1+sep2+2:])
   except:
      raise
   else:
      pass

   return(sign*(deg+(arcmin*5./3.+arcsec*5./180.)/100.))


# Convert RA (deg) to H.M.S:
def deg2HMS( RAin ):

   if(RAin<0):
      sign = -1
      ra   = -RAin
   else:
      sign = 1
      ra   = RAin

   h = int( ra/15. )
   ra -= h*15.
   m = int( ra*4.)
   ra -= m/4.
   s = ra*240.

   if(sign == -1):
      out = '-%02d:%02d:%06.3f'%(h,m,s)
   else: out = '+%02d:%02d:%06.3f'%(h,m,s)
   
   return out
   
# Convert Decl. (deg) to D.M.S:
def deg2DMS( Decin ):

   if(Decin<0):
      sign = -1
      dec  = -Decin
   else:
      sign = 1
      dec  = Decin

   d = int( dec )
   dec -= d
   dec *= 100.
   m = int( dec*3./5. )
   dec -= m*5./3.
   s = dec*180./5.

   if(sign == -1):
      out = '-%02d:%02d:%06.3f'%(d,m,s)
   else: out = '+%02d:%02d:%06.3f'%(d,m,s)

   return out

def remove_nan(input_list):
    output_list = []
    for i in input_list:
        if isnan(i):
            output_list.append(i)
    output_list = array(output_list)
    return output_list

def mean_nonan(input_list):
    output_list = []
    for i in input_list:
        if isnan(i):
            output_list.append(i)
    return mean(output_list)

def sigma_clipping(input_list,sigma):
    output_list = []
    for i in input_list:
        if abs(i-median(input_list)) < sigma * std(input_list):
            output_list.append(i)
    return array(output_list)

def sigma_clipping_twolists(input_list,input_list2,sigma):
    output_list = []
    output_list2 = []
    for i in range(len(input_list)):
        if isnan(input_list[i]):
            if abs(input_list[i]-median(remove_nan(input_list))) < sigma * std(remove_nan(input_list)):
                output_list.append(input_list[i])
                output_list2.append(input_list2[i])
    return array(output_list),array(output_list2)

def weighted_mean(input_list,input_weights):
    wgtmean = sum(input_list*input_weights) / sum(input_weights)
    return wgtmean


def find_cadence(input_lightcurve):
    diff = []
    for i in range(len(input_lightcurve)-1):
        diff.append(abs(input_lightcurve[i,0]-input_lightcurve[i+1,0]))
    if median(diff) < 0.005:
        cadence = "short"
    else:
        cadence = "long"
    return cadence
