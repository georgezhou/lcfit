#!/usr/bin/env python
# openfile.py this file is used for data io control
from types import *
import math
import os
import string
import numpy as np

#===============================================================
#=====open different type of files==============================
#===============================================================
def readcolumn(var,col,infile,datformat='float',div=None,fill=False):
	fin=open(infile,mode='r')
	data=fin.readline().split(div)
	if (datformat=='float'):
		while((not len(data)==0) and (not data == [''])):			
			if(not (data[0].startswith('#'))):
				try:
					var.append(float(data[col-1]))
				except (IndexError or ValueError):
					if(fill):
						var.append(None)
					else:
						raise 
			data=fin.readline().split(div)
	
	if (datformat=='str'):
		while((not len(data)==0) and (not data == [''])):			
			if(not (data[0].startswith('#'))):
				try:
					var.append(data[col-1])
				except IndexError:
					if(fill):
						var.append(None)
					else:
						raise IndexError

			data=fin.readline().split(div)
	fin.close()
	return 


def open_tab(time,flux,infile,coljd=0,colmag=0):
	if(infile.endswith('.tab')):
		if(coljd==0):
			coljd=1
			colmag=2
	if (infile.endswith('.pflt')):
		if(coljd==0):
			coljd=3
			colmag=4
        if (infile.endswith('.ltf')):
		if(coljd==0):
			coljd=5
			colmag=6
	if (infile.endswith('.tfa') or infile.endswith('.tfam')):
		if(coljd==0):
			coljd=5
			colmag=7
	if(coljd==0):
		return False
	else:	
		fin=open(infile,mode='r')
		name=os.path.splitext(infile)[0]
		header=fin.readline()
		data=fin.readline().split()
#		print len(data)
		if(colmag>len(data)):
			return False
		else:
			while(not len(data)==0):			
				if(not (data[0].startswith('#'))):
					if((not math.isnan(float(data[colmag-1]))) and (float(data[coljd-1])>964)):
						time.append(float(data[coljd-1]))
						flux.append(float(data[colmag-1]))
				data=fin.readline().split()
		print len(time),len(flux)
		fin.close()
		return True



def open_ref(infile,refchart,colrms=13):
	#open data in ref format
	if (infile.endswith('ref')):
		fin=open(infile,mode='r')
		line=fin.readline().split()

		while(not len(line)==0):
		#	print line
                        if(not (line[0].startswith('#'))):				
				info=curveinfo()
				info.name=line[0]
                                info.id=int(float(line[1]))
				info.type=line[2]
                                info.qn=int(line[3])
				info.sn=int(line[4])
				info.sgn=int(line[5])
				info.mod=int(line[6])
				info.chan=int(line[7])
				info.ra=float(line[8])
				info.dec=float(line[9])
				if (not (line[10] == 'NAN')):
					info.Kmag=float(line[10])
				else:
					info.Kmag='NAN'
				info.Kpmag=float(line[11])
				info.median=float(line[12])
				info.rms=float(line[colrms])
				refchart.append(info)	
			line=fin.readline().split()	
		fin.close()
		return True
	else:
		return False
	
#==============================================================
#==========output block========================================
#==============================================================
	#dump data to .tab files
def dump_tab(time,flux,id,ext='tab',path='./'):
#	from readin import *
	outfile=str(filecart(id,ext=ext,path=path))
        if not os.path.exists(path):
        	os.mkdir(path)
        fout=open(outfile,mode='w')
	try:
               	length=len(time)
                header1='#%s %s\n' % ('BJD','MAG')
       	        fout.write(header1)
                for i in range(length):
                        line='%s %s\n' % (time[i],flux[i]) 
			fout.write(line)
       	except IOError:
                print 'file write failed from %s' % (outfile)
      	finally:
                fout.close()
	return

def dump_ltab(time,flux,mom1,mom2,id,ext='ltab',path='./'):
#	from readin import *
	outfile=str(filecart(id,ext=ext,path=path))
        if not os.path.exists(path):
        	os.mkdir(path)
        fout=open(outfile,mode='w')
	try:
               	length=len(time)
                header1='#%s %s %s %s\n' % ('BJD','MAG','MOM_CENT1','MOM_CENT2')
       	        fout.write(header1)
                for i in range(length):
                        line='%.9f %.9f %.9f %.9f\n' % (time[i],flux[i],mom1[i],mom2[i]) 
			fout.write(line)
       	except IOError:
                print 'file write failed from %s' % (outfile)
      	finally:
                fout.close()
	return
def dump_pflt(time,flux,colx,coly,outfile,ext='pflt',path='./'):
    #dump data after prefilter
	if(not (outfile.endswith(ext))):
		outfile=os.path.splitext(outfile)+ext
       	if not os.path.exists(path):
                os.mkdir(path)
       	fout=open(outfile,mode='w')
        try:
                lengths=len(time)
		lengtho=len(colx)
       	        header1='#%s %s %s %s\n' % ('BJD','MAG','BJD1','MAG1')
               	fout.write(header1)
		length=min(lengths,lengtho)
                for i in range(length):
			line='%s %s %s %s\n' % (time[i],flux[i],colx[i],coly[i])
                        fout.write(line)
	#here assume length<=lengths
		for i in range(length,lengths):
			line='%s %s NAN NAN\n' % (time[i],flux[i])
                        fout.write(line)

	except IOError:
                print 'file write failed from %s' % (outfile)
        finally:
                fout.close()

def dump_ltf(name,time,flux,colw,colx,coly,colz,colu,qn,outfile='',ext='.ltf',path='./'):
       #dump data after long trend fiting
	
        if not os.path.exists(path):
        	os.mkdir(path)
	if outfile=='' :
		outfile=path+name+ext
       	fout=open(outfile,mode='w')
        try:
              	length=len(time)
		print length,len(flux),len(colw),len(colx),len(coly),len(colz),len(colu),len(qn)
	        header1='#%14s %10s %10s %10s %10s %10s %10s %10s\n' % ('BJD','RawMAG','RawMagerr','CMag','CMagerr','CHmag','CHmagdetrend','Quarter')
	        header2='#%14s %10s %10s %10s %10s %10s %10s %10s\n' % ('[1]','[2]','[3]','[4]','[5]','[6]','[7]','[8]')
                fout.write(header1)
                fout.write(header2)
		for i in range(length):
                        line='%-14.7f %-10.6f %-10.6f %-10.6f %-10.6f %-10.6f %-10.6f %2d\n' %(time[i],flux[i],colw[i],colx[i],coly[i],colz[i],colu[i],qn[i]) 
                        fout.write(line)

	except IOError:
                print 'file write failed from %s' % (outfile)
       	finally:
                fout.close()
	return

#===================================================================
#=================column operation==================================
#===================================================================

def replacecol(infile,coljd,colmag,jd,mag):
    	fin=open(infile,mode='r')
        filei=[]
        line=fin.readline().split()
        while(not len(line)==0):
                filei.append(line)
                line=fin.readline().split()
        fin.close()
        fout=open(infile,mode='w')
        length=len(jd)
 #       print len(filei),length
        for i in range(length):
                newline=filei[i]
		if(not newline[0].startswith('#')):
			newline[coljd-1]=str(jd[i-2])
			newline[colmag-1]=str(mag[i-2])
               	fout.write(' '.join(newline))
	return
	
def pendcol(infile,mag,colmag=-1,ext='#',outfile='',opath='',maskflag=True):
	fin=open(infile,mode='r')
        filei=[]
	head=[]
        line=fin.readline().split()
        while(not len(line)==0):
		if(line[0].startswith('#')):
			head.append(line)
		else:
                	filei.append(line)
                line=fin.readline().split()
        fin.close()
	if(outfile==''):
       		if((ext=='#')):
			outfile=infile
		elif(opath==''):
			outfile=os.path.splitext(infile)[0]+ext
	if not (opath==''):
		outfile=opath+os.path.basename(outfile)	
		
	fout=open(outfile,mode='w')
	for i in range(len(head)):
		newline=head[i]
		fout.write(' '.join(newline))
		fout.write(' %s\n' % ('Prefmag'))

	if(colmag < 0):	
		if (maskflag):
			indcmag=mag>0
			mag=mag[indcmag]
			filei=np.array(filei)
			filei=filei[indcmag]
			length=len(mag)
        		for i in range(length):
               			newline=filei[i]
            			fout.write(' '.join(newline))
                		fout.write(' %13.7f\n' % (mag[i]))
		else:
			length=len(mag)
	#		indcmag=mag<0
	#		mag[indcmag]=np.nan	
        		for i in range(length):
               			newline=filei[i]
            			fout.write(' '.join(newline))
                		fout.write(' %13.7f\n' % (mag[i]))
                		#fout.write(' %d\n' % (mag[i]))

	else:
		if (maskflag):
			indcmag=mag>0
			mag=mag[indcmag]
			filei=np.array(filei)
			filei=filei[indcmag]
			length=len(mag)
        		for i in range(length):
               			newline=filei[i]
				newline[colmag-1]='%13.7f' % mag[i]
            			fout.write(' '.join(newline))
				fout.write('\n')
		else:
			length=len(mag)
	#		indcmag=mag<0
	#		mag[indcmag]=np.nan	
        		for i in range(length):
               			newline=filei[i]
				newline[colmag-1]='%13.7f' % mag[i]
            			fout.write(' '.join(newline))
                		fout.write('\n')
	
       	fout.close()
	return

def readcol(infile,coljd,colmag):
	jd=[]
	mag=[]
	fin=open(infile,mode='r')
	line=fin.readline().split()
	while(not len(line)==0):
		if (not line[0].startswith('#')):
			if (not math.isnan(float(line[colmag-1]))):
				jd.append(float(line[coljd-1]))
				mag.append(float(line[colmag-1]))
		line=fin.readline().split()
	fin.close()
	return [jd,mag]

def countcol(infile):	
	fin=open(infile,mode='r')
	line=fin.readline().split()
	while(line[0].startswith('#')):
		line=fin.readline().split()
	return len(line)

#=====================================================================
if __name__=='__main__':
	print 'this is openfile.py'
#	infile='../reffile/kic.txt'
#	starchart=[]
	if(datain.open_tab(curve,infile)):
		for i in range(len(curve.id)):
			print 'the id is %f, BJD is %f' %(curve.id[i],curve.BJD[i])
	else:
		print 'tab file not found'	

