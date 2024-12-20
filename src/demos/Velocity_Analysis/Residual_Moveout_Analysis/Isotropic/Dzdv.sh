#! /bin/sh

# determine depth derivatives from ratios from two migration outputs

#set -x

 nx=5 dx=15 fx=800 noff=5 doff=200 off0=100 
	fz=0 dz=10 nz=251 
	nxw=2 nzw=2
suwind<kd.data.su key=cdp min=800 max=2000|sustrip>infile 
suwind<outfile1 key=cdp min=800 max=2000|sustrip>afile 

> dzfile \
dzdv <infile par=cig.par nx=$nx nz=$nz fx=$fx fz=$fz dx=$dx dz=$dz \
 afile=afile dfile=dfile  \
 off0=$off0 noff=$noff doff=$doff nxw=$nxw nzw=$nzw

echo
echo "Run Velpert.sh to obtain the value of dlambda."
echo "Add the value of dlambda to the original velocity to obtain new velocity"
      
exit 0
