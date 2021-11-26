#! /bin/sh

input=synthetic_shot.su
output=ghost+synthetic_shot.su

r=.3
dp=10
f1=90
f2=125
pmin=-400
pmax=2000
verbose=1
# add ghosts using the forward ghosting operator
suradon < $input choose=0 igopt=3 interoff=-262 offref=-3237 \
     pmin=$pmin pmax=$pmax dp=$dp f1=$f1 f2=$f2 cdpkey=ep anderson=0  |
sudeghost h=10 r=$r deghost=0 verbose=$verbose  | 
suradon choose=4 igopt=3 interoff=-262 offref=-3237 \
     pmin=$pmin pmax=$pmax dp=$dp f1=$f1 f2=$f2 cdpkey=ep anderson=0 > $output 

# remove ghost with the deghosting operator
input=ghost+synthetic_shot.su
output=deghosted.su
suradon < $input choose=0 igopt=3 interoff=-262 offref=-3237 \
     pmin=$pmin pmax=$pmax dp=$dp f1=$f1 f2=$f2 cdpkey=ep anderson=0  |
sudeghost h=10 r=$r deghost=1 |
suradon choose=4 igopt=3 interoff=-262 offref=-3237 \
     pmin=$pmin pmax=$pmax dp=$dp f1=$f1 f2=$f2 cdpkey=ep anderson=0  > $output

sudiff ghost+synthetic_shot.su synthetic_shot.su > ghost.su

sudiff ghost+synthetic_shot.su deghosted.su > ghost_estimate.su

sudiff synthetic_shot.su deghosted.su > residual.su

# plot results and differences
suxwigb < synthetic_shot.su interp=1 key=offset title="Synthetic shot gather" \
     label1="time (s)" label2="offset (m)" xbox=10 ybox=10 wbox=400 hbox=600 & 

suxwigb < ghost+synthetic_shot.su interp=1 key=offset title="Ghost + Synthetic shot "  clip=1.03 \
     label1="time (s)" label2="offset (m)" xbox=400 ybox=10 wbox=400 hbox=600 & 

suxwigb < ghost.su interp=1 key=offset title="Ghost " clip=1.03 \
     label1="time (s)" label2="offset (m)" xbox=800 ybox=10 wbox=400 hbox=600  & 

suxwigb < deghosted.su interp=1 key=offset title="Deghosted data " clip=1.03 \
     label1="time (s)" label2="offset (m)" xbox=10 ybox=200 wbox=400 hbox=600 & 

suxwigb < ghost_estimate.su interp=1 key=offset title="Ghost estimated " clip=1.03 \
     label1="time (s)" label2="offset (m)" xbox=400 ybox=200 wbox=400 hbox=600 & 
suxwigb < residual.su  interp=1 key=offset title="Residual= synthetic - deghosted " clip=1.03 \
     label1="time (s)" label2="offset (m)" xbox=800 ybox=200 wbox=400 hbox=600 & 
cat synthetic_shot.su ghost+synthetic_shot.su ghost.su deghosted.su ghost_estimate.su residual.su |  suxwigb interp=1 title="syn         ghost+syn         degh         ghost_est         residual" &
rm radontmp*
exit 0

