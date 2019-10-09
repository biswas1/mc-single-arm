#!/bin/bash 

cd  /w/hallc-scifs17exp/xem2/biswas/cross-section-code/mc-single-arm/worksim

for f in shms_*_carbon.rzdat

do 

h2root $f

done 
