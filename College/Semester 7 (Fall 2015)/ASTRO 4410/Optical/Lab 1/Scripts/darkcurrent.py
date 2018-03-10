#!/usr/bin/env python
from pylab import *
counts=[189.422294,189.43106,189.51489,189.8589,189.8461,190.45648,191.43513,193.12361,195.67702]
exposure=[0.1,1,5,10,20,40,80,160,320]
(m,b)=polyfit(exposure, counts,1)
print(b)
print(m)
yp=polyval([m,b],exposure)
plot(exposure,yp)
scatter(exposure,counts)
grid(True)
xlabel('Exposure time (s)')
ylabel('Average counts per pixel (DN)')
title('Average counts per pixel vs Exposure time')
show()