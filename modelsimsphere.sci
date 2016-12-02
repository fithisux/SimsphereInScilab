C_TO_K = 273.15
/*update by Yannis Konstas(start)*/               
EdA = poly([12, - 0.0569, - 0.000202, 0.00000825,  - 0.000000315] ,"xlat","coeff")                
EdB = poly([0.123,  - 0.000310, 0.0000008, 0.00000825,  0.000000499] ,"xlat","coeff") * %s 

edA = horner(EdA, xlat)
edB = horner(EdB, xlat)

temp = sin(PI * ( (dayOfTheYear + 10.0) / 365.0))
edN = (0.945 * ( edA + edB * temp^2 )) /24.0

deltaz = 250	//250 meter intervals for vertical spacing
tdew = ts - dep
temp = 2.49e6 / (461.51 * C_TO_K)
ew = 6.11 * exp(temp *(1  - C_TO_K ./ tdew))    //es in millibars
qs = (0.622 * ew) ./ ps                         //Specific humidity in g/kg
pot_temp = ts * (1000.0 ./ ps)^0.286




ext_ts = [0 ts]
tbar = (ext_ts(1:$-1) + ext_ts(2:$)) / 2            /* Average temperature in kelvin */
ext_ps = [1 ps]
thick = (287/9.8) * (tbar .* log(ext_ps(1:$-1) ./ ext_ps(2:$) ) //log(ps[i] / ps[i+1]
zls = cumsum(thick)

ext_pot_temp = [0 pot_temp]
gm = (ext_pot_temp(2:$) - ext_pot_temp(1:$-1)) ./ thick
ext_qs = [0 qs]
gmq = (ext_qs(1:$-1) - ext_qs(2:$)) ./ thick

precip_h20 = 
		/* This is where the sounding data was read in */
		/* Do some initial calculations on the soundings. */
		height = 0;

		for (j = 1; j <= nobs_ptq; j++) {
                        i = j - 1;
                        tdew = ts[j] - dep[j] + C_TO_K;	 /* Dew point in kelvin */
                        ew[j] = (6.11 * Math.exp((2.49e6 / 461.51) 
                                                * (1 / C_TO_K - 1 / tdew)));	 /* es in millibars */
                        qs[j] = (0.622 * ew[j] / ps[j]);	/* Specific humidity in g/kg */
                
                        pot_temp[j] = ((ts[j] + C_TO_K) * Math.pow((1000.0 / ps[j]), 0.286));	
                                         /* Theta in kelvin */
                
                        if (j > 1) {
                                tbar = ((ts[i] + ts[i + 1]) / 2) + C_TO_K;	/* Average temperature in kelvin */
                                thick = (287 * (tbar / 9.8) * Math.log(ps[i] / ps[i + 1]));	  
                                        /* Thickness in meters */
                                height = height + thick;	 
                                /* Height of pressure level above station */
                                zls[i + 1] = height;
                                gm[i] = (pot_temp[i + 1] - pot_temp[i]) / thick;		 			/* d(theta)/dZ */
                                gmq[i] = (qs[i] - qs[i + 1]) / thick;	 /* d(q)/dZ */
                                precip_h20 = (-0.622 / GRAV * 10.0 
                                                  * ((ew[j - 1] * ps[j] - ew[j] * ps[j - 1]) 
                                                  * Math.log(ps[j] / ps[j - 1]) 
                                                         / (ps[j] - ps[j - 1]) + ew[j] - ew[j - 1]));
                                sum_precip_h20 = sum_precip_h20 + precip_h20;
                                omega = sum_precip_h20;
                        } 
		} 
