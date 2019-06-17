/*
 * Copyright (C) 2015 Toby N. Carlson
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */

//param ps
//param ts
//param dep
//param nobs_wind
//param station_height
//param zh
//param ff0
//param dd0
//param strtim
//param timend
//param outtt
//param ugs
//param vgs
//param int dayOfTheYear
//param dual_ti
//param ti_a
//param ti_b
//param f
//param fsub
//param wmax
//param btemp
//param latitude
//param longitude
//param input_year
//param input_month
//param input_day
//param input_time
//param slope
//param aspect
//param albg
//param albf
//param frveg
//param xlai
//param cloud_flag
//param cld_fract
//param epsi
//param epsf

//constants we will use
C_TO_K = 273.15
GRAV = 9.78
FT_TO_METRES = 3.281
KTS_TO_METRES = 1.944
ROT_RATE_EARTH = 7.27E-5
SIGMA =  5.6521e-8

//constants we can tune
deltaz = 250	//250 meter intervals for vertical spacing

/*update by Yannis Konstas(start)*/               
EdA = poly([12, - 0.0569, - 0.000202, 0.00000825,  - 0.000000315] ,"xlat","coeff")                
EdB = poly([0.123,  - 0.000310, 0.0000008, 0.00000825,  0.000000499] ,"xlat","coeff") * %s 

edA = horner(EdA, xlat)
edB = horner(EdB, xlat)

temp = sin(PI * ( (dayOfTheYear + 10.0) / 365.0))
edN = (0.945 * ( edA + edB * temp^2 )) /24.0

fine_h = 50 + 50 * (0:5*ntrp)  //finer subdivision, each 50m from 50m (height of surface layer), onwards

ts = ts + C_TO_K //convert sounding temperature to Kelvin
tdew = ts - dep
temp_const = 2.49e6 / (461.51 * C_TO_K)
ew = 6.11 * exp(temp_cinst *(1  - C_TO_K ./ tdew))    //es in millibars
tbar = (ts(1:$-1) + ts(2:$)) / 2            /* Average temperature in kelvin */
delta_log_ps = log(ps(2:$) ./ ps(1:$-1) ) //log(ps[i+1] / ps[i]
thick = - (287/9.8) * (tbar .* delta_log_ps)
zls = cumsum(thick)

qs = (0.622 * ew) ./ ps                         //Specific humidity in g/kg
derivative_qs = splin(zls, qs)
q_fine = interp(fine_h, zls, qs, derivative_qs)
qd = q_fine(1:5:$)

pot_temp = ts * (1000.0 ./ ps)^0.286
derivative_pot_temp = splin(zls, pot_temp)
t_fine = interp(fine_h, zls, pot_temp, derivative_pot_temp)
td = t_fine(1:5:$)

gm = (pot_temp(2:$) - pot_temp(1:$-1)) ./ thick
gmq = (qs(1:$-1) - qs(2:$)) ./ thick
convo = ew(1:$-1) .* ps(2:$) + ps(1:$-1) .* ew(2:$)
delta_ps = ps(2:$) - ps(1:$-1)
delta_ew = ew(2:$) - ew(1:$-1)
ratio_ps = delta_log_ps ./ delta_ps
precip_h20 = -0.622 / GRAV * 10.0 * (convo .* ratio_ps + delta_ew)
omega = cumsum(precip_h20)($)

// Numpber of levels it's possible to interpolate (ntrp)
height_at_ntrp = min([height zh[nobs_wind]])
ntrp = floor((height_at_ntrp - 50.0) / deltaz)

//Winds - note the way that they're recorded by the weather service
zh = ((zh - station_height) * 1000) / FT_TO_METRES;
velocity = -(ff0 / KTS_TO_METRES) * exp( (%i * %pi * dd0) / 180 ) //velocity = vcomp + j ucomp

ucomp = imag(velocity)
derivative_ucomp = splin(zh, ucomp)
u_fine = interp(fine_h, zh, ucomp, derivative_ucomp)
ud = u_fine(1:5:$)

vcomp = real(velocity)
derivative_vcomp = splin(zh, vcomp)
vd = interp(fine_h, zh, vcomp, derivative_vcomp)
vd = v_fine(1:5:$)

//Calculate the lapse rate in the first layer and ew at screen level
//for use in the calc of screen level sat'n spec humidity

atemp = 50 * (ts(2) - ts(1)) / zls(2) + ts(1)
tscren = ts(1)	// Screen temperature

screen_level_ew= zeros(1,50)
screen_level_ew =  ew(1)

t= zeros(1,13))
// Store and initialise temperatures
t(1) = (atemp + 0.5)
otemp = tscren - 2
frveg = frveg / 100	//Vegetation %age as a fraction

oshum = (0.622 * ew(1) / ps(1))
ahum = qs(1)
old_ahum = qs(1)


//Changes 2/10/92
pres_50 = ps(1) * Math.exp(-9.8 * 50 / (287 * ts(1)))
pot_50 = (atemp * Math.pow((1000 / pres_50), 0.286));
aptemp = pot_50
o_pot_tmp = pot_50
tdif_50 = pot_50 - atemp
tdif_s = tdif_50 - 0.5

//Some basic calculations
k = floor(xlat)
renormalized_xlat = ((xlat - k) / 0.6 + k)
k = floor(xlong)
renormalized_xlong = ((xlong - k) / 0.6 + k)

//Convert lat to radians, calculate coliolis force
cf = 2 * ROT_RATE_EARTH * Math.sin( (%pi * renormalized_xlat) / 180)

no_rows = (floor(timend) - floor(strtim)) / floor(outtt) + 1

//interpolation for geostrophic winds
//Generate the daytime and nighttime vertical wind profiles of the
//geostrophic wind components at intervals of 250m from 50m (top of the
//surface layer) from the surface geostrophic winds.

v_geostrophic_daytime = zeros(1:20)	//per 250m
u_geostrophic_daytime = zeros(1:20)	//per 250m
v_geostrophic_nighttime = zeros(1:46)	//per 50m
u_geostrophic_nighttime = zeros(1:46)	//per 50m

v_geostrophic_daytime(1) = vgs  //initial value
v_geostrophic_daytime(5:ntrp) = vd(5:ntrp) //put values from vd from 1050m and above
v_geostrophic_daytime(ntrp:20) = v_geostrophic_daytime(ntrp) //make constant from ntrp and above
v_geostrophic_daytime(2:4)=interpln([1 5; v_geostrophic_daytime(1) v_geostrophic_daytime(5)],2:4) //interpolate linearly in between

u_geostrophic_daytime(1) = ugs  //initial value
u_geostrophic_daytime(5:ntrp) = ud(5:ntrp) //put values from ud from 1050m and above
u_geostrophic_daytime(ntrp:20) = u_geostrophic_daytime(ntrp) //make constant from ntrp and above
u_geostrophic_daytime(2:4)=interpln([1 5; u_geostrophic_daytime(1) u_geostrophic_daytimeeostrophic_daytime(5)],2:4) //interpolate linearly in between

v_geostrophic_nighttime(1)  = vgs //initial value
v_geostrophic_nighttime(21:46) = v_fine(21:46) //put values from v_fine from 1050m and above
v_geostrophic_nighttime(2:20)=interpln([1 21; v_geostrophic_nighttime(1) v_geostrophic_nighttime(21)],2:20) //interpolate linearly in between

u_geostrophic_nighttime(1)  = ugs //initial value
u_geostrophic_nighttime(21:46) = u_fine(21:46) //put values from u_fine from 1050m and above
u_geostrophic_nighttime(2:20)=interpln([1 21; u_geostrophic_nighttime(1) u_geostrophic_nighttime(21)],2:20) //interpolate linearly in between

// Calculates the geostrophic temperature advection based on the thermal
// wind equation and the vertical distribution of geostrophic wind.

//The depth of the layer over which we calculate the geostrophic
//temperature advection is 1000m, z[5]-z[1].
temp_const = cf * otemp / (GRAV * 1000)
dtdx = temp_const * (v_geostrophic_daytime(5) - v_geostrophic_daytime(1))
dtdy = temp_const * (u_geostrophic_daytime(5) - u_geostrophic_daytime(1))
advgt = v_geostrophic_daytime(3) * dtdy - u_geostrophic_daytime(3) * dtdx 

//Assume that the actual temperature change is one half that of
//geostrophic temperature advection.
advgt = advgt / 2

//Calculate the initial value of the wind at 50m... vel, bri
awind = Math.sqrt(ud(1)^2 + vd(1)^2)

// This section of the code takes the precipitable water (omega) and
// calculates the transmission coefficient for absorption using a lookup
// table. Interpolate linearly on omega. The lookup table contains 46-entry
// tables for omega from 0 to 5 in increments of 0.5.
// Note that the maximum value allowed for omega is 5.0.
// Code also copies the scattering (scatbl) and the backscattering (bsctbl)
// tables from the file into the common block.

omega = min(omega, 5)
i = floor(1 + omega * 2)

abstbl = LUT(1+(i-1)*46:1+i*46,1)'
scatbl = LUT(1+(i-1)*46:1+i*46,2)'
bsctbl = LUT(1+(i-1)*46:1+i*46,3)'


gt_a = zeros(1, 46)
gt_b = zeros(1, 46)
gt_c = zeros(1, 46)
if i > 1 then
    gt_a = LUT(1+(i-2)*46:1+(i-1)*46,1)'
    gt_b = LUT(1+(i-2)*46:1+(i-1)*46,2)'
    gt_c = LUT(1+(i-2)*46:1+(i-1)*46,3)'
else
    
//convex re-weighting
weight_factor = 2 * omega + 1 - i //the fractional part
abstbl = abstbl * (1-weight_factor) + gt_a * weight_factor

//Calc. lambda, kappa, and volumetric heat capacity of the ground (cg)
//See manual (HA!) for explanation of lambda and kappa.

//Gillies, 1/11/95: The thermal inertia (tp) is entered in SI (Wm-2K-1)
//units. Convert to CGS to be able to use the regression equation derived
//from Sellers. This regression equation should be redone with MKS units.
//Initiate dual ti option, tnc dajr marth 1996. */
		
tp = tp / 356.89466 // Conversion to cal cm-1 K-1 s-1/2

if dual_ti == 'y' then
    tp = ti_a * f + ti_b
end

lambdaPoly = poly([-0.00013, 0.0502, 1.21] ,"tp","coeff")
lambda = horner(lambdaPoly, tp)
kappa = (lambda / tp)^2
kappa = kappa / 10000  // Now convert to MKS; tp entered in CGS.
lambda = lambda * 418.68
cg = lambda / kappa;

//Find the best depth profile as a function of KAPPA.
//The following calcs are used in below. Creates a temp profile
//in the soil. Figure levels using scheme by Deardorff.
//del is set to numerically stable value. Initial temps at substrate
//levels (nlvls) calculated based on a linear interpolation between
//otemp and btemp

oma = 0.0000727
psoil_pass = 1
btemp = btemp + C_TO_K
sffk = 4
delta = 90
nlvl1 = 6
del = sqrt(sffk * delta * kappa) / (exp(1) - 1)
tt=interpln([1 6; otemp btemp],1:6)


//Here the substrate spacing (depth) is being calculated based on the
//scale depth (h=1+z/d).
z = del * (exp(0:5) - 1); xfun = 1 + z / del

//Correction of lambda to account for reduction in temp wave in 1st
//soil layer. Consult manual for explanation.
//There was an orphan "7678 continue" here in the fortran.
lambda  = 2 * lambda / (1 + exp( -z(2) * sqrt(oma / (2 * kappa))))

//Initialise w2g and wgg, the substrate water budget parameters. */
w2g = wmax * fsub; wgg = wmax * f

//This is the start of the diurnal loop (timend-strtim). Nominally 22 hours   
days_per_month = [31 28 31 30 31 30 31 31 30 31 30 31] //for non-leap
cumulative_days_per_month = cumsum(days_per_month)
ss = [0.967 0.971 0.982 0.999 1.016 1.029 1.034 1.030 1.018 1.002 0.985 0.973]

//Computes up and down longwave fluxes and the net radiation. 
input_month_plus_one = 1 + modulo(input_month, 12) // recycle from 1 when we find next month
//Compute the solar distance factor
solar_distance_factor = ss[input_month] + ((input_day - 1) / days_per_month[input_month]) * (ss[input_month_plus_one] - ss[input_month]);
//Set up the number of days in the year
day_in_year = days_per_month[input_month] - days_per_month[1] + input_day

//Amended to properly determine leap years
check_leap_year = 0
if modulo(input_year, 100) == 0 then
    check_leap_year = modulo(input_year, 400)
else
    check_leap_year = modulo(input_year, 4)
end

if  day_in_year >= 60 then    
    day_in_year = day_in_year + min(check_leap_year, 1) //add booleanized correction
end    
  
day_in_year_in_radians = (day_in_year * 2 * %pi) / 365.242
day_in_year_in_degrees = (day_in_year * 360) / 365.24

sig = 279.9348 + day_in_year_in_degrees
sig = sig - (0.079525 + 1.914827 * %i) * exp(%i * day_in_year_in_radians)
sig = sig - (0.00162 + 0.019938 * %i) * exp(2 *%i * day_in_year_in_radians)
sig = real(sig)

effdec = asin(0.39784988 * sin(sig*%pi / 180)) // Declination of the sun
                
//True solar noon
eqt = 12
eqt = eqt - (0.004289 + 0.12357 * %i) * exp(%i * day_in_year_in_radians)
eqt = eqt + (0.060783 - 0.153809 * %i) * exp(2 *%i * day_in_year_in_radians)
eqt = real(eqt) 

//Calculate the solar hour angle in radians
hrangl = (15 * (gmt - eqt) * %pi) / 180 - longitude
//Now we can finally compute the solar elevation sinus. We calculate solar elevation angle and azimuth for sloping 
//terrain when slope is non-zero. Elevation for northwest corner of square box is znw, etc. Grid size (called xmeshl)
//is in same units as znw.

truncated_slope = max(slope, 0); dipaz = aspect / 57.2958;  dipan = truncated_slope / 57.2958
noslope_solar_elevation_sinus = sin(effdec) * sin(latitude) + cos(effdec) * cos(latitude) * cos(hrangl)
solar_elevation_sinus = cos(dipan) * noslope_solar_elevation_sinus 
                      + sin(dipan) * (
                            cos(dipaz) * (sin(latitude) * noslope_solar_elevation_sinus - sin(effdec)) / cos(latitude)
                          + sin(dipaz) * cos(effdec) * sin(hrangl))                          
solar_elevation = asin(sinsolelslope)
solar_elevation = max(0.01, solar_elevation)

if ~isdef("albg", "local") then
    albg = 0.25 - 0.20 * wgg / wmax
end

if ~isdef("albf", "local") then
    albf = 0.032 / (0.1 + 0.1 * solar_elevation_sinus) + 0.1 * (1 - solar_elevation_sinus^2)
end

albedo_weight_factor = exp(-0.4 * xlai)
albdoe = ((1-albedo_weight_factor) * albf + albedo_weight_factor * albg) * frveg + (1 - frveg) * albg

//If the solar altitude is less than or equal to zero is night (swave = 0), so
//give up and go home, otherwise compute absolute optical air mass.
truncated_noslope_solar_elevation_sinus = max(noslope_solar_elevation_sinus, 0)
ps1 = ps(1)
noslope_solar_elevation = asin(truncated_noslope_solar_elevation_sinus)
rlpath = 1 / (noslope_solar_elevation_sinus + 0.15 * (180 * noslope_solar_elevation / %pi + 3.885)^ -1.253)
path = 0.001 * ps1 * rlpath
//transm calculates the solar transmission using the three-way
//lookup table produced in gettbl, we will use interpolation
truncated_rlpath = max(min(10, rlpath), 1), table_sample = 5 * (truncated_rlpath - 1) + 1
sample_index = floor(table_sample); next_sample_index = min(sample_index+1,46); sample_weight = table_sample -  sample_index
ftabs = (1-sample_weight) * abstbl(sample_index) + sample_weight * abstbl(next_sample_index,46); ftabs = (ps1 / 1013.25) * (ftabs - 1) + 1
ftscat = (1-sample_weight) * scatbl(sample_index) + sample_weight * scatbl(next_sample_index,46); ftscat = (ps1 / 1013.25) * (ftscat - 1) + 1
fbscat = (1-sample_weight) * bsctbl(sample_index) + sample_weight * bsctbl(next_sample_index)
tabs = ftabs; tscat = ftscat; bscat = fbscat;

//Set the path length for diffuse radiation equal to 1.7 
ftabs  = (abstbl(4) + abstbl(5))/2; ftabs = (ps1 / 1013.25) * (ftabs - 1) + 1
ftscat = (scatbl(4) + scatbl(5))/2; ftscat = (ps1 / 1013.25) * (ftscat - 1) + 1
fbscat = (bsctbl(4) + bsctbl(5))/2 
tabsd = ftabs; tscatd = ftscat; bscatd = fbscat;

// sheat >>> sunlight amount on horizontal plane outside atmosphere
// xser  >>> Diffuse shortware radiation at the ground.
// hi    >>> Direct shortwave radiation reaching the ground
// swave >>> Diffuse+direct = "Global"

skonst = (1.94 * 4.1868e4 / 60)
sheat = skonst * solar_elevation_sinus / solar_distance_factor
xser = bscatd * albdoe * (1 - tscatd) * tabsd  * truncated_noslope_solar_elevation_sinus
hi = (sheat * tabs * tscat) + (skonst / solar_distance_factor * tabs * (1 - tscat) * (1 - bscat)) * truncated_noslope_solar_elevation_sinus
swave = (hi * (1 - albdoe)) / (1 - xser)

if cloud_flag then            
        swave = swave * (1 - (0.7 * (0.01 * cld_fract))) //Impose cloud fraction; reduce swave accordingly.
end

//Calcluate the net radiation for the appropriate ground conditions.
//the decision is made on the fraction of vegetation.
//Code added 28/11/90 to solve the partial routine issue.
//Initialise and calculate the effective emissivity of the air and
//longwave down using a weighted average of surface and air temperatures.


aepsi = 0.7 + 0.17 * log10(omega)
if cloud_flag then            
       aepsi = aepsi + (1 - aepsi) * 0.8 * (0.01 * cld_fract)
end,
                
bareradiotemp = tscren - 2;  vegnradiotemp = tscren - 2

lwDown = aepsi * SIGMA * (t_fine(3) - tdif_s - 1.5)^4 //stefan boltzman law
lwUp = epsi * SIGMA * bareradiotemp^4 //stefan boltzman law
lwTg = SIGMA * tg^4
lwTf = SIGMA * tf^4

//Calculate incident solar flux at top of the canopy (sol)
sol = swave / (1 - albdoe);
rsg = sol * (1 - albedo_weight_factor) * (1 - albg) / (1 - albedo_weight_factor * albg * albf);
rsf = sol * (1 - albf) * albedo_weight_factor * (1 + albg * (1 - albedo_weight_factor) / (1 - albedo_weight_factor * albf * albg));
rlg = (1 - albedo_weight_factor) * (lwDown - lwTg -  epsf * albedo_weight_factor * (lwTg - lwTf); rlg = epsi * rlg
rlg = rlg / (1 - albedo_weight_factor * (1 - epsf) * (1 - epsi))
rlf = lwDown - lwTf + epsi * (lwTg - lwTf)  +  (1 - albedo_weight_factor) * (1 - epsi) * (lwDown - lwTf); rlf = albedo_weight_factor * epsf * rlf
rlf = rlf / (1 - albedo_weight_factor * (1 - epsf) * (1 - epsi))


//Note that the radiometric temperature at the surface is denoted
//by tzero and is the equivalent to otemp in the bare soil mode.
rnetg = rlg + rsg
rnetf = rlf + rsf
vegnnetradn= rnetg + rnetf
vegnshortwave = (rsg + rsf)
rnet = vegnnetradn
barenetradn = lwDown * epsi + swave - lwUp

rnet  = frveg * vegnnetradn   + (1 - frveg) * barenetradn
swave = frveg * vegnshortwave + (1 - frveg) * swave

//Resistance values in the transition and surface layers.
//entry to nighttime formulations (bri and mom) through this subroutine.                         
//Vel computes the Monin Obukhov length, the friction velocity, and
//the integral of heat diffusivity. *?
//Code altered 5th may 1992... Code transfrerred to bri.for
//alterations 16th july to account for different roughness lengths
//associated with partial vegetation calculations.
//zo roughness height, za top of surface layer. (50m.)
//zten - height at 10m, reflev -"Screen or anemometer height."

reflev = 2.0; zten = 10.0; ks = 0.0249951; kw = 0.0297681; cmh = 1; cmw = 1;

//line 2150

//The model assumes neutral conditions at the start of the run where
//heat=0. Therefore calc surface wind profile and resistances for the
//surface layer on the basis Math.log wind profile law.

/*
 * Vel computes the Monin Obukhov length, the friction velocity, and
 * the integral of heat diffusivity. *?
 * 
 * Code altered 5th may 1992... Code transfrerred to bri.for
 * alterations 16th july to account for different roughness lengths
 * associated with partial vegetation calculations.
 */
/* zo roughness height, za top of surface layer. (50m.) */
/* zten - height at 10m, reflev -"Screen or anemometer height." */
reflev = 2.0; zten = 10.0; ks = 0.0249951; kw = 0.0297681; cmh = 1; cmw = 1;


 /*
 * The model assumes neutral conditions at the start of the run where
 * heat=0. Therefore calc surface wind profile and resistances for the
 * surface layer on the basis Math.log wind profile law.
 */
if stabcriteria == 0 | ((stabcriteria == 1) && (heat <= 0.0)) then
    
                || ((stabcriteria == 1) && (heat <= 0.0))) {

        /* Neutral profile - rare */
        /* Allow for negative heat flux */
        ustar = _you_star(awind, za, zo, 0.0);	  /* friction velocity */
        uscrn = _wind(ustar, reflev, zo, 0.0);	  /* wind speeds */
        uten = _wind(ustar, zten, zo, 0.0);
        rzazo = _r_ohms(ustar, za, zo, 0.0);	  /* Resistances */
        rzascr = _r_ohms(ustar, za, reflev, 0.0);

        if (dual_regime == 1) {
                u_patch = _wind(ustar, obst_hgt, zo, 0.0);
                ustar_patch = _you_star(u_patch, obst_hgt, zo_patch, 0.0);
                rza_obst = _r_ohms(ustar, za, obst_hgt, 0.0);
                robst_patch = _r_ohms(ustar_patch, obst_hgt, zo_patch, 
                          0.0);
        } 

        /*
         * Potential, actual temp, specific humidity at 2m. Pass
         * to vegetation component.
         */
ptmp20 = aptemp + (heat * rzascr / (DENS * CP));
ta = ptmp20 - tdif_s;
qa = ahum + (evap * rzascr / (DENS * le));
                                
//line 4837                                
