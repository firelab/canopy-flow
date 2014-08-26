%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
%                          Canopy_Flow_Model.m
%
%    [RMRS Massman Version 1 = Uniform Foliage Distribution -- May 2012] 
%                      
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%__________________________________________________________________________
%
% INPUT
% 
% FAI = Foliage Area Index; 
% z0h = ratio of ground surface roughness length (z_0) to canopy height (h)
% Uph = ratio mean background pressure gradient to wind speed at canopy 
%       height (U_h) - This is a dimensionlees parameter: Wang's model 
% H   = height above canopy top in feet (or meters): Albini/Baughmann model 
% h   = canopy height in feet (or meters)
% UH  = wind speed at H in feet/s (or m/s)
%__________________________________________________________________________

FAI  = 1;
z0h  = 0.0025;
Uph  = 0;
H    = 20;
h    = 10;
UH   = 5;

%__________________________________________________________________________
%
% Set up vertical array    
%__________________________________________________________________________

delh = 0.0001;
z01  = z0h:delh:1;

%__________________________________________________________________________
%__________________________________________________________________________
%
% Wang (2012) - Massman (1997) - Albini & Baughmann (1979)  
%__________________________________________________________________________
%__________________________________________________________________________

%__________________________________________________________________________
%
% The Boundary Conditions detrmine C1 and C2 for Wang's (2012) model  
%__________________________________________________________________________

A    = 4.52*FAI + 0.62*FAI*FAI;
ZA0  = 2*sqrt(A*z0h);
ZA1  = 2*sqrt(A);
K0z0 = besselk(0,ZA0);
I0z0 = besseli(0,ZA0);
K0h1 = besselk(0,ZA1);
I0h1 = besseli(0,ZA1);
K0r  = K0h1/K0z0;
C1   = (1-Uph*(1-K0r))/(I0h1-I0z0*K0r);
C2   = -(Uph+C1*I0z0)/K0z0;

%__________________________________________________________________________
%
% Calculate Wang's (normalized) canopy wind speed profile = UzUh = U(z)/U_h
% Calculate the corresponding normalized wind shear profile = dUdz
% Calculate the corresponding normalized momentum flux profile = Ust2 
%
% Then Plot the results 
%__________________________________________________________________________

ZA01 = 2*sqrt(A*z01); 
UzUh = C1*besseli(0,ZA01) + C2*besselk(0,ZA01) + Uph;
dUun = 0.5*ZA01.*(C1*besseli(1,ZA01)-C2*besselk(1,ZA01))./z01;
dUdz = dUun/dUun(end);
Ust2 = z01.*dUdz;

figure(1);
plot(UzUh,z01);
figure(2);
plot(dUdz,z01);
figure(3);
plot(Ust2,z01);

%__________________________________________________________________________
%
% Use Massman (1994) to calculate the displacement height = dish = d/h 
%     and the surface roughness length = z0Ch = z_0C/h - To determine d/h 
%     it is necessary to integrate the momentum flux profile = Ust2 
%     over the canopy depth - 
% d/h and z_0/h are used with Albini and Baughmann (1979) to compute the 
%     wind speed ratios (which are the OUTPUT)
%
% NOTE: Up to this point the model results are fomulated in non-dimensional
%       terms - Therefore the only input parameters that have been used are
%       FAI, z0h, and Uph (all dimesionless) 
%__________________________________________________________________________

vkar = 0.4;
b1   = 0.35+vkar/log(z0h);
beta = 0.35-b1*exp(-4*FAI);

Int1 = delh*trapz(Ust2);
dish = 1-Int1;
z0Ch = Int1*exp(-vkar/beta);

%__________________________________________________________________________
%
% Use Albini and Baughmann (1979) and the input parameters H, h, and UH to 
%    calcualte the wind speed ratios UBUH and UCUH - 
%
% OUTPUT 
%
% UBUH = ratio of wind speed above the canopy (UBAR) to UH
% UCUH = ratio of canopy wind speed at the canopy displacement height (UC)
%        to UH
%
% NOTE: There are two ways of calculating UdUh - 
%       The model uses the approximate form, not the exact form -   
%       This allows for a better comparison with the variable foliage model       
%__________________________________________________________________________

dz0r = (H/h+Int1)/z0Ch;
Unrm = UH/log(dz0r);
UBAR = z0Ch*(h/H)*(dz0r*(UH-Unrm)-Unrm*(Int1/z0Ch)*(log(Int1/z0Ch)-1));
UBUH = UBAR/UH;

% -- ZAdh = 2*sqrt(A*dish); 
% -- UdUh = C1*besseli(0,ZAdh) + C2*besselk(0,ZAdh) + Uph;
UdUh = UzUh(round(round(1/delh)*dish));
UCUH = UdUh*log(Int1/z0Ch)/log(dz0r);