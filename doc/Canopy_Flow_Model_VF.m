%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
%                        Canopy_Flow_Model_VF.m
%
%  [RMRS Massman Version 2 = Variable Foliage Distribution -- June 2012]   
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
% INPUT Canopy Parameters
%
% NOTE: The canopy is assumed to be comprised of several layers, each 
%       with a fraction of the total of FAI - The foliage in each of 
%       these layers is assumed constant and each layer is bounded by 
%       an upper and a lower node - The physical height of each node 
%       is expressed as a fraction (zhfrac) of the total canopy height (h)
%
% Vertical array for foliage distribution [fraction of FAI] = folfrac
% Vertical array for height of canopy nodes [fraction of h] = zhfrac 
% NOTE: Nz = Nfol+1 (always) 
%       folfrac(end) = trunk space fraction of FAI 
%       zhfrac(1) = 1 (always) and zhfrac(end)  = z0h (always) 
%       zhfrac(1) and zhfrac(end) are the top and bottom canopy nodes
%__________________________________________________________________________

folfrac = [0.01 0.02 0.03 0.04 0.05 0.20 0.30 0.25 0.10 0.02];
%-a-folfrac = [0.22 0.20 0.22 0.20 0.16 0.02];
%-b-folfrac = [1 0.02];
Nfol    = length(folfrac);
zhfrac  = [1 0.98 0.96 0.94 0.92 0.90 0.75 0.45 0.40 0.35 z0h];
%-a-zhfrac  = [1 0.90 0.80 0.60 0.40 0.02 z0h];
%-b-zhfrac  = [1 0.40 z0h];
Nz      = length(zhfrac);

%__________________________________________________________________________
%
% Set up vertical array for estimating d/h by integration
% Also flip array so that the top of the canopy is the first node,  
%      which is opposite of the Uniform Foliage Wind Speed Model 
%__________________________________________________________________________

delh = 0.0001;
z01  = fliplr(z0h:delh:1);

%__________________________________________________________________________
%
% Calculate the Bessel function arrays at the nodes  
%__________________________________________________________________________

Aj  = (4.52 + 0.62*FAI*folfrac).*folfrac*FAI;

n   = 1:Nfol;
aru = 2*sqrt(Aj(n).*zhfrac(n));
arl = 2*sqrt(Aj(n).*zhfrac(n+1));
I0u = besseli(0,aru);
I0l = besseli(0,arl); 
K0u = besselk(0,aru);
K0l = besselk(0,arl); 
I1u = besseli(1,aru);
I1l = besseli(1,arl); 
K1u = besselk(1,aru);
K1l = besselk(1,arl); 

%__________________________________________________________________________
%__________________________________________________________________________
%
% Generalize Wang (2012) canopy wind model to variable foliage distribution 
%__________________________________________________________________________
%__________________________________________________________________________

%__________________________________________________________________________
%
% Initialize arrays for the vertical wind profile
%__________________________________________________________________________

xar = [1-Uph 0*(1:2*Nz-4) -Uph]';
Aar = zeros(2*Nz-2,2*Nz-2); 

%__________________________________________________________________________
%
% Populate array (Aar) for the vertical wind profile
%__________________________________________________________________________

Aar(1,1:2) = [I0u(1) K0u(1)];
for n      = 2:Nz-1;
  ql       = sqrt(Aj(n-1));
  qu       = sqrt(Aj(n));
  Aar(2*n-2,2*n-3:2*n) = [I0l(n-1)  K0l(n-1) -I0u(n) -K0u(n)];
  Aar(2*n-1,2*n-3:2*n) = [ql*I1l(n-1) -ql*K1l(n-1) -qu*I1u(n)  qu*K1u(n)];
end
Aar(end,end-1:end)     = [I0l(end) K0l(end)];

%__________________________________________________________________________
%
% Calculate the wind profile coefficients, C and D, for each canopy layer  
%__________________________________________________________________________

CD = Aar\xar;

%__________________________________________________________________________
%
% Calculate the wind profile and the normalized shear profile from the 
% C and D coefficients and the general form of Wang's analytical model
%
% But first calculate the integer depth (Ncan) of each foliage layer 
%__________________________________________________________________________

Ncan(1:Nz-1) = round(round(1/delh)*(zhfrac(1)-zhfrac(2:Nz)));

%__________________________________________________________________________
%
% Calculate Wang's (normalized) canopy wind speed profile = UzUh = U(z)/U_h
% Calculate the corresponding normalized wind shear profile = dUdz
% Calculate the corresponding normalized momentum flux profile = Ust2 
%
% Then Plot the results 
%__________________________________________________________________________

m       = 1:Ncan(1);
ZA(m)   = 2*sqrt(Aj(1)*z01(m));
UzUh(m) = CD(1)*besseli(0,ZA(m)) + CD(2)*besselk(0,ZA(m)) + Uph;
dUun(m) = ZA(m).*(CD(1)*besseli(1,ZA(m))-CD(2)*besselk(1,ZA(m)))./z01(m);
for n = 2:Nz-2
    m = Ncan(n-1)+1:Ncan(n);
    ZA(m)   = 2*sqrt(Aj(n)*z01(m));
    UzUh(m) = CD(2*n-1)*besseli(0,ZA(m))+CD(2*n)*besselk(0,ZA(m))+Uph;
    dUun(m) = ZA(m).*(CD(2*n-1)*besseli(1,ZA(m))-CD(2*n)*besselk(1,ZA(m)))./z01(m);
end
m       = Ncan(Nz-2)+1:Ncan(Nz-1)+1;
ZA(m)   = 2*sqrt(Aj(Nz-1)*z01(m));
UzUh(m) = CD(2*Nz-3)*besseli(0,ZA(m))+CD(2*Nz-2)*besselk(0,ZA(m))+Uph;
dUun(m) = ZA(m).*(CD(2*Nz-3)*besseli(1,ZA(m))-CD(2*Nz-2)*besselk(1,ZA(m)))./z01(m);

dUdz = dUun/dUun(1);

vkar = 0.4;
b1   = 0.35+vkar/log(z0h);
beta = 0.35-b1*exp(-4*FAI);
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

Int1 = delh*trapz(Ust2);
dish = 1-Int1;
z0Ch = Int1*exp(-vkar/beta);

%__________________________________________________________________________
%
% Use Albini and Baughmann (1979) and the input parameters H, h, and UH to 
%     calcualte the wind speed ratios UBUH and UCUH - 
%
% OUTPUT 
%
% UBUH = ratio of wind speed above the canopy (UBAR) to UH
% UCUH = ratio of canopy wind speed at the canopy displacement height (UC)
%        to UH   
%__________________________________________________________________________

dz0r = (H/h+Int1)/z0Ch;
Unrm = UH/log(dz0r);
UBAR = z0Ch*(h/H)*(dz0r*(UH-Unrm)-Unrm*(Int1/z0Ch)*(log(Int1/z0Ch)-1));
UBUH = UBAR/UH;

UdUh = UzUh(round(round(1/delh)*dish));
UCUH = UdUh*log(Int1/z0Ch)/log(dz0r);