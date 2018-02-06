%this script calculates the capture cross section and demarcation energies
%posited by Bube

%manually import attempt frequencies in units of Hz. Check lab notebook for
%details

nu0totHz = 1.3e12;
nu0totHz_err = .1e12;

nu00totHz = 1.5e7;
nu00totHz_err = .1e7;

%input constants used to calculate the capture cross section in cm^{-2}
m_0 = 9.10938356e-31; %electron mass in kg
m_h = .34*m_0; %hole effective mass
kBJ = 1.3806485e-23; %Boltzmann constant in joules
area_cm = 4e-2; %area of the device in cm^{-2}
area = 4e-6; %area of the device in m^{-2}

%calculate thermal velocuty cm/s
v_th = 1e2* sqrt((3*kBJ*300)/m_0);

%calculate particles per cm^3
pf = .6;
v_nc = (4/3)*pi*1.5e-9^3;
vfilm= area*.9e-7;
Nparticles = 1e6* (vfilm/v_nc) *pf; 

%assuming we have 2 states per QD (ask Cherie for a ref...), then
%calculated the density of valence band states N_{v} per cm^{-3}
Nv = 2*Nparticles; 

%calculate the hole capture cross-sections an error in cm^{2} for the traps detected
ccsHl = zeros(2,1);

ccsHl(1,:) = nu0totHz/(v_th*1e19);
ccsHl(2,:) = (nu0totHz_err/nu0totHz) .* ccsHl(1,:);

%calculate the critical electron capture cross section in cm^{2} assuming
%that E_{F} is pinned at the interface state and the bulk trap at E_{t} is
%at the demarcation energy separating effective recombination centers from
%charge trap centers
Edp = 0.32;
Efn = 1.0;
Efp = 0.19;
me = 0.3;
mh = 0.34;

p = 1e16;
n = 2e17; 

%ccsEl = ccsHl(1,1) * 1/(exp( (Edp - Efn - (3/2)*kB*300*log(mh/me))/(kB*300) ))  

ccsEl = (ccsHl*p)/n * exp(- (Edp - Efp)/(kB*300)) 