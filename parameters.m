close all
clear all

%% Tirzepatide parameters

%Define tirzepatide PK parameters taken from:
%https://www.ac cessdata.fda.gov/drugsatfda_docs/nda/2022/215866Orig1s000ClinPharmR.pdf
tirz.name="tirzepatide";
tirz.ka=0.0373; %absorption rate (1/hour)
tirz.Vc=2.47;
tirz.Vp=3.98;
tirz.cl=0.0329*(104.8/70)^(0.8);
tirz.q=0.126*(104.8/70)^(0.8);
tirz.F=0.8;
tirz.molar_mass=4814; % g/mol 
tirz.bw0=104.8;

% PD
tirz.kout=2.1644e-04;
tirz.emaxi=0.0259;
tirz.koutpl=3.6732e-05;
tirz.plbmax=0.3411;
tirz.tfac=0.0016;
tirz.emaxs=0.2761;
tirz.ec50i=45.6765;
tirz.ec50s=112.8899;
tirz.newss=0.9781;
tirz.kout_s_plcb=5.8012e-04;

tirz.n1=1;
tirz.n2=1;

%% Semaglutide parameters

semag=struct();
semag.name="semaglutide";

% PK
% Table S2A in Strathe et al 2023
semag.ka=0.0296; % absorption rate (1/hour)
semag.V=12.1; % volume of distribution (in L)
semag.ke=0.0488/semag.V; % ke=Cl/V, elimination rate is clearance divided by volume (1/hour)
semag.F=1;
semag.molar_mass=4114; % g/mol 

% PD
% Table S3 in Strathe et al 2023 (except emaxi and emaxs)
% the rates have been converted from 1/week to 1/hour
semag.kout= 0.0319/(7*24);
semag.kout_s_plcb= 0.0319/(7*24);
semag.emaxi=0.0419; % fitted for an average, demographic-blind patient
semag.koutpl= 0.00885/(7*24);
semag.plbmax=0.01*34.4;
semag.tfac=0.0794/(7*24);
semag.emaxs=0.2257; % fitted for an average, demographic-blind patient
semag.ec50i=48.0;
semag.ec50s=48.0;
semag.newss=0.990;
semag.bw0=104.8;
semag.weekly_dose=2.4;

semag.n1=1;
semag.n2=1;

%% save parameters

save('data_parameters.mat');
