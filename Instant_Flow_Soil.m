
% This program is used to solve the instantaneous flow under different soil moisture
% The default simulation is HA1. If you want to simulate HA2, you need to hide the programs under HA1 and display the programs under HA2

%--------------------------------
% Load the different soil water gradients between shallow and deep layers
load('H'); % load the numerical sequence of soil water potential, cm
H1=H; % shallow soil water potential, cm
H2=H; % deep soil water potential, cm
SWC1=h2swc(H1, SWCr, SWCs, a, n); % shallow soil water content, cm3 cm-3
SWC2=h2swc(H2, SWCr, SWCs, a, n); % deep soil water content, cm3 cm-3
%--------------------------------
% Select one of the following four to set soil characteristic parameters
SWCr=0.065; SWCs=0.41; a=0.075; n=1.89; Ks=106.1; l=-1; % Sandy loam
% SWCr=0.045; SWCs=0.43; a=0.145; n=2.68; Ks=712.8; l=-1; % Sand
% SWCr=0.078; SWCs=0.43; a=0.036; n=1.56; Ks=24.96; l=-1; % Loam
% SWCr=0.034; SWCs=0.46; a=0.016; n=1.37; Ks=6; l=-1; % SIlt
%--------------------------------
% Set soil and plant parameters
D1=30; % depth of shallow layer, cm
D2=370; % depth of deep layer, cm, cm
S=6000; % plant rooting area, cm2
N1=10000; % vessel number of shallow root
N2=50000; % vessel number of deep root
L_root1=200; % length of shallow root, cm
L_root2=400; % length of deep root, cm
L_stem=200; % length of stem, cm
L_branch=400; % length of branch, cm
L=50; % vessel length, cm
Droot=100*10^-4; % vessel diameter of root, cm (the value before "*10^-4" represents the unit of "μm")
Dstem=70*10^-4; %vessel diameter of stem, cm (the value before "*10^-4" represents the unit of "μm")
Dbranch=30*10^-4; %vessel diameter of branch, cm (the value before "*10^-4" represents the unit of "μm")
rpit=600*10^-2; % , MPa s cm-1 (the value before "*10^-2" represents the unit of "MPa s m-1")
Fc=0.2; % faction of contact area between adjacent vessels
Fp=0.5; % fraction of pit membrane area in pitted wall
RLD_s=0.4; % root length density of shallow layer, cm cm-3
RLD_d=0.2; % root length density of deep layer, cm cm-3
RSD_s=0.036; % root surface density of shallow layer, cm2 cm-3
RSD_d=0.018; % root surface density of deep layer, cm2 cm-3
l_s=S*D1*RLD_s; % fine root length of shallow layer, cm
l_d=S*D2*RLD_d; % fine root length of deep layer, cm
SAr_s=S*D1*RSD_s; % fine root surface area of shallow layer, cm2
SAr_d=S*D2*RSD_d; % fine root surface area of deep layer, cm2
ln_s=log((1/pi/RLD_s)^(0.5)/0.012); % the radius of the soil cylinder and the fine root of shallow layer
ln_d=log((1/pi/RLD_d)^(0.5)/0.012); % the radius of the soil cylinder and the fine root of deep layer
rr_fineroot=2*10^5; % area-specific radial resistance of fine root, MPa s cm-1
%--------------------------------
% Calculate hydraulic diameter in each organ
d_root1=(N1*Droot^2)^0.5; % hydraulic diameter of shallow root, cm
d_root2=(N2*Droot^2)^0.5; % hydraulic diameter of deep root, cm
d_stem1=(N1*Dstem^2)^0.5; % hydraulic diameter of stem in shallow root part, cm
d_stem2=(N2*Dstem^2)^0.5; % hydraulic diameter of stem in deep root part, cm
d_branch1=(N1*Dbranch^2)^0.5; % hydraulic diameter of branch in shallow root part, cm
d_branch2=(N2*Dbranch^2)^0.5; % % hydraulic diameter of branch in shallow root part, cm
%--------------------------------
% Calculate axial resistance per unit length and area in each organ
ra_root=(128*(1.002*10^-9)*L/(pi*Droot^4)+rpit/(Droot*L*Fc*Fp))*pi*Droot^2/4/L; % specific resistance of root, MPa s cm-2
ra_stem=(128*(1.002*10^-9)*L/(pi*Dstem^4)+rpit/(Dstem*L*Fc*Fp))*pi*Dstem^2/4/L; % specific resistance of stem, MPa s cm-2
ra_branch=(128*(1.002*10^-9)*L/(pi*Dbranch^4)+rpit/(Dbranch*L*Fc*Fp))*pi*Dbranch^2/4/L; % specific resistance of branch, MPa s cm-2



%--------------------------------
% If select HA1, please select the following programs and DO NOT secelt the programs under HA2
R.Ra_r1=ra_root*4*L_root1/pi/d_root1^2; % axial resistance of shallow root, MPa s cm-3
R.Ra_r2=ra_root*4*L_root2/pi/d_root2^2; % axial resistance of deep root, MPa s cm-3
R.Ra_r1_2=ra_root*8*L_root2/pi/d_root1^2; % axial resistance of deep root from shallow root, MPa s cm-3
R.Ra_s1=ra_stem*8*L_stem/pi/d_stem1^2; % axial resistance of stem from shallow root, MPa s cm-3
R.Ra_s2=ra_stem*4*L_stem/pi/d_stem2^2; % axial resistance of stem from deep root, MPa s cm-3
R.Ra_b1=ra_branch*8*L_branch/pi/d_branch1^2; % axial resistance of branch from shallow root, MPa s cm-3
R.Ra_b2=ra_branch*4*L_branch/pi/d_branch2^2; % axial resistance of branch from deep root, MPa s cm-3
R.Rc_r=rpit/(L_root1*d_root1*Fc*Fp); % circulferential resistance in shallow root, MPa s cm-3
R.Rc_r2=rpit/(L_root2*d_root1*Fc*Fp);% circulferential resistance in deep root, MPa s cm-3
R.Rc_s=rpit/(L_stem*d_stem1*Fc*Fp); % circulferential resistance in stem, MPa s cm-3
R.Rr_s=rr_fineroot/SAr_s; % radial resistance of shallow root, MPa s cm-3
R.Rr_d=rr_fineroot/SAr_d; % radial resistance of shallow root, MPa s cm-3
Q=Ins_Flow_Soil_HA1(SWC1, SWC2, H1, H2, R, l_s, l_d, ln_s, ln_d, SWCr, SWCs, n, Ks, l, D1, D2); % Solve the instantaneous flow in the network
Qc=Q.p1; % circumferential flow within shallow root
H_int_s=Q.H_int_shallow; % root surface water potential of shallow layer, cm
H_int_d=Q.H_int_deep; % root surface water potential of deep layer, cm

%--------------------------------
% If select HA2, please select the following programs and DO NOT secelt the programs under HA1
% R.Ra_r1=ra_root*4*L_root1/pi/d_root1^2; % axial resistance of shallow root, MPa s cm-3
% R.Ra_r2=ra_root*4*L_root2/pi/d_root2^2; % axial resistance of deep root, MPa s cm-3
% R.Ra_s1=ra_stem*4*L_stem/pi/d_stem1^2; % axial resistance of deep root from shallow root, MPa s cm-3
% R.Ra_s2=ra_stem*4*L_stem/pi/d_stem2^2; % axial resistance of stem from shallow root, MPa s cm-3
% R.Ra_b1=ra_branch*4*L_branch/pi/d_branch1^2; % axial resistance of branch from shallow root, MPa s cm-3
% R.Ra_b2=ra_branch*4*L_branch/pi/d_branch2^2; % axial resistance of branch from deep root, MPa s cm-3
% R.Rc_r=rpit/(L_root1*d_root1*Fc*Fp); % circulferential resistance in shallow root, MPa s cm-3
% R.Rc_s1=rpit/(L_stem*d_stem1*Fc*Fp); % circulferential resistance in stem (left), MPa s cm-3
% R.Rc_s2=2*rpit/(pi*L_stem*d_stem1*Fc*Fp); % circulferential resistance in stem (right), MPa s cm-3
% R.Rr_s=rr_fineroot/SAr_s; % radial resistance of shallow root, MPa s cm-3
% R.Rr_d=rr_fineroot/SAr_d; % radial resistance of shallow root, MPa s cm-3
% Q=Ins_Flow_Soil_HA2(SWC1, SWC2, H1, H2, R, l_s, l_d, ln_s, ln_d, SWCr, SWCs, n, Ks, l, D1, D2); % Solve the instantaneous flow in the network
% Qc=Q.p1; % circumferential flow within shallow root
% H_int_s=Q.H_int_shallow; % root surface water potential of shallow layer, cm
% H_int_d=Q.H_int_deep; % root surface water potential of deep layer, cm

%--------------------------------
% Draw the relationship figure between soil water potential and Qc
subplot(3,1,1);
fig=surf(H1, H2, Qc)
fig.LineStyle='none'
% Draw the relationship figure between soil water content and Qc
subplot(3,1,2);
fig=surf(SWC1, SWC2, Qc)
fig.LineStyle='none'
% Draw the relationship figure between root surface water potential and Qc
subplot(3,1,3);
fig=surf(H_int_s, H_int_d, Qc)
fig.LineStyle='none'