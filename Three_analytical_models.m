%%%%%%
% Three_analytical_models: this script can be used to reproduce
% figure 7
% Copyright (C) 2025  Yiwei Chen, Pingchuan Dong and Youheng Zhang

% This is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version. 

% Please cite us if you use our routine: 
% Yiwei Chen, Pingchuan Dong and Youheng Zhang (2025). 
% An accurate analytical model and application method for squirt flow in anisotropic fractured rocks
% 
clear,clc
%%
%%
%%
%% Analytical model: torus
%%
%%
%%
%% background moduli: (dry) torus in grains. VTI Num 
Cij(1,1) = 87230044269;
Cij(1,2) = 6180414279;
Cij(1,3) = 7068404997;
Cij(2,2) = 87230005255;
Cij(2,3) = 7068426541;
Cij(3,3) = 84040612409;
Cij(4,4) = 38788214253;
Cij(5,5) = 38788240021;
Cij(6,6) = 40421562888;
     for i=1:6
        for j=i:6
      Cij(j,i) = Cij(i,j);
        end
     end
%% moduli: (dry)torus+(dry)crack in grains. VTI Num 
Cij_cr(1,1) = 86959886906;
Cij_cr(1,2) = 5950757722;
Cij_cr(1,3) = 4849479475;
Cij_cr(2,2) = 86957367522;
Cij_cr(2,3) = 4850230042;
Cij_cr(3,3) = 61468382429;
Cij_cr(4,4) = 29247397555;
Cij_cr(5,5) = 29248310700;
Cij_cr(6,6) = 40399653594;
     for i=1:6
        for j=i:6
      Cij_cr(j,i) = Cij_cr(i,j);
        end
     end
%% Zn Zt calculation
Sij     = inv(Cij);
Sij_cr  = inv(Cij_cr);
S_diff  = Sij - Sij_cr;
Z_N     = -S_diff(3,3) ;
Z_T     = -S_diff(4,4);
%%
L                    = 0.2;    % cube length (m)
Volume_L             = L.^3;   % volume (m)
a                    = 0.0005  * 4; % thickness of the crack (cylinder) (m)
b                    = 0.1;    % length of the half crack (m)
Stiff_Pore_volume    = pi*(b/2.5*0.6).^2 * (2*pi*(b + b/2.5*0.6))/4; % torus volume
Compl_pore_Volume    = 1   *  pi*b^2*a/4; % crack volume 
Porosity_Volume      = Stiff_Pore_volume + Compl_pore_Volume; 
Porosity1            = Porosity_Volume./Volume_L;   % total porosity (fraction)
phi_d                = Stiff_Pore_volume./Volume_L; % stiff(torus) porosity (fraction)
phi_c                = Compl_pore_Volume./Volume_L; % crack porosity (fraction)

b                    = 0.1 ;    % length of the half crack (m)
aspect_ratio         = a/(b*2) ; % aspect ratio of the crack (cylinder)
%%
K_g                  = 36.0e9;                   % Grain bulk modulus (Pa)
K_f                  = 4.3e9;                    % fluid bulk modulus (Pa)
eta                  = 1.414;                    % fluid viscosity (Pa * s)
%% validation Zn
Z_Nsat = Z_N./(1 +  Z_N./(phi_c.*(1./(4.3*1e9) - 1./36e9)) );

Cij_cr_Zn = Cij; Sij_cr_Zn=inv(Cij_cr_Zn);
Sij_cr_Zn(3,3) = Sij_cr_Zn(3,3) +1*Z_Nsat + 0*Z_N;
Sij_cr_Zn(4,4) = Sij_cr_Zn(4,4) +Z_T;
Sij_cr_Zn(5,5) = Sij_cr_Zn(5,5) +Z_T;
Cij_cr_Zn_fin  = inv(Sij_cr_Zn);

Leak_Z_N    = 1*(Sij_cr_Zn(3,3) - Sij(3,3));
Z_N         = Z_N - 1*Leak_Z_N;
Sij_b       = inv(Cij); % compliance of background 
%% Numeric
f               = 10.^(-5:0.1:8);  % frequency
w               = 2*pi*f;
%% Physics
ka              = 1./aspect_ratio .*( - 3.*1i.*w.*eta./K_f).^(1./2);
CC              = 2*pi.^2 *(b + b/2.5*0.6).*(b/2.5*0.6).^2/K_f; 
bess1           = besselj(0,(ka));
bess2           = besselj(1,(ka));
K_fSTAR_Mur     = ((1 - 2.*CC.*K_f.*bess2./ ( ((ka.*CC.*K_f).*bess1) + 2*pi*b^2 *a.*bess2  )).*K_f);

ka_low          = 1./aspect_ratio .*( - 3.*1i.*1.*eta./K_f).^(1./2); % ka at omega=1Hz.
K_fSTAR_appr    = K_f * (-(ka_low .^ 2 .* CC ^ 3 .* (ka_low .^ 2 + 6) .* K_f ^ 3) / 0.48e2 -...
    (a .* b .^ 2 * pi * CC .^ 2 .* (ka_low .^ 4 + 24 .* ka_low .^ 2 - 192) .* K_f ^ 2) / 0.192e3 +...
    (2 * CC * K_f * a .^ 2 * b ^ 4 * pi ^ 2) + (a ^ 3 * b ^ 6 * pi ^ 3)) /...
    ((pi * b ^ 2 * a + CC * K_f) .^ 3); % LF asymptot

K_fSTAR_appr    = pi*b^2*a*K_f/(a*b^2*pi + CC*K_f);   % valid for big torus
const           = K_fSTAR_appr;  %15221957.6976328 +      84240.3977312137i       % stiffness of the fluid in the low-frequency limit
C_11            = K_f;                  % high-frequency limit for fluid bulk modulus
C_00            = const;                % low-frequency limit for fluid bulk modulus asymptote

wc              = (4*sqrt(3)*aspect_ratio^2*sqrt(const)*sqrt(K_f))/(3*eta); % omega corresponding to 1/Q max of K_f Stokes eq.
%Ksi             = (128*C_11*aspect_ratio^6*C_00^2)/(27*eta^3*wc^3); % shape parameter for branching function
%taub            = Ksi*3/4/aspect_ratio^2/C_11*eta *1;               % time scale parameter for branching function

%Kf_branch4      = (  (  C_11 - ( (C_11-0*C_00) )./(  1 -Ksi + Ksi.*( 1 + 1i.*w.*taub./(Ksi.^2) ).^0.5  )  )  ); % branching function for Fluid bulk modulus K_f
%Kf_branch       = complex( real(Kf_branch4), imag(Kf_branch4)); % Fluid bulk modulus as a function of frequency
%% SHAPE fctor
cg=C_00; cg1=C_11; cg3=aspect_ratio;cg7=eta; nn=0.45;  cg9=nn;wbr=wc*1.0;
Ksi = (0.3e1 / 0.4e1) * cg7 ./ cg1 ./ (cg3 ^ 2) .* exp((cg9 .* log(wbr) - log(sin(cg9 * pi ./ 0.2e1)) - log(cg1) -...
    0.4e1 .* log(cg3) - log(cg) + log(cg9) + 0.2e1 .* log(cg7) - 0.4e1 .* log(0.2e1) + 0.2e1 .* log(0.3e1) + log(wbr)) ./ (cg9 - 0.1e1));
taub0 = Ksi*3/4/aspect_ratio^2/C_11*eta *nn; 

cg=C_00; cg1=C_11; cg3=aspect_ratio; cg5=eta; cg7=nn; cg9=wc;
Ksi = -0.3e1 / 0.8e1 * (-cg1 + cg) * cg5 * exp((-log(sin(cg7 * pi / 0.2e1)) + 0.2e1 * log(cg1 - cg) + (cg7 + 0.1e1) * log(cg9) - log(cg7) - log(cg) - 0.3e1 * log(cg1) + 0.2e1 * log(cg5) - 0.4e1 * log(cg3) - 0.6e1 * log(0.2e1) + 0.2e1 * log(0.3e1)) / (cg7 - 0.1e1)) / cg1 ^ 2 / cg3 ^ 2 / cg7;

cg=C_00; cg1=C_11; cg3=K_f;cg5=aspect_ratio;cg9=eta;cg11=nn; cg13 = Ksi;
taub = - (  0.3e1 / 0.8e1 * cg9 * cg13 * (cg ^ 2 - 2 * cg * cg3 + cg3 ^ 2) / (cg3 ^ 2) / (cg5 ^ 2) / cg11 / (-cg1 + cg)  );
%tau_sh = -3*(-C_11 + C_00)*eta*Ksi./(8*C_11^2*aspect_ratio.^2*nn);
%%
Kf_branch4      = (  (  C_11 - ( (C_11-0*C_00) )./(  1 -Ksi + Ksi.*( 1 + 1i.*w.*taub./(Ksi.^2) ).^(nn)  )  )  ); %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% branching function for Fluid bulk modulus K_f
Kf_branch       = complex( real(Kf_branch4), imag(Kf_branch4)); % Fluid bulk modulus as a function of frequency
%%
Z_N_mf_eval     = Z_N.*phi_c.*( K_g - Kf_branch)./((K_g - Kf_branch ).*phi_c + Kf_branch.*K_g.*Z_N);
Z_N_crack       = complex( real(Z_N_mf_eval), imag(Z_N_mf_eval));

del_Sij_crack   = zeros(6,6,size(w,2)); 
Sij_mf_blue          = zeros(6,6,size(w,2)); 
Cij_final       = zeros(6,6,size(w,2)); 
for nn = 1:size(w,2) % add crack into the saturated model (saturated model = saturated torus in grains)
    del_Sij_crack(3,3,nn)      = Z_N_crack(1,nn)  +  Leak_Z_N;
    del_Sij_crack(4,4,nn)      = Z_T;
    del_Sij_crack(5,5,nn)      = Z_T;
    Sij_mf_blue(:,:,nn)        = Sij_b + (del_Sij_crack(:,:,nn));  
    Cij_final(:,:,nn)          = inv(Sij_mf_blue(:,:,nn));        % Stiffness matrix of the torus analytical model, saturated rock, HTI      
end

%% Fluid saturation of the mod frame using anisotropic Gassmann eq.
Cij_mf_black      = complex(real(Cij_final), imag(Cij_final));
K_star            = zeros(1,1,size(w,2));
for nn = 1:size(w,2)
    for i=1:3
        for j=1:3
            K_star(1,1,nn)  = K_star(1,1,nn) + Cij_mf_black(i,j,nn);
        end
    end
end
K_star      = K_star./9;
alf_i       = zeros(6,1,size(w,2));
for i=1:3
    for nn = 1:size(w,2)
    alf_i(i,1,nn)           = 1 - (Cij_mf_black(i,1,nn) + Cij_mf_black(i,2,nn) + Cij_mf_black(i,3,nn))/3./K_g;
    end
end
alf_j       = zeros(1,6,size(w,2));
for j=1:3
    for nn = 1:size(w,2)
    alf_j(1,j,nn)           = 1 - (Cij_mf_black(1,j,nn) + Cij_mf_black(2,j,nn) + Cij_mf_black(3,j,nn))/3./K_g;
    end
end
for nn = 1:size(w,2)
    M(nn)                   = K_g./((1 - K_star(nn)./K_g) - (phi_d+0*phi_c)*(1 - K_g./K_f));
end
for nn = 1:size(w,2)
    Cij_mf_black(1:6,1:6,nn)= Cij_mf_black(1:6,1:6,nn) + (alf_i(:,:,nn).*alf_j(:,:,nn)).*M(1,nn);
end
Cij_final = complex( real(Cij_mf_black), imag(Cij_mf_black)); % Final stiffness matrix, saturated rock, HTI
%%
%% END OF the torus analytical model
%%
figure(1),clf
Ceff_33_black    = squeeze( Cij_final(3,3,:)    ./10^9); 
plot(f,             real(Ceff_33_black)   ,'-','color',[0.921875000000000,0.363281250000000,0.230468750000000],...
    'LineWidth',5); hold on;

figure(2), clf
plot(f,            imag(Ceff_33_black)./real(Ceff_33_black),'-','color',[0.921875000000000,0.363281250000000,0.230468750000000], ...
    'LineWidth',5); hold on;
%% 

Numerics_C33_tor = [ ...
                 69630856332.146 +      7676985.03989344i; ...
                69632073263.9805 +      11942816.8053862i; ...
                69633659434.2473 +      19098597.7723095i; ...
                69635934627.6087 +      31259389.0635758i; ...
                69639586730.3029 +      51890427.3525852i; ...
                69645931257.2315 +      86665646.2291922i; ...
                69657510128.8093 +      144742588.532453i; ...
                69679077528.7969 +      240702187.032405i; ...
                69719582351.7284 +      397246377.204404i; ...
                 69794794363.774 +      649804295.475527i; ...
                69939028887.8668 +      1058675837.81359i; ...
                70254216910.6059 +      1709126176.76162i; ...
                70986177373.1914 +       2614610794.1676i; ...
                72425950627.9346 +      3443738928.87272i; ...
                  74302859878.12 +      3550990503.01102i; ...
                75786224387.3327 +      2989875489.48735i; ...
                76725839264.4886 +      2424251122.16542i; ...
                77476387249.3163 +      2033776965.95028i; ...
                78128123973.3895 +      1676268811.37579i; ...
                78637432453.7883 +      1353174563.65816i; ...
                79036589442.2584 +        1085801791.054i; ...
                79348131730.8276 +      867946347.555462i; ...
                79586890003.7773 +      697978913.601367i; ...
                79765636809.2162 +       579090091.26814i; ...
                79897049336.4423 +      521053089.511857i; ...
                79994982367.6358 +      543418922.852606i; ...
                80078228785.7404 +      682401711.215658i; ...
                 80182466447.797 +      1000144835.52898i; ...
                80391259058.1397 +      1588020513.54976i; ...
     ];
 
freqOUT_tor = [ ...
                             10; ...
               17.7827941003892; ...
               31.6227766016838; ...
               56.2341325190349; ...
                            100; ...
               177.827941003892; ...
               316.227766016838; ...
               562.341325190349; ...
                           1000; ...
               1778.27941003892; ...
               3162.27766016838; ...
               5623.41325190349; ...
                          10000; ...
               17782.7941003892; ...
               31622.7766016838; ...
               56234.1325190349; ...
                         100000; ...
               177827.941003892; ...
               316227.766016838; ...
               562341.325190349; ...
                        1000000; ...
               1778279.41003892; ...
               3162277.66016838; ...
               5623413.25190349; ...
                       10000000; ...
               17782794.1003892; ...
               31622776.6016838; ...
               56234132.5190349; ...
                      100000000; ...
     ];

SZ = 18;     
figure(1)%,clf
plot(freqOUT_tor(1:end-1),real(Numerics_C33_tor(1:end-1))./10^9,'s','color',[0.222656250000000,0.316406250000000,0.632812500000000],...
    ...
    'MarkerFaceColor','none','MarkerSize',SZ,...
    'LineWidth',2.9); hold on
figure(2)%,clf
plot(freqOUT_tor(1:end-1),imag(Numerics_C33_tor(1:end-1)) ./ real(Numerics_C33_tor(1:end-1)),'s','color',[0.222656250000000,0.316406250000000,0.632812500000000],...
    'MarkerFaceColor','none','MarkerSize',SZ,...
    'LineWidth',2.9); hold on

%%
%%
%%
%% The sphere analytical model: 1 sphere
%%
%%
%%
%% background moduli: (dry) sphere in grains. VTI Num 
Cij(1,1) = 87230044269*0 + 92*1e9;
Cij(1,2) = 6180414279 *0 + 6.66*1e9;
Cij(1,3) = 7068404997 *0 + 6.66*1e9;
Cij(2,2) = 87230005255*0 + 92*1e9;
Cij(2,3) = 7068426541 *0 + 6.66*1e9;
Cij(3,3) = 84040612409*0 + 92*1e9;
Cij(4,4) = 42788214253;
Cij(5,5) = 42788240021;
Cij(6,6) = 42421562888;
     for i=1:6
        for j=i:6
      Cij(j,i) = Cij(i,j);
        end
     end
%% moduli: (dry)sphere+(dry)crack in grains. VTI Num full 
Cij_cr(1,1) = 86959886906*0 + 91.98*1e9;
Cij_cr(1,2) = 5950757722*0 + 5.9*1e9;
Cij_cr(1,3) = 4849479475*0 + 6.62*1e9;
Cij_cr(2,2) = 86957367522*0 + 91.98*1e9;
Cij_cr(2,3) = 4850230042*0 + 6.62*1e9;
Cij_cr(3,3) = 61468382429*0 + 79.04*1e9;
Cij_cr(4,4) = 29247397555;
Cij_cr(5,5) = 29248310700;
Cij_cr(6,6) = 40399653594;
     for i=1:6
        for j=i:6
      Cij_cr(j,i) = Cij_cr(i,j);
        end
     end
%% Zn Zt calculation
Sij     = inv(Cij);
Sij_cr  = inv(Cij_cr);
S_diff  = Sij - Sij_cr;
Z_N     = -S_diff(3,3);
Z_T     = -S_diff(4,4);

Z_N_save = Z_N;
Z_T_save = Z_T;
%%
L                    = 0.2;    % cube length (m)
Volume_L             = 0.2*2 *0.2*2 *0.2;   % volume (m)
a                    = 0.0005  * 4; % thickness of the crack (cylinder) (m)
b                    = 0.1;    % length of the half crack (m)

Stiff_Pore_volume    = 4/3*pi*(b/2)^3  ;
Compl_pore_Volume    = pi*b^2*a; % crack volume 
Porosity_Volume      = Stiff_Pore_volume + Compl_pore_Volume; 
Porosity1            = Porosity_Volume./Volume_L;   % total porosity (fraction)
phi_d                = Stiff_Pore_volume./Volume_L; % stiff(torus) porosity (fraction)
phi_c                = 1*Compl_pore_Volume./Volume_L; % crack porosity (fraction)
b                    = 0.1  ;    % length of the half crack (m)

aspect_ratio         = a/(b*2 ) /2.0;
%%
K_g                  = 36.0e9;       % Grain bulk modulus (Pa)
K_f                  = 4.3e9;        % fluid bulk modulus (Pa)
eta                  = 1.414;        % fluid viscosity (Pa * s)
%% validation Zn
Z_Nsat = Z_N./(1 +  Z_N./(phi_c.*(1./(4.3*1e9) - 1./36e9)) );

Cij_cr_Zn = Cij; Sij_cr_Zn=inv(Cij_cr_Zn);
Sij_cr_Zn(3,3) = Sij_cr_Zn(3,3) +1*Z_Nsat + 0*Z_N;
Sij_cr_Zn(4,4) = Sij_cr_Zn(4,4) +Z_T;
Sij_cr_Zn(5,5) = Sij_cr_Zn(5,5) +Z_T;
Cij_cr_Zn_fin  = inv(Sij_cr_Zn);

Leak_Z_N    = 1*0.5*(Sij_cr_Zn(3,3) - Sij(3,3));
Z_N         = Z_N - 1*Leak_Z_N;
%%
Sij_b           = inv(Cij); % compliance of background 

%% Numeric
f               = 10.^(-5:0.1:8);  % frequency
w               = 2*pi*f;
%% Physics
Gw = 1i.*w*eta;
abL = 1/aspect_ratio .*sqrt(3*Gw ./ (C_11 + 4/3*Gw));
K_layr  = C_11 + 4/3.*Gw - ( C_11 - 2/3.*Gw).^2./(C_11 + 4/3*Gw) .* tanh(abL) ./abL;
Q_layr = imag(K_layr)./real(K_layr);
[max1,max2] = max(Q_layr);
Kf_branch   = K_layr;
%%
Kf_branch = Kf_branch.*1 + 0*4.3*1e9;
Z_N_mf_eval     = Z_N./(1 +  Z_N./(phi_c.*(1./Kf_branch - 1./36e9)) ); %1./K_g  K_starNew
Z_N_crack       = complex( real(Z_N_mf_eval), imag(Z_N_mf_eval));

del_Sij_crack   = zeros(6,6,size(w,2)); 
Sij_mf_blue     = zeros(6,6,size(w,2)); 
Cij_final       = zeros(6,6,size(w,2)); 
for nn = 1:size(w,2) % add crack into the saturated model (saturated model = saturated sphere in grains)
    del_Sij_crack(3,3,nn)      = Z_N_crack(1,nn)  +  Leak_Z_N;
    del_Sij_crack(4,4,nn)      = Z_T;
    del_Sij_crack(5,5,nn)      = Z_T;
    Sij_mf_blue(:,:,nn)        = Sij_b + (del_Sij_crack(:,:,nn));  
    Cij_final(:,:,nn)          = inv(Sij_mf_blue(:,:,nn));        % Stiffness matrix of the sphere analytical model, saturated rock, HTI      
end
%% Fluid saturation of the mod frame using anisotropic Gassmann eq.
Cij_mf_black      = complex(real(Cij_final), imag(Cij_final));
K_star            = zeros(1,1,size(w,2));
for nn = 1:size(w,2)
    for i=1:3
        for j=1:3
            K_star(1,1,nn)  = K_star(1,1,nn) + Cij_mf_black(i,j,nn);
        end
    end
end
K_star      = K_star./9;
alf_i       = zeros(6,1,size(w,2));
for i=1:3
    for nn = 1:size(w,2)
    alf_i(i,1,nn)           = 1 - (Cij_mf_black(i,1,nn) + Cij_mf_black(i,2,nn) + Cij_mf_black(i,3,nn))/3./K_g;
    end
end
alf_j       = zeros(1,6,size(w,2));
for j=1:3
    for nn = 1:size(w,2)
    alf_j(1,j,nn)           = 1 - (Cij_mf_black(1,j,nn) + Cij_mf_black(2,j,nn) + Cij_mf_black(3,j,nn))/3./K_g;
    end
end
for nn = 1:size(w,2)
    M(nn)                   = K_g./((1 - K_star(nn)./K_g) - (phi_d+0*phi_c)*(1 - K_g./K_f));
end
for nn = 1:size(w,2)
    Cij_mf_black(1:6,1:6,nn)= Cij_mf_black(1:6,1:6,nn) + (alf_i(:,:,nn).*alf_j(:,:,nn)).*M(1,nn);
end
Cij_final = complex( real(Cij_mf_black), imag(Cij_mf_black)); % Final stiffness matrix, saturated rock, HTI
%%
%% END OF the sphere analytical model
%%
figure(1)%,clf
Ceff_33_black    = squeeze( Cij_final(3,3,:)    ./10^9); 
plot(f,             real(Ceff_33_black)   ,'-','color',[0.988281250000000,0.722656250000000,0.417968750000000],...
    'LineWidth',5); hold on;

figure(2)%, clf
plot(f,            imag(Ceff_33_black)./real(Ceff_33_black)*1.0,'-','color',[0.988281250000000,0.722656250000000,0.417968750000000],...
    'LineWidth',5); hold on;
%% Numerical solution; 1 sphere
Numerics_C33_1ball = [ ...
                83504350913.4802 +      39995269.7538538i; ...
                83507383241.9635 +      68556577.2318278i; ...
                83513584399.8542 +      118071019.322557i; ...
                83527579429.6496 +      202916820.164106i; ...
                83560372952.4877 +      345443280.658381i; ...
                83637221719.4238 +      576591123.220246i; ...
                 83809806862.531 +      928408330.45429i; ...
                84168218040.1114 +      1412349789.06946i; ...
                84849288471.6974 +      1963444008.21717i; ...
                85935052920.3421 +      2318113835.04279i; ...
                87136632144.6463 +      2155509369.78648i; ...
                87984739558.8995 +      1614015325.8535i; ...
                88416006277.5426 +      1082868820.33454i; ...
                88626256873.3499 +      722466324.667097i; ...
                88759107058.9091 +      498774453.297965i; ...
                88856380871.2203 +      350583307.825787i; ...
                88926960708.6066 +      248996999.113108i; ...
                88977901033.7602 +      180139801.810242i; ...
                89016964796.8729 +      134524802.268338i; ...
                89048259535.2811 +      103222195.815278i; ...
                89072678325.5885 +      81838434.3480703i; ...
                89091446782.2818 +      68957551.584513i; ...
                89105650657.4134 +      64838705.8383689i; ...
                89116292979.0515 +      71923639.3430778i; ...
                89124634035.4752 +      95713053.1190183i; ...
                89132861608.4126 +      146598760.257286i; ...
                89146180446.5334 +      243146262.374364i; ...
                89179215976.5726 +      416679876.076046i; ...
                89273766164.6718 +      712501282.35376i; ...
     ];



freqOUT_1ball = [ ...
                             10; ...
               17.7827941003892; ...
               31.6227766016838; ...
               56.2341325190349; ...
                            100; ...
               177.827941003892; ...
               316.227766016838; ...
               562.341325190349; ...
                           1000; ...
               1778.27941003892; ...
               3162.27766016838; ...
               5623.41325190349; ...
                          10000; ...
               17782.7941003892; ...
               31622.7766016838; ...
               56234.1325190349; ...
                         100000; ...
               177827.941003892; ...
               316227.766016838; ...
               562341.325190349; ...
                        1000000; ...
               1778279.41003892; ...
               3162277.66016838; ...
               5623413.25190349; ...
                       10000000; ...
               17782794.1003892; ...
               31622776.6016838; ...
               56234132.5190349; ...
                      100000000; ...
     ];

figure(1)
plot(freqOUT_1ball(1:end-1),real(Numerics_C33_1ball(1:end-1))./10^9,'s','color',[0.921875000000000,0.363281250000000,0.230468750000000],...
    'MarkerFaceColor','none','MarkerSize',SZ,...
    'LineWidth',2.9); hold on


figure(2)
plot(freqOUT_1ball(1:end-1),imag(Numerics_C33_1ball(1:end-1)) ./ real(Numerics_C33_1ball(1:end-1)),...
    's','color',[0.921875000000000,0.363281250000000,0.230468750000000],...
    'MarkerFaceColor','none','MarkerSize',SZ,...
    'LineWidth',2.9); hold on


%%
%%
%%
%% Present analytical model: 2 cracks
%% 
%%
%%
%% background moduli: (dry) stiff crack in grains. VTI Num
Cij(1,1) = 87230044269*0 + 92380000000;
Cij(1,2) = 6180414279 *0 + 6524000000;
Cij(1,3) = 7068404997 *0 + 5994000000;
Cij(2,2) = 87230005255*0 + 92380000000;
Cij(2,3) = 7068426541 *0 + 5994000000;
Cij(3,3) = 84040612409*0 + 81084000000;
Cij(4,4) = 38904000000;
Cij(5,5) = 38904000000;
Cij(6,6) = 42928000000;
     for i=1:6
        for j=i:6
      Cij(j,i) = Cij(i,j);
        end
     end
%% moduli: (dry) stiff crack + (dry) compliant crack in grains. VTI Num
Cij_cr(1,1) = 86959886906*0 + 92164000000;
Cij_cr(1,2) = 5950757722*0 + 6460000000;
Cij_cr(1,3) = 4849479475*0 + 5060000000;
Cij_cr(2,2) = 86957367522*0 + 92164000000;
Cij_cr(2,3) = 4850230042*0 + 5060000000;
Cij_cr(3,3) = 61468382429*0 + 68086000000;
Cij_cr(4,4) = 34766000000;
Cij_cr(5,5) = 34766000000;
Cij_cr(6,6) = 42852000000;
     for i=1:6
        for j=i:6
      Cij_cr(j,i) = Cij_cr(i,j);
        end
     end
%% Zn Zt calculation
Sij     = inv(Cij);
Sij_cr  = inv(Cij_cr);
S_diff  = Sij - Sij_cr;
Z_N     = -S_diff(3,3);
Z_T     = -S_diff(4,4);

Z_N_save = Z_N;
Z_T_save = Z_T;
%%
L                    = 0.2;    % cube length (m)
Volume_L             = 0.2*2 *0.2*2 *0.2;   % volume (m)
a                    = 0.0005  * 4; % thickness of the crack (cylinder) (m)
b                    = 0.1;    % length of the half crack (m)

h                    =0.02; % thickness of the stiff pore (cylinder) (m)
Stiff_Pore_volume    = 1   *  pi*b^2*h/4;

Compl_pore_Volume    =   pi*b^2*a/4; % crack volume 
Porosity_Volume      = Stiff_Pore_volume + Compl_pore_Volume; 
Porosity1            = Porosity_Volume./Volume_L;   % total porosity (fraction)
phi_d                = Stiff_Pore_volume./Volume_L; % stiff(torus) porosity (fraction)
phi_c                = 1*Compl_pore_Volume./Volume_L; % crack porosity (fraction)

aspect_ratio         = a/(b*2 ) /2.0;
%%
K_g                  = 36.0e9;         % Grain bulk modulus (Pa)
K_f                  = 4.3e9;          % fluid bulk modulus (Pa)
eta                  = 1.414;          % fluid viscosity (Pa * s)
%% validation Zn
Z_Nsat = Z_N./(1 +  Z_N./(phi_c.*(1./(4.3*1e9) - 1./36e9)) );

Cij_cr_Zn = Cij; Sij_cr_Zn=inv(Cij_cr_Zn);
Sij_cr_Zn(3,3) = Sij_cr_Zn(3,3) +1*Z_Nsat + 0*Z_N;
Sij_cr_Zn(4,4) = Sij_cr_Zn(4,4) +Z_T;
Sij_cr_Zn(5,5) = Sij_cr_Zn(5,5) +Z_T;
Cij_cr_Zn_fin  = inv(Sij_cr_Zn);

Leak_Z_N    = 2*0.5*(Sij_cr_Zn(3,3) - Sij(3,3));
Z_N         = Z_N - 1*Leak_Z_N;
%%
Sij_b           = inv(Cij); % compliance of background 
%% Numeric
f               = 10.^(-5:0.1:8);  % frequency
w               = 2*pi*f;
%% Physics
Gw = 1i.*w*eta;
abL = 1/aspect_ratio .*sqrt(3*Gw ./ (K_f + 4/3*Gw));
K_layr  = K_f + 4/3.*Gw - ( K_f - 2/3.*Gw).^2./(K_f + 4/3*Gw) .* tanh(abL) ./abL;
Q_layr = imag(K_layr)./real(K_layr);
%[max1,max2] = max(Q_layr);

Kf_branch = K_layr;
%%
Kf_branch = Kf_branch.*1 + 0*4.3*1e9;
Z_N_mf_eval     = Z_N./(1 +  Z_N./(phi_c.*(1./Kf_branch - 1./36e9)) ); %1./K_g  K_starNew
Z_N_crack       = complex( real(Z_N_mf_eval), imag(Z_N_mf_eval));

del_Sij_crack   = zeros(6,6,size(w,2)); 
Sij_mf_blue          = zeros(6,6,size(w,2)); 
Cij_final       = zeros(6,6,size(w,2)); 
for nn = 1:size(w,2) % add crack into the saturated model (saturated model = saturated stiff crack in grains)
    del_Sij_crack(3,3,nn)      = Z_N_crack(1,nn)  +  Leak_Z_N;
    del_Sij_crack(4,4,nn)      = Z_T;
    del_Sij_crack(5,5,nn)      = Z_T;
    Sij_mf_blue(:,:,nn)        = Sij_b + (del_Sij_crack(:,:,nn));  
    Cij_final(:,:,nn)          = inv(Sij_mf_blue(:,:,nn));        % Stiffness matrix of the present analytical model, saturated rock, HTI      
end
%% Fluid saturation of the mod frame using anisotropic Gassmann eq.
Cij_mf_black      = complex(real(Cij_final), imag(Cij_final));
K_star            = zeros(1,1,size(w,2));
for nn = 1:size(w,2)
    for i=1:3
        for j=1:3
            K_star(1,1,nn)  = K_star(1,1,nn) + Cij_mf_black(i,j,nn);
        end
    end
end
K_star      = K_star./9;
alf_i       = zeros(6,1,size(w,2));
for i=1:3
    for nn = 1:size(w,2)
    alf_i(i,1,nn)           = 1 - (Cij_mf_black(i,1,nn) + Cij_mf_black(i,2,nn) + Cij_mf_black(i,3,nn))/3./K_g;
    end
end
alf_j       = zeros(1,6,size(w,2));
for j=1:3
    for nn = 1:size(w,2)
    alf_j(1,j,nn)           = 1 - (Cij_mf_black(1,j,nn) + Cij_mf_black(2,j,nn) + Cij_mf_black(3,j,nn))/3./K_g;
    end
end
for nn = 1:size(w,2)
    M(nn)                   = K_g./((1 - K_star(nn)./K_g) - (phi_d+0*phi_c)*(1 - K_g./K_f));
end
for nn = 1:size(w,2)
    Cij_mf_black(1:6,1:6,nn)= Cij_mf_black(1:6,1:6,nn) + (alf_i(:,:,nn).*alf_j(:,:,nn)).*M(1,nn);
end
Cij_final = complex( real(Cij_mf_black), imag(Cij_mf_black)); % Final stiffness matrix, saturated rock, HTI
%%
%% END OF Present analytical model: stiff crack
%%
figure(1)
Ceff_33_black    = squeeze( Cij_final(3,3,:)    ./10^9);
plot(f,             real(Ceff_33_black)   ,'-','color',[0.445312500000000,0.664062500000000,0.808593750000000],...
    'LineWidth',5); hold on;

figure(2)
plot(f,            imag(Ceff_33_black)./real(Ceff_33_black),'-','color',[0.445312500000000,0.664062500000000,0.808593750000000],...
    'LineWidth',5); hold on;
%% Numerical solution; stiff crack
Numerics_C33_1ball = [
    86027365873.1785+24725673.54048i ; ...
86035473475.3765+43969126.75122i ; ...
86037346546.7684+78189359.74179i ; ...
86047685677.4699+139042525.78788i ; ...
86063457567.5722+247256465.42949i ; ...
86135767564.2943+439248743.58585i ; ...
86244234344.3791+682893575.91972i ; ...
86424637457.8562+862425259.94922i ; ...
86845475735.1326+855564608.69649i ; ...
87201829345.5479+663911312.93157i ; ...
87343163436.5982+405693556.49625i ; ...
87435876946.2257+237948736.22685i ; ...
87502575675.4324+143428262.55747i ; ...
87542893475.2852+81488277.70182i ; ...
87592897468.6481+50963561.54472i ; ...
87647812359.1534+34221210.49545i ; ...
87648273469.3544+25030544.73711i ; ...
87648476893.359+18218116.07319i ; ...
87645865948.4463+13149527.72034i ; ...
87647832178.2786+9565209.47712i ; ...
87649046231.4706+7308063.06285i ; ...
87647895829.4269+6413813.53041i ; ...
87645742459.2876+5752720.86807i ; ...
87641878956.3809+6427409.48256i ; ...
87643241740.5587+8063309.87814i ; ...
87645632897.2199+12034330.69245i ; ...
87644533895.1762+20273506.62942i ; ...
87643547155.2263+33446735.96715i ; ...
87649964562.5792+49624359.33015i ; ...

];

freqOUT_1ball = [ ...
                             10; ...
               17.7827941003892; ...
               31.6227766016838; ...
               56.2341325190349; ...
                            100; ...
               177.827941003892; ...
               316.227766016838; ...
               562.341325190349; ...
                           1000; ...
               1778.27941003892; ...
               3162.27766016838; ...
               5623.41325190349; ...
                          10000; ...
               17782.7941003892; ...
               31622.7766016838; ...
               56234.1325190349; ...
                         100000; ...
               177827.941003892; ...
               316227.766016838; ...
               562341.325190349; ...
                        1000000; ...
               1778279.41003892; ...
               3162277.66016838; ...
               5623413.25190349; ...
                       10000000; ...
               17782794.1003892; ...
               31622776.6016838; ...
               56234132.5190349; ...
                      100000000; ...
     ];

figure(1)
plot(freqOUT_1ball(1:end-1),real(Numerics_C33_1ball(1:end-1))./10^9,'s','color',[0.656250000000000,0.011718750000000,0.148437500000000],...
    'MarkerFaceColor','none','MarkerSize',SZ,...
    'LineWidth',2.9); hold on
set(gca, 'XScale', 'log') ; grid on; 
xlim([1 10^6.7]);ylim([67 92]); 

 legend('Analytical model, Big torus','Numerical solution, Big torus',...
     'Analytical model, 1 sphere','Numerical solution, 1 sphere',...
          'Analytical model, 1 crack','Numerical solution, 1 crack','Location','west');
 
xlabel('  Frequency (Hz)','fontsize',12); ylabel(' Re(C33) (GPa)','fontsize',12);
x_range = [1e0 1e1 1e2 1e3 1e4 1e5 1e6];
set(gca,'xtick',x_range);

figure(2)
plot(freqOUT_1ball(1:end-1),imag(Numerics_C33_1ball(1:end-1)) ./ real(Numerics_C33_1ball(1:end-1)),...
    's','color',[0.656250000000000,0.011718750000000,0.148437500000000],...
    'MarkerFaceColor','none','MarkerSize',SZ,...
    'LineWidth',2.9); hold on
set(gca, 'XScale', 'log') ; grid on; 
set(gca, 'YScale', 'log') ; 
xlim([1 10^6.7]); ylim([5e-6 1e-1]);
    xlabel('  Frequency (Hz)','fontsize',12); ylabel('  1/Q','fontsize',12); hold on;
x_range = [1e-1 1e0 1e1 1e2 1e3 1e4 1e5 1e6];
set(gca,'xtick',x_range); 
%% The high-frequency asymptote
plot(freqOUT_1ball(17:end-0)*1.5,4.5*1./(freqOUT_1ball(17:end-0).^(4/10) ) ,'-','color',[0 0 0],...
    'LineWidth',1.5); hold on; 