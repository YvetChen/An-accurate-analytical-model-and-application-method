%%%%%%
% Thicknesses_of_analytical_models: this script can be used to reproduce
% figure 8
% Copyright (C) 2025  Yiwei Chen, Pingchuan Dong and Youheng Zhang

% This is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version. 

% Please cite us if you use our version: 
% Yiwei Chen, Pingchuan Dong and Youheng Zhang (2025). 
% An accurate analytical model and application method for squirt flow in anisotropic fractured rocks
% 
% This version is modified from the following script:
% ALKHIMENKOV Y, QUINTAL B. (2022)
% An accurate analytical model for squirt flow in anisotropic porous rocks Part 2: Complex geometry
%
clear,clc
%%
%%
%%
%% Present analytical model: hs=0.004 (thickness of the stiff crack)
%% 
%%
%%
%% background moduli: (dry) stiff crack in grains. VTI Num full only stiff
Cij(1,1) = 94224000000;
Cij(1,2) = 6592000000;
Cij(1,3) = 5876000000;
Cij(2,2) = 94224000000;
Cij(2,3) = 5876000000;
Cij(3,3) = 82136000000;
Cij(4,4) = 40214000000;
Cij(5,5) = 40214000000;
Cij(6,6) = 43816000000;
     for i=1:6
        for j=i:6
      Cij(j,i) = Cij(i,j);
        end
     end
%% moduli: (dry) stiff crack + (dry) compliant crack in grains. VTI Num full  
Cij_cr(1,1) = 93980000000;
Cij_cr(1,2) = 6516000000;
Cij_cr(1,3) = 4978000000;
Cij_cr(2,2) = 93980000000;
Cij_cr(2,3) = 4978000000;
Cij_cr(3,3) = 69016000000;
Cij_cr(4,4) = 35990000000;
Cij_cr(5,5) = 35990000000;
Cij_cr(6,6) = 43732000000;
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
Volume_L             = 4*L.^3;   % volume (m)
a                    = 0.0005  * 4; % thickness of the crack (cylinder) (m)
b                    = 0.1;    % length of the half crack (m)

h                    =0.004; % thickness of the stiff pore (cylinder) (m)
Stiff_Pore_volume    =1   *  pi*b^2*h/4;
Compl_pore_Volume    = 1   *  pi*b^2*a/4; % crack volume 
Porosity_Volume      = Stiff_Pore_volume + Compl_pore_Volume; 
Porosity1            = Porosity_Volume./Volume_L;   % total porosity (fraction)
phi_d                = Stiff_Pore_volume./Volume_L; % stiff porosity (fraction)
phi_c                = Compl_pore_Volume./Volume_L; % crack porosity (fraction)

aspect_ratio         =  a/(b*2 ) /2.0 ; % aspect ratio of the crack (cylinder)
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

const           = Compl_pore_Volume*K_f/(Compl_pore_Volume + Stiff_Pore_volume);  %15221957.6976328 +      84240.3977312137i       % stiffness of the fluid in the low-frequency limit
C_11            = K_f;                  % high-frequency limit for fluid bulk modulus
C_00            = const;                % low-frequency limit for fluid bulk modulus asymptote

%% Physics
Gw = 1i.*w*eta;
abL = 1/aspect_ratio .*sqrt(3*Gw ./ (C_11 + 4/3*Gw));
K_layr  = C_11 + 4/3.*Gw - ( C_11 - 2/3.*Gw).^2./(C_11 + 4/3*Gw) .* tanh(abL) ./abL;
Q_layr = imag(K_layr)./real(K_layr);
[max1,max2] = max(Q_layr);

Kf_branch = K_layr;
Kf_branch = Kf_branch.*1 + 0*4.3*1e9;
%%
Z_N_mf_eval     = Z_N./(1 +  Z_N./(phi_c.*(1./Kf_branch - 1./36e9)) ); %1./K_g  K_starNew
Z_N_crack       = complex( real(Z_N_mf_eval), imag(Z_N_mf_eval));

del_Sij_crack   = zeros(6,6,size(w,2)); 
Sij_mf_black     = zeros(6,6,size(w,2)); 
Cij_final       = zeros(6,6,size(w,2)); 
for nn = 1:size(w,2) % add crack into the saturated model (saturated model = saturated the stiff crack in grains)
    del_Sij_crack(3,3,nn)      = Z_N_crack(1,nn)  +  Leak_Z_N;
    del_Sij_crack(4,4,nn)      = Z_T;
    del_Sij_crack(5,5,nn)      = Z_T;
    Sij_mf_black(:,:,nn)        = Sij_b + (del_Sij_crack(:,:,nn));  
    Cij_final(:,:,nn)          = inv(Sij_mf_black(:,:,nn));        % Stiffness matrix of the present analytical model, saturated rock, HTI      
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
%% END OF the present analytical model: hs=0.004 (thickness of the stiff crack)
%%
figure(1),clf
Ceff_33_black    = squeeze( Cij_final(3,3,:)    ./10^9); 
plot(f,             real(Ceff_33_black)   ,'-','color',[0.222656250000000,0.316406250000000,0.632812500000000],...
    'LineWidth',5); hold on;

figure(2), clf
plot(f,            imag(Ceff_33_black)./real(Ceff_33_black),'-','color',[0.222656250000000,0.316406250000000,0.632812500000000], ...
    'LineWidth',5); hold on;


%%
%%
%%
%% Present analytical model: hs=0.005 (thickness of the stiff crack)
%% 
%%
%%
%% background moduli: (dry) stiff crack in grains. VTI Num full
Cij(1,1) = 87230044269*0 + 94124000000;
Cij(1,2) = 6180414279 *0 + 6588000000;
Cij(1,3) = 7068404997 *0 + 5886000000;
Cij(2,2) = 87230005255*0 + 94124000000;
Cij(2,3) = 7068426541 *0 + 5886000000;
Cij(3,3) = 84040612409*0 + 82056000000;
Cij(4,4) = 40116000000;
Cij(5,5) = 40116000000;
Cij(6,6) = 43768000000;
     for i=1:6
        for j=i:6
      Cij(j,i) = Cij(i,j);
        end
     end
%% moduli: (dry) stiff crack + (dry) compliant crack in grains. VTI Num full
Cij_cr(1,1) = 86959886906*0 + 93884000000;
Cij_cr(1,2) = 5950757722*0 + 6508000000;
Cij_cr(1,3) = 4849479475*0 + 4992000000;
Cij_cr(2,2) = 86957367522*0 + 93884000000;
Cij_cr(2,3) = 4850230042*0 + 4992000000;
Cij_cr(3,3) = 61468382429*0 + 69068000000;
Cij_cr(4,4) = 35952000000;
Cij_cr(5,5) = 35952000000;
Cij_cr(6,6) = 43688000000;
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
Volume_L             = 4*L.^3;   % volume (m)
a                    = 0.0005  * 4; % thickness of the crack (cylinder) (m)
b                    = 0.1;    % length of the half crack (m)
h                    =0.005; % thickness of the stiff pore (cylinder) (m)

Stiff_Pore_volume    = 1   *  pi*b^2*h/4;
Compl_pore_Volume    = 1   *  pi*b^2*a/4; % crack volume 
Porosity_Volume      = Stiff_Pore_volume + Compl_pore_Volume; 
Porosity1            = Porosity_Volume./Volume_L;   % total porosity (fraction)
phi_d                = Stiff_Pore_volume./Volume_L; % stiff porosity (fraction)
phi_c                = Compl_pore_Volume./Volume_L; % crack porosity (fraction)
aspect_ratio         = a/(b*2 ) /2.0; % aspect ratio of the crack (cylinder)
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

Kf_branch = K_layr;
%%
Kf_branch = Kf_branch.*1 + 0*4.3*1e9;
Z_N_mf_eval     = Z_N./(1 +  Z_N./(phi_c.*(1./Kf_branch - 1./36e9)) ); %1./K_g  K_starNew
Z_N_crack       = complex( real(Z_N_mf_eval), imag(Z_N_mf_eval));

del_Sij_crack   = zeros(6,6,size(w,2)); 
Sij_mf_blue          = zeros(6,6,size(w,2)); 
Cij_final       = zeros(6,6,size(w,2)); 
for nn = 1:size(w,2) % add crack into the saturated model (saturated model = saturated the stiff crack in grains)
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
%% END OF the present analytical model: hs=0.005 (thickness of the stiff crack)
%%
figure(1)
Ceff_33_black    = squeeze( Cij_final(3,3,:)    ./10^9);
plot(f,             real(Ceff_33_black)   ,'-','color',[0.445312500000000,0.664062500000000,0.808593750000000],...
    'LineWidth',5); hold on;

figure(2)
plot(f,            imag(Ceff_33_black)./real(Ceff_33_black),'-','color',[0.445312500000000,0.664062500000000,0.808593750000000],...
    'LineWidth',5); hold on;

%%
%%
%%
%% Present analytical model: hs=0.008 (thickness of the stiff crack)
%% 
%%
%%
%% background moduli: (dry) stiff crack in grains. VTI Num full
Cij(1,1) = 87230044269*0 + 93808000000;
Cij(1,2) = 6180414279 *0 + 6576000000;
Cij(1,3) = 7068404997 *0 + 5924000000;
Cij(2,2) = 87230005255*0 + 93808000000;
Cij(2,3) = 7068426541 *0 + 5924000000;
Cij(3,3) = 84040612409*0 + 81900000000;
Cij(4,4) = 39880000000;
Cij(5,5) = 39880000000;
Cij(6,6) = 43616000000;
     for i=1:6
        for j=i:6
      Cij(j,i) = Cij(i,j);
        end
     end
%% moduli: (dry) stiff crack + (dry) compliant crack in grains. VTI Num full
Cij_cr(1,1) = 86959886906*0 + 93980000000;
Cij_cr(1,2) = 5950757722*0 + 6516000000;
Cij_cr(1,3) = 4849479475*0 + 4978000000;
Cij_cr(2,2) = 86957367522*0 + 93980000000;
Cij_cr(2,3) = 4850230042*0 + 4978000000;
Cij_cr(3,3) = 61468382429*0 + 69016000000;
Cij_cr(4,4) = 35990000000;
Cij_cr(5,5) = 35990000000;
Cij_cr(6,6) = 43732000000;
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
Volume_L             = 4*L.^3;   % volume (m)
a                    = 0.0005  * 4; % thickness of the crack (cylinder) (m)
b                    = 0.1;    % length of the half crack (m)

h                    =0.008; % thickness of the stiff pore (cylinder) (m)
Stiff_Pore_volume    = 1   *  pi*b^2*h/4;
Compl_pore_Volume    = 1   *  pi*b^2*a/4; % crack volume 
Porosity_Volume      = Stiff_Pore_volume + Compl_pore_Volume; 
Porosity1            = Porosity_Volume./Volume_L;   % total porosity (fraction)
phi_d                = Stiff_Pore_volume./Volume_L; % stiff porosity (fraction)
phi_c                = Compl_pore_Volume./Volume_L; % crack porosity (fraction)
aspect_ratio         = a/(b*2 ) /2.0; % aspect ratio of the crack (cylinder)
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

Leak_Z_N    = 2*0.5*(Sij_cr_Zn(3,3) - Sij(3,3));
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
%[max1,max2] = max(Q_layr);

Kf_branch = K_layr;
%%
Kf_branch = Kf_branch.*1 + 0*4.3*1e9;
Z_N_mf_eval     = Z_N./(1 +  Z_N./(phi_c.*(1./Kf_branch - 1./36e9)) ); %1./K_g  K_starNew
Z_N_crack       = complex( real(Z_N_mf_eval), imag(Z_N_mf_eval));

del_Sij_crack   = zeros(6,6,size(w,2)); 
Sij_mf_blue          = zeros(6,6,size(w,2)); 
Cij_final       = zeros(6,6,size(w,2)); 
for nn = 1:size(w,2) % add crack into the saturated model (saturated model = saturated the stiff crack in grains)
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
%% END OF the present analytical model: hs=0.008 (thickness of the stiff crack)
%%
figure(1)
Ceff_33_black    = squeeze( Cij_final(3,3,:)    ./10^9);
plot(f,             real(Ceff_33_black)   ,'-','color',[0.789062500000000,0.906250000000000,0.945312500000000],...
    'LineWidth',5); hold on;

figure(2)
plot(f,            imag(Ceff_33_black)./real(Ceff_33_black),'-','color',[0.789062500000000,0.906250000000000,0.945312500000000],...
    'LineWidth',5); hold on;

%%
%%
%%
%% Present analytical model: hs=0.01 (thickness of the stiff crack)
%% 
%%
%%
%% background moduli: (dry) stiff crack in grains. VTI Num 
Cij(1,1) = 87230044269*0 + 93588000000;
Cij(1,2) = 6180414279 *0 + 6564000000;
Cij(1,3) = 7068404997 *0 + 5944000000;
Cij(2,2) = 87230005255*0 + 93588000000;
Cij(2,3) = 7068426541 *0 + 5944000000;
Cij(3,3) = 84040612409*0 + 81778000000;
Cij(4,4) = 39720000000;
Cij(5,5) = 39720000000;
Cij(6,6) = 43512000000;
     for i=1:6
        for j=i:6
      Cij(j,i) = Cij(i,j);
        end
     end
%% moduli: (dry) stiff crack +(dry) compliant crack in grains. VTI Num full 
Cij_cr(1,1) = 86959886906*0 + 93352000000;
Cij_cr(1,2) = 5950757722*0 + 6496000000;
Cij_cr(1,3) = 4849479475*0 + 5030000000;
Cij_cr(2,2) = 86957367522*0 + 93352000000;
Cij_cr(2,3) = 4850230042*0 + 5030000000;
Cij_cr(3,3) = 61468382429*0 + 68714000000;
Cij_cr(4,4) = 35518000000;
Cij_cr(5,5) = 35518000000;
Cij_cr(6,6) = 43428000000;
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

h                    =0.01; % thickness of the stiff pore (cylinder) (m)
Stiff_Pore_volume    = 1   *  pi*b^2*h/4;

Compl_pore_Volume    =   pi*b^2*a/4; % crack volume 
Porosity_Volume      = Stiff_Pore_volume + Compl_pore_Volume; 
Porosity1            = Porosity_Volume./Volume_L;   % total porosity (fraction)
phi_d                = Stiff_Pore_volume./Volume_L; % stiff porosity (fraction)
phi_c                = 1*Compl_pore_Volume./Volume_L; % crack porosity (fraction)

aspect_ratio         = a/(b*2 ) /2.0;
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
for nn = 1:size(w,2) % add crack into the saturated model (saturated model = saturated the stiff crack in grains)
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
%% END OF the present analytical model: hs=0.01 (thickness of the stiff crack)
%%
figure(1)%,clf
Ceff_33_black    = squeeze( Cij_final(3,3,:)    ./10^9); 
plot(f,             real(Ceff_33_black)   ,'-','color',[0.992187500000000,0.980468750000000,0.726562500000000],...
    'LineWidth',5); hold on;

figure(2)%, clf
plot(f,            imag(Ceff_33_black)./real(Ceff_33_black)*1.0,'-','color',[0.992187500000000,0.980468750000000,0.726562500000000],...
    'LineWidth',5); hold on;

%%
%%
%%
%% Present analytical model: hs=0.015 (thickness of the stiff crack)
%% 
%%
%%
%% background moduli: (dry) stiff crack in grains. VTI Num 
Cij(1,1) = 87230044269*0 + 93004000000;
Cij(1,2) = 6180414279 *0 + 6548000000;
Cij(1,3) = 7068404997 *0 + 5976000000;
Cij(2,2) = 87230005255*0 + 93004000000;
Cij(2,3) = 7068426541 *0 + 5976000000;
Cij(3,3) = 84040612409*0 + 81416000000;
Cij(4,4) = 39298000000;
Cij(5,5) = 39298000000;
Cij(6,6) = 43228000000;
     for i=1:6
        for j=i:6
      Cij(j,i) = Cij(i,j);
        end
     end
%% moduli: (dry) stiff crack +(dry) compliant crack in grains. VTI Num full 
Cij_cr(1,1) = 86959886906*0 + 92780000000;
Cij_cr(1,2) = 5950757722*0 + 6476000000;
Cij_cr(1,3) = 4849479475*0 + 5054000000;
Cij_cr(2,2) = 86957367522*0 + 92780000000;
Cij_cr(2,3) = 4850230042*0 + 5054000000;
Cij_cr(3,3) = 61468382429*0 + 68410000000;
Cij_cr(4,4) = 35144000000;
Cij_cr(5,5) = 35144000000;
Cij_cr(6,6) = 43152000000;
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

h                    =0.015; % thickness of the stiff pore (cylinder) (m)
Stiff_Pore_volume    = 1   *  pi*b^2*h/4;

Compl_pore_Volume    =   pi*b^2*a/4; % crack volume 
Porosity_Volume      = Stiff_Pore_volume + Compl_pore_Volume; 
Porosity1            = Porosity_Volume./Volume_L;   % total porosity (fraction)
phi_d                = Stiff_Pore_volume./Volume_L; % stiffporosity (fraction)
phi_c                = 1*Compl_pore_Volume./Volume_L; % crack porosity (fraction)

aspect_ratio         = a/(b*2 ) /2.0;
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
for nn = 1:size(w,2) % add crack into the saturated model (saturated model = saturated the stiff crack in grains)
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
%% END OF the present analytical model: hs=0.015 (thickness of the stiff crack)
%%
figure(1)%,clf
Ceff_33_black    = squeeze( Cij_final(3,3,:)    ./10^9); 
plot(f,             real(Ceff_33_black)   ,'-','color',[0.988281250000000,0.722656250000000,0.417968750000000],...
    'LineWidth',5); hold on;

figure(2)%, clf
plot(f,            imag(Ceff_33_black)./real(Ceff_33_black)*1.0,'-','color',[0.988281250000000,0.722656250000000,0.417968750000000],...
    'LineWidth',5); hold on;

%%
%%
%%
%% Present analytical model: hs=0.02 (thickness of the stiff crack)
%% 
%%
%%
%% background moduli: (dry) stiff crack in grains. VTI Num
Cij(1,1) = 87230044269*0 + 92380000000;
Cij(1,2) = 6180414279 *0 + 6524000000;
Cij(1,3) = 7068404997 *0 + 5994000000;
Cij(2,2) = 87230005255*0 + 92380000000;
Cij(2,3) = 7068426541 *0 + 5994000000;
Cij(3,3) = 84040612409*0 + 81120000000;
Cij(4,4) = 38850000000;
Cij(5,5) = 38850000000;
Cij(6,6) = 42928000000;
     for i=1:6
        for j=i:6
      Cij(j,i) = Cij(i,j);
        end
     end
%% moduli: (dry) stiff +(dry) compliant crack in grains. VTI Num full
Cij_cr(1,1) = 86959886906*0 + 92156000000;
Cij_cr(1,2) = 5950757722*0 + 6460000000;
Cij_cr(1,3) = 4849479475*0 + 5058000000;
Cij_cr(2,2) = 86957367522*0 + 92156000000;
Cij_cr(2,3) = 4850230042*0 + 5058000000;
Cij_cr(3,3) = 61468382429*0 + 68024000000;
Cij_cr(4,4) = 34730000000;
Cij_cr(5,5) = 34730000000;
Cij_cr(6,6) = 42848000000;
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
phi_d                = Stiff_Pore_volume./Volume_L; % stiff porosity (fraction)
phi_c                = 1*Compl_pore_Volume./Volume_L; % crack porosity (fraction)

aspect_ratio         = a/(b*2 ) /2.0;
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
for nn = 1:size(w,2) % add crack into the saturated model (saturated model = saturated the stiff crack in grains)
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
%% END OF the present analytical model: hs=0.02 (thickness of the stiff crack)
%%
figure(1)%,clf
Ceff_33_black    = squeeze( Cij_final(3,3,:)    ./10^9); 
plot(f,             real(Ceff_33_black)   ,'-','color',[0.921875000000000,0.363281250000000,0.230468750000000],...
    'LineWidth',5); hold on;

figure(2)%, clf
plot(f,            imag(Ceff_33_black)./real(Ceff_33_black)*1.0,'-','color',[0.921875000000000,0.363281250000000,0.230468750000000],...
    'LineWidth',5); hold on;
%% 
%% END
%% 


figure(1)
set(gca, 'XScale', 'log') ; grid on; 
 

xlim([0.5 10^6.5]); ylim([83.5 93.5]);


legend('hc=0.002, hs=0.004','hc=0.002, hs=0.005',...
     'hc=0.002, hs=0.008','hc=0.002, hs=0.01',...
          'hc=0.002, hs=0.015','hc=0.002, hs=0.02','Location','southeast');
  
 
xlabel('  Frequency (Hz)','fontsize',12); ylabel(' Re(C33) (GPa)','fontsize',12);
x_range = [1e0 1e1 1e2 1e3 1e4 1e5 1e6];
set(gca,'xtick',x_range);

figure(2)

set(gca, 'XScale', 'log') ; grid on; 
set(gca, 'YScale', 'log') ; 
xlim([0.5 10^6.5]); ylim([6e-7 2e-2]);
     
xlabel('  Frequency (Hz)','fontsize',12); ylabel('  1/Q','fontsize',12); hold on;
x_range = [1e0 1e1 1e2 1e3 1e4 1e5 1e6];
set(gca,'xtick',x_range); 
