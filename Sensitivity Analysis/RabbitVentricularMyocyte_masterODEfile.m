function dydt = RabbitVentricularMyocyte_masterODEfile(t,y,p)
% This function calls the ode files for EC coupling, CaM reactions, CaMKII
% phosphorylation module, and PKA phosphorylation module.

%% Collect params and ICs for each module

% CaMKII
camkii_exp = p(1); % 0, 1 or 6
% ISO
Ligtot = p(2); % [uM] ISO
% Protocol
prot_index = p(3);
cycleLength = p(4);     % [ms]
prot_input_par = p(5); % Input parameter for stimulation protocols
% Ca clamp
Ca_clamp = p(6); % 0 Ca-free, 1 Ca-clamp to initial value
% Na clamp
Na_clamp = p(7); % 0 Na-free, 1 Na-clamp to initial value
% Modules to use
flag_ECC = p(8);    % if 0, module clamped
flag_cam = p(9);    % if 0, module clamped
flag_CaMKII = p(10); % if 0, module clamped
flag_BAR = p(11);   % if 0, module clamped
% State variables for each module
ny_ECC = p(12);     % state vars in ECC module
ny_cam = p(13);     % state vars in each CaM module
ny_CaMKII = p(14);  % state vars in CaMKII module
ny_BAR = p(15);     % state vars in BAR module
% Sensitivity Analysis
par_SA = p(16:end);
%% Allocate ICs for each model

% Ca_j is y(36), Ca_sl is y(37), Ca_cytosol is y(38)
y_ecc = y(1:ny_ECC); 
y_camDyad = y(ny_ECC+1:ny_ECC+ny_cam);
y_camSL = y(ny_ECC+ny_cam+1:ny_ECC+2*ny_cam);
y_camCyt = y(ny_ECC+2*ny_cam+1:ny_ECC+3*ny_cam);
y_CaMKII = y(ny_ECC+3*ny_cam+1:ny_ECC+3*ny_cam+ny_CaMKII);
y_BAR = y(ny_ECC+3*ny_cam+ny_CaMKII+1:ny_ECC+3*ny_cam+ny_CaMKII+ny_BAR);
%% Parameters

K = y(35); %135; % [mM]
Mg = 1;  % [mM]

CKIIOE = 0+(camkii_exp>1);

CaMtotDyad = 418;           % [uM]
BtotDyad = 1.54/8.293e-4;   % [uM]
CaMKIItotDyad = 120*camkii_exp;        % [uM] 
CaNtotDyad = 3e-3/8.293e-4; % [uM] 
PP1totDyad = 96.5;          % [uM]
CaMtotSL = 5.65;            % [uM]
BtotSL = 24.2;              % [uM]
CaMKIItotSL = 120*8.293e-4*camkii_exp; % [uM]
CaNtotSL = 3e-3;            % [uM]
PP1totSL = 0.57;            % [uM]
CaMtotCyt = 5.65;           % [uM]
BtotCyt = 24.2;             % [uM]
CaMKIItotCyt = 120*8.293e-4*camkii_exp;% [uM]
CaNtotCyt = 3e-3;           % [uM] 
PP1totCyt = 0.57;           % [uM]

LCCtotDyad = 31.4*.9;       % [uM] - Total Dyadic [LCC] - (umol/l dyad)
LCCtotSL = 0.0846;          % [uM] - Total Subsarcolemmal [LCC] (umol/l sl)
RyRtot = 382.6;             % [uM] - Total RyR (in Dyad)
PP1_dyad = 95.7;            % [uM] - Total dyadic [PP1]
PP1_SL = 0.57;              % [uM] - Total Subsarcolemmal [PP1]
PP2A_dyad = 95.76;          % [uM] - Total dyadic PP2A
OA = 0;                     % [uM] - PP1/PP2A inhibitor Okadaic Acid
PLBtot = 38;                % [uM] - Total [PLB] in cytosolic units

% Parameters for BAR module
LCCtotBA = 0.025;           % [uM] - [umol/L cytosol]
RyRtotBA = 0.135;           % [uM] - [umol/L cytosol]
PLBtotBA = PLBtot;          % [uM] - [umol/L cytosol]
TnItotBA = 70;              % [uM] - [umol/L cytosol]
IKstotBA = 0.025;           % [uM] - [umol/L cytosol]
ICFTRtotBA = 0.025;         % [uM] - [umol/L cytosol]
PP1_PLBtot = 0.89;          % [uM] - [umol/L cytosol]
PLMtotBA = 48;              % [uM] - [umol/L cytosol] as in Yang & Saucerman (mouse) model
MyototBA = 70;              % [uM] - [umol/L cytosol] as TnI
IKrtotBA = 0.025;           % [uM] - [umol/L cytosol] as IKs
IClCatotBA = 0.025;         % [uM] - [umol/L cytosol] as ICFTR
ItototBA = 0.025;         % [uM] - [umol/L cytosol] as IKr
IK1totBA = 0.025;         % [uM] - [umol/L cytosol] as IKr
INatotBA = 0.025;         % [uM] - [umol/L cytosol] as IKr
%% Distribute parameters by module

% CaM module
CaDyad = y(36)*1e3; % from ECC model, *** Converting from [mM] to [uM] ***
compart_dyad = 2;
% ** NOTE: Btotdyad being sent to the dyad camODEfile is set to zero, but is used below for transfer between SL and dyad
pCaMDyad = [K, Mg, CaMtotDyad, 0, CaMKIItotDyad, CaNtotDyad, PP1totDyad, CaDyad, cycleLength, compart_dyad];
CaSL = y(37)*1e3; % from ECC model, *** Converting from [mM] to [uM] ***
compartSL = 1;
pCaMSL = [K, Mg, CaMtotSL, BtotSL, CaMKIItotSL, CaNtotSL, PP1totSL, CaSL, cycleLength, compartSL];
CaCyt = y(38)*1e3; % from ECC model, *** Converting from [mM] to [uM] ***
compartCyt = 0;
pCaMCyt = [K, Mg, CaMtotCyt, BtotCyt, CaMKIItotCyt, CaNtotCyt, PP1totCyt, CaCyt, cycleLength, compartCyt];

% CaMKII phosphorylation module 
CaMKIIact_Dyad = CaMKIItotDyad.*(y(ny_ECC+8)+y(ny_ECC+9)+y(ny_ECC+10)+y(ny_ECC+11)); % Multiply total by fraction
CaMKIIact_SL = CaMKIItotSL.*(y(ny_ECC+ny_cam+8)+y(ny_ECC+ny_cam+9)+y(ny_ECC+ny_cam+10)+y(ny_ECC+ny_cam+11));
%PP1_PLB_avail = y(ny_ECC+3*ny_cam+ny_CaMKII+25)./PP1_PLBtot + .0091;  % Active PP1 near PLB / total PP1 conc + basal value
%PP1_PLB_avail = 1 - y(ny_ECC+3*ny_cam+ny_CaMKII+25)/PP1_PLBtot + 0.0102290996035997; % NEW ODEs
    PP1_PLB_avail = 1 - y(ny_ECC+3*ny_cam+ny_CaMKII+25)/PP1_PLBtot + 0.0196035491719399; % PROVA
% PP1_PLB_avail is 100% (1) with NO ISO (Derived)
pCaMKII = [CaMKIIact_Dyad,LCCtotDyad,RyRtot,PP1_dyad,PP2A_dyad,OA,PLBtot,...
           CaMKIIact_SL,LCCtotSL,PP1_SL,...
           PP1_PLB_avail];

LCC_CKdyadp = y(ny_ECC+3*ny_cam+2)./LCCtotDyad;   % fractional CaMKII-dependent LCC dyad phosphorylation
RyR_CKp = y(ny_ECC+3*ny_cam+4)./RyRtot;           % fractional CaMKII-dependent RyR phosphorylation
PLB_CKp = y(ny_ECC+3*ny_cam+5)./PLBtot;           % fractional CaMKII-dependent PLB phosphorylation

% BAR (PKA phosphorylation) module
pBAR = [Ligtot,PP1_PLBtot,PLBtotBA,PLMtotBA,LCCtotBA,...
    RyRtotBA,TnItotBA,MyototBA,IKstotBA,IKrtotBA,...
    ICFTRtotBA,IClCatotBA,ItototBA,IK1totBA,INatotBA];
PLB_PKAn = (PLBtotBA - y(ny_ECC+3*ny_cam+ny_CaMKII+27))./PLBtotBA; % non-phosphorylated PLB targets
PLM_PKAp = y(ny_ECC+3*ny_cam+ny_CaMKII+28)./PLMtotBA;
LCCa_PKAp = y(ny_ECC+3*ny_cam+ny_CaMKII+29)./LCCtotBA;
LCCb_PKAp = y(ny_ECC+3*ny_cam+ny_CaMKII+30)./LCCtotBA;
RyR_PKAp = y(ny_ECC+3*ny_cam+ny_CaMKII+31)./RyRtotBA;
TnI_PKAp = y(ny_ECC+3*ny_cam+ny_CaMKII+32)./TnItotBA;
Myo_PKAp = y(ny_ECC+3*ny_cam+ny_CaMKII+33)./MyototBA;
IKs_PKAp = y(ny_ECC+3*ny_cam+ny_CaMKII+34)./IKstotBA;
IKr_PKAp = y(ny_ECC+3*ny_cam+ny_CaMKII+35)./IKrtotBA;
ICFTR_PKAp = y(ny_ECC+3*ny_cam+ny_CaMKII+36)./ICFTRtotBA;
IClCa_PKAp = y(ny_ECC+3*ny_cam+ny_CaMKII+37)./IClCatotBA;
Ito_PKAp = y(ny_ECC+3*ny_cam+ny_CaMKII+38)./ItototBA;
IK1_PKAp = y(ny_ECC+3*ny_cam+ny_CaMKII+39)./IK1totBA;
INa_PKAp = y(ny_ECC+3*ny_cam+ny_CaMKII+40)./INatotBA;

% ECC module
pECC = [prot_index,cycleLength,prot_input_par,CKIIOE,Ca_clamp,Na_clamp,...
    LCC_CKdyadp,RyR_CKp,PLB_CKp,... % CaMKII
    LCCa_PKAp,LCCb_PKAp,PLB_PKAn,RyR_PKAp,TnI_PKAp,... % PKA
    IKs_PKAp,ICFTR_PKAp,PLM_PKAp,Myo_PKAp,IKr_PKAp,IClCa_PKAp,...
    Ito_PKAp,IK1_PKAp,INa_PKAp,...
    par_SA];  
%% Solve dydt in each module

global JCaCyt JCaSL JCaDyad;

if flag_ECC == 0
    dydt_ecc = zeros(1,length(y_ecc))';
else
    dydt_ecc = RabbitVentricularMyocyte_eccODEfile(t,y_ecc,pECC);
end

if flag_cam == 0
    dydt_camDyad = zeros(1,length(y_camDyad))';
	dydt_camSL = zeros(1,length(y_camSL))';
	dydt_camCyt = zeros(1,length(y_camCyt))';
    JCaDyad = 0; JCaSL = 0; JCaCyt = 0;
else
    dydt_camDyad = RabbitVentricularMyocyte_camODEfile(t,y_camDyad,pCaMDyad);
    dydt_camSL = RabbitVentricularMyocyte_camODEfile(t,y_camSL,pCaMSL);
    dydt_camCyt = RabbitVentricularMyocyte_camODEfile(t,y_camCyt,pCaMCyt);
end

if flag_CaMKII == 0
    dydt_CaMKIIDyad = zeros(1,length(y_CaMKII))';
else
    dydt_CaMKIIDyad = RabbitVentricularMyocyte_camkiiODEfile(t,y_CaMKII,pCaMKII);
end

if flag_BAR == 0
    dydt_BAR = zeros(1,length(y_BAR))';
else
    dydt_BAR = RabbitVentricularMyocyte_barODEfile(t,y_BAR,pBAR);
end

% incorporate Ca buffering from CaM, convert JCaCyt from uM/msec to mM/msec
if Ca_clamp == 0
    dydt_ecc(36) = dydt_ecc(36) + 1e-3*JCaDyad;     % 0 with Ca-clamp
    dydt_ecc(37) = dydt_ecc(37) + 1e-3*JCaSL;       % 0 with Ca-clamp
    dydt_ecc(38) = dydt_ecc(38) + 1e-3*JCaCyt;      % 0 with Ca-clamp
end
    
% incorporate CaM diffusion between compartments
cellLength = 100;     % cell length [um]
cellRadius = 10.25;   % cell radius [um]
Vcell = pi*cellRadius^2*cellLength*1e-15;    % [L]
Vmyo = 0.65*Vcell; VSL = 0.02*Vcell; Vdyad = 0.0539*.01*Vcell; 

kDyadSL = 3.6363e-16;	% [L/msec]
kSLmyo = 8.587e-15;     % [L/msec]
k0Boff = 0.0014;        % [s^-1] 
k0Bon = k0Boff/0.2;     % [uM^-1 s^-1] kon = koff/Kd
k2Boff = k0Boff/100;    % [s^-1] 
k2Bon = k0Bon;          % [uM^-1 s^-1]
k4Boff = k2Boff;        % [s^-1]
k4Bon = k0Bon;          % [uM^-1 s^-1]
CaMtotDyad = sum(y_camDyad(1:6))+CaMKIItotDyad*sum(y_camDyad(7:10))+sum(y_camDyad(13:15));
Bdyad = BtotDyad - CaMtotDyad; % [uM dyad]
J_cam_dyadSL = 1e-3*(k0Boff*y_camDyad(1) - k0Bon*Bdyad*y_camSL(1)); % [uM/msec dyad]
J_ca2cam_dyadSL = 1e-3*(k2Boff*y_camDyad(2) - k2Bon*Bdyad*y_camSL(2)); % [uM/msec dyad]
J_ca4cam_dyadSL = 1e-3*(k2Boff*y_camDyad(3) - k4Bon*Bdyad*y_camSL(3)); % [uM/msec dyad]
J_cam_SLmyo = kSLmyo*(y_camSL(1)-y_camCyt(1)); % [umol/msec]
J_ca2cam_SLmyo = kSLmyo*(y_camSL(2)-y_camCyt(2)); % [umol/msec]
J_ca4cam_SLmyo = kSLmyo*(y_camSL(3)-y_camCyt(3)); % [umol/msec]
dydt_camDyad(1) = dydt_camDyad(1) - J_cam_dyadSL;
dydt_camDyad(2) = dydt_camDyad(2) - J_ca2cam_dyadSL;
dydt_camDyad(3) = dydt_camDyad(3) - J_ca4cam_dyadSL;
dydt_camSL(1) = dydt_camSL(1) + J_cam_dyadSL*Vdyad/VSL - J_cam_SLmyo/VSL;
dydt_camSL(2) = dydt_camSL(2) + J_ca2cam_dyadSL*Vdyad/VSL - J_ca2cam_SLmyo/VSL;
dydt_camSL(3) = dydt_camSL(3) + J_ca4cam_dyadSL*Vdyad/VSL - J_ca4cam_SLmyo/VSL;
dydt_camCyt(1) = dydt_camCyt(1) + J_cam_SLmyo/Vmyo;
dydt_camCyt(2) = dydt_camCyt(2) + J_ca2cam_SLmyo/Vmyo;
dydt_camCyt(3) = dydt_camCyt(3) + J_ca4cam_SLmyo/Vmyo;
%% Collect all dydt terms

dydt = [dydt_ecc; dydt_camDyad; dydt_camSL; dydt_camCyt; dydt_CaMKIIDyad; dydt_BAR];