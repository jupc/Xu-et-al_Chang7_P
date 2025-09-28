function dy = MBOX_equations(t,y)
%%%%%% MBOX 

%%%%%%% setup dy array  (FIX: must be 58)
dy = zeros(58,1);  

%%%%%%% set up global parameters
global stepnumber
global pars
global workingstate

%%%%%%%%%%%%% get variables from Y 
CO2_a  = y(1) ;
DIC_p  = y(2) ;
DIC_di = y(3) ;
DIC_s  = y(4) ;
DIC_h  = y(5) ;
DIC_d  = y(6) ;

ALK_p  = y(7) ;
ALK_di = y(8) ;
ALK_s  = y(9) ;
ALK_h  = y(10) ;
ALK_d  = y(11) ;

d13c_atm    = y(12) / y(1) ;
d13c_DIC_p  = y(13) / y(2) ;
d13c_DIC_di = y(14) / y(3) ;
d13c_DIC_s  = y(15) / y(4) ;
d13c_DIC_h  = y(16) / y(5) ;
d13c_DIC_d  = y(17) / y(6) ;

POC_p  = y(18);
POC_di = y(19);
POC_s  = y(20);
POC_h  = y(21);
POC_d  = y(22);

DP_p  = y(23);
DP_di = y(24);
DP_s  = y(25);
DP_h  = y(26);
DP_d  = y(27);

d13c_POC_p  = y(28) / y(18) ;
d13c_POC_di = y(29) / y(19) ;
d13c_POC_s  = y(30) / y(20) ;
d13c_POC_h  = y(31) / y(21) ;
d13c_POC_d  = y(32) / y(22) ;

O2_a  = y(33);
O2_p  = y(34);
O2_di = y(35);
O2_s  = y(36);
O2_h  = y(37);
O2_d  = y(38);

FeIII_p  = y(39);
FeIII_di = y(40);
FeIII_s  = y(41);
FeIII_h  = y(42);
FeIII_d  = y(43);

SO4_p  = y(44);
SO4_di = y(45);
SO4_s  = y(46);
SO4_h  = y(47);
SO4_d  = y(48);

FeII_p  = y(49);
FeII_di = y(50);
FeII_s  = y(51);
FeII_h  = y(52);
FeII_d  = y(53);

H2S_p  = y(54);
H2S_di = y(55);
H2S_s  = y(56);
H2S_h  = y(57);
H2S_d  = y(58);

Atmospheric_CO2_ppm = ( CO2_a / pars.CO2_a_0 ) * 280 ;

%%%% pCO2 in PAL
pCO2_a = ( CO2_a / pars.CO2_a_0 ) ;

%%%% pO2 in PAL
pO2_a = ( O2_a / pars.O2_a_0 ) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   兼容桥：从 t_T0/T0_cm 推导 win_*    %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(pars,'volc')
    % 支持 width_m 或 width_km（若只给 km，就转成 m，供后续体积计算）
    if isfield(pars.volc,'width_km') && ~isfield(pars.volc,'width_m')
        pars.volc.width_m = pars.volc.width_km * 1e3;
    end
    % 若未显式给 win_*，但给了 t_T0/T0_cm，则把 >0 的平台段转为窗口区间
    if ~isfield(pars.volc,'win_t1') && isfield(pars.volc,'t_T0') && isfield(pars.volc,'T0_cm') ...
            && ~isempty(pars.volc.t_T0) && numel(pars.volc.t_T0)==numel(pars.volc.T0_cm)
        tt = pars.volc.t_T0(:);
        vv = pars.volc.T0_cm(:);
        win_t1 = []; win_t2 = []; win_T0 = [];
        k = 1;
        while k < numel(tt)
            if vv(k) > 0
                t1 = tt(k); T = vv(k);
                j = k + 1;
                while j <= numel(tt) && vv(j) > 0
                    j = j + 1;
                end
                if j <= numel(tt), t2 = tt(j); else, t2 = tt(end); end
                win_t1(end+1,1) = t1; %#ok<AGROW>
                win_t2(end+1,1) = t2; %#ok<AGROW>
                win_T0(end+1,1) = T;  %#ok<AGROW>
                k = j;
            else
                k = k + 1;
            end
        end
        if ~isempty(win_t1)
            pars.volc.win_t1 = win_t1;
            pars.volc.win_t2 = win_t2;
            pars.volc.win_T0 = win_T0;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Flux calculations   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Thermohaline speed (Sv)
circ_TH_Sv = 10 ;
circ_coastal_Sv = 2 ;

%%%% Thermohaline speed (m3/yr)
circ_TH_m3_yr = ( circ_TH_Sv * 1e6 ) * 3.15e7 ; 
circ_coastal_m3_yr = ( circ_coastal_Sv * 1e6 ) * 3.15e7 ;

%%%% mixing coefficient (m3/yr)
mixcoeff_dip_m3_yr = pars.mixcoeff_dip * 1e6 * 3.15e7 ;
mixcoeff_ds_m3_yr  = pars.mixcoeff_ds * 1e6 * 3.15e7 ;
mixcoeff_dh_m3_yr  = pars.mixcoeff_dh * 1e6 * 3.15e7 ;

%%%% DIC concentration in mol/m3
DIC_conc_p  = DIC_p / pars.vol_p ;
DIC_conc_di = DIC_di / pars.vol_di ;
DIC_conc_s  = DIC_s / pars.vol_s ;
DIC_conc_h  = DIC_h / pars.vol_h ;
DIC_conc_d  = DIC_d / pars.vol_d ;

ALK_conc_p  = ALK_p / pars.vol_p ;
ALK_conc_di = ALK_di / pars.vol_di ;
ALK_conc_s  = ALK_s / pars.vol_s ;
ALK_conc_h  = ALK_h / pars.vol_h ;
ALK_conc_d  = ALK_d / pars.vol_d ;

POC_conc_p  = POC_p / pars.vol_p ;
POC_conc_di = POC_di / pars.vol_di ;
POC_conc_s  = POC_s / pars.vol_s ;
POC_conc_h  = POC_h / pars.vol_h ;
POC_conc_d  = POC_d / pars.vol_d ;

DP_conc_p  = DP_p / pars.vol_p ;
DP_conc_di = DP_di / pars.vol_di ;
DP_conc_s  = DP_s / pars.vol_s ;
DP_conc_h  = DP_h / pars.vol_h ;
DP_conc_d  = DP_d / pars.vol_d ;

O2_conc_p  = O2_p / pars.vol_p ;
O2_conc_di = O2_di / pars.vol_di ;
O2_conc_s  = O2_s / pars.vol_s ;
O2_conc_h  = O2_h / pars.vol_h ;
O2_conc_d  = O2_d / pars.vol_d ;

FeIII_conc_p  = FeIII_p/ pars.vol_p ;
FeIII_conc_di = FeIII_di/ pars.vol_di ;
FeIII_conc_s  = FeIII_s/ pars.vol_s ;
FeIII_conc_h  = FeIII_h/ pars.vol_h ;
FeIII_conc_d  = FeIII_d/ pars.vol_d ;

SO4_conc_p  = SO4_p/ pars.vol_p ;
SO4_conc_di = SO4_di/ pars.vol_di ;
SO4_conc_s  = SO4_s/ pars.vol_s ;
SO4_conc_h  = SO4_h/ pars.vol_h ;
SO4_conc_d  = SO4_d/ pars.vol_d ;

FeII_conc_p  = FeII_p/ pars.vol_p ;
FeII_conc_di = FeII_di/ pars.vol_di ;
FeII_conc_s  = FeII_s/ pars.vol_s ;
FeII_conc_h  = FeII_h/ pars.vol_h ;
FeII_conc_d  = FeII_d/ pars.vol_d ;

H2S_conc_p  = H2S_p/ pars.vol_p ;
H2S_conc_di = H2S_di/ pars.vol_di ;
H2S_conc_s  = H2S_s/ pars.vol_s ;
H2S_conc_h  = H2S_h/ pars.vol_h ;
H2S_conc_d  = H2S_d/ pars.vol_d ;

%%%% Transport fluxes in mol/yr
Tran_DIC_p_s  = DIC_conc_p * circ_coastal_m3_yr ;
Tran_DIC_s_h  = DIC_conc_s * circ_TH_m3_yr ;
Tran_DIC_h_d  = DIC_conc_h * circ_TH_m3_yr ;
Tran_DIC_d_di = DIC_conc_d * circ_coastal_m3_yr ;
Tran_DIC_di_p = DIC_conc_di * circ_coastal_m3_yr ;
Tran_DIC_d_s  = DIC_conc_d * ( circ_TH_m3_yr  - circ_coastal_m3_yr) ;

Tran_ALK_p_s  = ALK_conc_p * circ_coastal_m3_yr ;
Tran_ALK_s_h  = ALK_conc_s * circ_TH_m3_yr ;
Tran_ALK_h_d  = ALK_conc_h * circ_TH_m3_yr ;
Tran_ALK_d_di = ALK_conc_d * circ_coastal_m3_yr ;
Tran_ALK_di_p = ALK_conc_di * circ_coastal_m3_yr ;
Tran_ALK_d_s  = ALK_conc_d * ( circ_TH_m3_yr  - circ_coastal_m3_yr) ;

Tran_POC_p_s  = POC_conc_p * circ_coastal_m3_yr ;
Tran_POC_s_h  = POC_conc_s * circ_TH_m3_yr ;
Tran_POC_p_di = 0.103*POC_p ;       
Tran_POC_s_d  = 0.167*POC_s ;       
Tran_POC_h_d  = 0.131*POC_h ;       
Tran_POC_d_di = POC_conc_d * circ_coastal_m3_yr ;

Tran_DP_p_s  = DP_conc_p * circ_coastal_m3_yr ;
Tran_DP_s_h  = DP_conc_s * circ_TH_m3_yr ;
Tran_DP_h_d  = DP_conc_h * circ_TH_m3_yr ;
Tran_DP_d_di = DP_conc_d * circ_coastal_m3_yr ;
Tran_DP_di_p = DP_conc_di * circ_coastal_m3_yr ;
Tran_DP_d_s  = DP_conc_d * ( circ_TH_m3_yr  - circ_coastal_m3_yr) ;

Tran_O2_p_s  = O2_conc_p * circ_coastal_m3_yr ;
Tran_O2_s_h  = O2_conc_s * circ_TH_m3_yr ;
Tran_O2_h_d  = O2_conc_h * circ_TH_m3_yr ;
Tran_O2_d_di = O2_conc_d * circ_coastal_m3_yr ;
Tran_O2_di_p = O2_conc_di * circ_coastal_m3_yr ;
Tran_O2_d_s  = O2_conc_d * ( circ_TH_m3_yr  - circ_coastal_m3_yr) ;

Tran_FeIII_p_s  = FeIII_conc_p * circ_coastal_m3_yr ;
Tran_FeIII_s_h  = FeIII_conc_s * circ_TH_m3_yr ;
Tran_FeIII_h_d  = FeIII_conc_h * circ_TH_m3_yr ;
Tran_FeIII_d_di = FeIII_conc_d * circ_coastal_m3_yr ;
Tran_FeIII_di_p = FeIII_conc_di * circ_coastal_m3_yr ;
Tran_FeIII_d_s  = FeIII_conc_d * ( circ_TH_m3_yr  - circ_coastal_m3_yr) ;

Tran_SO4_p_s  = SO4_conc_p * circ_coastal_m3_yr ;
Tran_SO4_s_h  = SO4_conc_s * circ_TH_m3_yr ;
Tran_SO4_h_d  = SO4_conc_h * circ_TH_m3_yr ;
Tran_SO4_d_di = SO4_conc_d * circ_coastal_m3_yr ;
Tran_SO4_di_p = SO4_conc_di * circ_coastal_m3_yr ;
Tran_SO4_d_s  = SO4_conc_d * ( circ_TH_m3_yr  - circ_coastal_m3_yr) ;

Tran_FeII_p_s  = FeII_conc_p * circ_coastal_m3_yr ;
Tran_FeII_s_h  = FeII_conc_s * circ_TH_m3_yr ;
Tran_FeII_h_d  = FeII_conc_h * circ_TH_m3_yr ;
Tran_FeII_d_di = FeII_conc_d * circ_coastal_m3_yr ;
Tran_FeII_di_p = FeII_conc_di * circ_coastal_m3_yr ;
Tran_FeII_d_s  = FeII_conc_d * ( circ_TH_m3_yr  - circ_coastal_m3_yr) ;

Tran_H2S_p_s  = H2S_conc_p * circ_coastal_m3_yr ;
Tran_H2S_s_h  = H2S_conc_s * circ_TH_m3_yr ;
Tran_H2S_h_d  = H2S_conc_h * circ_TH_m3_yr ;
Tran_H2S_d_di = H2S_conc_d * circ_coastal_m3_yr ;
Tran_H2S_di_p = H2S_conc_di * circ_coastal_m3_yr ;
Tran_H2S_d_s  = H2S_conc_d * ( circ_TH_m3_yr  - circ_coastal_m3_yr) ;

%%%% Water mixing terms
Mix_DP_di_p =  mixcoeff_ds_m3_yr*(DP_conc_di - DP_conc_p) ;
Mix_DP_d_s  =  mixcoeff_ds_m3_yr*(DP_conc_d - DP_conc_s) ;
Mix_DP_d_h  =  mixcoeff_dh_m3_yr*(DP_conc_d - DP_conc_h) ;

Mix_DIC_di_p = mixcoeff_ds_m3_yr*(DIC_conc_di - DIC_conc_p);
Mix_DIC_d_s  = mixcoeff_ds_m3_yr*(DIC_conc_d - DIC_conc_s);
Mix_DIC_d_h  = mixcoeff_dh_m3_yr*(DIC_conc_d - DIC_conc_h) ;

Mix_ALK_di_p = mixcoeff_ds_m3_yr*(ALK_conc_di - ALK_conc_p);
Mix_ALK_d_s  = mixcoeff_ds_m3_yr*(ALK_conc_d - ALK_conc_s);
Mix_ALK_d_h  = mixcoeff_dh_m3_yr*(ALK_conc_d - ALK_conc_h) ;

Mix_O2_di_p = mixcoeff_ds_m3_yr*(O2_conc_di - O2_conc_p);
Mix_O2_d_s  = mixcoeff_ds_m3_yr*(O2_conc_d - O2_conc_s);
Mix_O2_d_h  = mixcoeff_dh_m3_yr*(O2_conc_d - O2_conc_h) ;

Mix_FeIII_di_p = mixcoeff_ds_m3_yr*(FeIII_conc_di - FeIII_conc_p);
Mix_FeIII_d_s  = mixcoeff_ds_m3_yr*(FeIII_conc_d - FeIII_conc_s);
Mix_FeIII_d_h  = mixcoeff_dh_m3_yr*(FeIII_conc_d - FeIII_conc_h) ;

Mix_SO4_di_p = mixcoeff_ds_m3_yr*(SO4_conc_di - SO4_conc_p);
Mix_SO4_d_s  = mixcoeff_ds_m3_yr*(SO4_conc_d - SO4_conc_s);
Mix_SO4_d_h  = mixcoeff_dh_m3_yr*(SO4_conc_d - SO4_conc_h) ;

Mix_FeII_di_p = mixcoeff_ds_m3_yr*(FeII_conc_di - FeII_conc_p);
Mix_FeII_d_s  = mixcoeff_ds_m3_yr*(FeII_conc_d - FeII_conc_s);
Mix_FeII_d_h  = mixcoeff_dh_m3_yr*(FeII_conc_d - FeII_conc_h) ;

Mix_H2S_di_p = mixcoeff_ds_m3_yr*(H2S_conc_di - H2S_conc_p);
Mix_H2S_d_s  = mixcoeff_ds_m3_yr*(H2S_conc_d - H2S_conc_s);
Mix_H2S_d_h  = mixcoeff_dh_m3_yr*(H2S_conc_d - H2S_conc_h) ;

%%%% Global average surface temperature etc.
Climate_Sensitivity = 3 ;
t_geol = 0 ;
GAST = 288 + Climate_Sensitivity * ( log( Atmospheric_CO2_ppm / 280 ) / log(2) )  - 7.4*(t_geol/-570)  ; 
T_s = 298 + (GAST - 288)*0.66  ;  
T_h = max( 275.5 + (GAST - 288) , 271 ) ;
T_d = max( 275.5 + (GAST - 288) , 271 ) ;
T_cont = GAST ;
T_p = T_s ;
T_di = T_s ;

%%%% Carbonate chemistry parameters
k_2  = 7.4e-10 ;
k_carb_p = 0.000575 + 0.000006 * ( T_p - 278 ) ;
KCO2_p = 0.035 + 0.0019 * ( T_p - 278 ) ;
k_carb_di = 0.000575 + 0.000006 * ( T_di - 278 ) ;
k_carb_s = 0.000575 + 0.000006 * ( T_s - 278 ) ;
KCO2_s = 0.035 + 0.0019 * ( T_s - 278 ) ;
k_carb_h = 0.000575 + 0.000006 * ( T_h - 278 ) ;
KCO2_h = 0.035 + 0.0019 * ( T_h - 278 ) ;
k_carb_d = 0.000575 + 0.000006 * ( T_d - 278 ) ;

%%%% Carbonate speciation proximal
HCO3_p = ( DIC_conc_p - ( DIC_conc_p^2 - ALK_conc_p * ( 2 * DIC_conc_p - ALK_conc_p ) * ( 1 - 4 * k_carb_p ) )^0.5  ) / ( 1 - 4 * k_carb_p ) ;
CO3_p = ( ALK_conc_p - HCO3_p ) / 2 ;
H_p =  k_2 * HCO3_p / CO3_p ;
pH_p = -1 * log10(H_p) ;
%%%% Carbonate speciation distal
HCO3_di = ( DIC_conc_di - ( DIC_conc_di^2 - ALK_conc_di * ( 2 * DIC_conc_di - ALK_conc_di ) * ( 1 - 4 * k_carb_di ) )^0.5  ) / ( 1 - 4 * k_carb_di ) ;
CO3_di = ( ALK_conc_di - HCO3_di ) / 2 ;
H_di =  k_2 * HCO3_di / CO3_di ;
pH_di = -1 * log10(H_di) ;
%%%% Carbonate speciation lowlat
HCO3_s = ( DIC_conc_s - ( DIC_conc_s^2 - ALK_conc_s * ( 2 * DIC_conc_s - ALK_conc_s ) * ( 1 - 4 * k_carb_s ) )^0.5  ) / ( 1 - 4 * k_carb_s ) ;
CO3_s = ( ALK_conc_s - HCO3_s ) / 2 ;
H_s =  k_2 * HCO3_s / CO3_s ;
pH_s = -1 * log10(H_s) ;
%%%% Carbonate speciation hilat
HCO3_h = ( DIC_conc_h - ( DIC_conc_h^2 - ALK_conc_h * ( 2 * DIC_conc_h - ALK_conc_h ) * ( 1 - 4 * k_carb_h ) )^0.5  ) / ( 1 - 4 * k_carb_h ) ;
CO3_h = ( ALK_conc_h - HCO3_h ) / 2 ;
H_h =  k_2 * HCO3_h / CO3_h ;
pH_h = -1 * log10(H_h) ;
%%%% Carbonate speciation deep
HCO3_d = ( DIC_conc_d - ( DIC_conc_d^2 - ALK_conc_d * ( 2 * DIC_conc_d - ALK_conc_d ) * ( 1 - 4 * k_carb_d ) )^0.5  ) / ( 1 - 4 * k_carb_d ) ;
CO3_d = ( ALK_conc_d - HCO3_d ) / 2 ;
H_d = k_2 * HCO3_d / CO3_d ;
pH_d = -1 * log10(H_d) ;

%%%% Air-sea exchange (mol/yr)  
pCO2_p =    KCO2_p * ( ( HCO3_p^2 ) / CO3_p ) ;
AirSea_p = 5e16 * 0.1 * 0.1 * ( pCO2_a - pCO2_p ) ; 
pCO2_s =    KCO2_s * ( ( HCO3_s^2 ) / CO3_s ) ;
AirSea_s = 5e16 * 0.765 * 0.1 * ( pCO2_a - pCO2_s ) ; 
pCO2_h =   KCO2_h * ( ( HCO3_h^2 ) / CO3_h ) ;
AirSea_h = 5e16 * 0.135 * 0.1 * ( pCO2_a - pCO2_h ) ; 

AirSea_ap = 5e16 * 0.1 * 0.1 * ( pCO2_a) ;
AirSea_pa = 5e16 * 0.1 * 0.1 * ( pCO2_p) ;
AirSea_as = 5e16 * 0.765 * 0.1 * ( pCO2_a) ;
AirSea_sa = 5e16 * 0.765 * 0.1 * ( pCO2_s) ;
AirSea_ah = 5e16 * 0.135 * 0.1 * ( pCO2_a ) ;
AirSea_ha = 5e16 * 0.135 * 0.1 * ( pCO2_h ) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Continental fluxes  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
silconst = 0.33 ;
carbconst = 0.9 ;
UPLIFT = 1 ;
CARB_AREA = 1 ;
BAS_AREA = 1 ;
GRAN_AREA = 1 ;
ORG_AREA = 1 ;
PG = 1 ;
f_biota = 1 ;

f_T_bas =  exp(0.0608*(T_cont-288)) * ( (1 + 0.038*(T_cont - 288))^0.65 ) ; 
f_T_gran =  exp(0.0724*(T_cont-288)) * ( (1 + 0.038*(T_cont - 288))^0.65 ) ; 
g_T = 1 + 0.087*(T_cont - 288) ;

basw = pars.k_basw * BAS_AREA * PG * f_biota * f_T_bas ;
granw = pars.k_granw * UPLIFT^silconst * GRAN_AREA * PG * f_biota * f_T_gran ;
silw = basw + granw ;
carbw = pars.k_carbw * CARB_AREA * UPLIFT^carbconst * PG * f_biota * g_T ;
oxidw = pars.k_oxidw*UPLIFT^silconst*ORG_AREA*((O2_a/pars.O2_a_0)^0.5)  ;

Eaw = 8.3*4.184 ; % kJ mol-1
phoswi = pars.k_phosw ; 
FeIIIw = pars.k_FeIIIw*silw/(pars.k_basw+pars.k_granw);
phosw = interp1([-10e6 0 1000 5e5 1e6 5e6 10e6],[phoswi phoswi phoswi phoswi phoswi phoswi phoswi],t_geol) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Degassing / Volcanic-ash coupling  %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DEGASS = 1 ;
Bforcing = 1 ;
ccdeg = pars.k_ccdeg*DEGASS*Bforcing ;
ocdeg = pars.k_ocdeg*DEGASS ;

CO2_input = 0 ;
lake_P_release_yr = 0;   % mol P / yr
lake_C_burial_yr  = 0;   % mol C / yr

if isfield(pars,'volc') && isfield(pars.volc,'win_t1') && ~isempty(pars.volc.win_t1) ...
        && strcmp(pars.volc.coupling,'input')

    tau    = max(pars.volc.tau_release_yr, 1); % yr
    lambda = 1; if isfield(pars.volc,'lambda'), lambda = pars.volc.lambda; end

    for k = 1:numel(pars.volc.win_t1)
        t1   = pars.volc.win_t1(k);
        t2   = pars.volc.win_t2(k);
        T0cm = pars.volc.win_T0(k);

        Qp = event_total_P_mol_from_T0cm(T0cm, pars); % 每次喷发的总 P（mol）
        Qc = event_total_C_mol_from_T0cm(T0cm, pars); % 每次喷发导致的总 C（mol）

        if t < t1
            contP = 0; contC = 0;
        elseif t <= t2
            contP = lambda * Qp * (1 - exp(-(t - t1)/tau));
            contC = lambda * Qc * (1 - exp(-(t - t1)/tau));
        else
            contP = lambda * Qp * ( exp(-(t - t2)/tau) - exp(-(t - t1)/tau) );
            contC = lambda * Qc * ( exp(-(t - t2)/tau) - exp(-(t - t1)/tau) );
        end

        lake_P_release_yr = lake_P_release_yr + contP;
        lake_C_burial_yr  = lake_C_burial_yr  + contC;
    end

    CO2_input = - lake_C_burial_yr;   % 以 1:1 抵消大气 CO2
end

d13c_CO2_input = -5 ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Burial / productivity / redox  %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% CaCO3 saturation proximal/distal/deep
Ca_conc_p = 1.397e19 / 1.35e18 ; 
ksp_p = 0.76 ; 
sat_p = ( Ca_conc_p * CO3_p ) / ksp_p ;
satpresent_p = 3 ;

Ca_conc_di = 1.397e19 / 1.35e18 ; 
ksp_di = 0.76 ; 
sat_di = ( Ca_conc_di * CO3_di ) / ksp_di ;
satpresent_di = 3 ;

Ca_conc_d = 1.397e19 / 1.35e18 ; 
ksp_d = 0.76 ; 
sat_d = ( Ca_conc_d * CO3_d ) / ksp_d ;
satpresent_d = 1 ;

deepfrac_mccb = 0.1 ;
distalfrac_mccb = 0.2 ;
sat_minus_1p  = max( sat_p - 1 , 0 ) ;
mccb_p  = pars.k_mccb * (1-deepfrac_mccb-distalfrac_mccb) * (1 / satpresent_p) * ( (sat_minus_1p)^1.7 ) ;
sat_minus_1di = max( sat_di - 1 , 0 ) ;
mccb_di = pars.k_mccb * distalfrac_mccb * (1 / satpresent_di) * ( (sat_minus_1di)^1.7 ) ;
sat_minus_1d  = max( sat_d - 1 , 0 ) ;
mccb_d  = pars.k_mccb * deepfrac_mccb * (1 / satpresent_d) * ( (sat_minus_1d)^1.7 ) ;

locb = pars.k_locb ; % 不再额外叠加湖泊 C，避免与 CO2_input 双计

%%%% Primary productivity
pripr_p = 9.142857e+17*min(DP_conc_p*DP_conc_p/(DP_conc_p+pars.KP),FeIII_conc_p/pars.Red_Fe_P*FeIII_conc_p/(FeIII_conc_p+pars.KFe));
pripr_s = 7.80e+20*min(DP_conc_s*DP_conc_s/(DP_conc_s+pars.KP),FeIII_conc_s/pars.Red_Fe_P*FeIII_conc_s/(FeIII_conc_s+pars.KFe));
pripr_h =  2.22e+18*min(DP_conc_h*DP_conc_h/(DP_conc_h+pars.KP),FeIII_conc_h/pars.Red_Fe_P*FeIII_conc_h/(FeIII_conc_h+pars.KFe));

%%%% oxygen level that is in equilibrium with atmosphere in mol/m3
O2_eq_ap = exp(-173.4292+249.6339*100/T_p+143.3483*log(T_p/100)-21.8492*T_p/100+35*(-0.033096+0.014259*T_p/100- 0.0017000*T_p/100*T_p/100))/22.4*O2_a / pars.O2_a_0 ;
O2_eq_as = exp(-173.4292+249.6339*100/T_s+143.3483*log(T_s/100)-21.8492*T_s/100+35*(-0.033096+0.014259*T_s/100- 0.0017000*T_s/100*T_s/100))/22.4*O2_a / pars.O2_a_0 ;
O2_eq_ah = exp(-173.4292+249.6339*100/T_h+143.3483*log(T_h/100)-21.8492*T_h/100+35*(-0.033096+0.014259*T_h/100- 0.0017000*T_h/100*T_h/100))/22.4*O2_a / pars.O2_a_0 ;

AirSeaO2_p = 3.572831e+14 * ( O2_eq_ap - O2_conc_p) ;
AirSeaO2_s = 1.476284e+15 * ( O2_eq_as - O2_conc_s) ;
AirSeaO2_h = 2.210798e+16 * ( O2_eq_ah - O2_conc_h) ;

if pripr_p < 620e12
    pripr_p_O = pripr_p/pars.Red_C_O;
else
    pripr_p_O = 0.8*(pripr_p-620e12)/pars.Red_C_O + 620e12/pars.Red_C_O;
end

%%%% OC Remineralization
AR_p = 0.736*POC_p          *O2_conc_p/(O2_conc_p+pars.KmO2);
AR_di = 0.240*POC_di        *O2_conc_di/(O2_conc_di+pars.KmO2);
AR_s = 0.967*POC_s          *O2_conc_s/(O2_conc_s+pars.KmO2);
AR_h = 0.836*POC_h          *O2_conc_h/(O2_conc_h+pars.KmO2);
AR_d = 0.0098*POC_d         *O2_conc_d/(O2_conc_d+pars.KmO2);

ironR_p  = 0.736*POC_p  *pars.KmO2/(O2_conc_p+pars.KmO2)  * FeIII_conc_p/(FeIII_conc_p+pars.KmFeIII);
ironR_di = 0.240*POC_di*pars.KmO2/(O2_conc_di+pars.KmO2) * FeIII_conc_di/(FeIII_conc_di+pars.KmFeIII);
ironR_s  = 0.967*POC_s *pars.KmO2/(O2_conc_s+pars.KmO2)  * FeIII_conc_s/(FeIII_conc_s+pars.KmFeIII);
ironR_h  = 0.836*POC_h *pars.KmO2/(O2_conc_h+pars.KmO2)  * FeIII_conc_h/(FeIII_conc_h+pars.KmFeIII);
ironR_d  = 0.0098*POC_d*pars.KmO2/(O2_conc_d+pars.KmO2)  * FeIII_conc_d/(FeIII_conc_d+pars.KmFeIII);

sulfR_p  = 0.736*POC_p  *0.075*pars.KmO2/(O2_conc_p+pars.KmO2)  * pars.KmFeIII/(FeIII_conc_p+pars.KmFeIII);
sulfR_di = 0.240*POC_di*0.075*pars.KmO2/(O2_conc_di+pars.KmO2) * pars.KmFeIII/(FeIII_conc_di+pars.KmFeIII);
sulfR_s  = 0.967*POC_s *0.075*pars.KmO2/(O2_conc_s+pars.KmO2)  * pars.KmFeIII/(FeIII_conc_s+pars.KmFeIII);
sulfR_h  = 0.836*POC_h *0.075*pars.KmO2/(O2_conc_h+pars.KmO2)  * pars.KmFeIII/(FeIII_conc_h+pars.KmFeIII);
sulfR_d  = 0.0098*POC_d*0.075*pars.KmO2/(O2_conc_d+pars.KmO2)  * pars.KmFeIII/(FeIII_conc_d+pars.KmFeIII);

remin_p  = AR_p + ironR_p + sulfR_p;
remin_di = AR_di + ironR_di + sulfR_di;
remin_s  = AR_s + ironR_s + sulfR_s;
remin_h  = AR_h + ironR_h + sulfR_h;
remin_d  = AR_d + ironR_d + sulfR_d;

pyF_p  = pars.kpy * max( (FeII_conc_p * H2S_conc_p/(pars.KspFeSaq+H2S_conc_p) - pars.STFeSaq),  0 ) * H2S_conc_p * pars.vol_p ;
pyF_di = pars.kpy * max( (FeII_conc_di * H2S_conc_di/(pars.KspFeSaq+H2S_conc_di) - pars.STFeSaq),  0 ) * H2S_conc_di * pars.vol_di ;
pyF_s  = pars.kpy * max( (FeII_conc_s * H2S_conc_s/(pars.KspFeSaq+H2S_conc_s)  - pars.STFeSaq),  0 ) * H2S_conc_s * pars.vol_s ;
pyF_h  = pars.kpy * max( (FeII_conc_h * H2S_conc_h/(pars.KspFeSaq+H2S_conc_h)  - pars.STFeSaq),  0 ) * H2S_conc_h * pars.vol_h ;
pyF_d  = pars.kpy * max( (FeII_conc_d * H2S_conc_d/(pars.KspFeSaq+H2S_conc_d)  - pars.STFeSaq),  0 ) * H2S_conc_d * pars.vol_d ;

ironO_p  = pars.kironO * FeII_conc_p * O2_conc_p * pars.vol_p ;
ironO_di = pars.kironO * FeII_conc_di * O2_conc_di * pars.vol_di ;
ironO_s  = pars.kironO * FeII_conc_s * O2_conc_s * pars.vol_s ;
ironO_h  = pars.kironO * FeII_conc_h * O2_conc_h * pars.vol_h ;
ironO_d  = pars.kironO * FeII_conc_d * O2_conc_d * pars.vol_d ;

sulfO_p  = pars.ksulfO * H2S_conc_p * O2_conc_p * pars.vol_p ;
sulfO_di = pars.ksulfO * H2S_conc_di * O2_conc_di * pars.vol_di ;
sulfO_s  = pars.ksulfO * H2S_conc_s * O2_conc_s * pars.vol_s ;
sulfO_h  = pars.ksulfO * H2S_conc_h * O2_conc_h * pars.vol_h ;
sulfO_d  = pars.ksulfO * H2S_conc_d * O2_conc_d * pars.vol_d ;

SironR_p  = pars.kSironR * H2S_conc_p * FeIII_conc_p * pars.vol_p ;
SironR_di = pars.kSironR * H2S_conc_di * FeIII_conc_di * pars.vol_di ;
SironR_s  = pars.kSironR * H2S_conc_s * FeIII_conc_s * pars.vol_s ;
SironR_h  = pars.kSironR * H2S_conc_h * FeIII_conc_h * pars.vol_h ;
SironR_d  = pars.kSironR * H2S_conc_d * FeIII_conc_d * pars.vol_d ;

omegaside_p  = FeII_conc_p * CO3_p/pars.Kspside;
omegaside_di = FeII_conc_di * CO3_di/pars.Kspside;
omegaside_s  = FeII_conc_s * CO3_s/pars.Kspside;
omegaside_h  = FeII_conc_h * CO3_h/pars.Kspside;
omegaside_d  = FeII_conc_d * CO3_d/pars.Kspside;

SideP_p  = pars.kSide * max((omegaside_p-1),0) * pars.vol_p ;
SideP_di = pars.kSide * max((omegaside_di-1),0) * pars.vol_di ;
SideP_s  = pars.kSide * max((omegaside_s-1),0) * pars.vol_s ;
SideP_h  = pars.kSide * max((omegaside_h-1),0) * pars.vol_h ;
SideP_d  = pars.kSide * max((omegaside_d-1),0) * pars.vol_d ;

water_sediment_p  = 0.0476*POC_p ; 
water_sediment_di = 0.04*POC_di ; 
water_sediment_d  = 1.79e-4*POC_d;
BE_p   = 0.15; BE_di  = 0.15; BE_d   = 0.1;

mocb_p  = water_sediment_p*BE_p ;
mocb_di = water_sediment_di*BE_di ;
mocb_d  = water_sediment_d*BE_d ;

DICbenthic_p  = water_sediment_p - mocb_p;
DICbenthic_di = water_sediment_di - mocb_di;
DICbenthic_d  = water_sediment_d - mocb_d;

fanoxic_p = 1/(1+exp(-12*(0.5*pripr_p/560e12 - O2_a/pars.O2_a_0)));
CPratio_p = 250 ;
POPburial_p = mocb_p/CPratio_p ;
PFe_p   = 11.8e9*DP_p/pars.DP_p_0;
Pauth_p = 24e9*(0.3+exp(-3.5+3.15*O2_p/pars.O2_p_0)) ;

fanoxic_di = 1/(1+exp(-12*(0.5*pripr_p/560e12 - O2_a/pars.O2_a_0)));
CPratio_di = 250 ;
POPburial_di = mocb_di/CPratio_di ;
PFe_di   = 3.93e9*DP_di/pars.DP_di_0;
Pauth_di = 8e9*(0.3+exp(-3.5+3.15*O2_di/pars.O2_di_0)) ;

fanoxic_d = 1/(1+exp(-12*(0.5*(pripr_s+pripr_h)/3600e12 - O2_a/pars.O2_a_0)));
CPratio_d = 250 ;
POPburial_d = mocb_d/CPratio_d ;
PFe_d   = 7e9*DP_d/pars.DP_d_0;
Pauth_d = 14e9*(0.3+exp(-3.5+3.15*O2_d/pars.O2_d_0));

phosbenthic_p  = water_sediment_p/pars.Red_C_P - POPburial_p - PFe_p - Pauth_p ;
phosbenthic_di = water_sediment_di/pars.Red_C_P - POPburial_di - PFe_di - Pauth_di ;
phosbenthic_d  = water_sediment_d/pars.Red_C_P - POPburial_d - PFe_d - Pauth_d;

oxygenbentic_p  = -(water_sediment_p - mocb_p)/pars.Red_C_O;
oxygenbentic_di = -(water_sediment_di - mocb_di)/pars.Red_C_O;
oxygenbentic_d  = -(water_sediment_d - mocb_d)/pars.Red_C_O;

FeIIIscavenging_p  = max(3.18e+18*(FeIII_conc_p - 0.58e-6),0);
FeIIIscavenging_di = max(1.54e+18*(FeIII_conc_di - 0.58e-6),0);
FeIIIscavenging_s  = max(3.63e+18*(FeIII_conc_s - 0.58e-6),0);
FeIIIscavenging_h  = max(3.63e+18*(FeIII_conc_h - 0.58e-6),0);
FeIIIscavenging_d  = max(3.63e+18*(FeIII_conc_d - 0.58e-6),0);

pyritew  = pars.k_pyritew;
sulfatew = pars.k_sulfatew;
pyriteb_p  = pars.k_pyriteb_p*SO4_p/pars.SO4_p_0*mocb_p/4.585e12;
pyriteb_di = pars.k_pyriteb_di*SO4_p/pars.SO4_p_0*mocb_di/1.52e12;
sulfateb   = pars.k_sulfateb*SO4_p/pars.SO4_p_0 ;

COX_p    = (water_sediment_p - mocb_p)*1000/365.25/(3.6e14*0.05);
COX_di   = (water_sediment_di - mocb_di)*1000/365.25/(3.6e14*0.05);
COX_d    = (water_sediment_d - mocb_d)*1000/365.25/(3.6e14*0.9);
O2BW_p   = O2_conc_p *1000;
O2BW_di  = O2_conc_di *1000;
O2BW_d   = O2_conc_d *1000;

FeIIbenthic_p  = 4.1566e+12*tanh(COX_p/O2BW_p);
FeIIbenthic_di = 4.6888e+12*tanh(COX_di/O2BW_di);
FeIIbenthic_d  = 1.0382e+14*tanh(COX_d/O2BW_d);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Reservoir ODEs   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% CO2_a
dy(1) = - AirSea_p - AirSea_s - AirSea_h + ccdeg + ocdeg + oxidw - locb - carbw - 2*silw + CO2_input ; %%%% CO2_input = CO2_volcano

%%%% DIC_p
dy(2) = AirSea_p + Tran_DIC_di_p - Tran_DIC_p_s + Mix_DIC_di_p + 2*carbw + 2*silw - mccb_p - pripr_p + remin_p + DICbenthic_p - SideP_p;

%%%% DIC_di
dy(3) =  Tran_DIC_d_di - Tran_DIC_di_p - Mix_DIC_di_p - mccb_di  + remin_di + DICbenthic_di - SideP_di;

%%%% DIC_s
dy(4) = AirSea_s + Tran_DIC_d_s + Tran_DIC_p_s + Mix_DIC_d_s - Tran_DIC_s_h - pripr_s + remin_s - SideP_s;

%%%% DIC_h
dy(5) = AirSea_h + Tran_DIC_s_h + Mix_DIC_d_h - Tran_DIC_h_d - pripr_h + remin_h - SideP_h;

%%%% DIC_d
dy(6) =  Tran_DIC_h_d - Tran_DIC_d_s - Tran_DIC_d_di - Mix_DIC_d_s - Mix_DIC_d_h  - mccb_d + remin_d + DICbenthic_d - SideP_d;

%%%% ALK_p
dy(7)  =  Tran_ALK_di_p - Tran_ALK_p_s + Mix_ALK_di_p + 2*carbw + 2*silw - 2*mccb_p - 2*SideP_p   + 0.25*pyF_p - 2*ironO_p - sulfO_p + 15*SironR_p  -4*ironR_p + 0.5*sulfR_p - 3*FeIIIscavenging_p;

%%%% ALK_di
dy(8)  =  Tran_ALK_d_di - Tran_ALK_di_p - Mix_ALK_di_p - 2*mccb_di - 2*SideP_di   + 0.25*pyF_di - 2*ironO_di - sulfO_di + 15*SironR_di  -4*ironR_di + 0.5*sulfR_di - 3*FeIIIscavenging_di;

%%%% ALK_s
dy(9)  =  Tran_ALK_d_s + Tran_ALK_p_s - Tran_ALK_s_h + Mix_ALK_d_s - 2*SideP_s   + 0.25*pyF_s - 2*ironO_s - sulfO_s + 15*SironR_s -4*ironR_s + 0.5*sulfR_s - 3*FeIIIscavenging_s;

%%%% ALK_h
dy(10) =  Tran_ALK_s_h - Tran_ALK_h_d + Mix_ALK_d_h - 2*SideP_h   + 0.25*pyF_h - 2*ironO_h - sulfO_h + 15*SironR_h -4*ironR_h + 0.5*sulfR_h - 3*FeIIIscavenging_h;

%%%% ALK_d
dy(11) =  Tran_ALK_h_d - Tran_ALK_d_s - Tran_ALK_d_di - Mix_ALK_d_s - Mix_ALK_d_h - 2*mccb_d - 2*SideP_d + 0.25*pyF_d - 2*ironO_d - sulfO_d + 15*SironR_d  -4*ironR_d + 0.5*sulfR_d - 3*FeIIIscavenging_d;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Carbon isotopes   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d13c_C = 0 ;
d13c_G = -27 ;
capdelB = 27 ;
capdelDIC_g = -8 ;

dy(12) = - AirSea_ap*d13c_atm - AirSea_as*d13c_atm - AirSea_ah*d13c_atm + AirSea_pa*(d13c_DIC_p+capdelDIC_g) + AirSea_sa*(d13c_DIC_s+capdelDIC_g) + AirSea_ha*(d13c_DIC_h+capdelDIC_g)  + ccdeg*d13c_C + ocdeg*d13c_G + oxidw*d13c_G - locb*(d13c_atm - capdelB) - carbw*d13c_atm  - 2*silw*d13c_atm + CO2_input*d13c_CO2_input;
dy(13) = AirSea_ap*d13c_atm - AirSea_pa*(d13c_DIC_p+capdelDIC_g)  + Tran_DIC_di_p*d13c_DIC_di - Tran_DIC_p_s*d13c_DIC_p + mixcoeff_dip_m3_yr*(DIC_conc_di *d13c_DIC_di - DIC_conc_p*d13c_DIC_p) + carbw*d13c_C + carbw*d13c_atm + 2*silw*d13c_atm - mccb_p*d13c_DIC_p - pripr_p*(d13c_DIC_p - capdelB) + remin_p*d13c_POC_p + DICbenthic_p*d13c_POC_p - SideP_p*d13c_DIC_p ;
dy(14) = Tran_DIC_d_di*d13c_DIC_d - Tran_DIC_di_p*d13c_DIC_di - mixcoeff_dip_m3_yr*(DIC_conc_di *d13c_DIC_di - DIC_conc_p*d13c_DIC_p)  - mccb_di*d13c_DIC_di + remin_di*d13c_POC_di + DICbenthic_di*d13c_POC_di - SideP_di*d13c_DIC_di ;
dy(15) = AirSea_as*d13c_atm - AirSea_sa*(d13c_DIC_s+capdelDIC_g) + Tran_DIC_p_s*d13c_DIC_p + Tran_DIC_d_s*d13c_DIC_d - Tran_DIC_s_h*d13c_DIC_s + mixcoeff_ds_m3_yr*(DIC_conc_d *d13c_DIC_d - DIC_conc_s*d13c_DIC_s)  - pripr_s*(d13c_DIC_s - capdelB) + remin_s*d13c_POC_s - SideP_s*d13c_DIC_s;
dy(16) = AirSea_ah*d13c_atm - AirSea_ha*(d13c_DIC_h+capdelDIC_g)  + Tran_DIC_s_h*d13c_DIC_s - Tran_DIC_h_d*d13c_DIC_h +  mixcoeff_dh_m3_yr*(DIC_conc_d *d13c_DIC_d - DIC_conc_h*d13c_DIC_h) - pripr_h*(d13c_DIC_h - capdelB) + remin_h*d13c_POC_h - SideP_h*d13c_DIC_h;
dy(17) =  Tran_DIC_h_d*d13c_DIC_h - Tran_DIC_d_s*d13c_DIC_d  - Tran_DIC_d_di*d13c_DIC_d - mixcoeff_ds_m3_yr*(DIC_conc_d *d13c_DIC_d - DIC_conc_s*d13c_DIC_s) -  mixcoeff_dh_m3_yr*(DIC_conc_d *d13c_DIC_d - DIC_conc_h*d13c_DIC_h) - mccb_d*d13c_DIC_d + (remin_d+ DICbenthic_d)*d13c_POC_d - SideP_d*d13c_DIC_d;

%%%% POC & DP etc.
dy(18) = pripr_p - remin_p - Tran_POC_p_s - water_sediment_p - Tran_POC_p_di;
dy(19) =  - remin_di + Tran_POC_d_di + Tran_POC_p_di - water_sediment_di;
dy(20) = pripr_s - remin_s + Tran_POC_p_s - Tran_POC_s_h - Tran_POC_s_d ;
dy(21) = pripr_h - remin_h + Tran_POC_s_h - Tran_POC_h_d ;
dy(22) = - remin_d + Tran_POC_s_d + Tran_POC_h_d - Tran_POC_d_di - water_sediment_d ;

dy(23) = phosbenthic_p + phosw - Tran_DP_p_s + Tran_DP_di_p + Mix_DP_di_p  - pripr_p/pars.Red_C_P + remin_p/pars.Red_C_P ;
dy(24) = phosbenthic_di + Tran_DP_d_di - Tran_DP_di_p - Mix_DP_di_p  + remin_di/pars.Red_C_P ;
dy(25) = Tran_DP_p_s - Tran_DP_s_h + Tran_DP_d_s + Mix_DP_d_s - pripr_s/pars.Red_C_P + remin_s/pars.Red_C_P ;
dy(26) = Tran_DP_s_h - Tran_DP_h_d + Mix_DP_d_h  - pripr_h/pars.Red_C_P + remin_h/pars.Red_C_P ;
dy(27) = phosbenthic_d + Tran_DP_h_d - Tran_DP_d_di - Tran_DP_d_s - Mix_DP_d_s - Mix_DP_d_h  + remin_d/pars.Red_C_P ;

%%%% d13c_POC
dy(28) = pripr_p*(d13c_DIC_p - 27) - remin_p*d13c_POC_p - Tran_POC_p_s*d13c_POC_p - water_sediment_p*d13c_POC_p - Tran_POC_p_di*d13c_POC_p;
dy(29) =  - remin_di*d13c_POC_di + Tran_POC_d_di*d13c_POC_d + Tran_POC_p_di*d13c_POC_p - water_sediment_di*d13c_POC_di;
dy(30) = pripr_s*(d13c_DIC_s - 27) - remin_s*d13c_POC_s + Tran_POC_p_s*d13c_POC_p - Tran_POC_s_h*d13c_POC_s - Tran_POC_s_d*d13c_POC_s ;
dy(31) = pripr_h*(d13c_DIC_h - 27) - remin_h*d13c_POC_h + Tran_POC_s_h*d13c_POC_s - Tran_POC_h_d*d13c_POC_h ;
dy(32) = - remin_d*d13c_POC_d + Tran_POC_s_d*d13c_POC_s + Tran_POC_h_d*d13c_POC_h - Tran_POC_d_di*d13c_POC_d - water_sediment_d*d13c_POC_d ;

%%%% O2
dy(33) = - ocdeg/pars.Red_C_O - oxidw/pars.Red_C_O + locb/pars.Red_C_O + mocb_p/pars.Red_C_O + mocb_di/pars.Red_C_O + mocb_d/pars.Red_C_O ;
dy(34) = oxygenbentic_p + AirSeaO2_p + Tran_O2_di_p - Tran_O2_p_s + Mix_O2_di_p  + pripr_p_O - AR_p/pars.Red_C_O - 0.25*ironO_p - 2*sulfO_p;
dy(35) = oxygenbentic_di + Tran_O2_d_di - Tran_O2_di_p - Mix_O2_di_p - AR_di/pars.Red_C_O - 0.25*ironO_di - 2*sulfO_di;
dy(36) = AirSeaO2_s + Tran_O2_p_s + Tran_O2_d_s  - Tran_O2_s_h + Mix_O2_d_s + pripr_s/pars.Red_C_O - AR_s/pars.Red_C_O - 0.25*ironO_s - 2*sulfO_s;
dy(37) = AirSeaO2_h + Tran_O2_s_h - Tran_O2_h_d  + Mix_O2_d_h + pripr_h/pars.Red_C_O - AR_h/pars.Red_C_O - 0.25*ironO_h - 2*sulfO_h;
dy(38) = oxygenbentic_d + Tran_O2_h_d - Tran_O2_d_di - Tran_O2_d_s - Mix_O2_d_s - Mix_O2_d_h - AR_d/pars.Red_C_O - 0.25*ironO_d - 2*sulfO_d ;

%%%% Fe/S species
dy(39) = FeIIIw + Tran_FeIII_di_p - Tran_FeIII_p_s + Mix_FeIII_di_p - pripr_p/pars.Red_C_Fe + remin_p/pars.Red_C_Fe - FeIIIscavenging_p - 4*ironR_p + ironO_p - 8*SironR_p;
dy(40) = Tran_FeIII_d_di - Tran_FeIII_di_p - Mix_FeIII_di_p + remin_di/pars.Red_C_Fe - FeIIIscavenging_di - 4*ironR_di + ironO_di - 8*SironR_di;
dy(41) = pars.FeIIIa_s + Tran_FeIII_p_s - Tran_FeIII_s_h + Tran_FeIII_d_s + Mix_FeIII_d_s - pripr_s/pars.Red_C_Fe + remin_s/pars.Red_C_Fe - FeIIIscavenging_s - 4*ironR_s  + ironO_s - 8*SironR_s;
dy(42) = pars.FeIIIa_h + Tran_FeIII_s_h - Tran_FeIII_h_d + Mix_FeIII_d_h - pripr_h/pars.Red_C_Fe + remin_h/pars.Red_C_Fe - FeIIIscavenging_h - 4*ironR_h  + ironO_h - 8*SironR_h;
dy(43) = - Tran_FeIII_d_di - Tran_FeIII_d_s + Tran_FeIII_h_d - Mix_FeIII_d_s - Mix_FeIII_d_h + remin_d/pars.Red_C_Fe - FeIIIscavenging_d - 4*ironR_d  + ironO_d - 8*SironR_d;

dy(44) = pyritew + sulfatew - pyriteb_p - sulfateb + Tran_SO4_di_p - Tran_SO4_p_s + Mix_SO4_di_p - 0.5*sulfR_p + sulfO_p + SironR_p - 0.25*pyF_p;
dy(45) = - pyriteb_di + Tran_SO4_d_di - Tran_SO4_di_p - Mix_SO4_di_p - 0.5*sulfR_di + sulfO_di + SironR_di - 0.25*pyF_di;
dy(46) = Tran_SO4_p_s - Tran_SO4_s_h + Tran_SO4_d_s + Mix_SO4_d_s - 0.5*sulfR_s + sulfO_s + SironR_s - 0.25*pyF_s;
dy(47) = Tran_SO4_s_h - Tran_SO4_h_d + Mix_SO4_d_h - 0.5*sulfR_h + sulfO_h + SironR_h - 0.25*pyF_h;
dy(48) = Tran_SO4_h_d - Tran_SO4_d_di - Tran_SO4_d_s - Mix_SO4_d_s - Mix_SO4_d_h - 0.5*sulfR_d + sulfO_d + SironR_d - 0.25*pyF_d;

dy(49) = FeIIbenthic_p + Tran_FeII_di_p - Tran_FeII_p_s + Mix_FeII_di_p + 4*ironR_p - pyF_p - ironO_p + 8*SironR_p - SideP_p;
dy(50) = FeIIbenthic_di + Tran_FeII_d_di - Tran_FeII_di_p - Mix_FeII_di_p + 4*ironR_di - pyF_di - ironO_di + 8*SironR_di - SideP_di;
dy(51) = Tran_FeII_p_s - Tran_FeII_s_h + Tran_FeII_d_s + Mix_FeII_d_s + 4*ironR_s - pyF_s - ironO_s + 8*SironR_s - SideP_s;
dy(52) = Tran_FeII_s_h - Tran_FeII_h_d + Mix_FeII_d_h + 4*ironR_h - pyF_h - ironO_h + 8*SironR_h - SideP_h;
dy(53) = FeIIbenthic_d + pars.FeIIhydro + Tran_FeII_h_d - Tran_FeII_d_di - Tran_FeII_d_s - Mix_FeII_d_s - Mix_FeII_d_h + 4*ironR_d - pyF_d - ironO_d + 8*SironR_d - SideP_d;

dy(54) = Tran_H2S_di_p - Tran_H2S_p_s + Mix_H2S_di_p + 0.5*sulfR_p - 1.75*pyF_p - sulfO_p - SironR_p;
dy(55) = Tran_H2S_d_di - Tran_H2S_di_p - Mix_H2S_di_p  + 0.5*sulfR_di - 1.75*pyF_di - sulfO_di - SironR_di;
dy(56) = Tran_H2S_p_s - Tran_H2S_s_h + Tran_H2S_d_s + Mix_H2S_d_s + 0.5*sulfR_s - 1.75*pyF_s - sulfO_s - SironR_s;
dy(57) = Tran_H2S_s_h - Tran_H2S_h_d + Mix_H2S_d_h + 0.5*sulfR_h - 1.75*pyF_h - sulfO_h - SironR_h;
dy(58) = Tran_H2S_h_d - Tran_H2S_d_di - Tran_H2S_d_s - Mix_H2S_d_s - Mix_H2S_d_h + 0.5*sulfR_d - 1.75*pyF_d - sulfO_d - SironR_d;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Save output as working   %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
workingstate.CO2_input(stepnumber,1)          = CO2_input ;
workingstate.d13c_CO2_input(stepnumber,1)     = d13c_CO2_input ;
workingstate.lake_P_release(stepnumber,1)     = lake_P_release_yr;  % mol/yr
workingstate.lake_C_burial(stepnumber,1)      = lake_C_burial_yr;   % mol/yr
workingstate.CO2_a(stepnumber,1)              = CO2_a ;
workingstate.DIC_p(stepnumber,1)              = DIC_p ;
workingstate.DIC_di(stepnumber,1)             = DIC_di ;
workingstate.DIC_s(stepnumber,1)              = DIC_s ;
workingstate.DIC_h(stepnumber,1)              = DIC_h ;
workingstate.DIC_d(stepnumber,1)              = DIC_d ;
workingstate.ALK_p(stepnumber,1) = ALK_p ;
workingstate.ALK_di(stepnumber,1) = ALK_di ;
workingstate.ALK_s(stepnumber,1) = ALK_s ;
workingstate.ALK_h(stepnumber,1) = ALK_h ;
workingstate.ALK_d(stepnumber,1) = ALK_d ;
workingstate.DIC_conc_p(stepnumber,1) = DIC_conc_p ;
workingstate.DIC_conc_di(stepnumber,1) = DIC_conc_di ;
workingstate.DIC_conc_s(stepnumber,1) = DIC_conc_s ;
workingstate.DIC_conc_h(stepnumber,1) = DIC_conc_h ;
workingstate.DIC_conc_d(stepnumber,1) = DIC_conc_d ;
workingstate.ALK_conc_p(stepnumber,1) = ALK_conc_p ;
workingstate.ALK_conc_di(stepnumber,1) = ALK_conc_di ;
workingstate.ALK_conc_s(stepnumber,1) = ALK_conc_s ;
workingstate.ALK_conc_h(stepnumber,1) = ALK_conc_h ;
workingstate.ALK_conc_d(stepnumber,1) = ALK_conc_d ;
workingstate.pH_p(stepnumber,1) = pH_p ;
workingstate.pH_di(stepnumber,1) = pH_di ;
workingstate.pH_s(stepnumber,1) = pH_s ;
workingstate.pH_h(stepnumber,1) = pH_h ;
workingstate.pH_d(stepnumber,1) = pH_d ;
workingstate.T_p(stepnumber,1) = T_p ;
workingstate.T_di(stepnumber,1) = T_di ;
workingstate.T_s(stepnumber,1) = T_s ;
workingstate.T_h(stepnumber,1) = T_h ;
workingstate.T_d(stepnumber,1) = T_d ;
workingstate.T_cont(stepnumber,1) = T_cont ;
workingstate.GAST(stepnumber,1) = GAST ;
workingstate.Atmospheric_CO2_ppm(stepnumber,1) = Atmospheric_CO2_ppm ;
workingstate.granw(stepnumber,1) = granw ;
workingstate.basw(stepnumber,1) = basw ;
workingstate.silw(stepnumber,1) = silw ;
workingstate.carbw(stepnumber,1) = carbw ;
workingstate.ccdeg(stepnumber,1) = ccdeg ;
workingstate.phosw(stepnumber,1) = phosw ;
workingstate.mccb_p(stepnumber,1) = mccb_p ;
workingstate.mccb_di(stepnumber,1) = mccb_di ;
workingstate.mccb_d(stepnumber,1) = mccb_d ;
workingstate.mocb_p(stepnumber,1) = mocb_p ;
workingstate.mocb_di(stepnumber,1) = mocb_di ;
workingstate.mocb_d(stepnumber,1) = mocb_d ;
workingstate.HCO3_p(stepnumber,1) = HCO3_p ;
workingstate.CO3_p(stepnumber,1) = CO3_p ;
workingstate.H_p(stepnumber,1) = H_p ;
workingstate.HCO3_di(stepnumber,1) = HCO3_di ;
workingstate.CO3_di(stepnumber,1) = CO3_di ;
workingstate.H_di(stepnumber,1) = H_di ;
workingstate.HCO3_s(stepnumber,1) = HCO3_s ;
workingstate.CO3_s(stepnumber,1) = CO3_s ;
workingstate.H_s(stepnumber,1) = H_s ;
workingstate.HCO3_h(stepnumber,1) = HCO3_h ;
workingstate.CO3_h(stepnumber,1) = CO3_h ;
workingstate.H_h(stepnumber,1) = H_h ;
workingstate.HCO3_d(stepnumber,1) = HCO3_d ;
workingstate.CO3_d(stepnumber,1) = CO3_d ;
workingstate.H_d(stepnumber,1) = H_d ;
workingstate.d13c_atm(stepnumber,1) = d13c_atm ;
workingstate.d13c_DIC_p(stepnumber,1) = d13c_DIC_p ;
workingstate.d13c_DIC_di(stepnumber,1) = d13c_DIC_di ;
workingstate.d13c_DIC_s(stepnumber,1) = d13c_DIC_s ;
workingstate.d13c_DIC_h(stepnumber,1) = d13c_DIC_h ;
workingstate.d13c_DIC_d(stepnumber,1) = d13c_DIC_d ;
workingstate.POC_p(stepnumber,1) = POC_p ;
workingstate.POC_di(stepnumber,1) = POC_di ;
workingstate.POC_s(stepnumber,1) = POC_s ;
workingstate.POC_h(stepnumber,1) = POC_h ;
workingstate.POC_d(stepnumber,1) = POC_d ;
workingstate.POC_conc_p(stepnumber,1) = POC_conc_p ;
workingstate.POC_conc_di(stepnumber,1) = POC_conc_di ;
workingstate.POC_conc_s(stepnumber,1) = POC_conc_s ;
workingstate.POC_conc_h(stepnumber,1) = POC_conc_h ;
workingstate.POC_conc_d(stepnumber,1) = POC_conc_d ;
workingstate.DP_p(stepnumber,1) = DP_p ;
workingstate.DP_di(stepnumber,1) = DP_di ;
workingstate.DP_s(stepnumber,1) = DP_s ;
workingstate.DP_h(stepnumber,1) = DP_h ;
workingstate.DP_d(stepnumber,1) = DP_d ;
workingstate.DP_conc_p(stepnumber,1) = DP_conc_p ;
workingstate.DP_conc_di(stepnumber,1) = DP_conc_di ;
workingstate.DP_conc_s(stepnumber,1) = DP_conc_s ;
workingstate.DP_conc_h(stepnumber,1) = DP_conc_h ;
workingstate.DP_conc_d(stepnumber,1) = DP_conc_d ;
workingstate.d13c_POC_p(stepnumber,1) = d13c_POC_p ;
workingstate.d13c_POC_di(stepnumber,1) = d13c_POC_di ;
workingstate.d13c_POC_s(stepnumber,1) = d13c_POC_s ;
workingstate.d13c_POC_h(stepnumber,1) = d13c_POC_h ;
workingstate.d13c_POC_d(stepnumber,1) = d13c_POC_d ;

workingstate.pO2_a(stepnumber,1) = pO2_a ;
workingstate.O2_conc_p(stepnumber,1) = O2_conc_p ;
workingstate.O2_conc_di(stepnumber,1) = O2_conc_di ;
workingstate.O2_conc_s(stepnumber,1) = O2_conc_s ;
workingstate.O2_conc_h(stepnumber,1) = O2_conc_h ;
workingstate.O2_conc_d(stepnumber,1) = O2_conc_d ;

workingstate.FeIII_conc_p(stepnumber,1) = FeIII_conc_p ;
workingstate.FeIII_conc_di(stepnumber,1) = FeIII_conc_di ;
workingstate.FeIII_conc_s(stepnumber,1) = FeIII_conc_s ;
workingstate.FeIII_conc_h(stepnumber,1) = FeIII_conc_h ;
workingstate.FeIII_conc_d(stepnumber,1) = FeIII_conc_d ;

workingstate.pripr_p(stepnumber,1) = pripr_p ;
workingstate.pripr_s(stepnumber,1) = pripr_s ;
workingstate.pripr_h(stepnumber,1) = pripr_h ;

workingstate.remin_p(stepnumber,1) = remin_p ;
workingstate.remin_di(stepnumber,1) = remin_di ;
workingstate.remin_s(stepnumber,1) = remin_s ;
workingstate.remin_h(stepnumber,1) = remin_h ;
workingstate.remin_d(stepnumber,1) = remin_d ;

workingstate.SO4_conc_p(stepnumber,1) = SO4_conc_p ;
workingstate.SO4_conc_di(stepnumber,1) = SO4_conc_di ;
workingstate.SO4_conc_s(stepnumber,1) = SO4_conc_s ;
workingstate.SO4_conc_h(stepnumber,1) = SO4_conc_h ;
workingstate.SO4_conc_d(stepnumber,1) = SO4_conc_d ;

workingstate.FeII_conc_p(stepnumber,1) = FeII_conc_p ;
workingstate.FeII_conc_di(stepnumber,1) = FeII_conc_di ;
workingstate.FeII_conc_s(stepnumber,1) = FeII_conc_s ;
workingstate.FeII_conc_h(stepnumber,1) = FeII_conc_h ;
workingstate.FeII_conc_d(stepnumber,1) = FeII_conc_d ;

workingstate.H2S_conc_p(stepnumber,1) = H2S_conc_p ;
workingstate.H2S_conc_di(stepnumber,1) = H2S_conc_di ;
workingstate.H2S_conc_s(stepnumber,1) = H2S_conc_s ;
workingstate.H2S_conc_h(stepnumber,1) = H2S_conc_h ;
workingstate.H2S_conc_d(stepnumber,1) = H2S_conc_d ;

workingstate.Tran_FeIII_p_s(stepnumber,1) = Tran_FeIII_p_s ;
workingstate.Tran_FeIII_s_h(stepnumber,1) = Tran_FeIII_s_h ;
workingstate.Tran_FeIII_d_s(stepnumber,1) = Tran_FeIII_d_s ;
workingstate.Mix_FeIII_d_s(stepnumber,1) = Mix_FeIII_d_s ;

workingstate.O2_a(stepnumber,1) = O2_a ;
workingstate.O2_p(stepnumber,1) = O2_p ;
workingstate.O2_di(stepnumber,1) = O2_di ;
workingstate.O2_s(stepnumber,1) = O2_s ;
workingstate.O2_h(stepnumber,1) = O2_h ;
workingstate.O2_d(stepnumber,1) = O2_d ;

workingstate.FeIII_p(stepnumber,1) = FeIII_p ;
workingstate.FeIII_di(stepnumber,1) = FeIII_di ;
workingstate.FeIII_s(stepnumber,1) = FeIII_s ;
workingstate.FeIII_h(stepnumber,1) = FeIII_h ;
workingstate.FeIII_d(stepnumber,1) = FeIII_d ;

workingstate.SO4_p(stepnumber,1) = SO4_p ;
workingstate.SO4_di(stepnumber,1) = SO4_di ;
workingstate.SO4_s(stepnumber,1) = SO4_s ;
workingstate.SO4_h(stepnumber,1) = SO4_h ;
workingstate.SO4_d(stepnumber,1) = SO4_d ;

workingstate.FeII_p(stepnumber,1) = FeII_p ;
workingstate.FeII_di(stepnumber,1) = FeII_di ;
workingstate.FeII_s(stepnumber,1) = FeII_s ;
workingstate.FeII_h(stepnumber,1) = FeII_h ;
workingstate.FeII_d(stepnumber,1) = FeII_d ;

workingstate.H2S_p(stepnumber,1) = H2S_p ;
workingstate.H2S_di(stepnumber,1) = H2S_di ;
workingstate.H2S_s(stepnumber,1) = H2S_s ;
workingstate.H2S_h(stepnumber,1) = H2S_h ;
workingstate.H2S_d(stepnumber,1) = H2S_d ;

workingstate.DICbenthic_p(stepnumber,1) = DICbenthic_p ;
workingstate.DICbenthic_di(stepnumber,1) = DICbenthic_di ;
workingstate.DICbenthic_d(stepnumber,1) = DICbenthic_d ;
workingstate.FeIIbenthic_p(stepnumber,1) = FeIIbenthic_p ;
workingstate.FeIIbenthic_d(stepnumber,1) = FeIIbenthic_d ;

coeff_p  =  62.2e9/tanh(COX_p/O2BW_p);
coeff_di =  30.5e9/tanh(COX_di/O2BW_di);
coeff_d  =  57.6e9/tanh(COX_d/O2BW_d);
workingstate.coeff_p(stepnumber,1)  = coeff_p;
workingstate.coeff_di(stepnumber,1) = coeff_di ;
workingstate.coeff_d(stepnumber,1)  = coeff_d ;

workingstate.time(stepnumber,1)     = t ;
workingstate.time_myr(stepnumber,1) = t / 1e6 ;
stepnumber = stepnumber + 1 ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%   Nested helper: event total C (mol)   %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MolC_total = event_total_C_mol_from_T0cm(T0_cm, P)
    % 支持 width_m 或 width_km
    if isfield(P.volc,'width_m')
        width_m = P.volc.width_m;
    elseif isfield(P.volc,'width_km')
        width_m = P.volc.width_km * 1e3;
    else
        width_m = 300e3;
    end

    kexp = P.volc.k_exp; % 1/km

    x1 = (P.volc.r_near_km(1):P.volc.dx_km:P.volc.r_near_km(2))'; if numel(x1)<2, x1 = P.volc.r_near_km(1:2)'; end
    x2 = (P.volc.r_far_km(1):P.volc.dx_km:P.volc.r_far_km(2))';   if numel(x2)<2, x2 = P.volc.r_far_km(1:2)'; end
    dx1_m = (x1(2)-x1(1)) * 1e3;
    dx2_m = (x2(2)-x2(1)) * 1e3;

    T0_m   = T0_cm / 100;  % cm→m
    T_near = T0_m * exp(-kexp * x1);
    T_far  = T0_m * exp(-kexp * x2);

    volume_m3  = (sum(T_near)*dx1_m + sum(T_far)*dx2_m) * width_m;
    volume_km3 = volume_m3 / 1e9;
    mass_kg    = volume_m3 * P.volc.ash_density;

    % P → C 有机埋藏；以及物理吸附
    P_mass_kg  = mass_kg * (P.volc.P_mean_ppm/1e6) * P.volc.P_release_frac;
    P_mol      = P_mass_kg / P.volc.MW_P;
    C_org_mol  = P_mol * P.volc.CP_molar * P.volc.BE_lake;

    C_ads_kg   = P.volc.C_adsorb_ton_per_km3 * 1e3 * volume_km3;  % kg
    C_ads_mol  = (C_ads_kg * 1e3) / 12;                           % g / (g/mol)

    MolC_total = C_org_mol + C_ads_mol;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%   Nested helper: event total P (mol)   %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MolP_total = event_total_P_mol_from_T0cm(T0_cm, P)
    if isfield(P.volc,'width_m')
        width_m = P.volc.width_m;
    elseif isfield(P.volc,'width_km')
        width_m = P.volc.width_km * 1e3;
    else
        width_m = 300e3;
    end

    kexp = P.volc.k_exp; % 1/km
    x1 = (P.volc.r_near_km(1):P.volc.dx_km:P.volc.r_near_km(2))'; if numel(x1)<2, x1 = P.volc.r_near_km(1:2)'; end
    x2 = (P.volc.r_far_km(1):P.volc.dx_km:P.volc.r_far_km(2))';   if numel(x2)<2, x2 = P.volc.r_far_km(1:2)'; end
    dx1_m = (x1(2)-x1(1)) * 1e3;
    dx2_m = (x2(2)-x2(1)) * 1e3;

    T0_m   = T0_cm / 100;  % cm→m
    T_near = T0_m * exp(-kexp * x1);
    T_far  = T0_m * exp(-kexp * x2);

    volume_m3  = (sum(T_near)*dx1_m + sum(T_far)*dx2_m) * width_m;
    mass_kg    = volume_m3 * P.volc.ash_density;

    P_mass_kg  = mass_kg * (P.volc.P_mean_ppm/1e6) * P.volc.P_release_frac;
    MolP_total = P_mass_kg / P.volc.MW_P;
end

end
