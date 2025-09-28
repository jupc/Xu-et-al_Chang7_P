%%%%%% MBOX 5 box ocean model Zhao et al 2023

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%8e12%8e12
%%%%%%%%%%%%%%%%%   Define parameters   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% set up global structures
clearvars -global
global stepnumber
global pars
global workingstate

%%%%%% water reservoir sizes in m3 (m=margins, s=surface, h= hi-lat, d=deep)
pars.vol_p  = 2.6e15 ;  %%%% approx volume of all shelves and slope to depth 100m, area pecentage 5%
pars.vol_di = 5.4e15 ;  %%%% approx volume of all shelves and slope in depth 100-1000m, area pecentage 5%
pars.vol_s  = 2.75e16 ; %%%% approx volume of suface water to depth 100m, area pecentage 76.5%
pars.vol_h  = 1.22e16 ; %%%% approx volume of hi-lat to depth 250m, area pecentage 13.5%
pars.vol_d  = 1.35e18 ;
pars.vol_ocean = pars.vol_p + pars.vol_di + pars.vol_s + pars.vol_h + pars.vol_d ;

%%%% mixing coefficient (Sv)
pars.mixcoeff_dip = 30.28;
pars.mixcoeff_ds  = 46.33;
pars.mixcoeff_dh  = 54.9;

%%%%%% inorganic carbon reservoirs in moles C
pars.CO2_a_0  = 5e16 ; %%%5e16
pars.DIC_p_0  = 5.2e15 ;
pars.DIC_di_0 = 1.08e16 ;
pars.DIC_s_0  = 5.37e16 ;
pars.DIC_h_0  = 2.71e16 ;
pars.DIC_d_0  = 3e18 ;
pars.ALK_p_0  = 5.2e15 ;
pars.ALK_di_0 = 1.08e16 ;
pars.ALK_s_0  = 5.37e16 ;
pars.ALK_h_0  = 2.71e16 ;
pars.ALK_d_0  = 3e18 ;

%%%%%% C isotope composition
pars.d13c_atm_0    = -7 ;
pars.d13c_DIC_p_0  = 0.1 ;
pars.d13c_DIC_di_0 = 0.1 ;
pars.d13c_DIC_s_0  = 0.1 ;
pars.d13c_DIC_h_0  = 0.1 ;
pars.d13c_DIC_d_0  = 0.1 ;

%%%%%% POC reservoirs in moles C
pars.POC_p_0  = 630e12 ;
pars.POC_di_0 = 250e12 ;
pars.POC_s_0  = 2329e12 ;
pars.POC_h_0  = 1084e12 ;
pars.POC_d_0  = 56000e12 ;

%%%%%% Dissolved phosphate in moles
pars.DP_p_0  = 1.82e12 ;
pars.DP_di_0 = 7.56e12 ;
pars.DP_s_0  = 0.55e12 ;
pars.DP_h_0  = 16.27e12 ;
pars.DP_d_0  = 2970e12 ;

pars.d13c_POC_p_0  = -26;
pars.d13c_POC_di_0 = -26;
pars.d13c_POC_s_0  = -26 ;
pars.d13c_POC_h_0  = -26 ;
pars.d13c_POC_d_0  = -26 ;

%%%%%% initial amount of O2 in moles C
pars.O2_a_0  = 3.7e19 ;
pars.O2_p_0  = 6.705e14 ;
pars.O2_di_0 = 8.964e14 ;
pars.O2_s_0  = 9.139e15 ;
pars.O2_h_0  = 4.02e15 ;
pars.O2_d_0  = 1.823e17 ;

%%%%%% initial amount of FeIII in moles
pars.FeIII_p_0  = 1.56e9 ;
pars.FeIII_di_0 = 3.24e9 ;
pars.FeIII_s_0  = 9.625e9 ;
pars.FeIII_h_0  = 3.66e9 ;
pars.FeIII_d_0  = 810e9 ;

%%%%%% initial amount of sulfate in moles
pars.SO4_p_0  = 7.28e16;
pars.SO4_di_0 = 1.512e17;
pars.SO4_s_0  = 7.7e17;
pars.SO4_h_0  = 3.416e17;
pars.SO4_d_0  = 3.78e19;

%%%%%% initial amount of FeII in moles
pars.FeII_p_0  = 0;
pars.FeII_di_0 = 0;
pars.FeII_s_0  = 0;
pars.FeII_h_0  = 0;
pars.FeII_d_0  = 0;

%%%%%% initial amount of H2S in moles
pars.H2S_p_0  = 0;
pars.H2S_di_0 = 0;
pars.H2S_s_0  = 0;
pars.H2S_h_0  = 0;
pars.H2S_d_0  = 0;

%%%%%% present day rates
pars.k_ccdeg = 13e12 ; %8e12
pars.k_carbw = 12e12 ; %12e12
pars.k_sfw = 0 ; %%% Seafloor weathering
pars.k_mccb = 20e12 ;
pars.k_silw = pars.k_mccb - pars.k_carbw ;
basfrac = 0.3 ;
pars.k_granw = pars.k_silw * (1-basfrac) ;
pars.k_basw = pars.k_silw * basfrac ;

%%%%%% organic C cycle
pars.k_ocdeg = 1.25e12 ;
pars.k_locb = 2.5e12 ;
pars.k_mocb = 7e12  ;
pars.k_oxidw = pars.k_mocb + pars.k_locb - pars.k_ocdeg ;

%%%%%% present P, Fe, pyrite and sulfate weathering rate
pars.k_phosw = 0.0967e12;
pars.k_FeIIIw = 2e9;
pars.k_pyritew = 1.85e12;
pars.k_sulfatew = 1.25e12;

%%%%%% present pyrite and sulfate burial rate
pars.k_pyriteb_p  = 1.4e12;
pars.k_pyriteb_di = 0.45e12;
pars.k_sulfateb   = 1.25e12;

%%%%%% Redfield ratio & biogeo params
pars.Red_C_P = 106;
pars.Red_C_N = 106/16;
pars.Red_C_O = 106/138;
pars.Red_C_Fe = 106*2000;
pars.Red_Fe_P = pars.Red_C_P/pars.Red_C_Fe;

pars.BE_p  = 6/40;
pars.BE_di = 6/40;
pars.BE_d  = 1/10;

pars.KP = 0.1e-3;
pars.KFe = 0.1e-6;
pars.KmO2 = 8e-3;
pars.KmFeIII = 10;

%%%%%% reaction rate constant m3 mol-1 yr-1
pars.kpy     = 0.3708/1e3*24*365.25;
pars.kironO  = 1.4e5;
pars.ksulfO  = 1.6e2;
pars.kSironR = 8;
pars.kSide   = 4e3; % mol m-3 yr-1

%%%%%% Ksp
pars.Kspside  = 10^(-8.4)*1e6;   % mol^2 m^-6
pars.KspFeSaq = 10^(-5.08)*1e3;  % mol m^-3
pars.STFeSaq  = 10^(-5.7)*1e3;   % mol m^-3

pars.FeIIIa_s  = 1.443e9;
pars.FeIIIa_h  = 0;
pars.FeIIhydro = 13.5e9;

%%%%%% Volcano–ash to CO2 drawdown (time window & geometry & chemistry)
pars.whenstart = -10e6;
pars.whenend   =  15e6;

% time series
pars.volc.t_T0  = [0, 1e6, 5e6, 7e6];
pars.volc.T0_cm = [0,  3,    3,    0];

% Geometric configuration
pars.volc.width_m   = 300e3;     % 
pars.volc.k_exp     = 0.01;      % 1/km
pars.volc.r_near_km = [0, 150];
pars.volc.r_far_km  = [150, 300];
pars.volc.dx_km     = 0.5;

% 
pars.volc.ash_density    = 1.4e3;                 % kg/m^3
pars.volc.P_mean_ppm     = (785+2032)/2;          % 785-2032
pars.volc.P_release_frac = 0.55;
pars.volc.MW_P           = 30.97e-3;              % kg/mol
pars.volc.C_adsorb_ton_per_km3 = 1.8e7;           % t/km^3 

% Lake area
pars.volc.lake_area_m2  = 300e3 * 300e3;
pars.volc.lake_depth_m  = 90;
pars.volc.CP_molar      = 400;     
pars.volc.BE_lake       = 0.30;

% time response
pars.volc.tau_release_yr = 1e3;
pars.volc.coupling       = 'input';
pars.volc.use_lake_trap_in_phosw = 1;

% 
tt = pars.volc.t_T0(:); vv = pars.volc.T0_cm(:);
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
        win_t1(end+1,1) = t1; %#ok<SAGROW>
        win_t2(end+1,1) = t2; %#ok<SAGROW>
        win_T0(end+1,1) = T;  %#ok<SAGROW>
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Initialise   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options = odeset('maxstep',1e3) ;
fprintf('Beginning run: \n')
stepnumber = 1 ;

%%%% model start state
pars.startstate(1) = pars.CO2_a_0 ;
pars.startstate(2) = pars.DIC_p_0 ;
pars.startstate(3) = pars.DIC_di_0 ;
pars.startstate(4) = pars.DIC_s_0 ;
pars.startstate(5) = pars.DIC_h_0 ;
pars.startstate(6) = pars.DIC_d_0 ;

pars.startstate(7)  = pars.ALK_p_0 ;
pars.startstate(8)  = pars.ALK_di_0 ;
pars.startstate(9)  = pars.ALK_s_0 ;
pars.startstate(10) = pars.ALK_h_0 ;
pars.startstate(11) = pars.ALK_d_0 ;

pars.startstate(12) = pars.CO2_a_0 * pars.d13c_atm_0 ;
pars.startstate(13) = pars.DIC_p_0 * pars.d13c_DIC_p_0 ;
pars.startstate(14) = pars.DIC_di_0 * pars.d13c_DIC_di_0 ;
pars.startstate(15) = pars.DIC_s_0 * pars.d13c_DIC_s_0 ;
pars.startstate(16) = pars.DIC_h_0 * pars.d13c_DIC_h_0 ;
pars.startstate(17) = pars.DIC_d_0 * pars.d13c_DIC_d_0 ;

pars.startstate(18) = pars.POC_p_0;
pars.startstate(19) = pars.POC_di_0;
pars.startstate(20) = pars.POC_s_0;
pars.startstate(21) = pars.POC_h_0;
pars.startstate(22) = pars.POC_d_0;

pars.startstate(23) = pars.DP_p_0;
pars.startstate(24) = pars.DP_di_0;
pars.startstate(25) = pars.DP_s_0;
pars.startstate(26) = pars.DP_h_0;
pars.startstate(27) = pars.DP_d_0;

pars.startstate(28) = pars.d13c_POC_p_0*pars.POC_p_0;
pars.startstate(29) = pars.d13c_POC_di_0*pars.POC_di_0;
pars.startstate(30) = pars.d13c_POC_s_0*pars.POC_s_0;
pars.startstate(31) = pars.d13c_POC_h_0*pars.POC_h_0;
pars.startstate(32) = pars.d13c_POC_d_0*pars.POC_d_0;

pars.startstate(33) = pars.O2_a_0;
pars.startstate(34) = pars.O2_p_0;
pars.startstate(35) = pars.O2_di_0;
pars.startstate(36) = pars.O2_s_0;
pars.startstate(37) = pars.O2_h_0;
pars.startstate(38) = pars.O2_d_0;

pars.startstate(39) = pars.FeIII_p_0;
pars.startstate(40) = pars.FeIII_di_0;
pars.startstate(41) = pars.FeIII_s_0;
pars.startstate(42) = pars.FeIII_h_0;
pars.startstate(43) = pars.FeIII_d_0;

pars.startstate(44) = pars.SO4_p_0;
pars.startstate(45) = pars.SO4_di_0;
pars.startstate(46) = pars.SO4_s_0;
pars.startstate(47) = pars.SO4_h_0;
pars.startstate(48) = pars.SO4_d_0;

pars.startstate(49) = pars.FeII_p_0;
pars.startstate(50) = pars.FeII_di_0;
pars.startstate(51) = pars.FeII_s_0;
pars.startstate(52) = pars.FeII_h_0;
pars.startstate(53) = pars.FeII_d_0;

pars.startstate(54) = pars.H2S_p_0;
pars.startstate(55) = pars.H2S_di_0;
pars.startstate(56) = pars.H2S_s_0;
pars.startstate(57) = pars.H2S_h_0;
pars.startstate(58) = pars.H2S_d_0;

%%%%%%% run the system 
[rawoutput.T,rawoutput.Y] = ode15s(@MBOX_equations,[pars.whenstart pars.whenend],pars.startstate,options);

tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Postprocessing   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pars.output_length = length(rawoutput.T) ;
fprintf('Integration finished \t') ; fprintf('Total steps: %d \t' , stepnumber ) ; fprintf('Output steps: %d \n' , pars.output_length ) 
toc

fprintf('assembling state vectors... \t')
tic
[sharedvals,trecords] = intersect(workingstate.time,rawoutput.T,'stable') ;

field_names = fieldnames(workingstate) ;
for numfields = 1:length(field_names)
    eval([' state.' char( field_names(numfields) ) ' = workingstate.' char( field_names(numfields) ) '(trecords) ; '])
end

fprintf('Done: ')
endtime = toc ;
fprintf('time (s): %d \n', endtime )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Plotting script   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
subplot(3,3,1); hold on; box on
plot(state.time_myr,state.CO2_input,'k')
xlim([-2 15]); ylim([0 62e12])
xlabel('Time (Ma)'); ylabel('CO_{2} input (mol/yr)')
title('a', 'FontSize', 12)

subplot(3,3,2); hold on; box on
plot(state.time_myr,state.phosw,'k')
xlim([-2 15])
xlabel('Time (Ma)'); ylabel('P weathering flux (mol/yr)')
title('b', 'FontSize', 12)

subplot(3,3,3); hold on; box on
plot(state.time_myr,state.Atmospheric_CO2_ppm,'k')
xlim([-2 15])
xlabel('Time (Ma)'); ylabel('Atm. CO_{2} (ppm)')
title('c', 'FontSize', 12)

subplot(3,3,4); hold on; box on
plot(state.time_myr,state.pH_p,'k')
plot(state.time_myr,state.pH_di,'b')
plot(state.time_myr,state.pH_s,'r')
plot(state.time_myr,state.pH_h,'c')
plot(state.time_myr,state.pH_d,'m')
xlim([-2 15])
xlabel('Time (Ma)'); ylabel('pH')
title('d', 'FontSize', 12)

subplot(3,3,5); hold on; box on
plot(state.time_myr,state.d13c_DIC_p+2.5,'k')
plot(state.time_myr,state.d13c_DIC_di+2.5,'b')
plot(state.time_myr,state.d13c_DIC_s+2.5,'r')
plot(state.time_myr,state.d13c_DIC_h+2.5,'c')
plot(state.time_myr,state.d13c_DIC_d+2.5,'m')
xlim([-2 15])
xlabel('Time (Ma)'); ylabel(['\delta^{13}C of DIC (', char(8240), ')'])
title('e', 'FontSize', 12)

subplot(3,3,6); hold on; box on
plot(state.time_myr,state.T_p,'k')
plot(state.time_myr,state.T_di,'b')
plot(state.time_myr,state.T_s,'r')
plot(state.time_myr,state.T_h,'c')
plot(state.time_myr,state.T_d,'m')
xlim([-2 15])
xlabel('Time (Ma)'); ylabel('T (^{0}C)')
title('f', 'FontSize', 12)

subplot(3,3,7); hold on; box on
plot(state.time_myr,state.d13c_POC_p,'k')
plot(state.time_myr,state.d13c_POC_di,'b')
plot(state.time_myr,state.d13c_POC_s,'r')
plot(state.time_myr,state.d13c_POC_h,'c')
plot(state.time_myr,state.d13c_POC_d,'m')
xlim([-2 15])
xlabel('Time (Ma)'); ylabel(['\delta^{13}C of POC (', char(8240), ')'])
title('g', 'FontSize', 12)

subplot(3,3,8); hold on; box on
plot(state.time_myr,state.O2_conc_p,'k')
plot(state.time_myr,state.O2_conc_di,'b')
plot(state.time_myr,state.O2_conc_s,'r')
plot(state.time_myr,state.O2_conc_h,'c')
plot(state.time_myr,state.O2_conc_d,'m')
xlim([-2 15])
xlabel('Time (Ma)'); ylabel('O_{2} conc (mol/m^{3})')
title('h', 'FontSize', 12)

subplot(3,3,9); hold on; box on
plot(state.time_myr,state.H2S_conc_p,'k')
plot(state.time_myr,state.H2S_conc_di,'b')
plot(state.time_myr,state.H2S_conc_s,'r')
plot(state.time_myr,state.H2S_conc_h,'c')
plot(state.time_myr,state.H2S_conc_d,'m')
xlim([-2 15])
xlabel('Time (Ma)'); ylabel('H_{2}S conc (mol/m^{3})')
title('i', 'FontSize', 12)

% === 新增：湖泊的磷释放量与有机碳埋藏量 ===
figure
subplot(2,1,1); hold on; box on
plot(state.time_myr, state.lake_P_release, 'k')
xlim([0 15])
xlabel('Time (Ma)'); ylabel('Lake P release (mol/yr)')
title('Lake P release')

subplot(2,1,2); hold on; box on
plot(state.time_myr, state.lake_C_burial, 'k')
xlim([0 15])
xlabel('Time (Ma)'); ylabel('Lake OC burial (mol/yr)')
title('Lake organic C burial')

