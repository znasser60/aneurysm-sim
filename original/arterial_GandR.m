% --------------- WORKSPACE CLEANUP -------------

clear
clc
close all

% ---------------- Notational notes -----------------------

% Variablke prefixes
%	c - constant
%	v - variable
%	sv - storage variable

% ----------------------------------------  Parameters --------------------------------------------------------- #

c_diam_tzero_mm = 2.9;                                               % physiological diameter at t=0 in mm;
																	 % middle cerebral artery

c_radius_tzero         = c_diam_tzero_mm * 10^(-3) / (2 * 1.3);      % unloaded radius at=0 in m
c_thickness_tzero    = c_radius_tzero/5;                             % thickness at t=0
c_pressure_sys       = 16000;                                        % systolic blood pressure in Pa

    % --------------------------------- Initial Stretch Considerations ----------------------------------------------- #
    
    c_lambda_z      	= 1.3;	% in vivo axial stretch
    c_lambda_elastin    = 1.3;  % homeostatic circumferential stretch and elastin stretch
    
    c_collagen_ratio_ad_me = 8; % mass ratio of adventitial to medial collagen
        
        % Collagen distribution in media
        
        c_att_min_me = 1.00001; % attachment (homeostatic) stretches
        c_att_mod_me =  1.01;   % minimum, mode, maximum
        c_att_max_me = 1.07;
        
        c_rec_max_me = c_lambda_elastin / c_att_min_me; % respective recruitment stretches
        c_rec_min_me = c_lambda_elastin / c_att_max_me;
        c_rec_mod_me = c_lambda_elastin / c_att_mod_me;
        
        v_a_me = c_rec_min_me; % abbreviations
        v_c_me = c_rec_mod_me;
        v_b_me = c_rec_max_me;
        
        % Collagen distribution in adventitia
        
        c_att_min_ad = 0.8;
        c_att_mod_ad =  0.9;
        c_att_max_ad = 0.99999;
        
        c_rec_max_ad = c_lambda_elastin / c_att_min_ad;
        c_rec_min_ad = c_lambda_elastin / c_att_max_ad;
        c_rec_mod_ad = c_lambda_elastin / c_att_mod_ad;
        
        v_a_ad = c_rec_min_ad;
        v_c_ad = c_rec_mod_ad;
        v_b_ad = c_rec_max_ad;
        
        % Print values for collagen distribution
        
        fprintf('Medial collagen\n');
        fprintf('rec_min = %d,\nrec_mod = %d,\nrec_max = %d,\n\n', c_rec_min_me,c_rec_mod_me,c_rec_max_me)
        fprintf('Adventitial collagen\n');
        fprintf('rec_min = %d,\nrec_mod = %d,\nrec_max = %d,\n\n', c_rec_min_ad,c_rec_mod_ad,c_rec_max_ad)
        
    
	% Muscle stretches: homeostatic stretch + mean and min for active response
    
    c_lambda_muscle  = 1.15; % VSMC attachment stretch
    c_rec_muscle     = c_lambda_elastin / c_lambda_muscle; 
    
    c_musc_mean = 1.1;
    c_musc_min = 0.4;
    c_vasodil_conc = 0.68;   % Concentration of vasodilators to vasoconstrictors at homeostasis
                             % see Humphrey's vasospasm paper (part II)

    c_ge_muscle         = (c_lambda_muscle^2 - 1.0) / 2.0;      % from stretch to green strain
  
    
        % ----------------------------------------- Material Parameters -------------------------------------------------- #
        
        % Assign load bearing proportions to each constituent
        % then solve force-balance equation for k_(.)
        % using systolic blood pressure for P
        % and attachment stretches for lambda's
        %
        % for example for muscle passive
        % c_load_borne_muscle_p * P = 
        % (H/R) * ( 1 / (lambda * lambda_z) ) * sigma_M^pass (lambda_M)
        % k_M is inside sigma_M
    
    c_load_borne_elastin   = 0.50;
    c_load_borne_muscle_p  = 0.20;
    c_load_borne_muscle_a  = 0.20;
    c_load_borne_collagen  = 1 - c_load_borne_elastin - c_load_borne_muscle_p - c_load_borne_muscle_a;
    
    c_common_factor = ( c_pressure_sys * c_radius_tzero * c_lambda_elastin^2 * c_lambda_z ) / c_thickness_tzero;

    c_k_elastin  = ( c_load_borne_elastin * c_common_factor ) / ...
                     ( c_lambda_elastin^2 * ( 1 - (1 / ( c_lambda_z^2 * c_lambda_elastin^4)) ) );

    c_k_collagen = ( c_load_borne_collagen * c_common_factor ) / ...
                     ( ( 2 * c_lambda_elastin / ( (v_b_me-v_a_me)*(v_c_me-v_a_me) ) ) * ( (v_a_me+c_lambda_elastin)*log(c_lambda_elastin/v_a_me) + 2*(v_a_me-c_lambda_elastin) ) ) ;

    c_k_muscle_p = (c_load_borne_muscle_p * c_common_factor ) / ...
                     ( c_lambda_muscle^2 * ( 1 - (1/  ( c_lambda_z^2 * c_lambda_muscle^4) ) )   );
                            
    c_k_muscle_a = (c_load_borne_muscle_a * c_common_factor ) / ...
         ( c_vasodil_conc * ( c_lambda_muscle * ( 1 - ...
        ( (c_musc_mean - c_lambda_muscle)/(c_musc_mean - c_musc_min) )^2 ) )  );
    
    % ------------------------------------------- Printing Values ------------------------------------------------ #
        
    printf('K_E = %d,\nK_C = %d,\nK_M_p = %d,\nk_M_a = %d\n\n', c_k_elastin, c_k_collagen, c_k_muscle_p, c_k_muscle_a);
         
         
   %% ------------- HEALTHY ARTERY -------------------------
   
    % ------------------------STRESS FUNCTIONS LOOP -----
         
         n=235;
         
         % Storage arrays
         
         sv_stretch_var = zeros(1,n);
         sv_stress_var_elastin = zeros(1,n);
         sv_stress_var_collagen = zeros(1,n);
         sv_stress_var_muscle_a = zeros(1,n);
         sv_stress_var_muscle_p = zeros(1,n);
         sv_stress_var_muscle_t = zeros(1,n);
         sv_stress_var_total = zeros(1,n);
         sv_pressure_var = zeros(1,n);
         sv_pressure_var_elastin = zeros(1,n);
         sv_pressure_var_collagen = zeros(1,n);
         sv_pressure_var_muscle = zeros(1,n);
         sv_pressure_var_muscle_a = zeros(1,n);
         sv_pressure_var_muscle_p = zeros(1,n);
         sv_pressure_var_collagen_me = zeros(1,n);
         sv_pressure_var_collagen_ad = zeros(1,n);
         
		 % Loop
         for i=1:n
             
             % Initialize stretch
             
             sv_stretch_var(i) = 0.55 + (i-1)*0.01;
             
             % Define stress functions
             
            v_lambda_collagen   = @(x) x / c_rec_collagen;
            v_lambda_muscle     = @(x) x / c_rec_muscle;
            v_m = @(x) (x / c_rec_muscle);

            v_ge_collagen       = @(x) (v_lambda_collagen(x)^2 - 1.0) / 2.0;
            v_ge_muscle         = @(x) (v_lambda_muscle(x)^2 - 1.0) / 2.0;
            v_ge                = @(x) ( x^2 - 1 ) / 2;

            v_sigma_elastin     = @(x) x^2 * c_k_elastin * ( 1 - (1 /  ( c_lambda_z^2 * x^4 ) ) );
            v_sigma_muscle_p    = @(x) v_lambda_muscle(x)^2 * c_k_muscle_p * ...
                ( 1 - (1 / (c_lambda_z^2 * v_lambda_muscle(x)^4)) );
            v_sigma_muscle_a    = @(x)   c_vasodil_conc * c_k_muscle_a * ...
                v_m(x)  * ( 1 - ( (c_musc_mean - v_m(x)) / (c_musc_mean - c_musc_min) )^2 );
            v_sigma_muscle_t    = @(x) v_sigma_muscle_a(x) + v_sigma_muscle_p(x) ;
            
            % Collagen Cauchy stresses in media
            
            v_gamma_me                 = c_k_collagen / ( (v_b_me - v_a_me) * (v_c_me - v_a_me)  );
            v_delta_me                     = c_k_collagen / ( (v_b_me - v_a_me) * (v_b_me - v_c_me) );
            
            v_sigma_collagen_me_0      = @(x) x * 0;
            v_sigma_collagen_me_ac     = @(x) x * v_gamma_me * 2 * ( (x + v_a_me) * log(x/v_a_me) + 2*(v_a_me - x) ) ;
            v_sigma_collagen_me_cb     = @(x) x * v_gamma_me * 2 * ( (x + v_a_me)*log(v_c_me/v_a_me) + v_a_me - v_c_me + ( (v_a_me - v_c_me) / v_c_me ) * x) ...
                - x * v_delta_me * 2 * ( (x + v_b_me)*log(x/v_c_me) + v_b_me + v_c_me - ( (v_b_me + v_c_me) / v_c_me ) * x );
            v_sigma_collagen_me_b      = @(x) x * v_gamma_me * 2 * ( (x + v_a_me)*log(v_c_me/v_a_me) + v_a_me - v_c_me + ( (v_a_me - v_c_me) / v_c_me ) * x) ...
                - x * v_delta_me * 2 * ( (x + v_b_me)*log(v_b_me/v_c_me) - v_b_me + v_c_me - ( (v_b_me - v_c_me) / v_c_me ) * x );
            
            v_sigma_collagen_me        =@(x) v_sigma_collagen_me_0(x).*(x<v_a_me)...
                + v_sigma_collagen_me_ac(x).*( x>=v_a_me & x<v_c_me)...
                + v_sigma_collagen_me_cb(x).*(x>=v_c_me & x<=v_b_me)...
                + v_sigma_collagen_me_b(x).*(x>v_b_me);

            % Collagen Cauchy stresses in adventitia
            
            v_gamma_ad                 =  c_collagen_ratio_ad_me * c_k_collagen / ( (v_b_ad - v_a_ad) * (v_c_ad - v_a_ad)  );
            v_delta_ad                 =  c_collagen_ratio_ad_me * c_k_collagen / ( (v_b_ad - v_a_ad) * (v_b_ad - v_c_ad) );
            
            v_sigma_collagen_ad_0      = @(x) x * 0;
            v_sigma_collagen_ad_ac     = @(x) x * v_gamma_ad * 2 * ( (x + v_a_ad) * log(x/v_a_ad) + 2*(v_a_ad - x) ) ;
            v_sigma_collagen_ad_cb     = @(x) x * v_gamma_ad * 2 * ( (x + v_a_ad)*log(v_c_ad/v_a_ad) + v_a_ad - v_c_ad + ( (v_a_ad - v_c_ad) / v_c_ad ) * x) ...
                - x * v_delta_ad * 2 * ( (x + v_b_ad)*log(x/v_c_ad) + v_b_ad + v_c_ad - ( (v_b_ad + v_c_ad) / v_c_ad ) * x );
            v_sigma_collagen_ad_b      = @(x) x * v_gamma_ad * 2 * ( (x + v_a_ad)*log(v_c_ad/v_a_ad) + v_a_ad - v_c_ad + ( (v_a_ad - v_c_ad) / v_c_ad ) * x) ...
                - x * v_delta_ad * 2 * ( (x + v_b_ad)*log(v_b_ad/v_c_ad) - v_b_ad + v_c_ad - ( (v_b_ad - v_c_ad) / v_c_ad ) * x );
            
            v_sigma_collagen_ad        = @(x) v_sigma_collagen_ad_0(x).*(x<v_a_ad)...
                + v_sigma_collagen_ad_ac(x).*( x>=v_a_ad & x<v_c_ad)...
                + v_sigma_collagen_ad_cb(x).*(x>=v_c_ad & x<=v_b_ad)...
                + v_sigma_collagen_ad_b(x).*(x>v_b_ad);

            v_sigma_collagen            = @(x) v_sigma_collagen_me(x) + v_sigma_collagen_ad(x);
                
            
            v_pres_prefactor    = @(x) c_thickness_tzero / ( c_radius_tzero * c_lambda_z * x^2 );
       
            v_pressure_ECM = @(x) v_pres_prefactor(x) * ( v_sigma_elastin(x) + v_sigma_collagen(x) + v_sigma_muscle_t(x) );
            v_pressure_EC = @(x) v_pres_prefactor(x) * ( v_sigma_elastin(x) + v_sigma_collagen(x) );
            v_pressure_EM = @(x) v_pres_prefactor(x) * ( v_sigma_elastin(x) + v_sigma_muscle_t(x) );
            v_pressure_E = @(x) v_pres_prefactor(x) * ( v_sigma_elastin(x) );
                      
            v_pressure_elastin      = @(x) v_pres_prefactor(x) * v_sigma_elastin(x);
            v_pressure_collagen     = @(x) v_pres_prefactor(x) * v_sigma_collagen(x);
            v_pressure_muscle       = @(x) v_pres_prefactor(x) * v_sigma_muscle_t(x);
            v_pressure_muscle_a     = @(x) v_pres_prefactor(x) * v_sigma_muscle_a(x);
            v_pressure_muscle_p     = @(x) v_pres_prefactor(x) * v_sigma_muscle_p(x);
            v_pressure_collagen_me  = @(x) v_pres_prefactor(x) * v_sigma_collagen_me(x);
            v_pressure_collagen_ad  = @(x) v_pres_prefactor(x) * v_sigma_collagen_ad(x);
            
             
            % Store results
            
            sv_stress_var_elastin(i) = v_sigma_elastin(sv_stretch_var(i)); 
            sv_stress_var_collagen(i) = v_sigma_collagen(sv_stretch_var(i));  
            sv_stress_var_muscle_a(i) = v_sigma_muscle_a(sv_stretch_var(i));  
            sv_stress_var_muscle_p(i) = v_sigma_muscle_p(sv_stretch_var(i)); 
            sv_stress_var_muscle_t(i) = v_sigma_muscle_t(sv_stretch_var(i)); 
            sv_stress_var_total(i) = sv_stress_var_elastin(i) + sv_stress_var_collagen(i) + sv_stress_var_muscle_t(i);
                
            sv_pressure_var_elastin(i) = v_pressure_elastin(sv_stretch_var(i));    
            
           
            sv_pressure_var_elastin(i)      = max( v_pressure_elastin(sv_stretch_var(i)) , 0 );   
            sv_pressure_var_collagen(i)     = v_pressure_collagen(sv_stretch_var(i)); %max( v_pressure_collagen(sv_stretch_var(i)) , 0 );    
            sv_pressure_var_collagen_me(i)  =v_pressure_collagen_me(sv_stretch_var(i));
            sv_pressure_var_collagen_ad(i)  =v_pressure_collagen_ad(sv_stretch_var(i));
            sv_pressure_var_muscle_p(i)     = max( v_pressure_muscle_p(sv_stretch_var(i)) , 0 ); 
            sv_pressure_var_muscle_a(i)     = max( v_pressure_muscle_a(sv_stretch_var(i)) , 0 ); 
            sv_pressure_var_muscle(i)       = sv_pressure_var_muscle_a(i) + sv_pressure_var_muscle_p(i); ...
            sv_pressure_var(i)    = sv_pressure_var_elastin(i) + sv_pressure_var_collagen_me(i) +...
                 sv_pressure_var_collagen_ad(i) + sv_pressure_var_muscle(i);
                     
            
         end
         
         
             %% -------------- PLOTS PRESSURE VS STRETCH
    
            n_zoom = 120;
    
    figure

    hold on
    plot(sv_stretch_var(1:n_zoom),sv_pressure_var(1:n_zoom)./(10^3),'LineWidth',2)
    plot(sv_stretch_var(33:n_zoom),sv_pressure_var_elastin(33:n_zoom)./(10^3), '--','LineWidth',2)
%   plot(sv_stretch_var(65:n_zoom),sv_pressure_var_collagen(65:n_zoom)./(10^3), '--','LineWidth',2)
	plot(sv_stretch_var(65:n_zoom),sv_pressure_var_collagen_me(65:n_zoom)./(10^3), '--','LineWidth',2)
	plot(sv_stretch_var(65:n_zoom),sv_pressure_var_collagen_ad(65:n_zoom)./(10^3), '--','LineWidth',2)
    plot(sv_stretch_var(45:n_zoom),sv_pressure_var_muscle_p(45:n_zoom)./(10^3),'--', 'LineWidth',2)
    plot(sv_stretch_var(1:n_zoom),sv_pressure_var_muscle_a(1:n_zoom)./(10^3), '--', 'LineWidth',2)
    plot(sv_stretch_var(1:n_zoom),sv_pressure_var(1:n_zoom)./(10^3),'LineWidth',2)
    hold off
    legend('Total','Elastin','Collagen media', 'Collagen adve','Muscle Passive', 'Muscle Active','Location','northwest')
    
    xlabel('Stretch')
    ylabel('Pressure (kPa)')
    set(gca, 'fontsize', 16)
    

    
   
        %% ---------  PRESSURE VS DIAMETER - CONSTITUENTS
    
    
            n_zoom = 125;
        
        sv_diam_var = 2 * c_radius_tzero * 10^3 * sv_stretch_var;
        
        sv_pressure_var_smc = sv_pressure_var_muscle_a + sv_pressure_var_muscle_p;
        sv_pressure_var_coll = sv_pressure_var_collagen_me + sv_pressure_var_collagen_ad;
        sv_pressure_var2 = sv_pressure_var_elastin + sv_pressure_var_coll + sv_pressure_var_smc;
        
    figure

    hold on
    plot(sv_diam_var(1:n_zoom),sv_pressure_var2(1:n_zoom)./(10^3),'k','LineWidth',6)
    plot(sv_diam_var(33:n_zoom),sv_pressure_var_elastin(33:n_zoom)./(10^3), 'k--','LineWidth',3)
    plot(sv_diam_var(67:n_zoom),sv_pressure_var_coll(67:n_zoom)./(10^3), 'k-.','LineWidth',3)
    plot(sv_diam_var(45:n_zoom),sv_pressure_var_muscle_p(45:n_zoom)./(10^3),'k+', 'LineWidth',2)
    plot(sv_diam_var(1:n_zoom),sv_pressure_var_muscle_a(1:n_zoom)./(10^3), 'k.', 'LineWidth',2)
    plot(sv_diam_var(1:n_zoom),sv_pressure_var2(1:n_zoom)./(10^3),'k','LineWidth',6)
    line([2.9 2.9],[0 16],'color','red','LineStyle','--','LineWidth',2)
    line([1 2.9],[16 16],'color','red','LineStyle','--','LineWidth',2)
    plot(2.9, 16, 'ro', 'LineWidth',3)
    hold off
    legend('Total','E','C','VSMCp', 'VSMCa', 'Location','northwest')
    
    xlabel('Diameter (mm)')
    ylabel('Pressure (kPa)')
    set(gca, 'fontsize', 24)
    ylim([0 60])