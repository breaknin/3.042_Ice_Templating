% Freeze Casting Simulation 
clear; clc; close all;

%% VARIABLES
total_minutes = 10;       
LN2_state = 'boiling';     

% Initial temperatures
T_trigger       = 0;    % Slurry added when rod top hits this temp
T_cu_start      = 20;     % Rod starts at room temp
T_slurry_start  = 20;    
T_plastic_start = 20;    
T_LN2           = -196;  
T_air           = 20;
T_freeze        = 0;     

% Heat Transfer Coefficients
h_air = 15; 
h_LN2 = 150; % High convection for boiling LN2

% Material properties
k_cu = 400; rho_cu = 8960; cp_cu = 385; 
k_pl = 0.37; rho_pl = 1420; cp_pl = 1500; 
k_air = 0.026; rho_air = 1.2; cp_air = 1005;

% Slurry properties 
k_sl_liq = 0.61; % approximating slurrying conductivity = water conductivity   
k_sl_solid = (2.2 + 1.3)/2;  % Approximate frozen slurry as (k_ice + k_silica / 2)
cp_sl_liq = 3100; 
cp_sl_solid = 1800; 
rho_sl = 1300; 
L_f = 334000 * 0.66; 
dT_m = 1.5; 

% Geometry
in2m = 0.0254; 
L_cu = 12 * in2m;           
L_bath = 4 * in2m;          % 4 inch bath immersion
R_cu = (0.75 * in2m)/2;     
H_sl = 0.25 * in2m;         

R_sl = R_cu;                  % Slurry radius now perfectly matches copper rod
mold_thickness = 0.2 * in2m; % Define a thin wall (e.g., 0.05 inches)
R_pl = R_sl + mold_thickness; % Outer mold radius is simply slurry + thickness

%% MESH & INITIALIZATION
dr = 0.0015; dz = 0.002; 
r = 0:dr:R_pl; z = 0:dz:(L_cu + H_sl);
[R, Z] = meshgrid(r, z); [nz, nr] = size(R);

% MatID: 0=Air/Environment, 1=Copper, 2=Plastic, 3=Slurry
MatID = zeros(nz, nr);
MatID(Z <= L_cu & R <= R_cu) = 1; 

% Temperature Initialization
T = ones(nz, nr) * T_air;
T(Z <= L_bath & R > R_cu) = T_LN2;  % Environment bath starts at -196C
T(MatID == 1) = T_cu_start;         % Entire Copper Rod starts at 20C

% Initial Property Arrays (Everything outside Cu is air initially)
K = ones(nz,nr) * k_air; RHO = ones(nz,nr) * rho_air; CP = ones(nz,nr) * cp_air;
K(MatID==1)=k_cu; RHO(MatID==1)=rho_cu; CP(MatID==1)=cp_cu;

dt = 0.1 * (min(dr,dz)^2) / (k_cu / (rho_cu * cp_cu)); 
n_steps = floor((total_minutes * 60) / dt);

slurry_deployed = false;
cu_top_idx = find(z >= L_cu, 1);
cu_edge_idx = find(r >= R_cu, 1);

%% VIDEO SETUP
video_name = 'FreezeCasting_Bath_Simulation-60C.mp4';
v = VideoWriter(video_name, 'MPEG-4');
v.FrameRate = 15;
open(v);

%% SIMULATION LOOP
fig = figure('Color', 'w', 'Position', [50 100 1600 500]);
r_vis = [-flip(r(2:end)), r]; [R_vis, Z_vis] = meshgrid(r_vis, z);

for s = 1:n_steps
    % CHECK FOR SLURRY ADDITION (Triggered by rod surface temp)
    if ~slurry_deployed && T(cu_top_idx, 1) <= T_trigger
        slurry_deployed = true;
        % Deploy slurry (3) and plastic mold (2)
        MatID(Z > L_cu & R <= R_sl) = 3;    
        MatID(Z > L_cu & R > R_sl & R <= R_pl) = 2; 
        
        T(MatID == 2) = T_plastic_start;
        T(MatID == 3) = T_slurry_start;
        
        K(MatID==2)=k_pl; RHO(MatID==2)=rho_pl; CP(MatID==2)=cp_pl;
        K(MatID==3)=k_sl_liq; RHO(MatID==3)=rho_sl; CP(MatID==3)=cp_sl_liq;
        fprintf('Slurry added at %.1f minutes.\n', (s*dt)/60);
    end

    % PHASE CHANGE (Slurry properties change: the ice is more conductive than the liquid)
    CP_eff = CP;
    if slurry_deployed
        sl_mask = (MatID == 3);
        f_liq = 0.5 * (1 + tanh((T(sl_mask) - T_freeze) / dT_m));
        K(sl_mask) = k_sl_solid + (k_sl_liq - k_sl_solid) * f_liq;
        base_cp = cp_sl_solid + (cp_sl_liq - cp_sl_solid) * f_liq;
        CP_eff(sl_mask) = base_cp + (L_f/(dT_m*sqrt(pi))) * exp(-((T(sl_mask)-T_freeze)/dT_m).^2);
    end

    % Internal temperature update
    T_old = T;
    rho_cp = RHO(2:end-1, 2:end-1) .* CP_eff(2:end-1, 2:end-1);
    r_mid = R(2:end-1, 2:end-1);
    
    K_up    = 2 * K(2:end-1, 2:end-1) .* K(3:end, 2:end-1) ./ (K(2:end-1, 2:end-1) + K(3:end, 2:end-1));
    K_down  = 2 * K(2:end-1, 2:end-1) .* K(1:end-2, 2:end-1) ./ (K(2:end-1, 2:end-1) + K(1:end-2, 2:end-1));
    K_right = 2 * K(2:end-1, 2:end-1) .* K(2:end-1, 3:end) ./ (K(2:end-1, 2:end-1) + K(2:end-1, 3:end));
    K_left  = 2 * K(2:end-1, 2:end-1) .* K(2:end-1, 1:end-2) ./ (K(2:end-1, 2:end-1) + K(2:end-1, 1:end-2));
    
    flux_z = (K_up .* (T_old(3:end, 2:end-1) - T_old(2:end-1, 2:end-1)) - ...
              K_down .* (T_old(2:end-1, 2:end-1) - T_old(1:end-2, 2:end-1))) / dz^2;
    r_right = r_mid + dr/2; r_left  = r_mid - dr/2;
    flux_r = (r_right .* K_right .* (T_old(2:end-1, 3:end) - T_old(2:end-1, 2:end-1)) - ...
              r_left .* K_left .* (T_old(2:end-1, 2:end-1) - T_old(2:end-1, 1:end-2))) ./ (r_mid * dr^2);
    
    T_calc = T_old(2:end-1, 2:end-1) + (dt ./ rho_cp) .* (flux_z + flux_r);
    T(2:end-1, 2:end-1) = T_calc;

    % Boundary conditions
    
    % Axis of symmetry
    T(:, 1) = T(:, 2); 
    
    % Enforce LN2 enviroment (prevents the bath from warming up)
    bath_env = (Z <= L_bath & R > R_cu);
    T(bath_env) = T_LN2;
    
    % Bottom face of rod in LN2
    T(1, 1:cu_edge_idx) = (K(1, 1:cu_edge_idx).*T(2, 1:cu_edge_idx) + h_LN2*dz*T_LN2) ./ (K(1, 1:cu_edge_idx) + h_LN2*dz);
    
    % Submerged copper side (Convection from LN2)
    sub_mask = (Z(:, cu_edge_idx) <= L_bath);
    T(sub_mask, cu_edge_idx) = (K(sub_mask, cu_edge_idx).*T(sub_mask, cu_edge_idx-1) + h_LN2*dr*T_LN2) ./ ...
                               (K(sub_mask, cu_edge_idx) + h_LN2*dr);
    
    % Exposed copper size (Convection from Air)
    exp_mask = (Z(:, cu_edge_idx) > L_bath & Z(:, cu_edge_idx) <= L_cu);
    T(exp_mask, cu_edge_idx) = (K(exp_mask, cu_edge_idx).*T(exp_mask, cu_edge_idx-1) + h_air*dr*T_air) ./ ...
                               (K(exp_mask, cu_edge_idx) + h_air*dr);
                           
    % Top edge of simulation (Air Convection)
    T(end, :) = (K(end,:).*T(end-1,:) + h_air*dz*T_air) ./ (K(end,:) + h_air*dz);
    
    % Far outer radius (Environment limit)
    T(:, end) = (K(:, end).*T(:, end-1) + h_air*dr*T_air) ./ (K(:, end) + h_air*dr);
    
    % Plastic mold outer radius (Convection from air)
    if slurry_deployed
        pl_edge_idx = find(r >= R_pl, 1);
        pl_mask = (Z(:, pl_edge_idx) > L_cu);
        T(pl_mask, pl_edge_idx) = (K(pl_mask, pl_edge_idx).*T(pl_mask, pl_edge_idx-1) + h_air*dr*T_air) ./ ...
                                  (K(pl_mask, pl_edge_idx) + h_air*dr);
    end

    % VISUALIZATION
    if mod(s, floor(3/dt)) == 0 || s == 1
        T_vis = [flip(T(:, 2:end), 2), T];
        MatID_mir = [flip(MatID(:, 2:end), 2), MatID];
        curr_time = s * dt;
        
        subplot(1, 3, 1);
        pcolor(R_vis*100, Z_vis*100, T_vis); shading interp; colormap(gca, jet); 
        colorbar; caxis([-196 25]); 
        title(sprintf('Full Setup\nTime: %.1f min', curr_time/60));
        ylabel('Height (cm)'); xlabel('Width (cm)'); hold on;
        % Draw bath surface line
        line([-R_pl R_pl]*100, [L_bath L_bath]*100, 'Color', 'w', 'LineStyle', '--');
        % Draw copper outline
        line([-R_cu R_cu R_cu -R_cu -R_cu]*100, [0 0 L_cu L_cu 0]*100, 'Color', 'k', 'LineWidth', 1.2);
        if slurry_deployed
            line([-R_sl R_sl R_sl -R_sl -R_sl]*100, [L_cu L_cu L_cu+H_sl L_cu+H_sl L_cu]*100, 'Color', 'w', 'LineWidth', 1.5);
            line([-R_pl -R_pl R_pl R_pl]*100, [L_cu+H_sl L_cu L_cu L_cu+H_sl]*100, 'Color', 'm', 'LineWidth', 1.2);
        end
        hold off;
        
        subplot(1, 3, 2);
        T_slurry_vis = T_vis; 
        if slurry_deployed
            T_slurry_vis(MatID_mir ~= 3) = NaN;
            pcolor(R_vis*100, Z_vis*100, T_slurry_vis); shading interp; colormap(gca, jet);
            hold on;
            contour(R_vis*100, Z_vis*100, T_slurry_vis, [T_freeze T_freeze], 'Color', 'w', 'LineWidth', 2, 'LineStyle', '-.');
            hold off;
            title('Slurry Zoom');
        else
            cla; 
            text(0, (L_cu+H_sl/2)*100, sprintf('Surface Temp: %.1f C\n.', T(cu_top_idx, 1)), ...
                'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');
            title('Upper Region (Air)');
        end
        colorbar; caxis([-25 25]);
        ylim([(L_cu*100)-1, (L_cu+H_sl)*100+1]); xlim([-R_sl*180, R_sl*180]);
        
        subplot(1, 3, 3);
        [~, dTdZ] = gradient(T, dr, dz);
        dTdZ_vis = [flip(dTdZ(:, 2:end), 2), dTdZ];
        if slurry_deployed
            dTdZ_vis(MatID_mir ~= 3) = NaN;
            pcolor(R_vis*100, Z_vis*100, dTdZ_vis); shading interp; colormap(gca, parula);
            title('Temperature Gradient');
        else
            cla; title('Temperature Gradient');
        end
        colorbar, caxis([0 15000]); 
        ylim([(L_cu*100)-1, (L_cu+H_sl)*100+1]); xlim([-R_sl*180, R_sl*180]);
        
        drawnow;
        frame = getframe(fig);
        writeVideo(v, frame);
    end
end
close(v);
fprintf('Simulation finished. File: %s\n', video_name);