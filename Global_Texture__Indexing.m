%% === Batch Texture Index, Global PFs, and Local Profiles ===
clear variables; close all; clc;

% --- 1. Define Input and Output Folders ---
dataFolder = 'Q:\EBSD\Master Thesis\20nm\2026_07_29_EBSD_Pt_20nm_600c_30min_ArO2\2026_07_29_EBSD_Pt_20nm_Aro2_angfiles\sample1_2_merged';
outRoot    = 'Q:\EBSD\Master Thesis\20nm\2026_07_29_EBSD_Pt_20nm_600c_30min_ArO2\sample1_2_merged\newTextureIndex\187 diff marker size';
if ~exist(outRoot, 'dir'); mkdir(outRoot); end

fileList = dir(fullfile(dataFolder, '*.ang'));
numFiles = length(fileList);
if numFiles == 0
    error('No .ang files found in the specified directory.');
end

% --- 2. MTEX Setup & Symmetry ---
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection', 'intoPlane');
setMTEXpref('defaultColorMap', parula); 

cs_pt = crystalSymmetry('m-3m', [3.92 3.92 3.92], 'mineral', 'Pt');  
ph1 = 'Pt';
fixed_halfwidth = 10 * degree; % LOCKED half-width for unbiased comparison

% --- Plotting Convention ---
my_convention = plottingConvention;
my_convention.east = xvector;        
my_convention.outOfScreen = -zvector; 

% --- 3. Initialize Tracking Variables ---
master_filenames = strings(numFiles, 1);
master_texture_indices = zeros(numFiles, 1);
mean_ori_list = cell(numFiles, 1);
area_list = cell(numFiles, 1);
master_grain_counts = zeros(numFiles, 1); 

%% --- START LOOP OVER ALL MAPS ---
for i = 1:numFiles
    filename = fullfile(fileList(i).folder, fileList(i).name);
    [~, baseName, ~] = fileparts(filename);
    fprintf('Processing Map %d of %d: %s\n', i, numFiles, baseName);
    master_filenames{i} = baseName;
    
    % Setup output directory for local files
    outDir = fullfile(outRoot, baseName); 
    if ~exist(outDir, 'dir'); mkdir(outDir); end
    
  % === Load Data & Gridify (Matches Original) ===
    ebsd = EBSD.load(filename, 'convertEuler2SpatialReferenceFrame', 'setting 2');
    ebsd = ebsd.gridify;
    
    % --- Manual Crystal Symmetries ---
    cs_pt = crystalSymmetry('m-3m', [3.92 3.92 3.92], 'mineral', 'Pt');  
    ph1 = 'Pt';
    ph2 = 'alumina';
    try
        cs_alumina = ebsd(ph2).CS; 
    catch
        cs_alumina = crystalSymmetry('-3m1', [4.76 4.76 12.99], 'mineral', 'alumina');
    end
    
    if isempty(ebsd(ph1))
        disp('  -> No Platinum found in this map. Skipping.');
        master_texture_indices(i) = NaN;
        continue;
    end
    
    % === Clean Data EXACTLY like original script ===
    ebsd_pt = ebsd(ph1);
    bad_pt_condition = ebsd_pt.ci < 0.2;
    bad_pt_pixels = ebsd_pt(bad_pt_condition);
    ebsd(bad_pt_pixels.id).phase = 0; % Set to notIndexed instead of deleting
    
    [grains, ebsd.grainId] = calcGrains(ebsd, 'angle', 10*degree, 'minPixel', 2);
    
    alumina_grains = grains(ph2);
    if ~isempty(alumina_grains)
        ebsd(alumina_grains(alumina_grains.numPixel < 20)) = [];
    end
    
    % Second calculation & smooth (Matches Original)
    [grains, ebsd.grainId] = calcGrains(ebsd, 'angle', 10*degree);
    grains_smoothed = smooth(grains, 3, 'moveTriplePoints', 'umbrella');
    
    pt_grains = grains_smoothed(ph1); % Extract final Pt grains
    master_grain_counts(i) = length(pt_grains);
    
    % Extract mean orientations and areas for this map
    map_mean_oris = pt_grains.meanOrientation;
    map_areas = pt_grains.area;

    % Calculate Local ODF and Store for global pool
    odf_local = calcDensity(map_mean_oris, 'weights', map_areas, 'halfwidth', fixed_halfwidth);
    mean_ori_list{i} = map_mean_oris(:);
    area_list{i} = map_areas(:);
    
    %% === POLE FIGURES & LINE PROFILES (From User Code) ===
    
    % --- 1. Scatter Plot Z (111) ---
    figure('Name', sprintf('%s_Pt_111_ScatterPlot', baseName), 'Visible', 'off');
    ipfKey = ipfColorKey(cs_pt);
    ipfKey.inversePoleFigureDirection = zvector; % Color based on Z
    colors_Z = ipfKey.orientation2color(map_mean_oris);
    plotPDF(map_mean_oris, colors_Z, Miller(1,1,1,cs_pt), 'how2plot', my_convention, 'all', 'MarkerSize', 3);
    
    % --- 2. Scatter Plot Y (110) ---
    figure('Name', sprintf('%s_Pt_110_ScatterPlot', baseName), 'Visible', 'off');
    ipfKey.inversePoleFigureDirection = yvector; % Recalculate color based on Y
    colors_Y = ipfKey.orientation2color(map_mean_oris);
    plotPDF(map_mean_oris, colors_Y, Miller(1,1,0,cs_pt), 'how2plot', my_convention, 'all', 'MarkerSize', 3);
    
    % --- 3. Texture Index & Contours ---
   current_tex_index = norm(odf_local)^2; % Updated to use newer textureindex() function
    master_texture_indices(i) = current_tex_index;
    fprintf('Calculated Texture Index for %s: %.3f\n', baseName, current_tex_index);
    
    figure('Name', sprintf('%s_Pt_111_110_ContourPlot_DifferentScales', baseName), 'Visible', 'off');
    plotPDF(odf_local, {Miller(1,1,1,cs_pt), Miller(1,1,0,cs_pt)}, 'how2plot', my_convention);
    cb = mtexColorbar;   
    set(cb, 'FontSize', 24);    
    
    figure('Name', sprintf('%s_Pt_111_110_ContourPlot_SameScale', baseName), 'Visible', 'off');
    plotPDF(odf_local, {Miller(1,1,1,cs_pt), Miller(1,1,0,cs_pt)}, 'how2plot', my_convention);
    setColorRange('equal');
    cb = mtexColorbar;   
    set(cb, 'FontSize', 24);    
    
    % % --- 4. Line Profiles ---
    % % Setup points on a great circle for evaluation
    % v = rotation.byAxisAngle(yvector, (0:1:180)*degree) * xvector;
    % ang = (-90:1:90)'; 
    % 
    % % Lineprofile through x (111)
    % pdf1 = calcPDF(odf_local, Miller(1,1,1,cs_pt));
    % values_111 = eval(pdf1, v);
    % T111 = table(ang, values_111(:), 'VariableNames', {'Angle_deg','Intensity'});
    % 
    % figure('Name', sprintf('%s_Pt_111_LineProfile_FWHM', baseName), 'Visible', 'off');
    % plot(pdf1, 'how2plot', my_convention);
    % nextAxis;
    % plot(v, values_111);
    % 
    % % Lineprofile through x (110)
    % pdf2 = calcPDF(odf_local, Miller(1,1,0,cs_pt));
    % values_110 = eval(pdf2, v);
    % T110 = table(ang, values_110(:), 'VariableNames', {'Angle_deg','Intensity'});
    % 
    % figure('Name', sprintf('%s_Pt_110_LineProfile', baseName), 'Visible', 'off');
    % plot(pdf2, 'how2plot', my_convention);
    % nextAxis;
    % plot(v, values_110);
    % 
    % % --- 5. Save Line Profile Text Data ---
    % file111 = fullfile(outDir, sprintf('%s_lineprofile_111.txt', baseName));
    % file110 = fullfile(outDir, sprintf('%s_lineprofile_110.txt', baseName));
    % writetable(T111, file111, 'Delimiter','\t');
    % writetable(T110, file110, 'Delimiter','\t');

    % --- 6. Save all local hidden figures ---
    figs = findall(groot, 'Type', 'figure');
    for k = 1:numel(figs)
        figName = get(figs(k), 'Name');
        if contains(figName, baseName)
            pngPath = fullfile(outDir, [figName '.png']);
            exportgraphics(figs(k), pngPath, 'Resolution', 600, 'ContentType', 'image');
        end
    end
    close all; % clear memory before next map
end

%% no bootstrap 
% %% --- GLOBAL ODF CALCULATION & EXCEL EXPORT ---
% fprintf('\n--- Processing Global Data ---\n');
% all_mean_oris = vertcat(mean_ori_list{:}); 
% all_areas = vertcat(area_list{:});
% 
% if isempty(all_mean_oris)
%     error('No valid platinum grains found across any maps.');
% end
% 
% fprintf('Total pooled grains across all maps: %d\n', length(all_mean_oris));
% 
% % Calculate True Global ODF
% global_odf = calcDensity(all_mean_oris, 'weights', all_areas, 'halfwidth', fixed_halfwidth);
% global_tex_index = norm(global_odf)^2;
% fprintf('TRUE GLOBAL TEXTURE INDEX (J): %.3f\n', global_tex_index);
% 
% % Prepare Export Table
% global_tex_column = NaN(numFiles, 1);
% global_tex_column(1) = global_tex_index; % Place global index in row 1 for clean formatting
% 
% tex_table = table(master_filenames, master_texture_indices, global_tex_column, ...
%     'VariableNames', {'Filename', 'Map_Texture_Index', 'Global_Texture_Index'});
% 
% excelOutName = fullfile(outRoot, 'Texture_Indices_Report.xlsx');
% writetable(tex_table, excelOutName);
% fprintf('Saved Texture Indices to: %s\n', excelOutName);

%% --- BOOTSTRAPPED GLOBAL ODF CALCULATION & EXCEL EXPORT ---
fprintf('\n--- Processing Bootstrapped Global Data ---\n');
all_mean_oris = vertcat(mean_ori_list{:}); 
all_areas = vertcat(area_list{:});

if isempty(all_mean_oris)
    error('No valid platinum grains found across any maps.');
end

total_grains = length(all_mean_oris);
fprintf('Total pooled grains across all maps: %d\n', total_grains);

% 1. Set the Limit and Iterations
N_min = 187;           % Set this to the exact grain count of your vacuum map
num_bootstraps = 100;  % 100 pulls is standard for statistical convergence

% Initialize global variables
global_tex_index = NaN;
std_dev_J = NaN;

if total_grains > N_min
    fprintf('Bootstrapping data down to %d grains (%d iterations)...\n', N_min, num_bootstraps);
    
    bootstrap_J_values = zeros(num_bootstraps, 1);
    
    % 2. The Resampling Loop
    for b = 1:num_bootstraps
        % Randomly select N_min unique indices from the total pool
        random_idx = randperm(total_grains, N_min);
        
        % Extract those specific orientations and their matching areas
        temp_oris = all_mean_oris(random_idx);
        temp_areas = all_areas(random_idx);
        
        % Calculate ODF and J-index for this specific random pull
        temp_odf = calcDensity(temp_oris, 'weights', temp_areas, 'halfwidth', fixed_halfwidth);
        bootstrap_J_values(b) = norm(temp_odf)^2;
    end
    
    % 3. The Final Result
    global_tex_index = mean(bootstrap_J_values);
    std_dev_J = std(bootstrap_J_values); % Calculates the error bar
    
    fprintf('STANDARDIZED (BOOTSTRAPPED) GLOBAL J-INDEX: %.3f ± %.3f\n', global_tex_index, std_dev_J);
else
    % If the sample (like Vacuum) is already at or below limited number of grains, calculate it normally once
    fprintf('Grain count is at or below the threshold. Calculating baseline J-index.\n');
    global_odf = calcDensity(all_mean_oris, 'weights', all_areas, 'halfwidth', fixed_halfwidth);
    global_tex_index = norm(global_odf)^2;
    std_dev_J = 0; % No standard deviation for a single baseline calculation
    fprintf('BASELINE GLOBAL J-INDEX: %.3f\n', global_tex_index);
end

% --- 4. PREPARE EXCEL EXPORT ---
% Create NaN columns so the global values only print in the very first row of the Excel file
global_tex_column = NaN(numFiles, 1);
global_tex_column(1) = global_tex_index; 

global_std_dev_column = NaN(numFiles, 1);
global_std_dev_column(1) = std_dev_J;

global_total_grains_column = NaN(numFiles, 1);
global_total_grains_column(1) = total_grains;

% Build the final table with the new Grain Count and Bootstrapped columns
tex_table = table(master_filenames, master_grain_counts, master_texture_indices, ...
    global_total_grains_column, global_tex_column, global_std_dev_column, ...
    'VariableNames', {'Filename', 'Map_Grain_Count', 'Map_Texture_Index', ...
    'Total_Pooled_Grains', 'Global_Texture_Index', 'Bootstrapped_Std_Dev'});

excelOutName = fullfile(outRoot, 'Texture_Indices_Report.xlsx');
writetable(tex_table, excelOutName);
fprintf('Saved Texture Indices and Grain Statistics to: %s\n', excelOutName);

% %% --- PLOT GLOBAL SCATTER POLE FIGURES (IPF COLORED) ---
% fprintf('\nGenerating Colored Global Scatter Pole Figures as Separate Images...\n');
% 
% % --- 1. Figure for {111} IPF-Z ---
% fig_scatter_Z = figure('Name', 'Global_Scatter_PF_111_Colored', 'Color', 'w');
% set(fig_scatter_Z, 'Position', [100, 650, 600, 500]); % Adjusted width for single plot
% 
% ipfKey_Z = ipfColorKey(cs_pt);
% ipfKey_Z.inversePoleFigureDirection = vector3d.Z;
% colors_Z = ipfKey_Z.orientation2color(all_mean_oris);
% plotPDF(all_mean_oris, colors_Z, Miller(1,1,1,cs_pt), ...
%     'how2plot', my_convention, 'MarkerSize', 3, 'MarkerEdgeAlpha', 0.8, 'points', 187);
% title('Global Scatter Plot (111)', 'FontSize', 14);
% 
% scatter_png_Z = fullfile(outRoot, 'Global_Scatter_PF_111_Colored.png');
% exportgraphics(fig_scatter_Z, scatter_png_Z, 'Resolution', 600, 'ContentType', 'image');
% fprintf('Saved {111} Scatter Pole Figure to: %s\n', scatter_png_Z);
% 
% % --- 2. Figure for {110} IPF-Y ---
% fig_scatter_Y = figure('Name', 'Global_Scatter_PF_110_Colored', 'Color', 'w');
% set(fig_scatter_Y, 'Position', [720, 650, 600, 500]); % Shifted start position to avoid screen overlap
% 
% ipfKey_Y = ipfColorKey(cs_pt);
% ipfKey_Y.inversePoleFigureDirection = vector3d.Y;
% colors_Y = ipfKey_Y.orientation2color(all_mean_oris);
% plotPDF(all_mean_oris, colors_Y, Miller(1,1,0,cs_pt), ...
%     'how2plot', my_convention, 'MarkerSize', 3, 'MarkerEdgeAlpha', 0.8, 'points', 187);
% title('Global Scatter Plot (110)', 'FontSize', 14);
% 
% scatter_png_Y = fullfile(outRoot, 'Global_Scatter_PF_110_Colored.png');
% exportgraphics(fig_scatter_Y, scatter_png_Y, 'Resolution', 600, 'ContentType', 'image');
% fprintf('Saved {110} Scatter Pole Figure to: %s\n', scatter_png_Y);

% %% --- PREPARE DATA FOR PLOTTING ---
% fprintf('\n--- Preparing Data for Plotting ---\n');
% 
% % 1. True Global ODF (All Data)
% true_global_odf = calcDensity(all_mean_oris, 'weights', all_areas, 'halfwidth', fixed_halfwidth);
% 
% % 2. Limited Subset (N_min) for visual comparison
% if total_grains > N_min
%     fprintf('Extracting a random subset of %d grains for limited plots...\n', N_min);
%     random_idx_plot = randperm(total_grains, N_min);
%     limited_oris = all_mean_oris(random_idx_plot);
%     limited_areas = all_areas(random_idx_plot);
% else
%     limited_oris = all_mean_oris;
%     limited_areas = all_areas;
% end
% limited_odf = calcDensity(limited_oris, 'weights', limited_areas, 'halfwidth', fixed_halfwidth);

%% --- PREPARE DATA FOR PLOTTING ---
fprintf('\n--- Preparing Data for Plotting ---\n');

% 1. True Global ODF (All Data)
true_global_odf = calcDensity(all_mean_oris, 'weights', all_areas, 'halfwidth', fixed_halfwidth);

% 2. Limited Subset (N_min) for visual comparison
if total_grains > N_min
    fprintf('Extracting a random subset of %d grains for limited plots...\n', N_min);
    random_idx_plot = randperm(total_grains, N_min);
    limited_oris = all_mean_oris(random_idx_plot);
    limited_areas = all_areas(random_idx_plot);
else
    limited_oris = all_mean_oris;
    limited_areas = all_areas;
end
limited_odf = calcDensity(limited_oris, 'weights', limited_areas, 'halfwidth', fixed_halfwidth);

%% --- 1. TRUE GLOBAL SCATTER PLOTS (ALL DATA) ---
fprintf('\nGenerating True Global Scatter Pole Figures (All Data)...\n');

% {111} IPF-Z
fig_scatter_Z_all = figure('Name', 'True_Global_Scatter_Plot_111', 'Color', 'w');
set(fig_scatter_Z_all, 'Position', [100, 100, 600, 500]);
ipfKey_Z = ipfColorKey(cs_pt);
ipfKey_Z.inversePoleFigureDirection = vector3d.Z;
colors_Z_all = ipfKey_Z.orientation2color(all_mean_oris);
plotPDF(all_mean_oris, colors_Z_all, Miller(1,1,1,cs_pt), 'how2plot', my_convention, 'MarkerSize', 2, 'MarkerEdgeAlpha', 0.8, 'all');
title(sprintf('True Global Scatter Plot (111) (N = %d)', total_grains), 'FontSize', 14);
exportgraphics(fig_scatter_Z_all, fullfile(outRoot, 'True_Global_Scatter_Plot_111.png'), 'Resolution', 600);

% {110} IPF-Y
fig_scatter_Y_all = figure('Name', 'True_Global_Scatter_Plot_110', 'Color', 'w');
set(fig_scatter_Y_all, 'Position', [720, 100, 600, 500]);
ipfKey_Y = ipfColorKey(cs_pt);
ipfKey_Y.inversePoleFigureDirection = vector3d.Y;
colors_Y_all = ipfKey_Y.orientation2color(all_mean_oris);
plotPDF(all_mean_oris, colors_Y_all, Miller(1,1,0,cs_pt), 'how2plot', my_convention, 'MarkerSize', 2, 'MarkerEdgeAlpha', 0.8, 'all');
title(sprintf('True Global Scatter {110} (N = %d)', total_grains), 'FontSize', 14);
exportgraphics(fig_scatter_Y_all, fullfile(outRoot, 'True_Global_Scatter_Plot_110.png'), 'Resolution', 600);

%% --- 2. LIMITED SCATTER PLOTS (N_MIN) ---
fprintf('\nGenerating Limited Scatter Pole Figures (N = %d)...\n', N_min);

% {111} IPF-Z
fig_scatter_Z_lim = figure('Name', 'Limited_Scatter_111', 'Color', 'w');
set(fig_scatter_Z_lim, 'Position', [100, 650, 600, 500]);
colors_Z_lim = ipfKey_Z.orientation2color(limited_oris);
plotPDF(limited_oris, colors_Z_lim, Miller(1,1,1,cs_pt), 'how2plot', my_convention, 'MarkerSize', 2, 'MarkerEdgeAlpha', 0.9, 'all'); 
title(sprintf('Limited Scatter Plot (111) (N = %d)', length(limited_oris)), 'FontSize', 14);
exportgraphics(fig_scatter_Z_lim, fullfile(outRoot, 'Limited_Scatter_111.png'), 'Resolution', 600);

% {110} IPF-Y
fig_scatter_Y_lim = figure('Name', 'Limited_Scatter_plot_110', 'Color', 'w');
set(fig_scatter_Y_lim, 'Position', [720, 650, 600, 500]);
colors_Y_lim = ipfKey_Y.orientation2color(limited_oris);
plotPDF(limited_oris, colors_Y_lim, Miller(1,1,0,cs_pt), 'how2plot', my_convention, 'MarkerSize', 2, 'MarkerEdgeAlpha', 0.9, 'all');
title(sprintf('Limited Scatter Plot (110) (N = %d)', length(limited_oris)), 'FontSize', 14);
exportgraphics(fig_scatter_Y_lim, fullfile(outRoot, 'Limited_Scatter_Plot_110.png'), 'Resolution', 600);

%% --- 3. TRUE GLOBAL CONTOUR POLE FIGURES (ALL DATA) ---
fprintf('\nGenerating True Global Contour Pole Figures...\n');
fig_pf_all = figure('Name', 'True_Global_Contour_PFs', 'Color', 'w');
set(fig_pf_all, 'Position', [100, 100, 1200, 500]);
plotPDF(true_global_odf, {Miller(1,1,1,cs_pt), Miller(1,1,0,cs_pt)}, 'how2plot', my_convention);
setColorRange('equal'); 
cb = mtexColorbar;   
set(cb, 'FontSize', 16); 
drawNow(gcm); 
axes_list = gcm().children; % <-- FIXED SYNTAX
title(axes_list(1), sprintf('True Global Contour Plot (111) (N = %d)', total_grains), 'FontSize', 16);
title(axes_list(2), sprintf('True Global Contour Plot (110) (N = %d)', total_grains), 'FontSize', 16);
exportgraphics(fig_pf_all, fullfile(outRoot, 'True_Global_Contours_111_110.png'), 'Resolution', 600, 'Padding', 'tight');

%% --- 4. LIMITED CONTOUR POLE FIGURES (N_MIN) ---
fprintf('\nGenerating Limited Contour Pole Figures...\n');
fig_pf_lim = figure('Name', 'Limited_Contour_PFs', 'Color', 'w');
set(fig_pf_lim, 'Position', [100, 650, 1200, 500]);
plotPDF(limited_odf, {Miller(1,1,1,cs_pt), Miller(1,1,0,cs_pt)}, 'how2plot', my_convention);
setColorRange('equal'); 
cb = mtexColorbar;   
set(cb, 'FontSize', 16); 
drawNow(gcm); 
axes_list_lim = gcm().children; % <-- FIXED SYNTAX
title(axes_list_lim(1), sprintf('Limited Contour {111} (N = %d)', length(limited_oris)), 'FontSize', 16);
title(axes_list_lim(2), sprintf('Limited Contour {110} (N = %d)', length(limited_oris)), 'FontSize', 16);
exportgraphics(fig_pf_lim, fullfile(outRoot, 'Limited_Contours_111_110.png'), 'Resolution', 600, 'Padding', 'tight');

close all;
disp('--- Script Complete ---');