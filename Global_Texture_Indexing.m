%% ========================================================================
%  Batch Texture Index (J-Index), Bootstrapped Global ODF, & Pole Figures
%  ========================================================================
%  Processes EBSD (.ang) maps of Pt thin films on sapphire substrates:
%    1. Cleans low-CI Pt pixels and small alumina inclusions, then reconstructs
%       and smooths Pt grains.
%    2. Computes area-weighted local ODFs, per-map J-indices, and local PFs
%       (plus optional great-circle line profiles).
%    3. Pools orientations globally and performs bootstrap resampling (N_min)
%       to compute an unbiased global J-index ± standard deviation
%       (with an optional non-bootstrapped baseline block).
%    4. Exports summary statistics (.xlsx) and global/subset pole figures.
%
%  Dependencies: MTEX Toolbox (https://mtex-toolbox.github.io/)
%  ========================================================================

clear variables; close all; clc;

%% --- 1. Configuration & Parameters ---
dataFolder = '';
outRoot    = '';

% Phase names & cleaning thresholds
ph1                 = 'Pt';
ph2                 = 'alumina';
ciThreshold         = 0.2;         % Minimum Confidence Index for Pt pixels
grainMisorientation = 10 * degree; % Grain boundary misorientation threshold
minGrainPixels      = 2;           % Minimum pixels per grain (initial pass)
minAluminaPixels    = 20;          % Minimum pixels to retain alumina grains
smoothIterations    = 3;           % Umbrella smoothing passes

% ODF & Bootstrapping settings
fixedHalfwidth      = 10 * degree; % Locked kernel half-width for unbiased comparison
N_min               = 187;         % Target grain count (matched to vacuum baseline)
numBootstraps       = 100;         % Number of random resampling iterations
exportDPI           = 600;         % Output figure resolution

if ~exist(outRoot, 'dir'); mkdir(outRoot); end

fileList = dir(fullfile(dataFolder, '*.ang'));
numFiles = numel(fileList);
if numFiles == 0
    error('No .ang files found in: %s', dataFolder);
end

%% --- 2. MTEX Setup, Symmetries, & Plotting Conventions ---
setMTEXpref('xAxisDirection', 'east');
setMTEXpref('zAxisDirection', 'intoPlane');
setMTEXpref('defaultColorMap', parula);

cs_pt = crystalSymmetry('m-3m', [3.92 3.92 3.92], 'mineral', 'Pt');

my_convention = plottingConvention;
my_convention.east = xvector;
my_convention.outOfScreen = -zvector;

ipfKey_Z = ipfColorKey(cs_pt);
ipfKey_Z.inversePoleFigureDirection = vector3d.Z;

ipfKey_Y = ipfColorKey(cs_pt);
ipfKey_Y.inversePoleFigureDirection = vector3d.Y;

%% --- 3. Initialize Tracking Arrays ---
master_filenames       = strings(numFiles, 1);
master_texture_indices = nan(numFiles, 1);
master_grain_counts    = zeros(numFiles, 1);
mean_ori_list          = cell(numFiles, 1);
area_list              = cell(numFiles, 1);

%% --- 4. Per-Map Grain Reconstruction, Local ODF, & Pole Figures ---
for i = 1:numFiles
    filename = fullfile(fileList(i).folder, fileList(i).name);
    [~, baseName, ~] = fileparts(filename);
    fprintf('Processing Map %d of %d: %s\n', i, numFiles, baseName);
    master_filenames(i) = string(baseName);

    outDir = fullfile(outRoot, baseName);
    if ~exist(outDir, 'dir'); mkdir(outDir); end

    % Load and gridify EBSD data
    ebsd = EBSD.load(filename, 'convertEuler2SpatialReferenceFrame', 'setting 2');
    ebsd = ebsd.gridify;

    if isempty(ebsd(ph1))
        fprintf('  -> No Platinum found in %s. Skipping.\n', baseName);
        continue;
    end

    % Filter low-CI Pt pixels by reassigning to notIndexed (phase = 0)
    ebsd_pt = ebsd(ph1);
    bad_pt_pixels = ebsd_pt(ebsd_pt.ci < ciThreshold);
    ebsd(bad_pt_pixels.id).phase = 0;

    % Initial grain reconstruction & removal of small alumina fragments
    [grains, ebsd.grainId] = calcGrains(ebsd, 'angle', grainMisorientation, 'minPixel', minGrainPixels);
    alumina_grains = grains(ph2);
    if ~isempty(alumina_grains)
        ebsd(alumina_grains(alumina_grains.numPixel < minAluminaPixels)) = [];
    end

    % Final grain reconstruction & smoothing
    [grains, ebsd.grainId] = calcGrains(ebsd, 'angle', grainMisorientation);
    grains_smoothed = smooth(grains, smoothIterations, 'moveTriplePoints', 'umbrella');

    pt_grains = grains_smoothed(ph1);
    master_grain_counts(i) = numel(pt_grains);

    map_mean_oris = pt_grains.meanOrientation;
    map_areas     = pt_grains.area;

    % Compute area-weighted local ODF and J-index
    odf_local = calcDensity(map_mean_oris, 'weights', map_areas, 'halfwidth', fixedHalfwidth);
    current_tex_index = norm(odf_local)^2;
    master_texture_indices(i) = current_tex_index;
    fprintf('  -> Grains: %d | Texture Index (J): %.3f\n', master_grain_counts(i), current_tex_index);

    mean_ori_list{i} = map_mean_oris(:);
    area_list{i}     = map_areas(:);

    % --- 4a. Local Scatter Pole Figures ---
    exportLocalScatterPF(map_mean_oris, ipfKey_Z, Miller(1,1,1,cs_pt), my_convention, ...
        fullfile(outDir, [baseName '_Pt_111_ScatterPlot.png']), exportDPI);

    exportLocalScatterPF(map_mean_oris, ipfKey_Y, Miller(1,1,0,cs_pt), my_convention, ...
        fullfile(outDir, [baseName '_Pt_110_ScatterPlot.png']), exportDPI);

    % --- 4b. Local Contour Pole Figures (Independent vs. Equal Color Scales) ---
    exportLocalContourPF(odf_local, cs_pt, my_convention, false, ...
        fullfile(outDir, [baseName '_Pt_111_110_ContourPlot_DifferentScales.png']), exportDPI);

    exportLocalContourPF(odf_local, cs_pt, my_convention, true, ...
        fullfile(outDir, [baseName '_Pt_111_110_ContourPlot_SameScale.png']), exportDPI);

    % --- 4c. [Optional] Great-Circle Line Profiles & Text Data Export ---
    % % Setup points on a great circle for evaluation
    % v = rotation.byAxisAngle(yvector, (0:1:180)*degree) * xvector;
    % ang = (-90:1:90)';
    %
    % % Line profile through x for {111}
    % pdf1 = calcPDF(odf_local, Miller(1,1,1,cs_pt));
    % values_111 = eval(pdf1, v);
    % T111 = table(ang, values_111(:), 'VariableNames', {'Angle_deg', 'Intensity'});
    %
    % fig_lp111 = figure('Name', sprintf('%s_Pt_111_LineProfile_FWHM', baseName), 'Visible', 'off');
    % plot(pdf1, 'how2plot', my_convention);
    % nextAxis;
    % plot(v, values_111);
    % exportgraphics(fig_lp111, fullfile(outDir, [baseName '_Pt_111_LineProfile_FWHM.png']), ...
    %     'Resolution', exportDPI, 'ContentType', 'image');
    % close(fig_lp111);
    %
    % % Line profile through x for {110}
    % pdf2 = calcPDF(odf_local, Miller(1,1,0,cs_pt));
    % values_110 = eval(pdf2, v);
    % T110 = table(ang, values_110(:), 'VariableNames', {'Angle_deg', 'Intensity'});
    %
    % fig_lp110 = figure('Name', sprintf('%s_Pt_110_LineProfile', baseName), 'Visible', 'off');
    % plot(pdf2, 'how2plot', my_convention);
    % nextAxis;
    % plot(v, values_110);
    % exportgraphics(fig_lp110, fullfile(outDir, [baseName '_Pt_110_LineProfile.png']), ...
    %     'Resolution', exportDPI, 'ContentType', 'image');
    % close(fig_lp110);
    %
    % % Save line profile tables
    % writetable(T111, fullfile(outDir, sprintf('%s_lineprofile_111.txt', baseName)), 'Delimiter', '\t');
    % writetable(T110, fullfile(outDir, sprintf('%s_lineprofile_110.txt', baseName)), 'Delimiter', '\t');
end

%% --- 5. [Optional] Unbootstrapped Global ODF & Excel Export ---
% fprintf('\n--- Processing Unbootstrapped Global Data ---\n');
% all_mean_oris = vertcat(mean_ori_list{:});
% all_areas     = vertcat(area_list{:});
%
% if isempty(all_mean_oris)
%     error('No valid platinum grains found across any maps.');
% end
%
% fprintf('Total pooled grains across all maps: %d\n', numel(all_mean_oris));
%
% global_odf       = calcDensity(all_mean_oris, 'weights', all_areas, 'halfwidth', fixedHalfwidth);
% global_tex_index = norm(global_odf)^2;
% fprintf('TRUE GLOBAL TEXTURE INDEX (J): %.3f\n', global_tex_index);
%
% global_tex_column    = nan(numFiles, 1);
% global_tex_column(1) = global_tex_index;
%
% tex_table = table(master_filenames, master_texture_indices, global_tex_column, ...
%     'VariableNames', {'Filename', 'Map_Texture_Index', 'Global_Texture_Index'});
%
% excelOutName = fullfile(outRoot, 'Texture_Indices_Report.xlsx');
% writetable(tex_table, excelOutName);
% fprintf('Saved Texture Indices to: %s\n', excelOutName);

%% --- 6. Bootstrapped Global ODF & Texture Index ---
fprintf('\n--- Processing Bootstrapped Global Data ---\n');
all_mean_oris = vertcat(mean_ori_list{:});
all_areas     = vertcat(area_list{:});

if isempty(all_mean_oris)
    error('No valid platinum grains found across any maps.');
end

total_grains = numel(all_mean_oris);
fprintf('Total pooled grains across all maps: %d\n', total_grains);

if total_grains > N_min
    fprintf('Bootstrapping down to %d grains (%d iterations)...\n', N_min, numBootstraps);
    bootstrap_J_values = zeros(numBootstraps, 1);

    for b = 1:numBootstraps
        random_idx = randperm(total_grains, N_min);
        temp_odf = calcDensity(all_mean_oris(random_idx), ...
            'weights', all_areas(random_idx), 'halfwidth', fixedHalfwidth);
        bootstrap_J_values(b) = norm(temp_odf)^2;
    end

    global_tex_index = mean(bootstrap_J_values);
    std_dev_J        = std(bootstrap_J_values);
    fprintf('Standardized (Bootstrapped) Global J-Index: %.3f ± %.3f\n', global_tex_index, std_dev_J);
else
    fprintf('Grain count <= %d. Calculating baseline J-index without resampling.\n', N_min);
    global_odf       = calcDensity(all_mean_oris, 'weights', all_areas, 'halfwidth', fixedHalfwidth);
    global_tex_index = norm(global_odf)^2;
    std_dev_J        = 0;
    fprintf('Baseline Global J-Index: %.3f\n', global_tex_index);
end

%% --- 7. Export Bootstrapped Summary Report to Excel ---
global_tex_column          = nan(numFiles, 1);
global_std_dev_column      = nan(numFiles, 1);
global_total_grains_column = nan(numFiles, 1);

global_tex_column(1)          = global_tex_index;
global_std_dev_column(1)      = std_dev_J;
global_total_grains_column(1) = total_grains;

tex_table = table(master_filenames, master_grain_counts, master_texture_indices, ...
    global_total_grains_column, global_tex_column, global_std_dev_column, ...
    'VariableNames', {'Filename', 'Map_Grain_Count', 'Map_Texture_Index', ...
    'Total_Pooled_Grains', 'Global_Texture_Index', 'Bootstrapped_Std_Dev'});

excelOutName = fullfile(outRoot, 'Texture_Indices_Report.xlsx');
writetable(tex_table, excelOutName);
fprintf('Saved Excel report to: %s\n', excelOutName);

%% --- 8. [Optional] Standalone Global Scatter PFs (Fixed Point Count = N_min) ---
% fprintf('\nGenerating Colored Global Scatter Pole Figures (points = %d)...\n', N_min);
%
% % {111} IPF-Z
% fig_scatter_Z = figure('Name', 'Global_Scatter_PF_111_Colored', 'Color', 'w', 'Position', [100, 650, 600, 500]);
% colors_Z = ipfKey_Z.orientation2color(all_mean_oris);
% plotPDF(all_mean_oris, colors_Z, Miller(1,1,1,cs_pt), ...
%     'how2plot', my_convention, 'MarkerSize', 3, 'MarkerEdgeAlpha', 0.8, 'points', N_min);
% title('Global Scatter Plot (111)', 'FontSize', 14);
% exportgraphics(fig_scatter_Z, fullfile(outRoot, 'Global_Scatter_PF_111_Colored.png'), ...
%     'Resolution', exportDPI, 'ContentType', 'image');
% close(fig_scatter_Z);
%
% % {110} IPF-Y
% fig_scatter_Y = figure('Name', 'Global_Scatter_PF_110_Colored', 'Color', 'w', 'Position', [720, 650, 600, 500]);
% colors_Y = ipfKey_Y.orientation2color(all_mean_oris);
% plotPDF(all_mean_oris, colors_Y, Miller(1,1,0,cs_pt), ...
%     'how2plot', my_convention, 'MarkerSize', 3, 'MarkerEdgeAlpha', 0.8, 'points', N_min);
% title('Global Scatter Plot (110)', 'FontSize', 14);
% exportgraphics(fig_scatter_Y, fullfile(outRoot, 'Global_Scatter_PF_110_Colored.png'), ...
%     'Resolution', exportDPI, 'ContentType', 'image');
% close(fig_scatter_Y);

%% --- 9. Prepare Global & Limited-Subset ODFs for Plotting ---
fprintf('\n--- Preparing Data for Global & Subset Plotting ---\n');

true_global_odf = calcDensity(all_mean_oris, 'weights', all_areas, 'halfwidth', fixedHalfwidth);

if total_grains > N_min
    fprintf('Extracting a random subset of %d grains for limited plots...\n', N_min);
    random_idx_plot = randperm(total_grains, N_min);
    limited_oris    = all_mean_oris(random_idx_plot);
    limited_areas   = all_areas(random_idx_plot);
else
    limited_oris    = all_mean_oris;
    limited_areas   = all_areas;
end
limited_odf = calcDensity(limited_oris, 'weights', limited_areas, 'halfwidth', fixedHalfwidth);

%% --- 10. Global & Limited-Subset Pole Figures ---
% True Global Scatter PFs (All Data)
exportGlobalScatterPF(all_mean_oris, ipfKey_Z, Miller(1,1,1,cs_pt), my_convention, ...
    sprintf('True Global Scatter Plot (111) (N = %d)', total_grains), ...
    0.8, fullfile(outRoot, 'True_Global_Scatter_Plot_111.png'), exportDPI);

exportGlobalScatterPF(all_mean_oris, ipfKey_Y, Miller(1,1,0,cs_pt), my_convention, ...
    sprintf('True Global Scatter {110} (N = %d)', total_grains), ...
    0.8, fullfile(outRoot, 'True_Global_Scatter_Plot_110.png'), exportDPI);

% Limited-Subset Scatter PFs (N_min)
exportGlobalScatterPF(limited_oris, ipfKey_Z, Miller(1,1,1,cs_pt), my_convention, ...
    sprintf('Limited Scatter Plot (111) (N = %d)', numel(limited_oris)), ...
    0.9, fullfile(outRoot, 'Limited_Scatter_111.png'), exportDPI);

exportGlobalScatterPF(limited_oris, ipfKey_Y, Miller(1,1,0,cs_pt), my_convention, ...
    sprintf('Limited Scatter Plot (110) (N = %d)', numel(limited_oris)), ...
    0.9, fullfile(outRoot, 'Limited_Scatter_Plot_110.png'), exportDPI);

% True Global & Limited-Subset Contour PFs
exportGlobalContourPF(true_global_odf, cs_pt, my_convention, ...
    sprintf('True Global Contour Plot (111) (N = %d)', total_grains), ...
    sprintf('True Global Contour Plot (110) (N = %d)', total_grains), ...
    fullfile(outRoot, 'True_Global_Contours_111_110.png'), exportDPI);

exportGlobalContourPF(limited_odf, cs_pt, my_convention, ...
    sprintf('Limited Contour {111} (N = %d)', numel(limited_oris)), ...
    sprintf('Limited Contour {110} (N = %d)', numel(limited_oris)), ...
    fullfile(outRoot, 'Limited_Contours_111_110.png'), exportDPI);

disp('--- Script Complete ---');

%% ========================================================================
%  Local Plotting Helper Functions
%  ========================================================================

function exportLocalScatterPF(oris, ipfKey, hkl, conv, outPath, dpi)
    fig = figure('Visible', 'off');
    colors = ipfKey.orientation2color(oris);
    plotPDF(oris, colors, hkl, 'how2plot', conv, 'all', 'MarkerSize', 3);
    exportgraphics(fig, outPath, 'Resolution', dpi, 'ContentType', 'image');
    close(fig);
end

function exportLocalContourPF(odf, cs, conv, equalScale, outPath, dpi)
    fig = figure('Visible', 'off');
    plotPDF(odf, {Miller(1,1,1,cs), Miller(1,1,0,cs)}, 'how2plot', conv);
    if equalScale
        setColorRange('equal');
    end
    cb = mtexColorbar;
    set(cb, 'FontSize', 24);
    exportgraphics(fig, outPath, 'Resolution', dpi, 'ContentType', 'image');
    close(fig);
end

function exportGlobalScatterPF(oris, ipfKey, hkl, conv, plotTitle, edgeAlpha, outPath, dpi)
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 600, 500]);
    colors = ipfKey.orientation2color(oris);
    plotPDF(oris, colors, hkl, 'how2plot', conv, ...
        'MarkerSize', 2, 'MarkerEdgeAlpha', edgeAlpha, 'all');
    title(plotTitle, 'FontSize', 14);
    exportgraphics(fig, outPath, 'Resolution', dpi);
    close(fig);
end

function exportGlobalContourPF(odf, cs, conv, title111, title110, outPath, dpi)
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 1200, 500]);
    plotPDF(odf, {Miller(1,1,1,cs), Miller(1,1,0,cs)}, 'how2plot', conv);
    setColorRange('equal');
    cb = mtexColorbar;
    set(cb, 'FontSize', 16);
    drawNow(gcm);
    axes_list = gcm().children;
    title(axes_list(1), title111, 'FontSize', 16);
    title(axes_list(2), title110, 'FontSize', 16);
    exportgraphics(fig, outPath, 'Resolution', dpi, 'Padding', 'tight');
    close(fig);
end
