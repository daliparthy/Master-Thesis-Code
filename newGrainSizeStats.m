%% === Master Batch Script: Batch Processing & Pooled Grain Statistics ===
clear variables; close all; clc;
cs_pt = crystalSymmetry('m-3m', [3.92 3.92 3.92], 'mineral', 'Pt');  

% --- 0. MTEX Preferences ---
setMTEXpref('xAxisDirection','east');
%setMTEXpref('yAxisDirection','south');
setMTEXpref('zAxisDirection', 'intoPlane');
% Added to ensure PF colors render nicely
setMTEXpref('defaultColorMap', parula); 

% Define the input and output directories
data_dir = 'Q:\EBSD\Master Thesis\5nm\30k\vacuum\ang files'; 
outDir = 'Q:\EBSD\Master Thesis\5nm\30k\vacuum\newGrainSizeStats';

if ~exist(outDir, 'dir'); mkdir(outDir); end

ang_files = dir(fullfile(data_dir, '*.ang'));

ph1 = 'Pt';
ph2 = 'Al2O3'; 

my_convention = plottingConvention;
my_convention.east = xvector;        
my_convention.outOfScreen = -zvector;

% Initialize empty arrays for pooling
pooled_diameters = [];
pooled_areas = [];
pooled_perimeters = [];
pooled_circularities = [];
pooled_map_names = strings(0, 1); 

%% --- 1. Map Processing Loop ---
for i = 1:length(ang_files)
    file_path = fullfile(data_dir, ang_files(i).name);
    [~, base_name, ~] = fileparts(ang_files(i).name);
    fprintf('Processing %s...\n', base_name);

    % === Plotting Convention ===
    my_convention = plottingConvention;
    my_convention.east = xvector;        
    my_convention.outOfScreen = -zvector; 
    
    % === Load Data & Gridify ===
    ebsd = EBSD.load(file_path, 'convertEuler2SpatialReferenceFrame', 'setting 2');
    ebsd = ebsd.gridify;

    % ---> CRITICAL FIX: OVERRIDE THE .ANG FILE SYMMETRY <---
    % Force the full cubic symmetry immediately so calcGrains, boundaries, 
    % and the IPF Color Key inherit the correct standard triangle.
    ebsd(ph1).CS = crystalSymmetry('m-3m', [3.92 3.92 3.92], 'mineral', 'Pt');

    phases_present = ebsd.mineralList;
    has_alumina = ismember(ph2, phases_present);
    
    % Extract the newly corrected symmetry for the color keys later
    cs_pt = ebsd(ph1).CS;
    if has_alumina
        cs_alumina = ebsd(ph2).CS;
    end

   %% 1a. Grain Calc & Filtering Logic
    ebsd_pt = ebsd(ph1);
    bad_pt_condition = ebsd_pt.ci < 0.2;
    bad_pt_pixels = ebsd_pt(bad_pt_condition);
    ebsd(bad_pt_pixels.id).phase = 0;
    
    [grains, ebsd.grainId] = calcGrains(ebsd, 'angle', 10*degree, 'minPixel', 2);
    
    if has_alumina
        alumina_grains = grains(ph2);
        if ~isempty(alumina_grains)
            % FIX 1: Set to phase 0 instead of deleting to preserve the spatial grid
            bad_alumina = alumina_grains(alumina_grains.numPixel < 20);
            ebsd(bad_alumina).phase = 0; 
        end
    end
    
    [grains, ebsd.grainId] = calcGrains(ebsd, 'angle', 10*degree);
    grains_smoothed = smooth(grains, 3, 'moveTriplePoints', 'umbrella');

    %% 1b. Boundary Calculations & Geometric Edge Filtering
    pt_grains = grains_smoothed(ph1);
    
   % FIX 2: Strict coordinate-based edge filtering (ignores internal holes)
    % Find the absolute geometric limits of the scan area (force scalars using (:))
    x_min = min(ebsd.x(:)); x_max = max(ebsd.x(:));
    y_min = min(ebsd.y(:)); y_max = max(ebsd.y(:));
    
    % Calculate tolerance based on average pixel step size
    step_size = sqrt((x_max - x_min) * (y_max - y_min) / length(ebsd));
    tol = 1.5 * step_size; % 1.5 pixels tolerance
    
    % Find any EBSD pixel that touches the absolute outer bounding box
    is_edge_pixel = (ebsd.x <= x_min + tol) | (ebsd.x >= x_max - tol) | ...
                    (ebsd.y <= y_min + tol) | (ebsd.y >= y_max - tol);
                    
    % Identify which grains own those specific border pixels
    edge_grain_ids = unique(ebsd(is_edge_pixel).grainId);
    
    % Filter out ONLY those outer geometric edge grains
    pt_clean = pt_grains(~ismember(pt_grains.id, edge_grain_ids));
   

    % Restrict boundary calculations to ONLY the clean grains
    gb_phase1       = pt_clean.boundary(ph1, ph1);
    gb_phase2       = pt_clean.boundary(ph1, 'notIndexed');
    gb_phase1_inner = pt_clean.innerBoundary;
    
    if has_alumina
        gb_phase3 = pt_clean.boundary(ph1, ph2);
    else
        gb_phase3 = []; 
    end
    %% 1c. Color Keys & IPF-Z Map with GBs
    if ~isempty(pt_clean)
        ph1ColorKey = ipfColorKey(cs_pt);
        ph1ColorKey.inversePoleFigureDirection = vector3d.Z;
        
        % Extract ONLY the EBSD pixels that belong to the non-edge grains
        ebsd_pt_clean = ebsd(pt_clean);
        
        % Calculate colors only for these clean pixels
        ph1Colors_Z = ph1ColorKey.orientation2color(ebsd_pt_clean.orientations);
    end
    
    if has_alumina
        ph2ColorKey = ipfColorKey(cs_alumina);
        ph2ColorKey.inversePoleFigureDirection = vector3d.Z;
        ph2Colors = ph2ColorKey.orientation2color(ebsd(ph2).orientations);
    end
    
    figure('Name','IPF-z_gb', 'Visible', 'off');
    
    % Plot the mapped area using the CLEANED pixel data
    plot(ebsd_pt_clean, ph1Colors_Z, 'figSize', 'huge', 'faceAlpha', 1, 'how2plot', my_convention);
    hold on;
    
    if has_alumina
        plot(ebsd(ph2), ph2Colors, 'faceAlpha', 1, 'how2plot', my_convention);
    end
    
    h1  = plot(gb_phase1, 'k', 'LineWidth', 1.0, 'LineJoin', 'round', 'Displayname', 'GBs > 10°', 'how2plot', my_convention);
    h11 = plot(gb_phase2, 'k', 'LineWidth', 1.0, 'LineJoin', 'round', 'how2plot', my_convention);
    if ~isempty(gb_phase3)
        h12 = plot(gb_phase3, 'k', 'LineWidth', 1.0, 'LineJoin', 'round', 'how2plot', my_convention);
    end
    
    h_leg = h1; 
    if isempty(gb_phase1_inner)
        disp('No inner boundaries found in this map.');
    else
        alpha_z = gb_phase1_inner.misorientation.angle / (5*degree);
        h2 = plot(gb_phase1_inner, 'LineWidth', 1.5, 'LineJoin', 'round', 'edgeAlpha', alpha_z, 'linecolor', 'grey', 'Displayname', 'GBs < 10°');
        h_leg = [h_leg, h2];
    end
    
    hold off;
    axis image off;
    legend(h_leg);
    saveas(gcf, fullfile(outDir, sprintf('%s_IPF_Z_GB.png', base_name)));
    close(gcf);
 
    %% 1d. Pool Metrics
    % Since pt_clean is already calculated above, we just extract the math
    valid_areas = pt_clean.area;
    valid_perims = pt_clean.perimeter;
    
    valid_diams = 2 * sqrt(valid_areas / pi);
    valid_circs = (4 * pi * valid_areas) ./ (valid_perims .^ 2);

    num_valid_grains = length(valid_areas);
    current_map_names = repmat(string(base_name), num_valid_grains, 1);

    pooled_areas = [pooled_areas; valid_areas];
    pooled_perimeters = [pooled_perimeters; valid_perims];
    pooled_diameters = [pooled_diameters; valid_diams];
    pooled_circularities = [pooled_circularities; valid_circs];
    pooled_map_names = [pooled_map_names; current_map_names];
end

%% --- 2. Export Individual Grain Data (Text File) ---
T = table((1:numel(pooled_areas))', pooled_map_names, pooled_areas, pooled_diameters, pooled_perimeters, pooled_circularities, ...
    'VariableNames', {'GrainID', 'Source_Map', 'Area (μm2)', 'Equivalent Circular Diameter (μm)', 'Perimeter (μm)', 'Circularity'});
outFile = fullfile(outDir, 'Pooled_grain_sizes_with_perimeter.txt');
writetable(T, outFile, 'Delimiter', ',');
fprintf('\nBatch processing complete. Exported %d pooled grains to:\n%s\n', length(pooled_diameters), outFile);

%% --- 3. Calculate Pooled Statistics ---
num_mean = mean(pooled_diameters);
num_std = std(pooled_diameters);

area_mean = sum(pooled_diameters .* pooled_areas) / sum(pooled_areas);
area_variance = sum(pooled_areas .* (pooled_diameters - area_mean).^2) / sum(pooled_areas);
area_std = sqrt(area_variance);

fprintf('\n--- Final Pooled Statistics ---\n');
fprintf('Number-Weighted Mean: %.3f ± %.3f µm\n', num_mean, num_std);
fprintf('Area-Weighted Mean: %.3f ± %.3f µm\n\n', area_mean, area_std);

%% --- 4. Export Summary Statistics to Excel ---
stats_file = fullfile(outDir, 'Pooled_Summary_Statistics.xlsx');
SummaryTable = table(length(pooled_diameters), num_mean, num_std, area_mean, area_std, ...
    'VariableNames', {'Total_Grains', 'Num_Weighted_Mean_um', 'Num_Weighted_Std_um', 'Area_Weighted_Mean_um', 'Area_Weighted_Std_um'});
writetable(SummaryTable, stats_file);

%% --- 5. Extract Histogram Data & Export to Excel ---
% Define FIXED bin edges to guarantee comparability across all samples
% Format: [Start_Size : Bin_Width : Max_Expected_Size]
fixed_edges = 0 : 0.018 : 0.6; 

% Use the fixed edges instead of a generic number of bins
[numberCounts, edges, binIdx] = histcounts(pooled_diameters, fixed_edges);
binCenters = edges(1:end-1) + diff(edges)/2;

numberFraction = numberCounts / sum(numberCounts);

% The number of bins is now determined by your fixed_edges array
numBins = length(binCenters);
weightedCounts = zeros(1, numBins);
for b = 1:numBins
    weightedCounts(b) = sum(pooled_areas(binIdx == b));
end
areaFraction = weightedCounts / sum(pooled_areas); 

col_BinCenter = binCenters(:);
col_NumCount  = numberCounts(:);
col_NumFrac   = numberFraction(:);
col_AreaFrac  = areaFraction(:);

hist_file = fullfile(outDir, 'Pooled_Histogram_Data_Only.xlsx');
HistTable = table(col_BinCenter, col_NumCount, col_NumFrac, col_AreaFrac, ...
    'VariableNames', {'Equivalent Circular Diameter (µm)', 'Number Count', 'Number Fraction', 'Area Fraction'});
writetable(HistTable, hist_file);

fprintf('Exported summary statistics to:\n%s\n', stats_file);
fprintf('Exported histogram plot data to:\n%s\n', hist_file);