%% === Master Batch Script: Batch Processing & Global CSL Length Fractions ===
clear variables; close all; clc;

% --- 0. MTEX Preferences ---
setMTEXpref('xAxisDirection','east');
%setMTEXpref('yAxisDirection','south');
setMTEXpref('zAxisDirection', 'intoPlane');
% Added to ensure PF colors render nicely
setMTEXpref('defaultColorMap', parula); 

% --- 1. Define Input and Output Folders ---
dataFolder = 'Q:\EBSD\Master Thesis\10nm_depRateOff\Ar O2\ang files';
outRoot    = 'Q:\EBSD\Master Thesis\10nm_depRateOff\Ar O2\2026_09_07';
fileList   = dir(fullfile(dataFolder, '*.ang'));

% --- 2. Global Variables for CSL Length Aggregation ---
maxSigma = 23;
sigmaValues = 3:2:maxSigma; % Starting at 3. Sigma 1 is LAGB (<10 deg), which are filtered out.
master_csl_lengths = zeros(length(sigmaValues), 1); 
master_total_gb_length = 0; 

numFiles   = length(fileList);
if numFiles == 0
    error('No .ang files found in the specified directory.');
end

% --- 3. Initialize Tracking Variables ---
master_filenames = strings(numFiles, 1);
master_texture_indices = zeros(numFiles, 1);
ori_list = cell(numFiles, 1); % To hold the raw data for pooling

% --- 4. Define Crystal Symmetry (for Platinum) ---
cs_pt = crystalSymmetry('m-3m', [3.92 3.92 3.92], 'mineral', 'Pt');  
ph1 = 'Pt';
%% --- START LOOP OVER ALL MAPS ---
for i = 1:length(fileList)
    filename = fullfile(fileList(i).folder, fileList(i).name);
    [~, baseName, ~] = fileparts(filename);
    fprintf('Processing Map %d of %d: %s\n', i, length(fileList), fileList(i).name);
    master_filenames{i} = baseName;
    outDir = fullfile(outRoot, baseName, 'IPF_z_and_y'); 
    if ~exist(outDir, 'dir'); mkdir(outDir); end
    
    cleanName = @(s) regexprep(strtrim(s),'[^\w\-]+','_');
    
    % === Plotting Convention ===
    my_convention = plottingConvention;
    my_convention.east = xvector;        
    my_convention.outOfScreen = -zvector; 
    
    % === Load Data & Gridify ===
    ebsd = EBSD.load(filename, 'convertEuler2SpatialReferenceFrame', 'setting 2');
    ebsd = ebsd.gridify;

    % disp(ebsd)
    % Adding (:) flattens the grid and forces MTEX to print the phase summary table
    ebsd(:)
    
    if isempty(ebsd(ph1))
        disp('  -> No Platinum found in this map. Skipping.');
        master_texture_indices(i) = NaN;
        continue;
    end

    % --- Manual Crystal Symmetries ---
    cs_pt = crystalSymmetry('m-3m', [3.92 3.92 3.92], 'mineral', 'Pt');  
    ph1 = 'Pt';
    ph2 = 'alumina';

    try
        cs_alumina = ebsd(ph2).CS; 
    catch
        cs_alumina = crystalSymmetry('-3m1', [4.76 4.76 12.99], 'mineral', 'alumina');
        disp('Notice: No alumina symmetry found in file. Using manual definition.');
    end
    
    has_alumina = length(ebsd(ph2)) > 0;
    if ~has_alumina
        disp('Notice: No alumina pixels found in this map. Proceeding with Pt only.');
    end
    
    % === Plot Confidence Index (CI) figure ===
    figure('Name','CI-Map', 'Visible', 'off'); 
    plot(ebsd,ebsd.ci,'how2plot', my_convention)
    mtexColorMap black2white
    setColorRange([0 1])
    mtexColorbar
    axis image off; % Call this BEFORE setting the colorbar position
    
    %% Grain Calc & Filtering
    ebsd_pt = ebsd('Pt');
    bad_pt_condition = ebsd_pt.ci < 0.2;
    bad_pt_pixels = ebsd_pt(bad_pt_condition);
    ebsd(bad_pt_pixels.id).phase = 0;
    
    [grains, ebsd.grainId] = calcGrains(ebsd, 'angle', 10*degree, 'minPixel', 2);
    
    pt_grains = grains(ph1);
    alumina_grains = grains(ph2);
    
    if ~isempty(alumina_grains)
        ebsd(alumina_grains(alumina_grains.numPixel < 20)) = [];
    end
    
    [grains, ebsd.grainId] = calcGrains(ebsd, 'angle', 10*degree);
    grains_smoothed = smooth(grains, 3, 'moveTriplePoints', 'umbrella');
    
    gb_phase1       = grains_smoothed.boundary(ph1, ph1);
    gb_phase2       = grains_smoothed.boundary(ph1, 'notIndexed');
    gb_phase3       = grains_smoothed.boundary(ph1, ph2);
    gb_phase1_inner = grains_smoothed(ph1).innerBoundary;
    
   %% CSL AGGREGATION (Background Global Tracking)
    if ~isempty(gb_phase1)
        master_total_gb_length = master_total_gb_length + sum(gb_phase1.segLength);
        
        for v = 1:length(sigmaValues)
            sig = sigmaValues(v);
            tol = (15 / sqrt(sig)) * degree;
            
            ideal_CSL = CSL(sig, ebsd(ph1).CS); 
            
            is_csl = gb_phase1.isTwinning(ideal_CSL, tol);
            
            csl_len = sum(gb_phase1(is_csl).segLength);
            master_csl_lengths(v) = master_csl_lengths(v) + csl_len;
        end
    end
    %% Color Keys (Z-Direction Initialized)
    if ~isempty(pt_grains)
        ph1ColorKey = ipfColorKey(cs_pt);
        ph1ColorKey.inversePoleFigureDirection = vector3d.Z;
        ph1Colors_Z = ph1ColorKey.orientation2color(ebsd(ph1).orientations);
    end
    if has_alumina
        ph2ColorKey = ipfColorKey(cs_alumina);
        ph2ColorKey.inversePoleFigureDirection = vector3d.Z;
        ph2Colors = ph2ColorKey.orientation2color(ebsd(ph2).orientations);
    end
    
    %% GB + IPF Map (Z-Direction)
    figure('Name','IPF-z_gb', 'Visible', 'off');
    plot(ebsd(ph1), ph1Colors_Z, 'figSize', 'huge', 'faceAlpha', 1, 'how2plot', my_convention);
    hold on;
    if has_alumina
        plot(ebsd(ph2), ph2Colors, 'faceAlpha', 1, 'how2plot', my_convention);
    end
    h1  = plot(gb_phase1, 'k', 'LineWidth', 1.0, 'LineJoin', 'round', 'Displayname', 'GBs > 10°', 'how2plot', my_convention);
    h11 = plot(gb_phase2, 'k', 'LineWidth', 1.0, 'LineJoin', 'round', 'how2plot', my_convention);
    if ~isempty(gb_phase3)
        h12 = plot(gb_phase3, 'k', 'LineWidth', 1.0, 'LineJoin', 'round', 'how2plot', my_convention);
    end
    if isempty(gb_phase1_inner)
        disp('No inner boundaries found in this map.');
    else
        % Safe to process or plot the boundaries
        alpha_z = gb_phase1_inner.misorientation.angle / (5*degree);
        h2 = plot(gb_phase1_inner, 'LineWidth', 1.5, 'LineJoin', 'round', 'edgeAlpha', alpha_z, 'linecolor', 'grey', 'Displayname', 'GBs < 10°');
    end
    
    hold off;
    axis image off;
    legend([h1 h2]);
    figure('Name','IPF color key', 'Visible', 'off');
    plot(ph1ColorKey);
    
    %% GB + IPF Map (Y-Direction)
    ph1ColorKey.inversePoleFigureDirection = vector3d.Y; 
    ph1Colors_Y = ph1ColorKey.orientation2color(ebsd(ph1).orientations);
    figure('Name','IPF-y_gb', 'Visible', 'off');
    plot(ebsd(ph1), ph1Colors_Y, 'figSize', 'huge', 'faceAlpha', 1, 'how2plot', my_convention);
    hold on;
    if has_alumina
        plot(ebsd(ph2), ph2Colors, 'faceAlpha', 1, 'how2plot', my_convention);
    end
    h1_y = plot(gb_phase1, 'k', 'LineWidth', 1.0, 'LineJoin', 'round','Displayname', 'GBs > 10°', 'how2plot', my_convention);
    h11_y = plot(gb_phase2, 'k', 'LineWidth', 1.0, 'LineJoin', 'round', 'how2plot', my_convention);
    if ~isempty(gb_phase3)
        h12_y = plot(gb_phase3, 'k', 'LineWidth', 1.0, 'LineJoin', 'round', 'how2plot', my_convention);
    end
    alpha_y = gb_phase1_inner.misorientation.angle / (5*degree);
    h2_y= plot(gb_phase1_inner, 'LineWidth', 1.5, 'LineJoin', 'round', 'edgeAlpha', alpha_y, 'linecolor', 'grey','Displayname', 'GBs < 10°');
    hold off;
    axis image off;
    legend([h1_y h2_y]);
    
    %% CSL Boundaries Map (on IPF-Y Map)
    CSL_Boundaries = [3 7 13 21];
    Boundary_color = ['r' 'g' 'y' 'c'];
    gB = grains_smoothed.boundary(ph1, ph1);
    
    figure('Name', 'CSL_Boundaries_Y', 'Visible', 'off');
    plot(ebsd(ph1), ph1Colors_Y, 'figSize', 'huge', 'faceAlpha', 1, 'how2plot', my_convention);
    hold on;
    if has_alumina
        plot(ebsd(ph2), ph2Colors, 'faceAlpha', 1, 'how2plot', my_convention);
    end
    plot(gb_phase1, 'k', 'LineWidth', 1.0, 'LineJoin', 'round', 'Displayname', 'GBs > 10°', 'how2plot', my_convention);
    plot(gb_phase2, 'k', 'LineWidth', 1.0, 'LineJoin', 'round', 'how2plot', my_convention);
    if ~isempty(gb_phase3)
        plot(gb_phase3, 'k', 'LineWidth', 1.0, 'LineJoin', 'round', 'how2plot', my_convention);
    end
    if ~isempty(gb_phase1_inner)
        plot(gb_phase1_inner, 'LineWidth', 1.5, 'edgeAlpha', alpha_y, 'linecolor', 'grey', 'Displayname', 'GBs < 10°');
    else
        plot(NaN, NaN, 'color', [0.5 0.5 0.5], 'DisplayName', 'GBs < 10°'); 
    end
    if ~isempty(gB)
        for v = 1:length(CSL_Boundaries)
            ideal_CSL = CSL(CSL_Boundaries(v), ebsd(ph1).CS);
            tol_visual = 15/sqrt(CSL_Boundaries(v))*degree;
            gB_CSL = gB(gB.isTwinning(ideal_CSL, tol_visual));
            if ~isempty(gB_CSL)
                plot(gB_CSL, 'lineColor', Boundary_color(v), 'linewidth', 2, 'LineJoin', 'round', ...
                    'DisplayName', ['CSL \Sigma' num2str(CSL_Boundaries(v))]);
            end
        end
    else
        disp('No Pt-Pt grain boundaries found in this map. Skipping CSL visual plot.');
    end
    hold off;
    axis image off;
    legend('Location', 'northeast');
    
    %% CSL Boundaries Map (on IPF-Z Map)
    CSL_Boundaries = [3 7 13 21];
    Boundary_color = ['r' 'g' 'y' 'c'];
    gB = grains_smoothed.boundary(ph1, ph1);
    figure('Name', 'CSL_Boundaries_Z', 'Visible', 'off');
    plot(ebsd(ph1), ph1Colors_Z, 'figSize', 'huge', 'faceAlpha', 1, 'how2plot', my_convention);
    hold on;
    if has_alumina
        plot(ebsd(ph2), ph2Colors, 'faceAlpha', 1, 'how2plot', my_convention);
    end
    plot(gb_phase1, 'k', 'LineWidth', 1.0, 'LineJoin', 'round', 'Displayname', 'GBs > 10°', 'how2plot', my_convention);
    plot(gb_phase2, 'k', 'LineWidth', 1.0, 'LineJoin', 'round', 'how2plot', my_convention);
    if ~isempty(gb_phase3)
        plot(gb_phase3, 'k', 'LineWidth', 1.0, 'LineJoin', 'round', 'how2plot', my_convention);
    end
    if ~isempty(gb_phase1_inner)
        plot(gb_phase1_inner, 'LineWidth', 1.5, 'edgeAlpha', alpha_y, 'linecolor', 'grey', 'Displayname', 'GBs < 10°');
    else
        plot(NaN, NaN, 'color', [0.5 0.5 0.5], 'DisplayName', 'GBs < 10°'); 
    end
    if ~isempty(gB)
        for v = 1:length(CSL_Boundaries)
            ideal_CSL = CSL(CSL_Boundaries(v), ebsd(ph1).CS);
            tol_visual = 15/sqrt(CSL_Boundaries(v))*degree;
            gB_CSL = gB(gB.isTwinning(ideal_CSL, tol_visual));
            if ~isempty(gB_CSL)
                plot(gB_CSL, 'lineColor', Boundary_color(v), 'linewidth', 2, 'LineJoin', 'round', ...
                    'DisplayName', ['CSL \Sigma' num2str(CSL_Boundaries(v))]);
            end
        end
    else
        disp('No Pt-Pt grain boundaries found in this map. Skipping CSL visual plot.');
    end
    hold off;
    axis image off;
    legend('Location', 'northeast');
    
    %% Compute Grain Perimeter and Circularity
    grainPerimeter = grains_smoothed('Pt').perimeter;  
    grainAreas = grains_smoothed('Pt').area;            
    grainDiam = 2 * sqrt(grainAreas / pi);
    grainCircularity = (4 * pi * grainAreas) ./ (grainPerimeter .^ 2);
    T = table((1:numel(grains_smoothed('Pt')))', grainAreas, grainDiam, grainPerimeter, grainCircularity, ...
        'VariableNames', {'GrainID','Area_um2','EqCircleDiam_um','Perimeter_um','Circularity'});
    outFile = fullfile(outDir, 'grain_sizes_with_perimeter.txt');
    writetable(T, outFile, 'Delimiter', ',');
    
    %% Grain Size Distributions & Unindexed Checks
    figure('Name', 'Unindexed Check', 'Visible', 'off');
    plot(ebsd('Pt'),'how2plot', my_convention, 'faceColor', [0.7 0.7 0.7]); 
    hold on;
    if has_alumina; plot(ebsd('alumina'), 'how2plot', my_convention,'faceColor', 'r'); end
    plot(ebsd('notIndexed'), 'how2plot', my_convention,'faceColor', 'y'); 
    hold off;
    axis image off;
    valid_idx = grainAreas > 0;         
    grainAreasValid = grainAreas(valid_idx);  
    grainDiamValid = grainDiam(valid_idx);    
    figure('Name', 'Area_Weighted_GSD', 'Visible', 'off');
    numBins = 40;
    [~, edges, binIdx] = histcounts(grainDiamValid, numBins);
    weightedCounts = zeros(1, numBins);
    for b = 1:numBins
        weightedCounts(b) = sum(grainAreasValid(binIdx == b));
    end
    areaFraction = weightedCounts / sum(grainAreasValid); 
    binCenters = edges(1:end-1) + diff(edges)/2;
    bar(binCenters, areaFraction, 'FaceColor', '#D95319', 'BarWidth', 1);
    xlabel('Equivalent Circular Diameter (\mu m)', 'FontSize', 14);
    ylabel('Area Fraction', 'FontSize', 14);
    title('Area-Weighted Grain Size Distribution', 'FontSize', 16);
    figure('Name', 'Number_Weighted_GSD', 'Visible', 'off');
    histogram(grainDiamValid, numBins, 'FaceColor', '#0072BD');
    xlabel('Equivalent Circular Diameter (\mu m)', 'FontSize', 14);
    ylabel('Number of Grains', 'FontSize', 14);
    title('Number-Weighted Grain Size Distribution', 'FontSize', 16);
    
  %% === POLE FIGURES & LINE PROFILES ===
    if ~isempty(ebsd(ph1))
        
        % --- 1. Scatter Plot Z (111) ---
        figure('Name', 'Pt_111_ScatterPlot', 'Visible', 'off');
        ipfKey = ipfColorKey(cs_pt);
        ipfKey.inversePoleFigureDirection = zvector; % Color based on Z
        colors_Z = ipfKey.orientation2color(ebsd(ph1).orientations);
        plotPDF(ebsd(ph1).orientations, colors_Z, Miller(1,1,1,cs_pt), 'all', 'MarkerSize', 3);
        
        % --- 2. Scatter Plot Y (110) ---
        figure('Name', 'Pt_110_ScatterPlot', 'Visible', 'off');
        ipfKey.inversePoleFigureDirection = yvector; % Recalculate color based on Y
        colors_Y = ipfKey.orientation2color(ebsd(ph1).orientations);
        plotPDF(ebsd(ph1).orientations, colors_Y, Miller(1,1,0,cs_pt), 'all', 'MarkerSize', 3);
        
       % --- 3. ODF Calculation & Contours ---
        odf = calcDensity(ebsd(ph1).orientations);
        
        % --- NEW: Calculate and Store the Texture Index (J-Index) ---
        current_tex_index = norm(odf)^2;
        master_texture_indices(i) = current_tex_index;
        fprintf('Calculated Texture Index for %s: %.3f\n', baseName, current_tex_index);
        
        % ---> ADD THIS LINE HERE: Store the CLEANED orientations for the global pool
        ori_list{i} = ebsd(ph1).orientations(:);

        figure('Name', 'Pt_111_110_ContourPlot_DifferentScales', 'Visible', 'off');
        plotPDF(odf, {Miller(1,1,1,cs_pt), Miller(1,1,0,cs_pt)});
        cb = mtexColorbar;   
        set(cb, 'FontSize', 24);    
        
        figure('Name', 'Pt_111_110_ContourPlot_SameScale', 'Visible', 'off');
        plotPDF(odf, {Miller(1,1,1,cs_pt), Miller(1,1,0,cs_pt)});
        setColorRange('equal');
        cb = mtexColorbar;   
        set(cb, 'FontSize', 24);    
        
        % --- 4. Line Profiles ---
        % Setup points on a great circle for evaluation
        v = rotation.byAxisAngle(yvector, (0:1:180)*degree) * xvector;
        ang = (-90:1:90)'; 
        
        % Lineprofile through x (111)
        pdf1 = calcPDF(odf, Miller(1,1,1,odf.CS));
        values_111 = eval(pdf1, v);
        T111 = table(ang, values_111(:), 'VariableNames', {'Angle_deg','Intensity'});
        
        figure('Name', 'Pt_111_LineProfile_FWHM', 'Visible', 'off');
        plot(pdf1);
        nextAxis;
        plot(v, values_111);
        
        % Lineprofile through x (110)
        pdf2 = calcPDF(odf, Miller(1,1,0,odf.CS));
        values_110 = eval(pdf2, v);
        T110 = table(ang, values_110(:), 'VariableNames', {'Angle_deg','Intensity'});
        
        figure('Name', 'Pt_110_LineProfile', 'Visible', 'off');
        plot(pdf2);
        nextAxis;
        plot(v, values_110);
        
        % --- 5. Save Line Profile Text Data ---
        file111 = fullfile(outDir, sprintf('%s_lineprofile_111.txt', baseName));
        file110 = fullfile(outDir, sprintf('%s_lineprofile_110.txt', baseName));
        writetable(T111, file111, 'Delimiter','\t');
        writetable(T110, file110, 'Delimiter','\t');
    end

    %% Save all opened figures
    figs = findall(groot, 'Type', 'figure');
    if ~isempty(figs)
        [~, idx] = sort([figs.Number]);
        figs = figs(idx);
        
        for k = 1:numel(figs)
            f = figs(k);
            nm = get(f, 'Name');
            if isempty(nm); nm = sprintf('Figure_%02d', f.Number); end
            
            baseFile = sprintf('%04d_%s_%02d_%s', i, baseName, k, cleanName(nm));
            pngPath = fullfile(outDir, [baseFile '.png']);
            figPath = fullfile(outDir, [baseFile '.fig']);
            
            % Force exact layout sizing for robust saving (especially for PFs)
            set(f, 'PaperPositionMode', 'auto'); 
            f.Position = [100, 100, 1000, 800]; 
            
            exportgraphics(f, pngPath, 'Resolution', 600, 'ContentType', 'image', 'Padding', 'tight');
            savefig(f, figPath);
        end
        
        close all; 
    else
        fprintf('No figures generated for Map %d: %s. Bypassing save.\n', i, baseName);
    end
end
fprintf('--- All Maps Processed and Saved! ---\n');

%% --- POOL DATA AND CALCULATE GLOBAL ODF ---
fprintf('\nConcatenating all orientations vertically...\n');

% vertcat safely stacks the columns on top of each other regardless of length
all_orientations = vertcat(ori_list{:}); 

fprintf('Total pooled data points: %d\n', length(all_orientations));
fprintf('Calculating Global ODF...\n');

% Calculate the single true global ODF and its index
global_odf = calcDensity(all_orientations);
global_tex_index = norm(global_odf)^2;

fprintf('  -> TRUE GLOBAL TEXTURE INDEX (J): %.3f\n\n', global_tex_index);


%% === Final Step: Calculate & Plot Global CSL Length Fractions ===
if master_total_gb_length > 0
    combinedLengthFractions = (master_csl_lengths / master_total_gb_length) * 100; 
    
    % Set to 'on' so the final summary pops up when the script finishes
    figure('Name', 'Combined_CSL_Length_Fraction_Distribution', 'Visible', 'on');
    bar(sigmaValues, combinedLengthFractions, 1, 'FaceColor', '#0072BD', 'EdgeColor', 'k');
    title('Combined CSL Boundaries Length Fraction', 'FontSize', 14);
    xlabel('\Sigma Value', 'FontSize', 12); 
    ylabel('Length Fraction (%)', 'FontSize', 12);
    xticks(sigmaValues);
    xticklabels(string(sigmaValues));
    xlim([0 maxSigma+2]); 
    box on; 
    grid on; 
    set(gca, 'TickDir', 'out', 'FontSize', 10, 'YGrid', 'on', 'XGrid', 'off');
    
    globalPlotOutName = 'Combined_CSL_LengthFraction_BarChart.png';
    exportgraphics(gcf, fullfile(outRoot, globalPlotOutName), 'Resolution', 600, 'ContentType', 'image');
    fprintf('Saved Global CSL Length Fraction Plot to: %s\n', fullfile(outRoot, globalPlotOutName));

    %% --- Final Step: Export to CSV ---
    % Create a column of empty values (NaN writes as a blank cell in CSV)
    global_tex_column = NaN(numFiles, 1);

    % Place the global index into the very first cell of that column
    global_tex_column(1) = global_tex_index;

    % Create the final table
    tex_table = table(master_filenames, master_texture_indices, global_tex_column, ...
        'VariableNames', {'Filename', 'Individual_Texture_Index', 'Global_Texture_Index'});

    csvOutName = fullfile(outRoot, 'Batch_and_Global_Texture_Indices.csv');
    writetable(tex_table, csvOutName);
    fprintf('Saved Texture Indices to: %s\n', csvOutName);

%% === Export Global CSL Data for OriginPro ===
    % Transpose sigmaValues to match the column orientation of combinedLengthFractions
    sigma_col = sigmaValues(:);
    fraction_col = combinedLengthFractions(:);
    % Create a table with clear headers for Origin
    csl_table = table(sigma_col, fraction_col, ...
        'VariableNames', {'Σ Value', 'Length Fraction (%)'});
    % Define the output file name and path in the root output folder
    csvOutName = fullfile(outRoot, 'Global_CSL_Length_Fractions.csv');
    % Write the table to a CSV file
    writetable(csl_table, csvOutName);
    fprintf('Saved Global CSL Data for OriginPro to: %s\n', csvOutName);
end

