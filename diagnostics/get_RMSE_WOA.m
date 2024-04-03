% calculate RMSEs along trajectories compares to WOA data (DIN)

%% load model output 
close all
clear all


modTag = 'optimisedParams_AH_wPOM_10_pmax_b_-0,08_Qmin_QC_a_0.035_theta_4.2_Gmax_a_11_m2_5e-2_atLargeSC_rDOM_0.04_icelessArcSummer_adjustedVin_fullRunForPaper';

%% add paths 
addpath(genpath('~/Documents/microARC model')) % location of my work with the model
addpath(genpath('~/microARC/')) % location of model


% set directories
Directories.resultsDir  = ['~/Documents/microARC model/2nd paper/mod output/' modTag '/'];
Directories.plotDir  = ['~/Documents/microARC model/2nd paper/mod output/' modTag '/plots/'];


% load model output
summer = load([Directories.resultsDir  'summer_output_' modTag '.mat']);
autumn = load([Directories.resultsDir  'autumn_output_' modTag '.mat']);

%%

% WOA DIN profiles

% plot trajectory means of WOA data extracted along the model trajectories
%  for model validation


% _________________________________________________________________________

% 4 x 2 tiled layout
% summer arctic model   | summer arctic WOA
% summer atlantic model | summer atlantic WOA
% autumn arctic model   | autumn arctic WOA
% autumn atlantic model | autumn atlantic WOA
% 

% load WOA data
WOAtable.summer = readtable('~/Documents/microARC model/Sat&WOA trajectories/monthlyTrajNSummer.csv'); % mu mol kg^-1
WOAtable.autumn = readtable('~/Documents/microARC model/Sat&WOA trajectories/monthlyTrajNAutumn.csv');

% Seawater Potential Density: 1025 kg/m^3  see https://www.nodc.noaa.gov/OC5/WOD/wod18-notes.html
WOAtable.summer.mean_N = WOAtable.summer.mean_N * 1025 / 1000; % convert mumol/kg to mmol/m3
WOAtable.autumn.mean_N = WOAtable.autumn.mean_N * 1025 / 1000;

% reassign model output
alldat.summer = summer;
alldat.autumn = autumn; 


seasons = [{'summer'}, {'autumn'}];
watermasses = [{'Arctic'}, {'Atlantic'}];


% clabelMod = 'DIN (mmol N m^{-3})';
% % clabelObs = 'DIN (µmol N kg^{-1})';
% clabelObs = 'DIN (mmol N m^{-3})';
% ul = 16;
% ll = 1;
% Ticks = logspace(log10(ll),log10(ul), 5); 
% Ticks = [0.3, Ticks]
% TickLabels = round(Ticks, 2, 'significant');
% ll = 0.3
% 
% t = tiledlayout(4,2)
% t.TileSpacing = 'compact';
% t.Padding = 'compact';
% title(t, 'DIN concentration, Model vs World Ocean Atlas')


% initialise result table
outtable = table();


for j = 1:length(seasons)

    season = seasons{j};
    
    % restructure WOAtable to array: [depth time traj]
    WOA = WOAtable.(season);
    Traj = unique(WOA.Traj);
    Month = unique(WOA.month);
    Depth = unique(WOA.depth);
   
    nTraj = length(Traj);
    nMonth = length(Month);
    nDepth = length(Depth);

    data_dates = [min(alldat.(season).Data.scalar.Date) max(alldat.(season).Data.scalar.Date)];
    
    WOA = reshape(WOA.mean_N, [nDepth nMonth nTraj]);
    
    % only keep until model depth
    maxDepth_mod = alldat.(season).FixedParams.Htot; 
    md_index = find(Depth == maxDepth_mod);
    
    Depth = Depth(1:md_index);
    WOA = WOA(1:md_index, :, :);

    trajwm_tab = unique(WOAtable.(season)(:,1:2));
    
    for i = 1:length(watermasses)
        wm = watermasses{i}; 
        wmindex = strcmp(alldat.(season).Forc.waterMass, wm);
        
        % pick Model and WOA Output 

        % WOA
        % WOA includes also the excluded ARC trajs, therefore wm index has
        % wrong dimensions. get wmindex_woa

        wmindex_woa = contains(trajwm_tab.wmTraj, wm);
        WOA_ = WOA(:,:,wmindex_woa);

        % model
        mod = squeeze(alldat.(season).out.N(:, :, :, wmindex));

        % interpolate mod to WOA depths
        depthsMOD = abs(alldat.(season).FixedParams.z);
        modINT = interpolate_mod(mod, depthsMOD, Depth);
        
        % get monthly means of mod 
        mod_months = month(datetime(alldat.(season).Forc.t(:,1),'ConvertFrom','datenum'));
        
        mod_monthav = nan([size(WOA_,1:2), size(modINT,3)]);
        for m = 1:length(unique(mod_months)) 
            mo_index = find(mod_months == m)
            mod_monthav(:,m,:) = mean(modINT(:,mo_index, :), 2, 'omitnan')
        end
        mod = mod_monthav;
        clear mod_monthav, mod_months

        % average both mod and WOA over trajectories
        mod = squeeze(mean(mod, 3, "omitnan"));
        WOA_ = squeeze(mean(WOA_, 3, "omitnan"));


        % calc RMSE
        RMSE = sqrt(mean((mod - WOA_).^2, "all", "omitnan"));

        mean_woa = mean(WOA_, "all");
        CV = RMSE / mean_woa;

        % assign to outtable
        out.row_label = {[season, '_', wm]};
        out.RMSE = RMSE;
        out.CV = CV; 

        out = struct2table(out);
        
        outtable = [outtable; out];
        clear out RMSE CV WOA_

    end

end

%%
% save 
writetable(outtable, [Directories.plotDir 'modelWOAFit.csv'])


close all
clearvars -except autumn summer Directories modTag pSetName