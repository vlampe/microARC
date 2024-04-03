 %% load model output 
close all
clear all


modTag = 'optimisedParams_AH_wPOM_10_pmax_b_-0,08_Qmin_QC_a_0.035_theta_4.2_Gmax_a_11_m2_5e-2_atLargeSC_rDOM_0.04_icelessArcSummer_adjustedVin_fullRunForPaper';

%% add paths 
addpath(genpath('~/Documents/microARC model')) % location of my work with the model
addpath(genpath('~/microARC/')) % location of model


% set directories
Directories.resultsDir  = ['~/Documents/microARC model/2nd paper/mod output/' modTag '/'];
Directories.plotDir  = ['~/Documents/microARC model/2nd paper/mod output/' modTag '/plots/publicationPlots/'];


% load model output
summer = load([Directories.resultsDir  'summer_output_' modTag '.mat']);
% structure that contains model input and output for summer set-up: 
% Data: observational data used for model optimisation
% FixedParams: constants, info about depth layers, size classes, useful indeces...
% Forc: Forcing data input as time series along model trajectories
% Params: variable model parameters, equal in summer and autumn set-up
% auxVars: more model output, derived quantities or rates... (eg cell densities)
% modData: model output in the same format as Data, for convinient comparisions (but there's probably an error in the size/sizeFull substruct)
% out: model output for N (DIN), P (Phy biomass), Z (zoo biomass), OM (Detritus)
% v0: initial model conditions

% common format:
% size(out.P) = [nSizeClass nDepth nVar nDays nTrajs]  with Var = {C, N, Chl}
% size(out.OM) = [nOMType nDepth nVar nDays nTrajs] with OMType = {DOM, POM}

autumn = load([Directories.resultsDir  'autumn_output_' modTag '.mat']);

%% RSMEs and CVs of DIN, Chl, POC, PON

% splitting by watermass is not really useful in autumn, as there are no
% stations in obs data that can be fully considered artic.


% this could be possibly used if water mass splitting is nessescary 
% compare "my" assignment of regime (warm/cold) to aidans classification 
t = table(autumn.Data.sizeFull.Label, autumn.Data.sizeFull.regime, autumn.Data.sizeFull.Event);
unique(t)

t = table(summer.Data.sizeFull.Label, summer.Data.sizeFull.regime, summer.Data.sizeFull.Event, 'VariableNames', {'Label', 'regime', 'Event'});
t = unique(t)
% is pretty good
% in summer, stations classified as mixed by aidan are "cold"
% in autumn, aidan finds no arctic, but I id'ed event 18/MSM77-52-01 as
% actic

%%
clearvars -except autumn summer Directories modTag pSetName
close all
%%

% PON, POC
plotVars = [{'DIN'}, {'Chl a'}, {'POC'}, {'PON'}];

seasons = [{'summer'}, {'autumn'}];
alldat.summer = summer;
alldat.autumn = autumn;
watermasses = [{'Arctic'}, {'Atlantic'}];

timeFrame = 0 ; % include simulations +- 7 days around sampling of observations

% initialise result table
outtable = table();


for s = 1:length(seasons)
    season = seasons{s};
    
    for wm =1:length(watermasses)

        watermass = watermasses{wm};
        
        wm_mod_index = strcmp(alldat.(season).Forc.waterMass, watermass);
        switch watermass
            case 'Arctic' 
                regime = 'cold';
            case 'Atlantic'
                regime = 'warm';
        end
        
        % get look-up table to find identifiers (Labels) for this
        % watermass in obs data
        data_lookup = table(alldat.(season).Data.sizeFull.Label, ...
            alldat.(season).Data.sizeFull.regime, ...
            alldat.(season).Data.sizeFull.Event, ...
            'VariableNames', {'Label', 'regime', 'Event'});
        data_lookup = unique(data_lookup);
        data_lookup.waterMass = alldat.(season).Data.sizeFull.waterMass; 

        data_labels_wm = data_lookup.Label(strcmp(data_lookup.regime, regime));
        wm_dat_index = contains(alldat.(season).Data.scalar.Label, data_labels_wm); 

        % extract all data (no subsetting to wm and var yet)
        dat = alldat.(season).Data.scalar.Value;
    
        depth = alldat.(season).Data.scalar.Depth;
    
        % pick model equivalents with a wider time frame
        samplingDates = unique(alldat.(season).Data.scalar.Yearday); 
        timeindexMOD = samplingDates(1)-timeFrame:samplingDates(end)+timeFrame; 
        
        for pv = 1:length(plotVars)
    
    
            plotVar = plotVars{pv};
            switch plotVar
                case 'DIN'
                    varindexDAT = strcmp(alldat.(season).Data.scalar.Variable, 'N');
                    elemVar = 'N';
                    ylab = 'DIN concentration (mmol m^{-3})';
    
                case 'Chl a'
                    varindexDAT = strcmp(alldat.(season).Data.scalar.Variable, 'chl_a');
                    elemVar = 'Chl';
                    ylab = 'Chl concentration (mg m^{-3})';
    
                case 'POC'
                    varindexDAT = strcmp(alldat.(season).Data.scalar.Variable, 'POC');
                    elemVar = 'C';
                    ylab = 'POC concentration (mmol m^{-3})';
                case 'PON'
                    varindexDAT = strcmp(alldat.(season).Data.scalar.Variable, 'PON');
                    ylab = 'PON concentration (mmol m^{-3})';
                    elemVar = 'N';
            end
    
    
            % pick obs data to plot (according to season, watermass, variable)
            plotdat = dat(varindexDAT&wm_dat_index);
            depthDAT = depth(varindexDAT&wm_dat_index);
    
            data_tab = table(plotdat, depthDAT, 'VariableNames',{'data', 'depth'});
            data_tab = groupsummary(data_tab, 'depth', 'median'); % get median per depth
            
            % pick model output 
            if strcmp(plotVar, 'POC')|strcmp(plotVar, 'PON') 
                 
                % detOM
                POM_index = alldat.(season).FixedParams.POM_index;
                pomvarindexMOD = strcmp(alldat.(season).FixedParams.OM_nut, elemVar);
                detOM = squeeze(alldat.(season).out.OM(POM_index, :, pomvarindexMOD, timeindexMOD, wm_mod_index));
                %size(detOM)
        
                % phyOM
                phyvarindexMOD = strcmp(alldat.(season).FixedParams.PP_nut, elemVar); 
                phyOM = alldat.(season).out.P(:,:,phyvarindexMOD, timeindexMOD, wm_mod_index); 
                phyOM = squeeze(sum(phyOM, 1)); % sum over size classes
                %size(phyOM)
        
                % zooOM
                zoovarindexMOD = strcmp(alldat.(season).FixedParams.ZP_nut, elemVar); 
                zooOM = alldat.(season).out.Z(:,:,zoovarindexMOD, timeindexMOD, wm_mod_index); 
                zooOM = squeeze(sum(zooOM, 1)); % sum over size classes
                %size(zooOM)
        
                % add all OM sources
                plotmod = detOM + phyOM + zooOM;
                %size(plotmod)
            elseif strcmp(plotVar, 'DIN')
                plotmod = squeeze(alldat.(season).out.N(:,:, timeindexMOD, wm_mod_index)); 
            elseif strcmp(plotVar, 'Chl a')
                phyvarindexMOD = strcmp(alldat.(season).FixedParams.PP_nut, 'Chl'); 
                plotmod = alldat.(season).out.P(:,:,phyvarindexMOD, timeindexMOD, wm_mod_index); 
                plotmod = squeeze(sum(plotmod,1)); % sum over size class
            end
    
            % repair: if plotdat is empty: skip to next iteration (plotting
            % of only mod results not possible here, no depthDAt vector
            % available for interpolation
            if isempty(plotdat)
                
                RMSE = NaN;
                CV = NaN;
    
            else

                % interpolate mod to obs depths
                depthsMOD = abs(alldat.(season).FixedParams.z);
                udepthDAT = unique(depthDAT);
                plotmodINT = NaN([length(udepthDAT), size(plotmod,2), size(plotmod,3)]);
        
                plotmodINT = interpolate_mod(plotmod, depthsMOD, udepthDAT);
        

                % get RMSE
                median_obs_i = data_tab.median_data;
                RMSE = sqrt(mean((plotmodINT - median_obs_i).^2,"all", "omitnan")); % has same unit as original data

                % get CV
                CV = RMSE / mean(data_tab.median_data); % unitless
                
            end

            % assign to outtable
            out.row_label = {[season, '_', watermass]};
            out.var = {plotVar};
            out.RMSE = RMSE;
            out.CV = CV; 

            out = struct2table(out);

            outtable = [outtable; out];
            clear out RMSE CV
        end % end wm loop
    end % end var loop
end % end season loop

clearvars -except autumn summer Directories modTag pSetName outtable

%%
% then spread table, 4x(1+2x4) 
outtable = unstack(outtable, ["RMSE", "CV"], ["var"])


%% calculate IQD for biovolume
% for phy nd zoo
% mean over trajs
% compare trajs with s-spectra

% size spectra
% abundance and biovolume
concTypes = [{'abundance'}, {'biovolume'}]; % 
alldat.summer = summer; 
alldat.autumn = autumn;

seasons = [{'summer'}, {'autumn'}];
waterMasses = [{'Arctic'}, {'Atlantic'}];

iDepth = [1:5]; 





IQD_table = table()

for ct = 1:length(concTypes)
    concType = concTypes{ct};
    
    switch concType
        case 'abundance'
            modVar = 'cellDensity';
            obsVar = 'cellDensity';
            miny = 1e3; 
            totunit = '(cells / m^3)';
            yunit = '(cells / m^3 log_{10}(ESD / 1 µm)^{-1})'; 
            ytext = 'cell conc. density';
        case 'biovolume'
            modVar = 'biovolume';
            obsVar = 'BioVolDensity';
            miny = 1e6; 
            totunit = '(µm^3 / m^3)';
            yunit = '(µm^3 / m^3 log_{10}(ESD / 1 µm)^{-1})';
            ytext = 'biovol. conc. density'; 
    end
    
    for s=1:length(seasons)
        season = seasons{s}; 

        switch season
            case 'summer'
               % iTime = [197:208];
                iTime = [190:215];
            case 'autumn'
               % iTime = [259:277]; 
                iTime = [252:284]; 
        end

        iPhyto = alldat.(season).FixedParams.phytoplankton;
        iZoo  = alldat.(season).FixedParams.zooplankton;
        zwidth = alldat.(season).FixedParams.zwidth;

        for w=1:length(waterMasses)
            wm = waterMasses{w}; 
            iTraj = strcmp(alldat.(season).Forc.waterMass, wm);

            % get observation data for comparison
            if strcmp(season, 'summer')&strcmp(wm, 'Atlantic')
                scenario = 'S1';
            elseif strcmp(season, 'summer')&strcmp(wm, 'Arctic')
                scenario = 'S2';
            elseif strcmp(season, 'autumn')&strcmp(wm, 'Atlantic')
                scenario = 'S3';
            elseif strcmp(season, 'autumn')&strcmp(wm, 'Arctic')
                scenario = 'S4';
            end
            
            % pick oberservational data (already biovolume densities, (µm^3 / m^3 log_{10}(ESD / 1 µm)^{-1}))
            scenarioindex = strcmp(alldat.(season).Data.size.scenario, scenario);
            obsESDs = alldat.(season).Data.size.ESD(scenarioindex & strcmp(alldat.(season).Data.size.trophicLevel, 'heterotroph')); 
            dlog10ESDobs = diff(log10(obsESDs(1:2)));

            obsPhy = alldat.(season).Data.size.(obsVar)(scenarioindex & strcmp(alldat.(season).Data.size.trophicLevel, 'autotroph'));
            % obsPhySE = alldat.(season).Data.size.(strcat(obsVar, 'SE'))(scenarioindex & strcmp(alldat.(season).Data.size.trophicLevel, 'autotroph'));
            obsPhyConc = sum(obsPhy .* dlog10ESDobs); 
            obsPhy_norm = obsPhy ./ obsPhyConc; % normalise

            obsZoo = alldat.(season).Data.size.(obsVar)(scenarioindex & strcmp(alldat.(season).Data.size.trophicLevel, 'heterotroph'));
            % obsZooSE = alldat.(season).Data.size.(strcat(obsVar, 'SE'))(scenarioindex & strcmp(alldat.(season).Data.size.trophicLevel, 'heterotroph'));
            obsZooConc = sum(obsZoo .* dlog10ESDobs); 
            obsZoo_norm = obsZoo ./ obsZooConc; % normalise



            % pick model output

            % phy
            modPhy = alldat.(season).auxVars.(modVar)(iPhyto, iDepth, iTime, iTraj); % abundance or biovol per size class
            % weighted average over depth
            modPhy = squeeze(sum(modPhy .* zwidth(iDepth)', 2)./ sum(zwidth(iDepth))); 
            % conversion to biovol density
            ESDs = alldat.(season).FixedParams.PPdia_intervals;
            dlog10ESDmod = diff(log10(ESDs(1:2))); % step size


            modPhy = modPhy ./ dlog10ESDmod; % conversion
            
            % interpolate to obs ESDs
            modPhy_interp = NaN([length(obsESDs), size(modPhy, 2:3)]);
            for i =1:size(modPhy,2)
                for j = 1:size(modPhy,3)

                    modPhy_rep = repelem(modPhy(:,i,j), 2);
        
                    modESDs_rep = [ESDs(1); repelem(ESDs(2:end-1),2); ESDs(end)];
                    log_modESDs_rep = log10(modESDs_rep);
                    log_modESDs_rep(2:2:end) = log_modESDs_rep(2:2:end) - dlog10ESDobs;
        
                    modPhy_interp(:,i,j) = interp1(log_modESDs_rep, modPhy_rep, log10(obsESDs));
                    
                end
            end

            % plot(obsESDs, modPhy_interp(:,:,1))
            % replace Nans with small non-zero number
            modPhy_interp = fillmissing(modPhy_interp, "constant", 1e-16);

            % normalise 
            modPhy_norm = modPhy_interp ./ sum(modPhy_interp .* dlog10ESDobs, 1);



            % repeat for zoo
            modZoo = alldat.(season).auxVars.(modVar)(iZoo, iDepth, iTime, iTraj); % abundance or biovol per size class
            % weighted average over depth
            modZoo = squeeze(sum(modZoo .* zwidth(iDepth)', 2)./ sum(zwidth(iDepth))); 
            % conversion to biovol density
            ESDs = alldat.(season).FixedParams.ZPdia_intervals;
            dlog10ESDmod = diff(log10(ESDs(1:2))); % step size


            modZoo = modZoo ./ dlog10ESDmod; % conversion
            
            % interpolate to obs ESDs
            modZoo_interp = NaN([length(obsESDs), size(modZoo, 2:3)]);
            for i =1:size(modZoo,2)
                for j = 1:size(modZoo,3)

                    modZoo_rep = repelem(modZoo(:,i,j), 2);
        
                    modESDs_rep = [ESDs(1); repelem(ESDs(2:end-1),2); ESDs(end)];
                    log_modESDs_rep = log10(modESDs_rep);
                    log_modESDs_rep(2:2:end) = log_modESDs_rep(2:2:end) - dlog10ESDobs;
        
                    modZoo_interp(:,i,j) = interp1(log_modESDs_rep, modZoo_rep, log10(obsESDs));
                    
                end
            end

            % plot(obsESDs, modPhy_interp(:,:,1))
            % replace Nans with small non-zero number
            modZoo_interp = fillmissing(modZoo_interp, "constant", 1e-16);

            % normalise 
            modZoo_norm = modZoo_interp ./ sum(modZoo_interp .* dlog10ESDobs, 1); % sepctrum / total concentration; tot conc = integral of spectrum

            
            % get IQDs, one phy one zoo
            IQDs.phy = NaN([size(modPhy_norm, 2:3)]);
            IQDs.zoo = NaN([size(modZoo_norm, 2:3)]);
            for i = 1:size(modPhy_norm, 2)
                for j = 1:size(modPhy_norm, 3)
                    
                    IQDs.phy(i, j) = sum( (cumsum(obsPhy_norm) - cumsum(modPhy_norm(:,i,j))).^2)*dlog10ESDobs;
                    IQDs.zoo(i, j) = sum( (cumsum(obsZoo_norm) - cumsum(modZoo_norm(:,i,j))).^2)*dlog10ESDobs;
                    
                   %  plot(1:450, cumsum(obsPhy_norm)*dlog10ESDobs)

                end
            end

            % get mean of IQDs and save
            out.row_label = {[season, '_', wm]};
            out.var = {concType};
            out.IQD_phy = mean(IQDs.phy, 'all');
            out.IQD_zoo = mean(IQDs.zoo, 'all');

            

            IQD_table = [IQD_table; struct2table(out)];

            % plot(obsESDs, modZoo_norm(:,i,j))
            % hold on 
            % plot(obsESDs, obsZoo_norm)
            % hold off
            % 
            % figure
            % plot(obsESDs, cumsum(modZoo_norm(:,i,j))*dlog10ESDobs)
            % hold on 
            % plot(obsESDs, cumsum(obsZoo_norm)*dlog10ESDobs)
            % hold off


            
        end
    end
end

clearvars -except autumn summer Directories modTag pSetName outtable IQD_table

% spread table
IQD_table = unstack(IQD_table, ["IQD_phy", "IQD_zoo"], "var")

%num_to_exp_string(IQD_table{:, 2:end},3)


%%


% % convert IQD columns to scientific notation and strings
% IQD_str_table = varfun(@(x) num_to_exp_string(x, 3), IQD_table, 'InputVariables',{'IQD_phy_abundance', 'IQD_phy_biovolume', 'IQD_zoo_abundance', 'IQD_zoo_biovolume'}) % apply num_to_exp_str (self made function) to selected columns
% IQD_str_table.row_label = IQD_table.row_label;
% IQD_str_table.Properties.VariableNames = regexprep(IQD_str_table.Properties.VariableNames, 'Fun_', '');
%% join and format results table

% then join with IQD table


results = join(outtable, IQD_table, "Keys","row_label") % or IQD_str_table

% make pretty (round, string, e-4 notation...
% filter for RMSE and CV columns
results{:,contains(results.Properties.VariableNames, "RMSE"|"CV")} = ...
    round(results{:,contains(results.Properties.VariableNames, "RMSE"|"CV")},4)

% sort 
order = [{'row_label'}, {'RMSE_DIN'}, {'CV_DIN'}, {'RMSE_ChlA'}, {'CV_ChlA'}, ...
    {'RMSE_POC'}, {'CV_POC'}, {'RMSE_PON'}, {'CV_PON'}, {'IQD_phy_biovolume'}, ...
    {'IQD_zoo_biovolume'},{'IQD_phy_abundance'}, {'IQD_zoo_abundance'}]

results = results(:,order)


% round and convert to stings
%fun = @(x) sprintf('%0.4f', x);
%results{:, 2:end} = cellfun(fun, num2cell(round(results{:, 2:end},4)) , 'UniformOutput',0);
results{:, 2:end} = round(results{:, 2:end},4)

% string(results{:,contains(results.Properties.VariableNames, "RMSE"|"CV")})
% save 
writetable(results, [Directories.plotDir 'modelDataFit.csv'])
clearvars -except autumn summer Directories modTag pSetName outtable IQD_table results


%% rearrage table (swap rows and cols)


results = readtable([Directories.plotDir 'modelDataFit.csv'])

rows2vars(results)
