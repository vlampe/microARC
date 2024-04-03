%% create table with integrated production and export 

close all
clear all
clearvars -except autumn summer Directories modTag pSetName



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


%% 

% get poductivity

watermasses = [{'Arctic'}, {'Atlantic'}]
seasons = [{'summer'}, {'autumn'}]


alldat.summer = summer;
alldat.autumn = autumn; 

depths = [5, 8, 9];

for d = 1:length(depths)
    maxdepth = depths(d);
    for j = 1:length(seasons)
    
        season = seasons{j};
        
        z = alldat.(season).FixedParams.z;
        zEdges = alldat.(season).FixedParams.zw;
        zWidth = alldat.(season).FixedParams.zwidth; 
        
        Phy_index = alldat.(season).FixedParams.phytoplankton;
        C_index = strcmp(alldat.(season).FixedParams.PP_nut, 'C');
        N_index = strcmp(alldat.(season).FixedParams.PP_nut, 'N');
        
        for i = 1:length(watermasses)
            wm = watermasses{i}; 
    
            % pick data
            wm_index = strcmp(alldat.(season).Forc.waterMass, wm);
              
            % all phy C uptake 
            modC = alldat.(season).auxVars.uptake(Phy_index, :, C_index, :, wm_index);
            % sum over size class (1st dim)
            modC = squeeze(sum(modC, 1, 'omitnan')); % mmol C /m3 / d
            % integrate over depth until maxdepth
            modC = modC(1:maxdepth, :, :) .* zWidth(1:maxdepth);
            modC = squeeze(sum(modC, 1, 'omitnan'));
    
            % all phy C uptake 
            modN = alldat.(season).auxVars.uptake(Phy_index, :, N_index, :, wm_index);
            % sum over size class (1st dim)
            modN = squeeze(sum(modN, 1, 'omitnan')); % mmol C /m3 / d
            % integrate over depth until maxdepth
            modN = modN(1:maxdepth, :, :) .* zWidth(1:maxdepth);
            modN = squeeze(sum(modN, 1, 'omitnan'));
    
            % export for further analyses
            subsetC = ['C_' season '_' wm '_' num2str(abs(round(alldat.summer.FixedParams.zw(maxdepth+1)))) 'm'];
            subsetN = ['N_' season '_' wm '_' num2str(abs(round(alldat.summer.FixedParams.zw(maxdepth+1)))) 'm'];
            productivity.(subsetC) = modC; 
            productivity.(subsetN) = modN; 
        end
    end
end

clearvars -except autumn summer Directories modTag pSetName export productivity
%% get export

watermasses = [{'Arctic'}, {'Atlantic'}]
seasons = [{'summer'}, {'autumn'}]


alldat.summer = summer;
alldat.autumn = autumn; 

depths = [5, 8, 9];

for d = 1:length(depths)
    maxdepth = depths(d);

    for j = 1:length(seasons)
    
        season = seasons{j};
        
        C_index = strcmp(alldat.(season).FixedParams.PP_nut, 'C'); % ist eh immer 1
        N_index = strcmp(alldat.(season).FixedParams.PP_nut, 'N'); % ist eh immer 1
        zwidth = alldat.(season).FixedParams.zwidth(maxdepth); 
        
        for i = 1:length(watermasses)
            wm = watermasses{i}; 
    
            % pick data for C
            wm_index = strcmp(alldat.(season).Forc.waterMass, wm);
            
            %Phy
            concPhy = alldat.(season).out.P(:, maxdepth, C_index, :, wm_index);  
            concPhy = squeeze(sum(concPhy, 1, 'omitnan')); % sum over size classes
            %Zoo
            concZoo = alldat.(season).out.Z(:, maxdepth, C_index, :, wm_index);  
            concZoo = squeeze(sum(concZoo, 1, 'omitnan')); % sum over size classes
            %Plankton
            concP = concPhy + concZoo; 
            %OM
            concOM = squeeze(alldat.(season).out.OM(:, maxdepth, C_index, :, wm_index)); 
            
            % extract sinking speeds from FixedParams
            sinksPlankton = alldat.(season).Params.wp(maxdepth, 1); % m/d: careful, if wp is not 0 or different for each size class, then per-sizeclass export has to be calculated before adding 
            sinksOM = [alldat.(season).Params.wDOM1; alldat.(season).Params.wPOM1];
    
            exOM = concOM .* sinksOM;
            exOM = squeeze(sum(exOM, 1)); % add POM and DOM export (DOM is 0 anyways)
    
            exP = concP .* sinksPlankton; 
    
             % get changes in conc due to diffusion
            diffP = alldat.(season).auxVars.B_diffuse(:, maxdepth, C_index, :, wm_index);  % all plankton
            diffP = squeeze(sum(diffP, 1, 'omitnan')) .* zwidth;  % sum over size classes (Phy+Zoo) and muliply with depth layer width
    
            diffOM = alldat.(season).auxVars.OM_diffuse(:, maxdepth, C_index, :, wm_index); 
            diffOM = squeeze(sum(diffOM, 1, 'omitnan')) .* zwidth; % sum over types (DOM+POM) and multiply with depth layer width
    
            extot = exOM + exP; % total export of C in mmol m^-2 d^-1 % for all trajs
            % extot = exOM' + exP; % for single trajs
            extot2 = extot + diffOM + diffP; % add changes due to diffusion to changes due to sinking
    
    
            subsetC = ['C_' season '_' wm '_' num2str(abs(round(alldat.summer.FixedParams.zw(maxdepth+1)))) 'm'];
            export.(subsetC) = extot2; 

            % ___
            % pick data for N
            wm_index = strcmp(alldat.(season).Forc.waterMass, wm);
            
            %Phy
            concPhy = alldat.(season).out.P(:, maxdepth, N_index, :, wm_index);  
            concPhy = squeeze(sum(concPhy, 1, 'omitnan')); % sum over size classes
            %Zoo
            concZoo = alldat.(season).out.Z(:, maxdepth, N_index, :, wm_index);  
            concZoo = squeeze(sum(concZoo, 1, 'omitnan')); % sum over size classes
            %Plankton
            concP = concPhy + concZoo; 
            %OM
            concOM = squeeze(alldat.(season).out.OM(:, maxdepth, N_index, :, wm_index)); 
            
            % extract sinking speeds from FixedParams
            sinksPlankton = alldat.(season).Params.wp(maxdepth, 1); % m/d: careful, if wp is not 0 or different for each size class, then per-sizeclass export has to be calculated before adding 
            sinksOM = [alldat.(season).Params.wDOM1; alldat.(season).Params.wPOM1];
    
            exOM = concOM .* sinksOM;
            exOM = squeeze(sum(exOM, 1)); % add POM and DOM export (DOM is 0 anyways)
    
            exP = concP .* sinksPlankton; 
    
             % get changes in conc due to diffusion
            diffP = alldat.(season).auxVars.B_diffuse(:, maxdepth, N_index, :, wm_index);  % all plankton
            diffP = squeeze(sum(diffP, 1, 'omitnan')) .* zwidth;  % sum over size classes (Phy+Zoo) and muliply with depth layer width
    
            diffOM = alldat.(season).auxVars.OM_diffuse(:, maxdepth, N_index, :, wm_index); 
            diffOM = squeeze(sum(diffOM, 1, 'omitnan')) .* zwidth; % sum over types (DOM+POM) and multiply with depth layer width
    
            extot = exOM + exP; % total export of C in mmol m^-2 d^-1 % for all trajs
            % extot = exOM' + exP; % for single trajs
            extot2 = extot + diffOM + diffP; % add changes due to diffusion to changes due to sinking
    
    
            subsetN = ['N_' season '_' wm '_' num2str(abs(round(alldat.summer.FixedParams.zw(maxdepth+1)))) 'm'];
            export.(subsetN) = extot2; 

        end
    end
end

clearvars -except autumn summer Directories modTag pSetName export productivity

%%
% integrate both

% integrate production over time

% stored in struct productivity
% unit is mmol N / m2 / d


exnames = fieldnames(productivity); 
for i = 1:length(exnames)
    prod = productivity.(exnames{i});
    ex = export.(exnames{i}); 

    % pe = ex ./ prod;
    
    % whole model period (year)
    prod_int = sum(prod, 1); 
    prod_int_avg = mean(prod_int, "omitnan"); 
    prod_int_sd = std(prod_int, "omitnan"); 

    prod_int_year.(exnames{i}) = prod_int_avg;
    prod_int_year_sd.(exnames{i}) = prod_int_sd;

    ex_int = sum(ex, 1, "omitnan");
    ex_int_avg = mean(ex_int, "omitnan");
    ex_int_sd = std(ex_int, "omitnan");

    ex_int_year.(exnames{i}) = ex_int_avg;
    ex_int_year_sd.(exnames{i}) = ex_int_sd;
   
    % get pe from time-integrated values, then average and sd them
    pe_traj = ex_int ./ prod_int; 
    pe_year = mean(pe_traj);
    pe_sd_year = std(pe_traj);

    % average pe over time and trajectories
    pe_av_year.(exnames{i}) = pe_year;
    pe_year_sd.(exnames{i}) = pe_sd_year;
     
    
    % pe_av_year.(exnames{i}) = mean(pe, 'all'); % leads to very weird numbers, since pe is veeeeery high outside of the productive period (what we saw earlier in the pe time series)
    % pe_year_sd.(exnames{i}) = std(pe, [], 'all'); % its better to calculate the pe from integrated prod and ex (see Serra-Pompei) 
    % 

%    clear prod_int prod_int_avg prod_int_sd ex_int ex_int_avg ex_int_sd 
    
    % Feb - Aug
    if contains(exnames{i}, 'summer')
        datetimes = datetime(summer.Forc.t(:,1),'ConvertFrom','datenum');
    elseif contains(exnames{i}, 'autumn')
        datetimes = datetime(autumn.Forc.t(:,1),'ConvertFrom','datenum');
    end
    selectMonths = any(month(datetimes) == 4:9, 2);

    prod_int = sum(prod(selectMonths,:), 1); 
    prod_int_avg = mean(prod_int);

    prod_int_AprSep.(exnames{i}) = prod_int_avg;

    ex_int = sum(ex(selectMonths,:), 1); 
    ex_int_avg = mean(ex_int);

    ex_int_AprSep.(exnames{i}) = ex_int_avg;
    clear prod_int prod_int_avg selectMonths datetimes ex_int ex_int_avg
end


% all units still in mmol m2 per time
clearvars -except autumn summer Directories modTag pSetName export ...
    productivity prod_int_AprSep ex_int_AprSep ...
    ex_int_year ex_int_year_sd prod_int_year prod_int_year_sd pe_av_year pe_year_sd

%% construct table, per 6 month productive season

t1 = struct2table(prod_int_AprSep, 'RowNames', {'production'})
t1 = rows2vars(t1)

t2 = struct2table(ex_int_AprSep, 'RowNames', {'export'})
t2 = rows2vars(t2)

% join prod and ex tables
t = join(t1, t2, "Keys","OriginalVariableNames")
t.split = split(t.OriginalVariableNames, '_')
t.var = t.split(:,1)
t.depth = t.split(:,4)
t.subset = strcat(t.split(:,2), '_', t.split(:,3))
t = removevars(t, "split")

t = removevars(t, "OriginalVariableNames")
t = unstack(t, ["production", "export"], ["var"])
t = movevars(t, ["subset", "depth", "production_C", "production_N", "export_C", "export_N" ])

t.pe_C = t.export_C ./ t.production_C
t.pe_N = t.export_N ./ t.production_N


% sort table
desiredOrder = {'summer_Arctic', 'summer_Atlantic', 'autumn_Arctic', 'autumn_Atlantic'}
t.subset = categorical(t.subset, desiredOrder, 'ordinal', true)
t = sortrows(t, "subset")

% convert to g C /m2 / 6mo and g N /m2 / 6mo 
% convert to g C /m2/d
t.production_C = t.production_C .* (12.011 * 1e-3)
t.export_C = t.export_C .* (12.011 * 1e-3)
t.production_N = t.production_N .* (14.007 * 1e-3)
t.export_N = t.export_N .* (14.007 * 1e-3)


% round and convert to stings
t{:, 3:end} = round(t{:, 3:end},2)
t{:, 3:end} = string(t{:, 3:end})

% export table
writetable(t, [Directories.plotDir 'prodution_export.csv'])



%% construct table, per 12 months (year)

t1 = struct2table(prod_int_year, 'RowNames', {'production'})
t1 = rows2vars(t1)

t11 = struct2table(prod_int_year_sd, 'RowNames', {'production_sd'})
t11 = rows2vars(t11)

t2 = struct2table(ex_int_year, 'RowNames', {'export'})
t2 = rows2vars(t2)

t22 = struct2table(ex_int_year_sd, 'RowNames', {'export_sd'})
t22 = rows2vars(t22)

t3 = struct2table(pe_av_year, 'RowNames', {'pe'})
t3 = rows2vars(t3)

t33 = struct2table(pe_year_sd, 'RowNames',{'pe_sd'})
t33 = rows2vars(t33)

% join prod and ex tables
t = join(t1, t11, "Keys","OriginalVariableNames")
t = join(t, t2, "Keys","OriginalVariableNames")
t = join(t, t22, "Keys","OriginalVariableNames")
t = join(t, t3, "Keys","OriginalVariableNames")
t = join(t, t33, "Keys","OriginalVariableNames")

t.split = split(t.OriginalVariableNames, '_')
t.var = t.split(:,1)
t.depth = t.split(:,4)
t.subset = strcat(t.split(:,2), '_', t.split(:,3))
t = removevars(t, "split")

t = removevars(t, "OriginalVariableNames")
t = unstack(t, ["production", "production_sd", "export", "export_sd", "pe", "pe_sd"], ["var"])
t = movevars(t, ["subset", "depth", "production_C", "production_sd_C", "production_N", ...
    "production_sd_N", "export_C", "export_sd_C" "export_N", "export_sd_N", ...
    "pe_C", "pe_sd_C", "pe_N", "pe_sd_N"])

% t.pe_C = t.export_C ./ t.production_C
% t.pe_N = t.export_N ./ t.production_N


% % sort table
% desiredOrder = {'summer_Arctic', 'summer_Atlantic', 'autumn_Arctic', 'autumn_Atlantic'}
% t.subset = categorical(t.subset, desiredOrder, 'ordinal', true)
% t = sortrows(t, "subset")

% convert to g C /m2 / 6mo and g N /m2 / 6mo 
% convert to g C /m2/d
t.production_C = t.production_C .* (12.011 * 1e-3)
t.production_sd_C = t.production_sd_C .* (12.011 * 1e-3)
t.export_C = t.export_C .* (12.011 * 1e-3)
t.export_sd_C = t.export_sd_C .* (12.011 * 1e-3)
t.production_N = t.production_N .* (14.007 * 1e-3)
t.production_sd_N = t.production_sd_N .* (14.007 * 1e-3)
t.export_N = t.export_N .* (14.007 * 1e-3)
t.export_sd_N = t.export_sd_N .* (14.007 * 1e-3)

% round and convert to stings
t{:, 3:end} = round(t{:, 3:end},2)
t{:, 3:end} = string(t{:, 3:end})

% export table
writetable(t, [Directories.plotDir 'prodution_export_year.csv'])



%% get average production and export during sampling times 

fn = fields(export)


for f=1:length(fn)
    f_  = fn{f};

    if contains(f_, "summer")
        data_dates = [min(summer.Data.scalar.Date) max(summer.Data.scalar.Date)];
    else
        data_dates = [min(autumn.Data.scalar.Date) max(autumn.Data.scalar.Date)];
    end

    % dates to indices
    data_index = day(data_dates, "dayofyear"); 

    % extract prod and ex
    prod = productivity.(f_)(data_index(1):data_index(2), :);
    ex = export.(f_)(data_index(1):data_index(2), :);

    % % convert to mg m^-3 d^-1
    % if contains(f_, "C_")
    %     prod = prod .* (12.011 * 1e-3);
    %     ex = ex .* (12.011 * 1e-3);
    % else
    %     prod = prod .* (14.007 * 1e-3);
    %     ex = ex .* (14.007 * 1e-3);
    % end

    disp([f_ ': mean prod: ' num2str(mean(prod, "all", "omitnan")) ' sd prod: ' num2str(std(prod, [], "all", "omitnan")) ...
        ': mean ex: ' num2str(mean(ex, "all", "omitnan")) ' sd ex: ' num2str(std(ex, [], "all", "omitnan"))])
end


%% get distances to sampling region

close all
clearvars -except autumn summer Directories modTag pSetName

alldat.summer = summer;
alldat.autumn = autumn; 

seasons = {"summer", "autumn"}
watermasses = [{'Arctic'}, {'Atlantic'}];

for s=1:length(seasons)
    season = seasons{s}

    if strcmp(season, "summer")
        data_dates = [min(summer.Data.scalar.Date) max(summer.Data.scalar.Date)];
    else
        data_dates = [min(autumn.Data.scalar.Date) max(autumn.Data.scalar.Date)];
    end


       for w=1:length(watermasses)
    
        watermass = watermasses{w}
        wm_filter = strcmp(alldat.(season).Forc.waterMass, watermass)


        dist_days = days(datetime(alldat.(season).Forc.t(:, wm_filter), 'ConvertFrom','datenum') - data_dates(1))
    
        % central point in HG
        HG_center_lat = 79
        HG_center_lon = 5

        % distance to central hausgarten
        dist_deg = distance(alldat.(season).Forc.y(:, wm_filter), alldat.(season).Forc.x(:, wm_filter), HG_center_lat, HG_center_lon)
         
        % % distance to coordinates of this traj on sampling day 1 (may not be
        % % in HG)   WORKS ONLY FOR ATLANTIC???? 
        % date_filter = datetime(alldat.(season).Forc.t(:, wm_filter), 'ConvertFrom','datenum') == data_dates(1);
        % dist_deg = distance(alldat.(season).Forc.y(:, wm_filter), ...
        %     alldat.(season).Forc.x(:, wm_filter), ...
        %     repmat(alldat.(season).Forc.y(date_filter)' ,size(dist_days,1), 1), ...
        %     repmat(alldat.(season).Forc.x(date_filter)' ,size(dist_days,1), 1));
        
        
        % repmat(alldat.(season).Forc.y(date_filter)' ,size(dist_days,1), 1)
        
        dist_km = deg2km(dist_deg)
        km_max = max(dist_km,[], 'all')

 
        
        figure
        patch([data_dates fliplr(data_dates)], [0*[1 1] (km_max+km_max*0.1)*[1 1]], 'k', 'facecolor', [0.8 0.8 0.8], ...
            'edgecolor',[0.8 0.8 0.8],  'facealpha',0.5,'edgealpha',0.5)
        hold on
        for t=1:size(dist_days, 2)
            scatter(datetime(alldat.(season).Forc.t(:,t), 'ConvertFrom','datenum'), dist_km(:,t),[], dist_days(:,t), 'filled')
        end
        ylabel("distance to HG center (km)")
        yline(25)   
        title([season " - " watermass])
        cb = colorbar()
        cb.Label.String = "distance to first sampling date (days)"
       end
end


wm_filter = strcmp(alldat.(season).Forc.waterMass, "Arctic")
dist_deg = distance(alldat.(season).Forc.y(1,wm_filter), alldat.(season).Forc.x(1,wm_filter), alldat.(season).Forc.y(end,wm_filter), alldat.(season).Forc.x(end,wm_filter))
dist_km = deg2km(dist_deg)
max(dist_km, [], "all")
