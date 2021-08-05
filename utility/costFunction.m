function [cost, costComponents, costFunctionChoices] = costFunction(varargin)

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Labels for cost function options -- code for each is below.
% To use a different cost function simply write its label here then include
% another case below -- it must be fully code-able from Data and modData
costFunctionChoices = { ...
    'LSS', ...
    'RMS', ...
    'syntheticLikelihood_ScalarNormal_SizeSpectraLogNormal_logisticNormal', ...
    'syntheticLikelihood_ScalarNormal_SizeSpectraLogNormalDirichlet', ...
    'syntheticLikelihood_ScalarNormalShape_SizeSpectraLogNormalDirichlet', ...
    'N_LN-Dir_groupWaterOrigin', ...
    'Hellinger_groupWaterOrigin', ...
    'Hellinger2_groupWaterOrigin', ...
    'Hellinger3_groupWaterOrigin', ...
    'Hellinger_MVN_groupWaterOrigin', ...
    'IQD_Hellinger_groupWaterOrigin', ...
    'meanCDFdist_Hellinger', ...
    'smoothCDFdist_Hellinger', ...
    'LeastAbsErr_Hellinger', ...
    'RMSsmooth_Hellinger', ...
    'RMS_Hellinger', ...
    'RMS_Hellinger2', ...
    'meanCDFdist_HellingerFullSpectrum', ...
    'meanCDFdist_HellingerFullSpectrum_averagedEventsDepths'
    };

% Evaluating cost requires passing name-value pairs for 'label', 'Data' and
% 'modData' as optional varargin values. The 'label' specifies which cost
% function to use; 'Data' are measurements; 'modData' are the modelled
% equivalents.
extractVarargin(varargin)
calculateCost = exist('Data', 'var') & exist('modData', 'var') & exist('label', 'var');

if ~calculateCost
    
    cost = nan;
    costComponents = nan;
    
else
    
    if ~ismember(label, costFunctionChoices)
        error('Cost function name must match one of the available options in costFunction.m.')
    end
    
    % Data types used to fit the model
    Vars = Data.scalar.obsInCostFunction;
    VarsSize = Data.size.obsInCostFunction;

    switch label
        
        case 'LSS' % Least sum of squares
            % Scalar data
            for i = 1:length(Vars)
                varLabel = Vars{i};
                yobs = Data.scalar.scaled_Value(strcmp(Data.scalar.Variable, varLabel));
                ymod = modData.scalar.scaled_Value(strcmp(modData.scalar.Variable, varLabel),:);
                squaredError = (yobs - ymod) .^ 2;
                costComponents.(varLabel) = sum(squaredError(:)) / numel(squaredError);
            end
            % Size spectra data
            for i = 1:length(VarsSize)
                varLabel = VarsSize{i};
                ind0 = strcmp(Data.size.dataBinned.Variable, varLabel);
                % autotrophs
                ind = ind0 & strcmp(Data.size.dataBinned.trophicLevel, 'autotroph');
                yobs = Data.size.dataBinned.scaled_Value(ind);
                if isfield(modData.size, 'scaled_Value')
                    ymod = modData.size.scaled_Value(ind,:);
                else
                    ymod = modData.size.Value(ind,:);
                    ymod = (ymod - mean(ymod)) ./ std(ymod);
                end
                squaredError = (yobs - ymod) .^ 2;
                costComponents.([varLabel '_autotroph']) = sum(squaredError(:)) / numel(squaredError);
                % heterotrophs
                ind = ind0 & strcmp(Data.size.dataBinned.trophicLevel, 'heterotroph');
                yobs = Data.size.dataBinned.scaled_Value(ind);
                if isfield(modData.size, 'scaled_Value')
                    ymod = modData.size.scaled_Value(ind,:);
                else
                    ymod = modData.size.Value(ind,:);
                    ymod = (ymod - mean(ymod)) ./ std(ymod);
                end
                squaredError = (yobs - ymod) .^ 2;
                costComponents.([varLabel '_heterotroph']) = sum(squaredError(:)) / numel(squaredError);
            end            
            fields = fieldnames(costComponents);
            x = struct2array(costComponents);
            scalarCosts = x(contains(fields, Vars));
            sizeCosts = x(contains(fields, VarsSize));
            cost = mean(scalarCosts) + mean(sizeCosts);
            
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~            
            
        case 'RMS' % Root mean square
            % Scalar data
            for i = 1:length(Vars)
                varLabel = Vars{i};
                yobs = Data.scalar.scaled_Value(strcmp(Data.scalar.Variable, varLabel));
                ymod = modData.scalar.scaled_Value(strcmp(modData.scalar.Variable, varLabel),:);
                absError = sqrt((yobs - ymod) .^ 2);
                costComponents.(varLabel) = sum(absError(:)) / numel(absError);
            end
            % Size spectra data
            for i = 1:length(VarsSize)
                varLabel = VarsSize{i};
                ind0 = strcmp(Data.size.dataBinned.Variable, varLabel);
                % autotrophs
                ind = ind0 & strcmp(Data.size.dataBinned.trophicLevel, 'autotroph');
                yobs = Data.size.dataBinned.scaled_Value(ind);
                if isfield(modData.size, 'scaled_Value')
                    ymod = modData.size.scaled_Value(ind,:);
                else
                    ymod = modData.size.Value(ind,:);
                    ymod = (ymod - mean(ymod)) ./ std(ymod);
                end
                absError = sqrt((yobs - ymod) .^ 2);
                costComponents.([varLabel '_autotroph']) = sum(absError(:)) / numel(absError);
                % heterotrophs
                ind = ind0 & strcmp(Data.size.dataBinned.trophicLevel, 'heterotroph');
                yobs = Data.size.dataBinned.scaled_Value(ind);
                if isfield(modData.size, 'scaled_Value')
                    ymod = modData.size.scaled_Value(ind,:);
                else
                    ymod = modData.size.Value(ind,:);
                    ymod = (ymod - mean(ymod)) ./ std(ymod);
                end
                absError = sqrt((yobs - ymod) .^ 2);
                costComponents.([varLabel '_heterotroph']) = sum(absError(:)) / numel(absError);
                costComponents.(varLabel) = costComponents.([varLabel '_autotroph']) + costComponents.([varLabel '_heterotroph']);
            end
            fields = fieldnames(costComponents);
            x = struct2array(costComponents);
            scalarCosts = x(contains(fields, Vars));
            sizeCosts = x(contains(fields, VarsSize));
            cost = mean(scalarCosts) + mean(sizeCosts);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~            

        case 'Hellinger2_groupWaterOrigin'
            % Calculate non-parameteric Hellinger distances for the scalar
            % (nutrient) data and for size spectra vectors of relative
            % abundance.
            % Use a different, but comparable, metric for the total abundance
            % data contained in the size spectra. This is necessary because
            % these data do not form distributions so Hellinger distance is not
            % applicable. Like the Hellinger distance, the metric takes values
            % in the in the interval [0,1], with 0 represetning a perfect fit
            % and 1 being an extremely poor fit.
            % This method makes zero assumptions about variability -- unlike
            % the other methods which use variability across particle
            % trajectories to create distributions. The Hellinger distances are
            % calculated numerically instead of using formulas based on
            % assumptions of Gaussian distributions. The scalar data are
            % continuous and approximately Gaussian whereas the size spectra
            % data are discrete, so Hellinger distance calculations are
            % different for each...
            
            % Scalar data
            for i = 1:length(Vars)
                varLabel = Vars{i};
                yobs = Data.scalar.scaled_Value(strcmp(Data.scalar.Variable, varLabel));
                ymod = modData.scalar.scaled_Value(strcmp(modData.scalar.Variable, varLabel),:);
                n = size(ymod, 1);
                m = size(ymod, 2);
                % As the standardised data are approximately normally distributed,
                % we could assume that identically transformed model outputs
                % are also normally distributed, then find their means and
                % variances to calculate Hellinger distances analytically.
                % However, assuming normality of transformed model output is
                % inappropriate because it could be any shape...
                % Thus, derive cdfs and pdfs directly without assuming that
                % transformed model outputs follow normal distributions.
                % To compare observed and modelled distributions, the pdfs will
                % need to be evaluated across identical domains.
                cdf = (1:n)' ./ n;
                yobsc = sort(yobs); % sort observations
                ymodc = sort(ymod);
                % remove any duplicate values (more likely in model output than
                % in data), retaining the largest probability values
                keep = ~[diff(yobsc) == 0; false];
                cdfobs = cdf(keep);
                yobsc = yobsc(keep);
                % store cdfmod as cell-array because vector lengths may differ
                % after removing duplicates
                keep = ~[diff(ymodc) == 0; false(1, m)];
                cdfmod = cell(1, m);
                ymodc_ = cell(1, m);
                for ij = 1:m
                    cdfmod{:,ij} = cdf(keep(:,ij));
                    ymodc_{:,ij} = ymodc(keep(:,ij),ij);
                end
                % define regular grid across measurement space -- define number
                % of grid nodes, ng, relative to number of observations, n.
                vrange = [min([yobsc(:); ymodc(:)]), max([yobsc(:); ymodc(:)])];
                nr = 1; % choice of grid resolution influences cost values... is there a better way to do this, or a principled way to choose grid resolution?
%                 ng = 9; 
                ng = ceil(nr * n);
                grid = nan(ng, 1);
                grid(2:end) = linspace(vrange(1), vrange(2), ng-1)';
                sp = diff(grid(2:3));
                grid(1) = grid(2) - sp;
                % interpolate cdfs and derive pdfs
                cdfobsi = interp1(yobsc, cdfobs, grid, 'linear', 'extrap');
                cdfobsi(cdfobsi < min(cdfobs)) = 0;
                cdfobsi(cdfobsi > 1) = 1;
                pdfobsi = diff(cdfobsi);
                cdfmodi = nan(ng, m);
                pdfmodi = nan(ng-1, m);
                for ij = 1:m
                    cdfmodi(:,ij) = interp1(ymodc_{ij}, cdfmod{ij}, grid, 'linear', 'extrap');
                    cdfmodi(cdfmodi(:,ij) < min(cdfmod{ij}), ij) = 0;
                    cdfmodi(cdfmodi(:,ij) > 1, ij) = 1;
                    pdfmodi(:,ij) = diff(cdfmodi(:,ij));
                end
                hellingerDistance = (1 - sum((pdfobsi .* pdfmodi) .^ 0.5)) .^ 0.5;
                costComponents.(varLabel) = mean(hellingerDistance); % average over trajectory selections
            end
            
            % Vector (size) data
            groupedByWaterOrigin = isfield(Data.size.dataBinned, 'waterMass');
            a = log(3) / log(2); % steepness of cost metric for total abundance
            for i = 1:length(VarsSize)
                varLabel = VarsSize{i};
                ind0 = strcmp(Data.size.dataBinned.Variable, varLabel);
                if groupedByWaterOrigin
                    waterMasses = unique(Data.size.dataBinned.waterMass);
                else
                    waterMasses = {[]};
                end
                for w = 1:length(waterMasses)
                    wm = waterMasses{w};
                    ind1 = ind0;
                    if groupedByWaterOrigin
                        ind1 = ind0 & strcmp(Data.size.dataBinned.waterMass, wm);
                    end
                    if groupedByWaterOrigin
                        label_autoRel = [varLabel '_' wm '_autotroph_Rel'];
                        label_autoTot = [varLabel '_' wm '_autotroph_Tot'];
                        label_heteroRel = [varLabel '_' wm '_heterotroph_Rel'];
                        label_heteroTot = [varLabel '_' wm '_heterotroph_Tot'];
                    else
                        label_autoRel = [varLabel '_autotroph_Rel'];
                        label_autoTot = [varLabel '_autotroph_Tot'];
                        label_heteroRel = [varLabel '_heterotroph_Rel'];
                        label_heteroTot = [varLabel '_heterotroph_Tot'];
                    end
                    
                    % autotrophs
                    ind = ind1 & strcmp(Data.size.dataBinned.trophicLevel, 'autotroph');
                    yobs = Data.size.dataBinned.Value(ind);
                    ymod = modData.size.Value(ind,:);
                    nsize = size(ymod, 1);
                    ymod = reshape(ymod(~isnan(ymod)), nsize, []); % remove NaNs ymod contains when fitting to Arctic AND Atlantic data
                    nsample = size(ymod, 2);
                    yobsTot = sum(yobs);
                    ymodTot = sum(ymod);
                    % Cost metric for relative abundance (vector)
                    cdfobs = cumsum(yobs) ./ yobsTot;
                    pdfobs = diff([0; cdfobs]);
                    cdfmod = cumsum(ymod) ./ ymodTot;
                    pdfmod = diff([zeros(1, nsample); cdfmod]);
                    hellingerDistance = (1 - sum((pdfobs .* pdfmod) .^ 0.5)) .^ 0.5;
                    costComponents.(label_autoRel) = mean(hellingerDistance); % average over trajectory selections
                    % Cost metric for total abundance (scalar)
                    u = abs(log(ymodTot / yobsTot));
                    u = exp(-a .* u);
                    z = (1 - u) ./ (1 + u);
                    costComponents.(label_autoTot) = mean(z); % average over trajectory selections
                    
                    % heterotrophs
                    ind = ind1 & strcmp(Data.size.dataBinned.trophicLevel, 'heterotroph');
                    yobs = Data.size.dataBinned.Value(ind);
                    ymod = modData.size.Value(ind,:);
                    nsize = size(ymod, 1);
                    ymod = reshape(ymod(~isnan(ymod)), nsize, []); % remove NaNs ymod contains when fitting to Arctic AND Atlantic data
                    nsample = size(ymod, 2);
                    yobsTot = sum(yobs);
                    ymodTot = sum(ymod);
                    % Cost metric for relative abundance (vector)
                    cdfobs = cumsum(yobs) ./ yobsTot;
                    pdfobs = diff([0; cdfobs]);
                    cdfmod = cumsum(ymod) ./ ymodTot;
                    pdfmod = diff([zeros(1, nsample); cdfmod]);
                    hellingerDistance = (1 - sum((pdfobs .* pdfmod) .^ 0.5)) .^ 0.5;
                    costComponents.(label_heteroRel) = mean(hellingerDistance);
                    % Cost metric for total abundance (scalar)
                    u = abs(log(ymodTot / yobsTot));
                    u = exp(-a .* u);
                    z = (1 - u) ./ (1 + u);
                    costComponents.(label_heteroTot) = mean(z); % average over trajectory selections
                end
            end
            
            % Within each size data group, weight relative abundance-at-size
            % relative to total abundance.
            weight_relVsTot = 3; % weighting factor of relative vs total abundance
            costSize = zeros(2,length(waterMasses)); % store weighted costs for all size data groups
            for i = 1:length(waterMasses)
                wm = waterMasses{i};
                if groupedByWaterOrigin
                    cta = costComponents.([varLabel '_' wm '_autotroph_Tot']);
                    cra = costComponents.([varLabel '_' wm '_autotroph_Rel']);
                    cth = costComponents.([varLabel '_' wm '_heterotroph_Tot']);
                    crh = costComponents.([varLabel '_' wm '_heterotroph_Rel']);
                else
                    cta = costComponents.([varLabel '_autotroph_Tot']);
                    cra = costComponents.([varLabel '_autotroph_Rel']);
                    cth = costComponents.([varLabel '_heterotroph_Tot']);
                    crh = costComponents.([varLabel '_heterotroph_Rel']);
                end
                costSize(1,i) = mean(2 .* [1, weight_relVsTot] ./ (weight_relVsTot+1) .* [cta, cra]);
                costSize(2,i) = mean(2 .* [1, weight_relVsTot] ./ (weight_relVsTot+1) .* [cth, crh]);
            end
            
            costNutrient = zeros(1,length(Vars));
            for i = 1:length(Vars)
                costNutrient(i) = costComponents.(Vars{i});
            end
            
            % Average the cost across data-types (nutrient & size)
            cost = [mean(costNutrient), mean(costSize(:))];
            % Assign size vs nutrient weighting
            weight_sizeVsNutrient = 1; % weighting factor of size vs nutrient data
            weights = 2 .* [1, weight_sizeVsNutrient] ./ (weight_sizeVsNutrient+1);
            cost = weights .* cost;
            
            cost = mean(cost); % finally, average cost over nutrient and size data components
            
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~            
            
        case 'Hellinger3_groupWaterOrigin'
            
            % Scalar data
            for i = 1:length(Vars)
                varLabel = Vars{i};
                yobs = Data.scalar.scaled_Value(strcmp(Data.scalar.Variable, varLabel));
                ymod = modData.scalar.scaled_Value(strcmp(modData.scalar.Variable, varLabel),:);
                n = size(ymod, 1);
                m = size(ymod, 2);
                yobsc = sort(yobs); % sort observations
                ymodc = sort(ymod);
                cdf = (1:n)' ./ n;
                % remove any duplicate values (more likely in model output than
                % in data), retaining the largest probability values
                keep = ~[diff(yobsc) == 0; false];
                cdfobs = cdf(keep);
                yobsc = yobsc(keep);
                % store cdfmod as cell-array because vector lengths may differ
                % after removing duplicates
                keep = ~[diff(ymodc) == 0; false(1, m)];
                cdfmod = cell(1, m);
                ymodc_ = cell(1, m);
                for ij = 1:m
                    cdfmod{:,ij} = cdf(keep(:,ij));
                    ymodc_{:,ij} = ymodc(keep(:,ij),ij);
                end
                
                % interpolate data-CDF over modelled query points
                cdfobsi = cell(1,m);
                for ij = 1:m
                    cdfobsi{ij} = interp1(yobsc, cdfobs, ymodc_{ij}, 'linear', 'extrap');
                    cdfobsi{ij}(cdfobsi{ij} < 0) = 0;
                    cdfobsi{ij}(cdfobsi{ij} > 1) = 1;
                end
                % distances between each modelled CDF value and equivalent
                % data-CDF value, and use average distance as the metric
                cdfDist = cell(1,m);
                avDist = nan(1,m);
                for ij = 1:m
                    cdfDist{ij} = abs(cdfobsi{ij} - cdfmod{ij});
                    avDist(ij) = mean(cdfDist{ij});
                end
                costComponents.(varLabel) = mean(avDist); % average over trajectory selections
            end
                
                
            % Vector (size) data
            groupedByWaterOrigin = isfield(Data.size.dataBinned, 'waterMass');
            a = log(3) / log(2); % steepness of cost metric for total abundance
            for i = 1:length(VarsSize)
                varLabel = VarsSize{i};
                ind0 = strcmp(Data.size.dataBinned.Variable, varLabel);
                if groupedByWaterOrigin
                    waterMasses = unique(Data.size.dataBinned.waterMass);
                else
                    waterMasses = {[]};
                end
                for w = 1:length(waterMasses)
                    wm = waterMasses{w};
                    ind1 = ind0;
                    if groupedByWaterOrigin
                        ind1 = ind0 & strcmp(Data.size.dataBinned.waterMass, wm);
                    end
                    if groupedByWaterOrigin
                        label_autoRel = [varLabel '_' wm '_autotroph_Rel'];
                        label_autoTot = [varLabel '_' wm '_autotroph_Tot'];
                        label_heteroRel = [varLabel '_' wm '_heterotroph_Rel'];
                        label_heteroTot = [varLabel '_' wm '_heterotroph_Tot'];
                    else
                        label_autoRel = [varLabel '_autotroph_Rel'];
                        label_autoTot = [varLabel '_autotroph_Tot'];
                        label_heteroRel = [varLabel '_heterotroph_Rel'];
                        label_heteroTot = [varLabel '_heterotroph_Tot'];
                    end
                    
                    % autotrophs
                    ind = ind1 & strcmp(Data.size.dataBinned.trophicLevel, 'autotroph');
                    yobs = Data.size.dataBinned.Value(ind);
                    ymod = modData.size.Value(ind,:);
                    nsize = size(ymod, 1);
                    ymod = reshape(ymod(~isnan(ymod)), nsize, []); % remove NaNs ymod contains when fitting to Arctic AND Atlantic data
                    nsample = size(ymod, 2);
                    yobsTot = sum(yobs);
                    ymodTot = sum(ymod);
                    % Cost metric for relative abundance (vector)
                    cdfobs = cumsum(yobs) ./ yobsTot;
                    pdfobs = diff([0; cdfobs]);
                    cdfmod = cumsum(ymod) ./ ymodTot;
                    pdfmod = diff([zeros(1, nsample); cdfmod]);
                    hellingerDistance = (1 - sum((pdfobs .* pdfmod) .^ 0.5)) .^ 0.5;
                    costComponents.(label_autoRel) = mean(hellingerDistance); % average over trajectory selections
                    % Cost metric for total abundance (scalar)
                    u = abs(log(ymodTot / yobsTot));
                    u = exp(-a .* u);
                    z = (1 - u) ./ (1 + u);
                    costComponents.(label_autoTot) = mean(z); % average over trajectory selections
                    
                    % heterotrophs
                    ind = ind1 & strcmp(Data.size.dataBinned.trophicLevel, 'heterotroph');
                    yobs = Data.size.dataBinned.Value(ind);
                    ymod = modData.size.Value(ind,:);
                    nsize = size(ymod, 1);
                    ymod = reshape(ymod(~isnan(ymod)), nsize, []); % remove NaNs ymod contains when fitting to Arctic AND Atlantic data
                    nsample = size(ymod, 2);
                    yobsTot = sum(yobs);
                    ymodTot = sum(ymod);
                    % Cost metric for relative abundance (vector)
                    cdfobs = cumsum(yobs) ./ yobsTot;
                    pdfobs = diff([0; cdfobs]);
                    cdfmod = cumsum(ymod) ./ ymodTot;
                    pdfmod = diff([zeros(1, nsample); cdfmod]);
                    hellingerDistance = (1 - sum((pdfobs .* pdfmod) .^ 0.5)) .^ 0.5;
                    costComponents.(label_heteroRel) = mean(hellingerDistance);
                    % Cost metric for total abundance (scalar)
                    u = abs(log(ymodTot / yobsTot));
                    u = exp(-a .* u);
                    z = (1 - u) ./ (1 + u);
                    costComponents.(label_heteroTot) = mean(z); % average over trajectory selections
                end
            end
            
            % Within each size data group, weight relative abundance-at-size
            % relative to total abundance.
            weight_relVsTot = 3; % weighting factor of relative vs total abundance
            costSize = zeros(2,length(waterMasses)); % store weighted costs for all size data groups
            for i = 1:length(waterMasses)
                wm = waterMasses{i};
                if groupedByWaterOrigin
                    cta = costComponents.([varLabel '_' wm '_autotroph_Tot']);
                    cra = costComponents.([varLabel '_' wm '_autotroph_Rel']);
                    cth = costComponents.([varLabel '_' wm '_heterotroph_Tot']);
                    crh = costComponents.([varLabel '_' wm '_heterotroph_Rel']);
                else
                    cta = costComponents.([varLabel '_autotroph_Tot']);
                    cra = costComponents.([varLabel '_autotroph_Rel']);
                    cth = costComponents.([varLabel '_heterotroph_Tot']);
                    crh = costComponents.([varLabel '_heterotroph_Rel']);
                end
                costSize(1,i) = mean(2 .* [1, weight_relVsTot] ./ (weight_relVsTot+1) .* [cta, cra]);
                costSize(2,i) = mean(2 .* [1, weight_relVsTot] ./ (weight_relVsTot+1) .* [cth, crh]);
            end
                        
            costNutrient = zeros(1,length(Vars));
            for i = 1:length(Vars)
                costNutrient(i) = costComponents.(Vars{i});
            end
            
            % Average the cost across data-types (nutrient & size)
            cost = [mean(costNutrient), mean(costSize(:))];
            % Assign size vs nutrient weighting
            weight_sizeVsNutrient = 1; % weighting factor of size vs nutrient data
            weights = 2 .* [1, weight_sizeVsNutrient] ./ (weight_sizeVsNutrient+1);
            cost = weights .* cost;
            
            cost = mean(cost); % finally, average cost over nutrient and size data components

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        case 'IQD_Hellinger_groupWaterOrigin'
            
            % Scalar data
            for i = 1:length(Vars)
                varLabel = Vars{i};
                yobs = Data.scalar.scaled_Value(strcmp(Data.scalar.Variable, varLabel));
                ymod = modData.scalar.scaled_Value(strcmp(modData.scalar.Variable, varLabel),:);
                n = size(ymod, 1);
                m = size(ymod, 2);
                yobsc = sort(yobs); % sort observations
                ymodc = sort(ymod);
                cdf = (1:n)' ./ n;
                % remove any duplicate values (more likely in model output than
                % in data), retaining the largest probability values
                keep = ~[diff(yobsc) == 0; false];
                cdfobs = cdf(keep);
                yobsc = yobsc(keep);
                % store cdfmod as cell-array because vector lengths may differ
                % after removing duplicates
                keep = ~[diff(ymodc) == 0; false(1, m)];
                cdfmod = cell(1, m);
                ymodc_ = cell(1, m);
                for ij = 1:m
                    cdfmod{:,ij} = cdf(keep(:,ij));
                    ymodc_{:,ij} = ymodc(keep(:,ij),ij);
                end
                
                % interpolate modelled CDFs over query points defined by
                % the data
                cdfmodi = cell(1,m);
                for ij =1:m
                    cdfmodi{ij} = interp1(ymodc_{ij}, cdfmod{ij}, yobsc, 'linear', 'extrap');
                    cdfmodi{ij}(cdfmodi{ij} < 0) = 0;
                    cdfmodi{ij}(cdfmodi{ij} > 1) = 1;
                end
                
                % distances between each data CDF value and equivalent 
                % modelled CDF value -- use average distance as the metric
                cdfDist = cell(1,m);
                avDist = nan(1,m);
                for ij = 1:m
                    cdfDist{ij} = (cdfobs - cdfmodi{ij}) .^ 2;
                    avDist(ij) = mean(cdfDist{ij}) .^ 0.5;
                end
                costComponents.(varLabel) = mean(avDist); % average over trajectory selections
            end    
                
            % Vector (size) data
            groupedByWaterOrigin = isfield(Data.size.dataBinned, 'waterMass');
            a = log(3) / log(2); % steepness of cost metric for total abundance
            for i = 1:length(VarsSize)
                varLabel = VarsSize{i};
                ind0 = strcmp(Data.size.dataBinned.Variable, varLabel);
                if groupedByWaterOrigin
                    waterMasses = unique(Data.size.dataBinned.waterMass);
                else
                    waterMasses = {[]};
                end
                for w = 1:length(waterMasses)
                    wm = waterMasses{w};
                    ind1 = ind0;
                    if groupedByWaterOrigin
                        ind1 = ind0 & strcmp(Data.size.dataBinned.waterMass, wm);
                    end
                    if groupedByWaterOrigin
                        label_autoRel = [varLabel '_' wm '_autotroph_Rel'];
                        label_autoTot = [varLabel '_' wm '_autotroph_Tot'];
                        label_heteroRel = [varLabel '_' wm '_heterotroph_Rel'];
                        label_heteroTot = [varLabel '_' wm '_heterotroph_Tot'];
                    else
                        label_autoRel = [varLabel '_autotroph_Rel'];
                        label_autoTot = [varLabel '_autotroph_Tot'];
                        label_heteroRel = [varLabel '_heterotroph_Rel'];
                        label_heteroTot = [varLabel '_heterotroph_Tot'];
                    end
                    
                    % autotrophs
                    ind = ind1 & strcmp(Data.size.dataBinned.trophicLevel, 'autotroph');
                    yobs = Data.size.dataBinned.Value(ind);
                    ymod = modData.size.Value(ind,:);
                    nsize = size(ymod, 1);
                    ymod = reshape(ymod(~isnan(ymod)), nsize, []); % remove NaNs ymod contains when fitting to Arctic AND Atlantic data
                    nsample = size(ymod, 2);
                    yobsTot = sum(yobs);
                    ymodTot = sum(ymod);
                    % Cost metric for relative abundance (vector)
                    cdfobs = cumsum(yobs) ./ yobsTot;
                    pdfobs = diff([0; cdfobs]);
                    cdfmod = cumsum(ymod) ./ ymodTot;
                    pdfmod = diff([zeros(1, nsample); cdfmod]);
                    hellingerDistance = (1 - sum((pdfobs .* pdfmod) .^ 0.5)) .^ 0.5;
                    costComponents.(label_autoRel) = mean(hellingerDistance); % average over trajectory selections
                    % Cost metric for total abundance (scalar)
                    u = abs(log(ymodTot / yobsTot));
                    u = exp(-a .* u);
                    z = (1 - u) ./ (1 + u);
                    costComponents.(label_autoTot) = mean(z); % average over trajectory selections
                    
                    % heterotrophs
                    ind = ind1 & strcmp(Data.size.dataBinned.trophicLevel, 'heterotroph');
                    yobs = Data.size.dataBinned.Value(ind);
                    ymod = modData.size.Value(ind,:);
                    nsize = size(ymod, 1);
                    ymod = reshape(ymod(~isnan(ymod)), nsize, []); % remove NaNs ymod contains when fitting to Arctic AND Atlantic data
                    nsample = size(ymod, 2);
                    yobsTot = sum(yobs);
                    ymodTot = sum(ymod);
                    % Cost metric for relative abundance (vector)
                    cdfobs = cumsum(yobs) ./ yobsTot;
                    pdfobs = diff([0; cdfobs]);
                    cdfmod = cumsum(ymod) ./ ymodTot;
                    pdfmod = diff([zeros(1, nsample); cdfmod]);
                    hellingerDistance = (1 - sum((pdfobs .* pdfmod) .^ 0.5)) .^ 0.5;
                    costComponents.(label_heteroRel) = mean(hellingerDistance);
                    % Cost metric for total abundance (scalar)
                    u = abs(log(ymodTot / yobsTot));
                    u = exp(-a .* u);
                    z = (1 - u) ./ (1 + u);
                    costComponents.(label_heteroTot) = mean(z); % average over trajectory selections
                end
            end
            
            % Within each size data group, weight relative abundance-at-size
            % relative to total abundance.
            weight_relVsTot = 3; % weighting factor of relative vs total abundance
            costSize = zeros(2,length(waterMasses)); % store weighted costs for all size data groups
            for i = 1:length(waterMasses)
                wm = waterMasses{i};
                if groupedByWaterOrigin
                    cta = costComponents.([varLabel '_' wm '_autotroph_Tot']);
                    cra = costComponents.([varLabel '_' wm '_autotroph_Rel']);
                    cth = costComponents.([varLabel '_' wm '_heterotroph_Tot']);
                    crh = costComponents.([varLabel '_' wm '_heterotroph_Rel']);
                else
                    cta = costComponents.([varLabel '_autotroph_Tot']);
                    cra = costComponents.([varLabel '_autotroph_Rel']);
                    cth = costComponents.([varLabel '_heterotroph_Tot']);
                    crh = costComponents.([varLabel '_heterotroph_Rel']);
                end
                costSize(1,i) = mean(2 .* [1, weight_relVsTot] ./ (weight_relVsTot+1) .* [cta, cra]);
                costSize(2,i) = mean(2 .* [1, weight_relVsTot] ./ (weight_relVsTot+1) .* [cth, crh]);
            end
                        
            costNutrient = zeros(1,length(Vars));
            for i = 1:length(Vars)
                costNutrient(i) = costComponents.(Vars{i});
            end
            
            % Average the cost across data-types (nutrient & size)
            cost = [mean(costNutrient), mean(costSize(:))];
            % Assign size vs nutrient weighting
            weight_sizeVsNutrient = 1; % weighting factor of size vs nutrient data
            weights = 2 .* [1, weight_sizeVsNutrient] ./ (weight_sizeVsNutrient+1);
            cost = weights .* cost;
            
            cost = mean(cost); % finally, average cost over nutrient and size data components

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        case 'meanCDFdist_Hellinger'
            
            % This method uses [0,1] metrics for all data types.
            % The scalar data are fitted by reference to their CDFs in a
            % manner that explicitly accounts for each data point rather
            % than fitting to a distribution.
            % The relative abundance-at-size data are fitted with Hellinger
            % distances.
            % The total abundance from the size data are fitted with my
            % custom metric (relies on a ratio rather like Shannon's entropy)
            
            % Scalar data
            for i = 1:length(Vars)
                varLabel = Vars{i};
                yobs = Data.scalar.scaled_Value(strcmp(Data.scalar.Variable, varLabel));
                ymod = modData.scalar.scaled_Value(strcmp(modData.scalar.Variable, varLabel),:);
                n = size(ymod, 1);
%                 m = size(ymod, 2);
                cdf = (1:n)' ./ n;
                [yobs_sort, o] = sort(yobs); % sort observations
                % put modelled output in same order -- these values will
                % almost certainly NOT be in ascending order...
                ymod_sort = ymod(o,:);
%                 figure(1)
%                 plot(yobs_sort, cdf)
%                 hold on
%                 plot(ymod_sort, cdf)
%                 hold off
                % find distances between each modelled data point and the
                % empirical data CDF
                ymodi = interp1(yobs_sort, cdf, ymod_sort, 'linear', 'extrap');
                ymodi = min(1, max(0, ymodi));
                CDFdist = abs(cdf - ymodi);
                avDist = mean(CDFdist); % average over data points
                costComponents.(varLabel) = mean(avDist); % average over trajectory selections
            end    
                
            % Vector (size) data
            groupedByWaterOrigin = isfield(Data.size.dataBinned, 'waterMass');
            a = log(3) / log(2); % steepness of cost metric for total abundance
            for i = 1:length(VarsSize)
                varLabel = VarsSize{i};
                ind0 = strcmp(Data.size.dataBinned.Variable, varLabel);
                if groupedByWaterOrigin
                    waterMasses = unique(Data.size.dataBinned.waterMass);
                else
                    waterMasses = {[]};
                end
                for w = 1:length(waterMasses)
                    wm = waterMasses{w};
                    ind1 = ind0;
                    if groupedByWaterOrigin
                        ind1 = ind0 & strcmp(Data.size.dataBinned.waterMass, wm);
                    end
                    if groupedByWaterOrigin
                        label_autoRel = [varLabel '_' wm '_autotroph_Rel'];
                        label_autoTot = [varLabel '_' wm '_autotroph_Tot'];
                        label_heteroRel = [varLabel '_' wm '_heterotroph_Rel'];
                        label_heteroTot = [varLabel '_' wm '_heterotroph_Tot'];
                    else
                        label_autoRel = [varLabel '_autotroph_Rel'];
                        label_autoTot = [varLabel '_autotroph_Tot'];
                        label_heteroRel = [varLabel '_heterotroph_Rel'];
                        label_heteroTot = [varLabel '_heterotroph_Tot'];
                    end
                    
                    % autotrophs
                    ind = ind1 & strcmp(Data.size.dataBinned.trophicLevel, 'autotroph');
                    yobs = Data.size.dataBinned.Value(ind);
                    ymod = modData.size.Value(ind,:);
                    nsize = size(ymod, 1);
                    ymod = reshape(ymod(~isnan(ymod)), nsize, []); % remove NaNs ymod contains when fitting to Arctic AND Atlantic data
                    nsample = size(ymod, 2);
                    yobsTot = sum(yobs);
                    ymodTot = sum(ymod);
                    % Cost metric for relative abundance (vector)
                    cdfobs = cumsum(yobs) ./ yobsTot;
                    pdfobs = diff([0; cdfobs]);
                    cdfmod = cumsum(ymod) ./ ymodTot;
                    pdfmod = diff([zeros(1, nsample); cdfmod]);
                    hellingerDistance = (1 - sum((pdfobs .* pdfmod) .^ 0.5)) .^ 0.5;
                    costComponents.(label_autoRel) = mean(hellingerDistance); % average over trajectory selections
                    % Cost metric for total abundance (scalar)
                    u = abs(log(ymodTot / yobsTot));
                    u = exp(-a .* u);
                    z = (1 - u) ./ (1 + u);
                    costComponents.(label_autoTot) = mean(z); % average over trajectory selections
                    
                    % heterotrophs
                    ind = ind1 & strcmp(Data.size.dataBinned.trophicLevel, 'heterotroph');
                    yobs = Data.size.dataBinned.Value(ind);
                    ymod = modData.size.Value(ind,:);
                    nsize = size(ymod, 1);
                    ymod = reshape(ymod(~isnan(ymod)), nsize, []); % remove NaNs ymod contains when fitting to Arctic AND Atlantic data
                    nsample = size(ymod, 2);
                    yobsTot = sum(yobs);
                    ymodTot = sum(ymod);
                    % Cost metric for relative abundance (vector)
                    cdfobs = cumsum(yobs) ./ yobsTot;
                    pdfobs = diff([0; cdfobs]);
                    cdfmod = cumsum(ymod) ./ ymodTot;
                    pdfmod = diff([zeros(1, nsample); cdfmod]);
                    hellingerDistance = (1 - sum((pdfobs .* pdfmod) .^ 0.5)) .^ 0.5;
                    costComponents.(label_heteroRel) = mean(hellingerDistance);
                    % Cost metric for total abundance (scalar)
                    u = abs(log(ymodTot / yobsTot));
                    u = exp(-a .* u);
                    z = (1 - u) ./ (1 + u);
                    costComponents.(label_heteroTot) = mean(z); % average over trajectory selections
                end
            end
            
%             % Within each size data group, weight relative abundance-at-size
%             % relative to total abundance.
%             weight_relVsTot = 3; % weighting factor of relative vs total abundance
%             costSize = zeros(2,length(waterMasses)); % store weighted costs for all size data groups
%             for i = 1:length(waterMasses)
%                 wm = waterMasses{i};
%                 if groupedByWaterOrigin
%                     cta = costComponents.([varLabel '_' wm '_autotroph_Tot']);
%                     cra = costComponents.([varLabel '_' wm '_autotroph_Rel']);
%                     cth = costComponents.([varLabel '_' wm '_heterotroph_Tot']);
%                     crh = costComponents.([varLabel '_' wm '_heterotroph_Rel']);
%                 else
%                     cta = costComponents.([varLabel '_autotroph_Tot']);
%                     cra = costComponents.([varLabel '_autotroph_Rel']);
%                     cth = costComponents.([varLabel '_heterotroph_Tot']);
%                     crh = costComponents.([varLabel '_heterotroph_Rel']);
%                 end
%                 costSize(1,i) = mean(2 .* [1, weight_relVsTot] ./ (weight_relVsTot+1) .* [cta, cra]);
%                 costSize(2,i) = mean(2 .* [1, weight_relVsTot] ./ (weight_relVsTot+1) .* [cth, crh]);
%             end
                        
            
            % Omit the total abundance info from the size data --
            % simplified compared to the commented code above...
            costSize = zeros(2,length(waterMasses)); % store weighted costs for all size data groups
            for i = 1:length(waterMasses)
                wm = waterMasses{i};
                if groupedByWaterOrigin
                    cra = costComponents.([varLabel '_' wm '_autotroph_Rel']);
                    crh = costComponents.([varLabel '_' wm '_heterotroph_Rel']);
                else
                    cra = costComponents.([varLabel '_autotroph_Rel']);
                    crh = costComponents.([varLabel '_heterotroph_Rel']);
                end
                costSize(1,i) = cra;
                costSize(2,i) = crh;
            end

            
            costNutrient = zeros(1,length(Vars));
            for i = 1:length(Vars)
                costNutrient(i) = costComponents.(Vars{i});
            end
            
            % Average the cost across data-types (nutrient & size)
            cost = [mean(costNutrient), mean(costSize(:))];
            % Assign size vs nutrient weighting
            weight_sizeVsNutrient = 1; % weighting factor of size vs nutrient data
            weights = 2 .* [1, weight_sizeVsNutrient] ./ (weight_sizeVsNutrient+1);
            cost = weights .* cost;
            
            cost = mean(cost); % finally, average cost over nutrient and size data components
%             cost = cost(2); % try fitting only to the size data
            

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        case 'smoothCDFdist_Hellinger'
            
            % Scalar data
            for i = 1:length(Vars)
                varLabel = Vars{i};
                yobs = Data.scalar.scaled_Value(strcmp(Data.scalar.Variable, varLabel));
                ymod = modData.scalar.scaled_Value(strcmp(modData.scalar.Variable, varLabel),:);
                n = size(ymod, 1);
                
                % Sort observations
                [yobs_sort, o] = sort(yobs);
                CDFobs = [yobs_sort, ...
                    (1:n)' ./ n]; % Empirical CDF of standardised data
                
%                 figure
%                 plot(CDFobs(:,1), CDFobs(:,2)) % empirical CDF of data
                
%                 % Smooth the CDF
%                 smoothFactor = 0.15;
%                 CDFobs(:,2) = smooth(CDFobs(:,1), CDFobs(:,2), ... 
%                     smoothFactor * n, 'loess');
%                 % Tidy up boundaries
%                 ind = find(CDFobs(:,2) >= 1, 1);
%                 CDFobs(ind:end,2) = 1;
%                 ind = find(CDFobs(:,2) == min(CDFobs(:,2)));
%                 CDFobs(1:ind,2) = min(CDFobs(:,2));
%                 CDFobs(:,2) = max(0, CDFobs(:,2));
%                 plot(CDFobs(:,1), CDFobs(:,2)) % smoothed CDF of standardised data
                
%                 PDFobs = [CDFobs(2:end,1), diff(CDFobs(:,2))];
%                 plot(PDFobs(:,1), PDFobs(:,2))
                
                % Reorder modelled values to match sorted data
                ymod_sort = ymod(o,:);
                CDFmod = [ymod_sort, CDFobs(:,2)]; % align modelled values with the smoothed CDF of the data
                
%                 hold on
%                 
%                 scatter(CDFmod(:,1), CDFmod(:,2))                
%                 for ij = 1:n
%                     plot([CDFobs(ij,1), CDFmod(ij,1)], [CDFobs(ij,2), CDFmod(ij,2)], 'Color', [0.85, 0.85, 0.85])
%                 end
                
                % Variability in model outputs may result in a
                % noisy/non-smooth cost function.
                % Smoothing model outputs before comparing to data may help
                % by improving the signal:noise ratio... fitting to pattern
                % in data rather than each individual data point.
                
                % Try fitting by aligning the CDFs of data and model
                
                CDFmod_ = [sort(CDFmod(:,1)), (1:n)' ./ n];
                ind = diff(CDFmod_(:,1)) == 0;
                CDFmod_r = CDFmod_;
                CDFmod_r(ind,:) = [];
                
%                 plot(CDFmod_r(:,1), CDFmod_r(:,2))

                % Evaluate model CDF value for each data point
                matchCDF = interp1(CDFmod_r(:,1), CDFmod_r(:,2), CDFobs(:,1), ...
                    'linear', 'extrap');
                matchCDF = min(1, max(0, matchCDF));
                % Difference between observed and modelled CDFs
                CDFdist = abs(CDFobs(:,2) - matchCDF);
                
%                 figure
%                 plot(CDFobs(:,1), CDFobs(:,2))
%                 hold on
%                 plot(CDFmod_(:,1), CDFmod_(:,2))
%                 for ij = 1:n
%                     plot([CDFobs(ij,1), CDFobs(ij,1)], ...
%                         [CDFobs(ij,2), matchCDF(ij)], 'Color', [0.85, 0.85, 0.85])
%                 end
%                 
                % Average over data points and trajectory selections
                avFun = @mean;
%                 avFun = @geomean;
                avDist = avFun(CDFdist);  % average over data points
                costComponents.(varLabel) = avFun(avDist); % average over trajectory selections
            end
                
            % Vector (size) data
            groupedByWaterOrigin = isfield(Data.size.dataBinned, 'waterMass');
            a = log(3) / log(2); % steepness of cost metric for total abundance
            for i = 1:length(VarsSize)
                varLabel = VarsSize{i};
                ind0 = strcmp(Data.size.dataBinned.Variable, varLabel);
                if groupedByWaterOrigin
                    waterMasses = unique(Data.size.dataBinned.waterMass);
                else
                    waterMasses = {[]};
                end
                for w = 1:length(waterMasses)
                    wm = waterMasses{w};
                    ind1 = ind0;
                    if groupedByWaterOrigin
                        ind1 = ind0 & strcmp(Data.size.dataBinned.waterMass, wm);
                    end
                    if groupedByWaterOrigin
                        label_autoRel = [varLabel '_' wm '_autotroph_Rel'];
                        label_autoTot = [varLabel '_' wm '_autotroph_Tot'];
                        label_heteroRel = [varLabel '_' wm '_heterotroph_Rel'];
                        label_heteroTot = [varLabel '_' wm '_heterotroph_Tot'];
                    else
                        label_autoRel = [varLabel '_autotroph_Rel'];
                        label_autoTot = [varLabel '_autotroph_Tot'];
                        label_heteroRel = [varLabel '_heterotroph_Rel'];
                        label_heteroTot = [varLabel '_heterotroph_Tot'];
                    end
                    
                    % autotrophs
                    ind = ind1 & strcmp(Data.size.dataBinned.trophicLevel, 'autotroph');
                    yobs = Data.size.dataBinned.Value(ind);
                    ymod = modData.size.Value(ind,:);
                    nsize = size(ymod, 1);
                    ymod = reshape(ymod(~isnan(ymod)), nsize, []); % remove NaNs ymod contains when fitting to Arctic AND Atlantic data
                    nsample = size(ymod, 2);
                    yobsTot = sum(yobs);
                    ymodTot = sum(ymod);
                    % Cost metric for relative abundance (vector)
                    cdfobs = cumsum(yobs) ./ yobsTot;
                    pdfobs = diff([0; cdfobs]);
                    cdfmod = cumsum(ymod) ./ ymodTot;
                    pdfmod = diff([zeros(1, nsample); cdfmod]);
                    hellingerDistance = (1 - sum((pdfobs .* pdfmod) .^ 0.5)) .^ 0.5;
                    costComponents.(label_autoRel) = mean(hellingerDistance); % average over trajectory selections
                    % Cost metric for total abundance (scalar)
                    u = abs(log(ymodTot / yobsTot));
                    u = exp(-a .* u);
                    z = (1 - u) ./ (1 + u);
                    costComponents.(label_autoTot) = mean(z); % average over trajectory selections
                    
                    % heterotrophs
                    ind = ind1 & strcmp(Data.size.dataBinned.trophicLevel, 'heterotroph');
                    yobs = Data.size.dataBinned.Value(ind);
                    ymod = modData.size.Value(ind,:);
                    nsize = size(ymod, 1);
                    ymod = reshape(ymod(~isnan(ymod)), nsize, []); % remove NaNs ymod contains when fitting to Arctic AND Atlantic data
                    nsample = size(ymod, 2);
                    yobsTot = sum(yobs);
                    ymodTot = sum(ymod);
                    % Cost metric for relative abundance (vector)
                    cdfobs = cumsum(yobs) ./ yobsTot;
                    pdfobs = diff([0; cdfobs]);
                    cdfmod = cumsum(ymod) ./ ymodTot;
                    pdfmod = diff([zeros(1, nsample); cdfmod]);
                    hellingerDistance = (1 - sum((pdfobs .* pdfmod) .^ 0.5)) .^ 0.5;
                    costComponents.(label_heteroRel) = mean(hellingerDistance);
                    % Cost metric for total abundance (scalar)
                    u = abs(log(ymodTot / yobsTot));
                    u = exp(-a .* u);
                    z = (1 - u) ./ (1 + u);
                    costComponents.(label_heteroTot) = mean(z); % average over trajectory selections
                end
            end
            
%             % Within each size data group, weight relative abundance-at-size
%             % relative to total abundance.
%             weight_relVsTot = 3; % weighting factor of relative vs total abundance
%             costSize = zeros(2,length(waterMasses)); % store weighted costs for all size data groups
%             for i = 1:length(waterMasses)
%                 wm = waterMasses{i};
%                 if groupedByWaterOrigin
%                     cta = costComponents.([varLabel '_' wm '_autotroph_Tot']);
%                     cra = costComponents.([varLabel '_' wm '_autotroph_Rel']);
%                     cth = costComponents.([varLabel '_' wm '_heterotroph_Tot']);
%                     crh = costComponents.([varLabel '_' wm '_heterotroph_Rel']);
%                 else
%                     cta = costComponents.([varLabel '_autotroph_Tot']);
%                     cra = costComponents.([varLabel '_autotroph_Rel']);
%                     cth = costComponents.([varLabel '_heterotroph_Tot']);
%                     crh = costComponents.([varLabel '_heterotroph_Rel']);
%                 end
%                 costSize(1,i) = mean(2 .* [1, weight_relVsTot] ./ (weight_relVsTot+1) .* [cta, cra]);
%                 costSize(2,i) = mean(2 .* [1, weight_relVsTot] ./ (weight_relVsTot+1) .* [cth, crh]);
%             end
                        
            
            % Omit the total abundance info from the size data --
            % simplified compared to the commented code above...
            costSize = zeros(2,length(waterMasses)); % store weighted costs for all size data groups
            for i = 1:length(waterMasses)
                wm = waterMasses{i};
                if groupedByWaterOrigin
                    cra = costComponents.([varLabel '_' wm '_autotroph_Rel']);
                    crh = costComponents.([varLabel '_' wm '_heterotroph_Rel']);
                else
                    cra = costComponents.([varLabel '_autotroph_Rel']);
                    crh = costComponents.([varLabel '_heterotroph_Rel']);
                end
                costSize(1,i) = cra;
                costSize(2,i) = crh;
            end

            
            costNutrient = zeros(1,length(Vars));
            for i = 1:length(Vars)
                costNutrient(i) = costComponents.(Vars{i});
            end
            
            % Average the cost across data-types (nutrient & size)
            cost = [mean(costNutrient), mean(costSize(:))];
            % Assign size vs nutrient weighting
            weight_sizeVsNutrient = 1; % weighting factor of size vs nutrient data
            weights = 2 .* [1, weight_sizeVsNutrient] ./ (weight_sizeVsNutrient+1);
            cost = weights .* cost;
            
            cost = mean(cost); % finally, average cost over nutrient and size data components
%             cost = cost(2); % try fitting only to the size data


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        case 'smooth2CDFdist_Hellinger'
            
            % There are problems emerging from using CDF values to
            % calculate the scalar cost values. When modelled values are
            % either much larger or smaller than the data range, the cost
            % becomes insensitive... Probably better to use differences
            % between data and model, rather than converting to CDF
            % values... even though it might be difficult to scale for
            % compatibility with the size-data Hellinger distances.
            
            % Scalar data
            for i = 1:length(Vars)
                varLabel = Vars{i};
                yobs = Data.scalar.scaled_Value(strcmp(Data.scalar.Variable, varLabel));
                ymod = modData.scalar.scaled_Value(strcmp(modData.scalar.Variable, varLabel),:);
                n = size(ymod, 1);
                
                % Sort observations
                [yobs_sort, o] = sort(yobs);
                CDFobs = [yobs_sort, ...
                    (1:n)' ./ n]; % Empirical CDF of standardised data
                
                figure
                plot(CDFobs(:,1), CDFobs(:,2)) % empirical CDF of data
                                
                % Reorder modelled values to match sorted data
                ymod_sort = ymod(o,:);
                CDFmod = [ymod_sort, CDFobs(:,2)]; % align modelled values with the smoothed CDF of the data
                
                hold on
                
                scatter(CDFmod(:,1), CDFmod(:,2))
                for ij = 1:n
                    plot([CDFobs(ij,1), CDFmod(ij,1)], [CDFobs(ij,2), CDFmod(ij,2)], 'Color', [0.85, 0.85, 0.85])
                end
                
                % Variability in model outputs may result in a
                % noisy/non-smooth cost function.
                % Smoothing model outputs before comparing to data may help
                % by improving the signal:noise ratio... fitting to pattern
                % in data rather than each individual data point.
                
                smoothFactor = 0.35 * n;
                xx = smooth(CDFmod(:,2), CDFmod(:,1), smoothFactor, 'loess');                

                
                
                
                
                % For each data point, find the smoothed model equivalent,
                % and evaluate the data CDF at those points.
                matchCDF = interp1(CDFobs(:,1), CDFobs(:,2), xx, 'linear', 'extrap');
                matchCDF = min(1, max(0, matchCDF));
                % Difference between observed and modelled CDFs
                CDFdist = abs(CDFobs(:,2) - matchCDF);
                                
                
                figure
                plot(CDFobs(:,1), CDFobs(:,2))
                hold on
                plot(xx, CDFmod(:,2))
                for ij = 1:n
                    plot([xx(ij), xx(ij)], ...
                        [CDFmod(ij,2), matchCDF(ij)], 'Color', [0.85, 0.85, 0.85])
                end
                
                % Average over data points and trajectory selections
%                 avFun = @mean;
                avFun = @geomean;
                avDist = avFun(CDFdist);  % average over data points
                costComponents.(varLabel) = avFun(avDist); % average over trajectory selections
            end
            
            % Vector (size) data
            groupedByWaterOrigin = isfield(Data.size.dataBinned, 'waterMass');
            a = log(3) / log(2); % steepness of cost metric for total abundance
            for i = 1:length(VarsSize)
                varLabel = VarsSize{i};
                ind0 = strcmp(Data.size.dataBinned.Variable, varLabel);
                if groupedByWaterOrigin
                    waterMasses = unique(Data.size.dataBinned.waterMass);
                else
                    waterMasses = {[]};
                end
                for w = 1:length(waterMasses)
                    wm = waterMasses{w};
                    ind1 = ind0;
                    if groupedByWaterOrigin
                        ind1 = ind0 & strcmp(Data.size.dataBinned.waterMass, wm);
                    end
                    if groupedByWaterOrigin
                        label_autoRel = [varLabel '_' wm '_autotroph_Rel'];
                        label_autoTot = [varLabel '_' wm '_autotroph_Tot'];
                        label_heteroRel = [varLabel '_' wm '_heterotroph_Rel'];
                        label_heteroTot = [varLabel '_' wm '_heterotroph_Tot'];
                    else
                        label_autoRel = [varLabel '_autotroph_Rel'];
                        label_autoTot = [varLabel '_autotroph_Tot'];
                        label_heteroRel = [varLabel '_heterotroph_Rel'];
                        label_heteroTot = [varLabel '_heterotroph_Tot'];
                    end
                    
                    % autotrophs
                    ind = ind1 & strcmp(Data.size.dataBinned.trophicLevel, 'autotroph');
                    yobs = Data.size.dataBinned.Value(ind);
                    ymod = modData.size.Value(ind,:);
                    nsize = size(ymod, 1);
                    ymod = reshape(ymod(~isnan(ymod)), nsize, []); % remove NaNs ymod contains when fitting to Arctic AND Atlantic data
                    nsample = size(ymod, 2);
                    yobsTot = sum(yobs);
                    ymodTot = sum(ymod);
                    % Cost metric for relative abundance (vector)
                    cdfobs = cumsum(yobs) ./ yobsTot;
                    pdfobs = diff([0; cdfobs]);
                    cdfmod = cumsum(ymod) ./ ymodTot;
                    pdfmod = diff([zeros(1, nsample); cdfmod]);
                    hellingerDistance = (1 - sum((pdfobs .* pdfmod) .^ 0.5)) .^ 0.5;
                    costComponents.(label_autoRel) = mean(hellingerDistance); % average over trajectory selections
                    % Cost metric for total abundance (scalar)
                    u = abs(log(ymodTot / yobsTot));
                    u = exp(-a .* u);
                    z = (1 - u) ./ (1 + u);
                    costComponents.(label_autoTot) = mean(z); % average over trajectory selections
                    
                    % heterotrophs
                    ind = ind1 & strcmp(Data.size.dataBinned.trophicLevel, 'heterotroph');
                    yobs = Data.size.dataBinned.Value(ind);
                    ymod = modData.size.Value(ind,:);
                    nsize = size(ymod, 1);
                    ymod = reshape(ymod(~isnan(ymod)), nsize, []); % remove NaNs ymod contains when fitting to Arctic AND Atlantic data
                    nsample = size(ymod, 2);
                    yobsTot = sum(yobs);
                    ymodTot = sum(ymod);
                    % Cost metric for relative abundance (vector)
                    cdfobs = cumsum(yobs) ./ yobsTot;
                    pdfobs = diff([0; cdfobs]);
                    cdfmod = cumsum(ymod) ./ ymodTot;
                    pdfmod = diff([zeros(1, nsample); cdfmod]);
                    hellingerDistance = (1 - sum((pdfobs .* pdfmod) .^ 0.5)) .^ 0.5;
                    costComponents.(label_heteroRel) = mean(hellingerDistance);
                    % Cost metric for total abundance (scalar)
                    u = abs(log(ymodTot / yobsTot));
                    u = exp(-a .* u);
                    z = (1 - u) ./ (1 + u);
                    costComponents.(label_heteroTot) = mean(z); % average over trajectory selections
                end
            end
            
            %             % Within each size data group, weight relative abundance-at-size
            %             % relative to total abundance.
            %             weight_relVsTot = 3; % weighting factor of relative vs total abundance
            %             costSize = zeros(2,length(waterMasses)); % store weighted costs for all size data groups
            %             for i = 1:length(waterMasses)
            %                 wm = waterMasses{i};
            %                 if groupedByWaterOrigin
            %                     cta = costComponents.([varLabel '_' wm '_autotroph_Tot']);
            %                     cra = costComponents.([varLabel '_' wm '_autotroph_Rel']);
            %                     cth = costComponents.([varLabel '_' wm '_heterotroph_Tot']);
            %                     crh = costComponents.([varLabel '_' wm '_heterotroph_Rel']);
            %                 else
            %                     cta = costComponents.([varLabel '_autotroph_Tot']);
            %                     cra = costComponents.([varLabel '_autotroph_Rel']);
            %                     cth = costComponents.([varLabel '_heterotroph_Tot']);
            %                     crh = costComponents.([varLabel '_heterotroph_Rel']);
            %                 end
            %                 costSize(1,i) = mean(2 .* [1, weight_relVsTot] ./ (weight_relVsTot+1) .* [cta, cra]);
            %                 costSize(2,i) = mean(2 .* [1, weight_relVsTot] ./ (weight_relVsTot+1) .* [cth, crh]);
            %             end
            
            
            % Omit the total abundance info from the size data --
            % simplified compared to the commented code above...
            costSize = zeros(2,length(waterMasses)); % store weighted costs for all size data groups
            for i = 1:length(waterMasses)
                wm = waterMasses{i};
                if groupedByWaterOrigin
                    cra = costComponents.([varLabel '_' wm '_autotroph_Rel']);
                    crh = costComponents.([varLabel '_' wm '_heterotroph_Rel']);
                else
                    cra = costComponents.([varLabel '_autotroph_Rel']);
                    crh = costComponents.([varLabel '_heterotroph_Rel']);
                end
                costSize(1,i) = cra;
                costSize(2,i) = crh;
            end
            
            
            costNutrient = zeros(1,length(Vars));
            for i = 1:length(Vars)
                costNutrient(i) = costComponents.(Vars{i});
            end
            
            % Average the cost across data-types (nutrient & size)
            cost = [mean(costNutrient), mean(costSize(:))];
            % Assign size vs nutrient weighting
            weight_sizeVsNutrient = 1; % weighting factor of size vs nutrient data
            weights = 2 .* [1, weight_sizeVsNutrient] ./ (weight_sizeVsNutrient+1);
            cost = weights .* cost;
            
            cost = mean(cost); % finally, average cost over nutrient and size data components
            %             cost = cost(2); % try fitting only to the size data


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        case 'RMS_Hellinger'
            
            % Use root-mean-squared (absolute) errors for scalar/nutrient
            % data. Fit to smoothed RMS values for robustness against
            % overfitting to data points the model cannot replicate
            % (usually the deepest samples).
            
            % Scalar data
            for i = 1:length(Vars)
                varLabel = Vars{i};
                yobs = Data.scalar.scaled_Value(strcmp(Data.scalar.Variable, varLabel));
                ymod = modData.scalar.scaled_Value(strcmp(modData.scalar.Variable, varLabel),:);
                n = size(ymod, 1);
                
                % Sort observations
                [yobs_sort, o] = sort(yobs);
                CDFobs = [yobs_sort, ...
                    (1:n)' ./ n]; % Empirical CDF of standardised data
                
%                 figure
%                 plot(CDFobs(:,1), CDFobs(:,2)) % empirical CDF of data
                                
                % Reorder modelled values to match sorted data
                ymod_sort = ymod(o,:);
                CDFmod = [ymod_sort, CDFobs(:,2)]; % align modelled values with the smoothed CDF of the data
                
%                 hold on
%                 
%                 scatter(CDFmod(:,1), CDFmod(:,2))
%                 for ij = 1:n
%                     plot([CDFobs(ij,1), CDFmod(ij,1)], [CDFobs(ij,2), CDFmod(ij,2)], 'Color', [0.85, 0.85, 0.85])
%                 end
                
                
%                 plot(xx, CDFmod(:,2))
                
                
%                 err2 = (CDFobs(:,1) - CDFmod(:,1)) .^ 2; % squared error
                err = abs(CDFobs(:,1) - CDFmod(:,1)); % absolute error
                
%                 avFun = @mean;
                avFun = @geomean; % Errors have skewed distributions => use geometric mean for robustness against overfitting data points that the model cannot reproduce
                
%                 RMS = (avFun(err2)) .^ 0.5; % average absolute error
                RMS = avFun(err); % average absolute error
                RMS = RMS ./ range(CDFobs(:,1)); % scale to get values more in line with Hellinger distance values... this is ad-hoc method, could be improved...
                
                costComponents.(varLabel) = avFun(RMS); % average over trajectory selections

            end
            
            % Vector (size) data
            groupedByWaterOrigin = isfield(Data.size.dataBinned, 'waterMass');
            a = log(3) / log(2); % steepness of cost metric for total abundance
            for i = 1:length(VarsSize)
                varLabel = VarsSize{i};
                ind0 = strcmp(Data.size.dataBinned.Variable, varLabel);
                if groupedByWaterOrigin
                    waterMasses = unique(Data.size.dataBinned.waterMass);
                else
                    waterMasses = {[]};
                end
                for w = 1:length(waterMasses)
                    wm = waterMasses{w};
                    ind1 = ind0;
                    if groupedByWaterOrigin
                        ind1 = ind0 & strcmp(Data.size.dataBinned.waterMass, wm);
                    end
                    if groupedByWaterOrigin
                        label_autoRel = [varLabel '_' wm '_autotroph_Rel'];
                        label_autoTot = [varLabel '_' wm '_autotroph_Tot'];
                        label_heteroRel = [varLabel '_' wm '_heterotroph_Rel'];
                        label_heteroTot = [varLabel '_' wm '_heterotroph_Tot'];
                    else
                        label_autoRel = [varLabel '_autotroph_Rel'];
                        label_autoTot = [varLabel '_autotroph_Tot'];
                        label_heteroRel = [varLabel '_heterotroph_Rel'];
                        label_heteroTot = [varLabel '_heterotroph_Tot'];
                    end
                    
                    % autotrophs
                    ind = ind1 & strcmp(Data.size.dataBinned.trophicLevel, 'autotroph');
                    yobs = Data.size.dataBinned.Value(ind);
                    ymod = modData.size.Value(ind,:);
                    nsize = size(ymod, 1);
                    ymod = reshape(ymod(~isnan(ymod)), nsize, []); % remove NaNs ymod contains when fitting to Arctic AND Atlantic data
                    nsample = size(ymod, 2);
                    yobsTot = sum(yobs);
                    ymodTot = sum(ymod);
                    % Cost metric for relative abundance (vector)
                    cdfobs = cumsum(yobs) ./ yobsTot;
                    pdfobs = diff([0; cdfobs]);
                    cdfmod = cumsum(ymod) ./ ymodTot;
                    pdfmod = diff([zeros(1, nsample); cdfmod]);
                    hellingerDistance = (1 - sum((pdfobs .* pdfmod) .^ 0.5)) .^ 0.5;
                    costComponents.(label_autoRel) = mean(hellingerDistance); % average over trajectory selections
                    % Cost metric for total abundance (scalar)
                    u = abs(log(ymodTot / yobsTot));
                    u = exp(-a .* u);
                    z = (1 - u) ./ (1 + u);
                    costComponents.(label_autoTot) = mean(z); % average over trajectory selections
                    
                    % heterotrophs
                    ind = ind1 & strcmp(Data.size.dataBinned.trophicLevel, 'heterotroph');
                    yobs = Data.size.dataBinned.Value(ind);
                    ymod = modData.size.Value(ind,:);
                    nsize = size(ymod, 1);
                    ymod = reshape(ymod(~isnan(ymod)), nsize, []); % remove NaNs ymod contains when fitting to Arctic AND Atlantic data
                    nsample = size(ymod, 2);
                    yobsTot = sum(yobs);
                    ymodTot = sum(ymod);
                    % Cost metric for relative abundance (vector)
                    cdfobs = cumsum(yobs) ./ yobsTot;
                    pdfobs = diff([0; cdfobs]);
                    cdfmod = cumsum(ymod) ./ ymodTot;
                    pdfmod = diff([zeros(1, nsample); cdfmod]);
                    hellingerDistance = (1 - sum((pdfobs .* pdfmod) .^ 0.5)) .^ 0.5;
                    costComponents.(label_heteroRel) = mean(hellingerDistance);
                    % Cost metric for total abundance (scalar)
                    u = abs(log(ymodTot / yobsTot));
                    u = exp(-a .* u);
                    z = (1 - u) ./ (1 + u);
                    costComponents.(label_heteroTot) = mean(z); % average over trajectory selections
                end
            end
            
            
            % Within each size data group, weight relative abundance-at-size
            % relative to total abundance.
            weight_relVsTot = 3; % weighting factor of relative vs total abundance (relative abundance assumed more reliable)
            costSize = zeros(2,length(waterMasses)); % store weighted costs for all size data groups
            for i = 1:length(waterMasses)
                wm = waterMasses{i};
                if groupedByWaterOrigin
                    cta = costComponents.([varLabel '_' wm '_autotroph_Tot']);
                    cra = costComponents.([varLabel '_' wm '_autotroph_Rel']);
                    cth = costComponents.([varLabel '_' wm '_heterotroph_Tot']);
                    crh = costComponents.([varLabel '_' wm '_heterotroph_Rel']);
                else
                    cta = costComponents.([varLabel '_autotroph_Tot']);
                    cra = costComponents.([varLabel '_autotroph_Rel']);
                    cth = costComponents.([varLabel '_heterotroph_Tot']);
                    crh = costComponents.([varLabel '_heterotroph_Rel']);
                end
                costSize(1,i) = mean(2 .* [1, weight_relVsTot] ./ (weight_relVsTot+1) .* [cta, cra]);
                costSize(2,i) = mean(2 .* [1, weight_relVsTot] ./ (weight_relVsTot+1) .* [cth, crh]);
            end
            
            costNutrient = zeros(1,length(Vars));
            for i = 1:length(Vars)
                costNutrient(i) = costComponents.(Vars{i});
            end
            
            % Treat PON & POC data as a single type, POM, thereby
            % downweighting their combined cost contribution
            POMi = ismember(Vars, {'PON','POC'});
            costPOM = mean(costNutrient(POMi));
            costNutrient(POMi) = [];
            costNutrient = [costNutrient, costPOM];
            
            % Average the cost across data-types (nutrient & size)
            cost = [mean(costNutrient), mean(costSize(:))];
            % Assign size vs nutrient weighting
            weight_sizeVsNutrient = 1; % weighting factor of size vs nutrient data
            weights = 2 .* [1, weight_sizeVsNutrient] ./ (weight_sizeVsNutrient+1);
            cost = weights .* cost;
            
            cost = mean(cost); % finally, average cost over nutrient and size data components
            %             cost = cost(2); % try fitting only to the size data

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        case 'RMS_Hellinger2'
            
            % As above, but for the total abundance in size data fit the
            % ratio that's zooplankton. This means that the chlorophyll
            % data should inform the phytoplankton abundance while the
            % zooplankton abundance is constrained relative to
            % phytoplankton.
            
            % Scalar data
            for i = 1:length(Vars)
                varLabel = Vars{i};
                yobs = Data.scalar.scaled_Value(strcmp(Data.scalar.Variable, varLabel));
                ymod = modData.scalar.scaled_Value(strcmp(modData.scalar.Variable, varLabel),:);
                n = size(ymod, 1);
                
                % Sort observations
                [yobs_sort, o] = sort(yobs);
                CDFobs = [yobs_sort, ...
                    (1:n)' ./ n]; % Empirical CDF of standardised data
                
%                 figure
%                 plot(CDFobs(:,1), CDFobs(:,2)) % empirical CDF of data
                                
                % Reorder modelled values to match sorted data
                ymod_sort = ymod(o,:);
                CDFmod = [ymod_sort, CDFobs(:,2)]; % align modelled values with the smoothed CDF of the data
                
%                 hold on
%                 
%                 scatter(CDFmod(:,1), CDFmod(:,2))
%                 for ij = 1:n
%                     plot([CDFobs(ij,1), CDFmod(ij,1)], [CDFobs(ij,2), CDFmod(ij,2)], 'Color', [0.85, 0.85, 0.85])
%                 end
                
                
%                 plot(xx, CDFmod(:,2))
                
                
%                 err2 = (CDFobs(:,1) - CDFmod(:,1)) .^ 2; % squared error
                err = abs(CDFobs(:,1) - CDFmod(:,1)); % absolute error
                
%                 avFun = @mean;
                avFun = @geomean; % Errors have skewed distributions => use geometric mean for robustness against overfitting data points that the model cannot reproduce
                
%                 RMS = (avFun(err2)) .^ 0.5; % average absolute error
                RMS = avFun(err); % average absolute error
                RMS = RMS ./ range(CDFobs(:,1)); % scale to get values more in line with Hellinger distance values... this is ad-hoc method, could be improved...
                
                costComponents.(varLabel) = avFun(RMS); % average over trajectory selections

            end
            
            % Vector (size) data
            groupedByWaterOrigin = isfield(Data.size.dataBinned, 'waterMass');
            for i = 1:length(VarsSize)
                varLabel = VarsSize{i};
                ind0 = strcmp(Data.size.dataBinned.Variable, varLabel);
                if groupedByWaterOrigin
                    waterMasses = unique(Data.size.dataBinned.waterMass);
                else
                    waterMasses = {[]};
                end
                for w = 1:length(waterMasses)
                    wm = waterMasses{w};
                    ind1 = ind0;
                    if groupedByWaterOrigin
                        ind1 = ind0 & strcmp(Data.size.dataBinned.waterMass, wm);
                    end
                    if groupedByWaterOrigin
                        label_autoRel = [varLabel '_' wm '_autotroph_Rel'];
                        label_autoTot = [varLabel '_' wm '_autotroph_Tot'];
                        label_heteroRel = [varLabel '_' wm '_heterotroph_Rel'];
                        label_heteroTot = [varLabel '_' wm '_heterotroph_Tot'];
                    else
                        label_autoRel = [varLabel '_autotroph_Rel'];
                        label_autoTot = [varLabel '_autotroph_Tot'];
                        label_heteroRel = [varLabel '_heterotroph_Rel'];
                        label_heteroTot = [varLabel '_heterotroph_Tot'];
                    end
                    
                    % autotrophs
                    ind = ind1 & strcmp(Data.size.dataBinned.trophicLevel, 'autotroph');
                    yobs = Data.size.dataBinned.Value(ind);
                    ymod = modData.size.Value(ind,:);
                    nsize = size(ymod, 1);
                    ymod = reshape(ymod(~isnan(ymod)), nsize, []); % remove NaNs ymod contains when fitting to Arctic AND Atlantic data
                    nsample = size(ymod, 2);
                    yobsTotP = sum(yobs);
                    ymodTotP = sum(ymod);
                    % Cost metric for relative abundance (vector)
                    cdfobs = cumsum(yobs) ./ yobsTotP;
                    pdfobs = diff([0; cdfobs]);
                    cdfmod = cumsum(ymod) ./ ymodTotP;
                    pdfmod = diff([zeros(1, nsample); cdfmod]);
                    hellingerDistance = (1 - sum((pdfobs .* pdfmod) .^ 0.5)) .^ 0.5;
                    costComponents.(label_autoRel) = mean(hellingerDistance); % average over trajectory selections
%                     % Cost metric for total abundance (scalar)
%                     u = abs(log(ymodTot / yobsTot));
%                     u = exp(-a .* u);
%                     z = (1 - u) ./ (1 + u);
%                     costComponents.(label_autoTot) = mean(z); % average over trajectory selections
                    
                    % heterotrophs
                    ind = ind1 & strcmp(Data.size.dataBinned.trophicLevel, 'heterotroph');
                    yobs = Data.size.dataBinned.Value(ind);
                    ymod = modData.size.Value(ind,:);
                    nsize = size(ymod, 1);
                    ymod = reshape(ymod(~isnan(ymod)), nsize, []); % remove NaNs ymod contains when fitting to Arctic AND Atlantic data
                    nsample = size(ymod, 2);
                    yobsTotZ = sum(yobs);
                    ymodTotZ = sum(ymod);
                    % Cost metric for relative abundance (vector)
                    cdfobs = cumsum(yobs) ./ yobsTotZ;
                    pdfobs = diff([0; cdfobs]);
                    cdfmod = cumsum(ymod) ./ ymodTotZ;
                    pdfmod = diff([zeros(1, nsample); cdfmod]);
                    hellingerDistance = (1 - sum((pdfobs .* pdfmod) .^ 0.5)) .^ 0.5;
                    costComponents.(label_heteroRel) = mean(hellingerDistance);
%                     % Cost metric for total abundance (scalar)
%                     u = abs(log(ymodTot / yobsTot));
%                     u = exp(-a .* u);
%                     z = (1 - u) ./ (1 + u);
%                     costComponents.(label_heteroTot) = mean(z); % average over trajectory selections

                    % Cost value for heterotroph total biovolume
                    yobsTot = yobsTotP + yobsTotZ;
                    ymodTot = ymodTotP + ymodTotZ;
                    cost_totZ = abs((yobsTotZ ./ yobsTot) - (ymodTotZ ./ ymodTot)); % absolute difference of probabilities (bounded in (0,1))
                    
                    costComponents.(label_autoTot) = nan;
                    costComponents.(label_heteroTot) = mean(cost_totZ);
                    
                end
            end
            
            
            % Within each size data group, weight relative abundance-at-size
            % relative to total abundance.
            weight_relVsTot = 3; % weighting factor of relative vs total abundance (relative abundance assumed more reliable)
            costSize = zeros(2,length(waterMasses)); % store weighted costs for all size data groups
            for i = 1:length(waterMasses)
                wm = waterMasses{i};
                if groupedByWaterOrigin
%                     cta = costComponents.([varLabel '_' wm '_autotroph_Tot']);
                    cra = costComponents.([varLabel '_' wm '_autotroph_Rel']);
                    cth = costComponents.([varLabel '_' wm '_heterotroph_Tot']);
                    crh = costComponents.([varLabel '_' wm '_heterotroph_Rel']);
                else
%                     cta = costComponents.([varLabel '_autotroph_Tot']);
                    cra = costComponents.([varLabel '_autotroph_Rel']);
                    cth = costComponents.([varLabel '_heterotroph_Tot']);
                    crh = costComponents.([varLabel '_heterotroph_Rel']);
                end
%                 costSize(1,i) = mean(2 .* [1, weight_relVsTot] ./ (weight_relVsTot+1) .* [cta, cra]);
                costSize(1,i) = cra;
                costSize(2,i) = mean(2 .* [1, weight_relVsTot] ./ (weight_relVsTot+1) .* [cth, crh]);
            end
            
            costNutrient = zeros(1,length(Vars));
            for i = 1:length(Vars)
                costNutrient(i) = costComponents.(Vars{i});
            end
            
            % Treat PON & POC data as a single type, POM, thereby
            % downweighting their combined cost contribution
            POMi = ismember(Vars, {'PON','POC'});
            costPOM = mean(costNutrient(POMi));
            costNutrient(POMi) = [];
            costNutrient = [costNutrient, costPOM];
            
            % Average the cost across data-types (nutrient & size)
            cost = [mean(costNutrient), mean(costSize(:))];
            % Assign size vs nutrient weighting
            weight_sizeVsNutrient = 1; % weighting factor of size vs nutrient data
            weights = 2 .* [1, weight_sizeVsNutrient] ./ (weight_sizeVsNutrient+1);
            cost = weights .* cost;
            
            cost = mean(cost); % finally, average cost over nutrient and size data components
            %             cost = cost(2); % try fitting only to the size data


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        case 'RMSsmooth_Hellinger'
            
            % Use root-mean-squared (absolute) errors for scalar/nutrient
            % data. Fit to smoothed RMS values for robustness against
            % overfitting to data points the model cannot replicate
            % (usually the deepest samples).
            
            % Scalar data
            smoothFactor = 0.35; % smoothFactor multiplies number of data points to give loess smoothing parameter
            for i = 1:length(Vars)
                varLabel = Vars{i};
                yobs = Data.scalar.scaled_Value(strcmp(Data.scalar.Variable, varLabel));
                ymod = modData.scalar.scaled_Value(strcmp(modData.scalar.Variable, varLabel),:);
                n = size(ymod, 1);
                
                % Sort observations
                [yobs_sort, o] = sort(yobs);
                CDFobs = [yobs_sort, ...
                    (1:n)' ./ n]; % Empirical CDF of standardised data
                
%                 figure
%                 plot(CDFobs(:,1), CDFobs(:,2)) % empirical CDF of data
                                
                % Reorder modelled values to match sorted data
                ymod_sort = ymod(o,:);
                CDFmod = [ymod_sort, CDFobs(:,2)]; % align modelled values with the smoothed CDF of the data
                
%                 hold on
%                 
%                 scatter(CDFmod(:,1), CDFmod(:,2))
%                 for ij = 1:n
%                     plot([CDFobs(ij,1), CDFmod(ij,1)], [CDFobs(ij,2), CDFmod(ij,2)], 'Color', [0.85, 0.85, 0.85])
%                 end
                
                % Variability in model outputs may result in a
                % noisy/non-smooth cost function.
                % Smoothing model outputs before comparing to data may help
                % by improving the signal:noise ratio... fitting to pattern
                % in data rather than each individual data point.
                % Smoothing will also reduce effects of overfitting to data
                % points the model cannot replicate.
                
                smoothF = smoothFactor * n;
                xx = smooth(CDFmod(:,2), CDFmod(:,1), smoothF, 'loess');
                
%                 plot(xx, CDFmod(:,2))
                
                
%                 err2 = (CDFobs(:,1) - CDFmod(:,1)) .^ 2; % squared error
                err2 = (CDFobs(:,1) - xx) .^ 2; % squared error, smoothed
                
                avFun = @mean;
%                 avFun = @geomean;
                
                RMS = (avFun(err2)) .^ 0.5; % average absolute error
                RMS = RMS ./ range(CDFobs(:,1)); % scale to get values more in line with Hellinger distance values... this is ad-hoc method, could be improved...
                
                costComponents.(varLabel) = avFun(RMS); % average over trajectory selections

            end
            
            % Vector (size) data
            groupedByWaterOrigin = isfield(Data.size.dataBinned, 'waterMass');
            a = log(3) / log(2); % steepness of cost metric for total abundance
            for i = 1:length(VarsSize)
                varLabel = VarsSize{i};
                ind0 = strcmp(Data.size.dataBinned.Variable, varLabel);
                if groupedByWaterOrigin
                    waterMasses = unique(Data.size.dataBinned.waterMass);
                else
                    waterMasses = {[]};
                end
                for w = 1:length(waterMasses)
                    wm = waterMasses{w};
                    ind1 = ind0;
                    if groupedByWaterOrigin
                        ind1 = ind0 & strcmp(Data.size.dataBinned.waterMass, wm);
                    end
                    if groupedByWaterOrigin
                        label_autoRel = [varLabel '_' wm '_autotroph_Rel'];
                        label_autoTot = [varLabel '_' wm '_autotroph_Tot'];
                        label_heteroRel = [varLabel '_' wm '_heterotroph_Rel'];
                        label_heteroTot = [varLabel '_' wm '_heterotroph_Tot'];
                    else
                        label_autoRel = [varLabel '_autotroph_Rel'];
                        label_autoTot = [varLabel '_autotroph_Tot'];
                        label_heteroRel = [varLabel '_heterotroph_Rel'];
                        label_heteroTot = [varLabel '_heterotroph_Tot'];
                    end
                    
                    % autotrophs
                    ind = ind1 & strcmp(Data.size.dataBinned.trophicLevel, 'autotroph');
                    yobs = Data.size.dataBinned.Value(ind);
                    ymod = modData.size.Value(ind,:);
                    nsize = size(ymod, 1);
                    ymod = reshape(ymod(~isnan(ymod)), nsize, []); % remove NaNs ymod contains when fitting to Arctic AND Atlantic data
                    nsample = size(ymod, 2);
                    yobsTot = sum(yobs);
                    ymodTot = sum(ymod);
                    % Cost metric for relative abundance (vector)
                    cdfobs = cumsum(yobs) ./ yobsTot;
                    pdfobs = diff([0; cdfobs]);
                    cdfmod = cumsum(ymod) ./ ymodTot;
                    pdfmod = diff([zeros(1, nsample); cdfmod]);
                    hellingerDistance = (1 - sum((pdfobs .* pdfmod) .^ 0.5)) .^ 0.5;
                    costComponents.(label_autoRel) = mean(hellingerDistance); % average over trajectory selections
                    % Cost metric for total abundance (scalar)
                    u = abs(log(ymodTot / yobsTot));
                    u = exp(-a .* u);
                    z = (1 - u) ./ (1 + u);
                    costComponents.(label_autoTot) = mean(z); % average over trajectory selections
                    
                    % heterotrophs
                    ind = ind1 & strcmp(Data.size.dataBinned.trophicLevel, 'heterotroph');
                    yobs = Data.size.dataBinned.Value(ind);
                    ymod = modData.size.Value(ind,:);
                    nsize = size(ymod, 1);
                    ymod = reshape(ymod(~isnan(ymod)), nsize, []); % remove NaNs ymod contains when fitting to Arctic AND Atlantic data
                    nsample = size(ymod, 2);
                    yobsTot = sum(yobs);
                    ymodTot = sum(ymod);
                    % Cost metric for relative abundance (vector)
                    cdfobs = cumsum(yobs) ./ yobsTot;
                    pdfobs = diff([0; cdfobs]);
                    cdfmod = cumsum(ymod) ./ ymodTot;
                    pdfmod = diff([zeros(1, nsample); cdfmod]);
                    hellingerDistance = (1 - sum((pdfobs .* pdfmod) .^ 0.5)) .^ 0.5;
                    costComponents.(label_heteroRel) = mean(hellingerDistance);
                    % Cost metric for total abundance (scalar)
                    u = abs(log(ymodTot / yobsTot));
                    u = exp(-a .* u);
                    z = (1 - u) ./ (1 + u);
                    costComponents.(label_heteroTot) = mean(z); % average over trajectory selections
                end
            end
            
            
            % Within each size data group, weight relative abundance-at-size
            % relative to total abundance.
            weight_relVsTot = 3; % weighting factor of relative vs total abundance (relative abundance assumed more reliable)
            costSize = zeros(2,length(waterMasses)); % store weighted costs for all size data groups
            for i = 1:length(waterMasses)
                wm = waterMasses{i};
                if groupedByWaterOrigin
                    cta = costComponents.([varLabel '_' wm '_autotroph_Tot']);
                    cra = costComponents.([varLabel '_' wm '_autotroph_Rel']);
                    cth = costComponents.([varLabel '_' wm '_heterotroph_Tot']);
                    crh = costComponents.([varLabel '_' wm '_heterotroph_Rel']);
                else
                    cta = costComponents.([varLabel '_autotroph_Tot']);
                    cra = costComponents.([varLabel '_autotroph_Rel']);
                    cth = costComponents.([varLabel '_heterotroph_Tot']);
                    crh = costComponents.([varLabel '_heterotroph_Rel']);
                end
                costSize(1,i) = mean(2 .* [1, weight_relVsTot] ./ (weight_relVsTot+1) .* [cta, cra]);
                costSize(2,i) = mean(2 .* [1, weight_relVsTot] ./ (weight_relVsTot+1) .* [cth, crh]);
            end
            
            costNutrient = zeros(1,length(Vars));
            for i = 1:length(Vars)
                costNutrient(i) = costComponents.(Vars{i});
            end
            
            % Treat PON & POC data as a single type, POM, thereby
            % downweighting their combined cost contribution
            POMi = ismember(Vars, {'PON','POC'});
            costPOM = mean(costNutrient(POMi));
            costNutrient(POMi) = [];
            costNutrient = [costNutrient, costPOM];
            
            % Average the cost across data-types (nutrient & size)
            cost = [mean(costNutrient), mean(costSize(:))];
            % Assign size vs nutrient weighting
            weight_sizeVsNutrient = 1; % weighting factor of size vs nutrient data
            weights = 2 .* [1, weight_sizeVsNutrient] ./ (weight_sizeVsNutrient+1);
            cost = weights .* cost;
            
            cost = mean(cost); % finally, average cost over nutrient and size data components
            %             cost = cost(2); % try fitting only to the size data


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        case 'LeastAbsErr_Hellinger'
            
            % There are problems emerging from using CDF values to
            % calculate the scalar cost values. When modelled values are
            % either much larger or smaller than the data range, the cost
            % becomes insensitive... Probably better to use differences
            % between data and model, rather than converting to CDF
            % values... even though it might be difficult to scale for
            % compatibility with the size-data Hellinger distances.
            
            % Scalar data
            for i = 1:length(Vars)
                varLabel = Vars{i};
                yobs = Data.scalar.scaled_Value(strcmp(Data.scalar.Variable, varLabel));
                ymod = modData.scalar.scaled_Value(strcmp(modData.scalar.Variable, varLabel),:);
                n = size(ymod, 1);
                
                % Sort observations
                [yobs_sort, o] = sort(yobs);
                CDFobs = [yobs_sort, ...
                    (1:n)' ./ n]; % Empirical CDF of standardised data
                
%                 figure
%                 plot(CDFobs(:,1), CDFobs(:,2)) % empirical CDF of data
                                
                % Reorder modelled values to match sorted data
                ymod_sort = ymod(o,:);
                CDFmod = [ymod_sort, CDFobs(:,2)]; % align modelled values with the smoothed CDF of the data
                
%                 hold on
%                 
%                 scatter(CDFmod(:,1), CDFmod(:,2))
%                 for ij = 1:n
%                     plot([CDFobs(ij,1), CDFmod(ij,1)], [CDFobs(ij,2), CDFmod(ij,2)], 'Color', [0.85, 0.85, 0.85])
%                 end
                
%                 smoothFactor = 0.35 * n;
%                 xx = smooth(CDFmod(:,2), CDFmod(:,1), smoothFactor, 'loess');                
%                 plot(xx, CDFmod(:,2))

                absErr = abs(CDFobs(:,1) - CDFmod(:,1));
%                 absErr = abs(CDFobs(:,1) - xx);
                
                absErr = absErr ./ range(CDFobs(:,1)); % scale by data range to make cost values more compatible with size data cost

                % Average over data points and trajectory selections
%                 avFun = @mean;
                avFun = @geomean;
%                 avFun = @median;
                avDist = avFun(absErr);  % average over data points
                avFun = @mean;
                costComponents.(varLabel) = avFun(avDist); % average over trajectory selections
            end
            
            % Vector (size) data
            groupedByWaterOrigin = isfield(Data.size.dataBinned, 'waterMass');
            a = log(3) / log(2); % steepness of cost metric for total abundance
            for i = 1:length(VarsSize)
                varLabel = VarsSize{i};
                ind0 = strcmp(Data.size.dataBinned.Variable, varLabel);
                if groupedByWaterOrigin
                    waterMasses = unique(Data.size.dataBinned.waterMass);
                else
                    waterMasses = {[]};
                end
                for w = 1:length(waterMasses)
                    wm = waterMasses{w};
                    ind1 = ind0;
                    if groupedByWaterOrigin
                        ind1 = ind0 & strcmp(Data.size.dataBinned.waterMass, wm);
                    end
                    if groupedByWaterOrigin
                        label_autoRel = [varLabel '_' wm '_autotroph_Rel'];
                        label_autoTot = [varLabel '_' wm '_autotroph_Tot'];
                        label_heteroRel = [varLabel '_' wm '_heterotroph_Rel'];
                        label_heteroTot = [varLabel '_' wm '_heterotroph_Tot'];
                    else
                        label_autoRel = [varLabel '_autotroph_Rel'];
                        label_autoTot = [varLabel '_autotroph_Tot'];
                        label_heteroRel = [varLabel '_heterotroph_Rel'];
                        label_heteroTot = [varLabel '_heterotroph_Tot'];
                    end
                    
                    % autotrophs
                    ind = ind1 & strcmp(Data.size.dataBinned.trophicLevel, 'autotroph');
                    yobs = Data.size.dataBinned.Value(ind);
                    ymod = modData.size.Value(ind,:);
                    nsize = size(ymod, 1);
                    ymod = reshape(ymod(~isnan(ymod)), nsize, []); % remove NaNs ymod contains when fitting to Arctic AND Atlantic data
                    nsample = size(ymod, 2);
                    yobsTot = sum(yobs);
                    ymodTot = sum(ymod);
                    % Cost metric for relative abundance (vector)
                    cdfobs = cumsum(yobs) ./ yobsTot;
                    pdfobs = diff([0; cdfobs]);
                    cdfmod = cumsum(ymod) ./ ymodTot;
                    pdfmod = diff([zeros(1, nsample); cdfmod]);
                    hellingerDistance = (1 - sum((pdfobs .* pdfmod) .^ 0.5)) .^ 0.5;
                    costComponents.(label_autoRel) = mean(hellingerDistance); % average over trajectory selections
                    % Cost metric for total abundance (scalar)
                    u = abs(log(ymodTot / yobsTot));
                    u = exp(-a .* u);
                    z = (1 - u) ./ (1 + u);
                    costComponents.(label_autoTot) = mean(z); % average over trajectory selections
                    
                    % heterotrophs
                    ind = ind1 & strcmp(Data.size.dataBinned.trophicLevel, 'heterotroph');
                    yobs = Data.size.dataBinned.Value(ind);
                    ymod = modData.size.Value(ind,:);
                    nsize = size(ymod, 1);
                    ymod = reshape(ymod(~isnan(ymod)), nsize, []); % remove NaNs ymod contains when fitting to Arctic AND Atlantic data
                    nsample = size(ymod, 2);
                    yobsTot = sum(yobs);
                    ymodTot = sum(ymod);
                    % Cost metric for relative abundance (vector)
                    cdfobs = cumsum(yobs) ./ yobsTot;
                    pdfobs = diff([0; cdfobs]);
                    cdfmod = cumsum(ymod) ./ ymodTot;
                    pdfmod = diff([zeros(1, nsample); cdfmod]);
                    hellingerDistance = (1 - sum((pdfobs .* pdfmod) .^ 0.5)) .^ 0.5;
                    costComponents.(label_heteroRel) = mean(hellingerDistance);
                    % Cost metric for total abundance (scalar)
                    u = abs(log(ymodTot / yobsTot));
                    u = exp(-a .* u);
                    z = (1 - u) ./ (1 + u);
                    costComponents.(label_heteroTot) = mean(z); % average over trajectory selections
                end
            end
            
            %             % Within each size data group, weight relative abundance-at-size
            %             % relative to total abundance.
            %             weight_relVsTot = 3; % weighting factor of relative vs total abundance
            %             costSize = zeros(2,length(waterMasses)); % store weighted costs for all size data groups
            %             for i = 1:length(waterMasses)
            %                 wm = waterMasses{i};
            %                 if groupedByWaterOrigin
            %                     cta = costComponents.([varLabel '_' wm '_autotroph_Tot']);
            %                     cra = costComponents.([varLabel '_' wm '_autotroph_Rel']);
            %                     cth = costComponents.([varLabel '_' wm '_heterotroph_Tot']);
            %                     crh = costComponents.([varLabel '_' wm '_heterotroph_Rel']);
            %                 else
            %                     cta = costComponents.([varLabel '_autotroph_Tot']);
            %                     cra = costComponents.([varLabel '_autotroph_Rel']);
            %                     cth = costComponents.([varLabel '_heterotroph_Tot']);
            %                     crh = costComponents.([varLabel '_heterotroph_Rel']);
            %                 end
            %                 costSize(1,i) = mean(2 .* [1, weight_relVsTot] ./ (weight_relVsTot+1) .* [cta, cra]);
            %                 costSize(2,i) = mean(2 .* [1, weight_relVsTot] ./ (weight_relVsTot+1) .* [cth, crh]);
            %             end
            
            
            % Omit the total abundance info from the size data --
            % simplified compared to the commented code above...
            costSize = zeros(2,length(waterMasses)); % store weighted costs for all size data groups
            for i = 1:length(waterMasses)
                wm = waterMasses{i};
                if groupedByWaterOrigin
                    cra = costComponents.([varLabel '_' wm '_autotroph_Rel']);
                    crh = costComponents.([varLabel '_' wm '_heterotroph_Rel']);
                else
                    cra = costComponents.([varLabel '_autotroph_Rel']);
                    crh = costComponents.([varLabel '_heterotroph_Rel']);
                end
                costSize(1,i) = cra;
                costSize(2,i) = crh;
            end
            
            
            costNutrient = zeros(1,length(Vars));
            for i = 1:length(Vars)
                costNutrient(i) = costComponents.(Vars{i});
            end
            
            % Average the cost across data-types (nutrient & size)
            cost = [mean(costNutrient), mean(costSize(:))];
            % Assign size vs nutrient weighting
            weight_sizeVsNutrient = 1; % weighting factor of size vs nutrient data
            weights = 2 .* [1, weight_sizeVsNutrient] ./ (weight_sizeVsNutrient+1);
            cost = weights .* cost;
            
            cost = mean(cost); % finally, average cost over nutrient and size data components
            %             cost = cost(2); % try fitting only to the size data

            
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        case 'meanCDFdist_HellingerFullSpectrum'
            
            % As above, but use model output representing the full size
            % spectra rather than binned (integrated) data.
            % In fact, it is potentially far better than the above because
            % it is now much easier to separately fit all the individual
            % size spectra rather than averaging over sampling events and
            % depths... of course, we could still take averages...
            % I don't know which will be the more useful approach. Fitting
            % to each spectra separately is quite likely to produce poor
            % fits, but the average may still be reasonable. Fitting to the
            % different depth layers could also prove very beneficial by
            % much more strongly informing dependence of size structure
            % upon depth.
            
            % Build the size data cost section to separately fit all size
            % spectra -- keep the returned cost components as autotrophs &
            % heterotrophs in Arctic and/or Atlantic waters => average the
            % costs over sampling events and depths within this function.
            
            
            % Scalar data
            for i = 1:length(Vars)
                varLabel = Vars{i};
                yobs = Data.scalar.scaled_Value(strcmp(Data.scalar.Variable, varLabel));
                ymod = modData.scalar.scaled_Value(strcmp(modData.scalar.Variable, varLabel),:);
                n = size(ymod, 1);
                %                 m = size(ymod, 2);
                cdf = (1:n)' ./ n;
                [yobs_sort, o] = sort(yobs); % sort observations
                % put modelled output in same order -- these values will
                % almost certainly NOT be in ascending order...
                ymod_sort = ymod(o,:);
                %                 figure(1)
                %                 plot(yobs_sort, cdf)
                %                 hold on
                %                 plot(ymod_sort, cdf)
                %                 hold off
                % find distances between each modelled data point and the
                % empirical data CDF
                ymodi = interp1(yobs_sort, cdf, ymod_sort, 'linear', 'extrap');
                ymodi = min(1, max(0, ymodi));
                CDFdist = abs(cdf - ymodi);
                avDist = mean(CDFdist); % average over data points
                costComponents.(varLabel) = mean(avDist); % average over trajectory selections
            end
            
            % Vector (size) data
            trophicLevels = unique(Data.size.trophicLevel);
            groupedByWaterOrigin = isfield(Data.size, 'waterMass');
            if groupedByWaterOrigin
                waterMasses = unique(Data.size.waterMass);
            else
                waterMasses = {[]};
            end
            % The VarsSize labels are wrong => redefine
            switch Data.size.obsInCostFunction{1}
                case 'BioVol'; VarsSize = {'BioVolDensity'};
                case 'CellConc'; VarsSize = {'cellDensity'};
            end
            allEvents = unique(Data.size.Event, 'stable');
            nVars = length(VarsSize);
            nWaterMasses = length(waterMasses);
            nTrophicLevels = length(trophicLevels);
            HellingerDistances = nan(nVars, nWaterMasses, nTrophicLevels);
            for i = 1:length(VarsSize)
                varLabel = VarsSize{i};
                for w = 1:length(waterMasses)
                    wm = waterMasses{w};
                    ind0 = ismember(Data.size.Event, allEvents(strcmp(Data.size.waterMass, wm)));
                    Events = unique(Data.size.Event(ind0));
                    counter = 0;
                    for ie = 1:length(Events)
                        event = Events(ie);
                        ind1 = ind0 & Data.size.Event == event;
                        depths = unique(Data.size.Depth(ind1));
                        for id = 1:length(depths)
                            counter = counter + 1;
                            depth = depths(id);
                            ind2 = ind1 & Data.size.Depth == depth;
                            for it = 1:length(trophicLevels)
                                trophicLevel = trophicLevels{it};
                                ind3 = ind2 & strcmp(Data.size.trophicLevel, trophicLevel);
                                yobs = Data.size.(varLabel)(ind3);
                                ymod = modData.size.(varLabel)(ind3);
                                % Fit only to relative abundance-at-size =>
                                % normalise the data & modelled output
                                yobs = yobs ./ sum(yobs);
                                ymod = ymod ./ sum(ymod);
                                HellingerDistance = (1 - sum((yobs .* ymod) .^ 0.5)) .^ 0.5;
                                HellingerDistances(i,w,it,counter) = HellingerDistance;
                            end
                        end
                    end
                end
            end
            
            % Average Hellinger distances over sample events & depths
            HellingerDistances_ = mean(HellingerDistances, ndims(HellingerDistances));
            % Store in costComponents
            for i = 1:length(VarsSize)
                varLabel = VarsSize{i};
                for it = 1:length(trophicLevels)
                    trophicLevel = trophicLevels{it};
                    if groupedByWaterOrigin
                        for w = 1:length(waterMasses)
                            waterMass = waterMasses{w};
                            label = [varLabel '_' waterMass '_' trophicLevel];
                            costComponents.(label) = HellingerDistances_(i,w,it);
                        end
                    else
                        label = [varLabel '_' trophicLevel];
                        costComponents.(label) = HellingerDistances_(i,1,it);
                    end
                end
            end
            
            % Calculate total cost -- a scalar value
            fields = fieldnames(costComponents);
            costComponents_ = struct2array(costComponents);
            costNutrient = costComponents_(ismember(fields, Vars));
            costSize = costComponents_(~ismember(fields, Vars));
            cost_ = [mean(costNutrient), mean(costSize)];
            cost = mean(cost_);
            
            
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        case 'meanCDFdist_HellingerFullSpectrum_averagedEventsDepths'
            
            % As above, but average over all events and depths before
            % comparing model to data
            
            
            % Scalar data
            for i = 1:length(Vars)
                varLabel = Vars{i};
                yobs = Data.scalar.scaled_Value(strcmp(Data.scalar.Variable, varLabel));
                ymod = modData.scalar.scaled_Value(strcmp(modData.scalar.Variable, varLabel),:);
                n = size(ymod, 1);
                %                 m = size(ymod, 2);
                cdf = (1:n)' ./ n;
                [yobs_sort, o] = sort(yobs); % sort observations
                ymod_sort = ymod(o,:);
                % find distances between each modelled data point and the
                % empirical data CDF
                ymodi = interp1(yobs_sort, cdf, ymod_sort, 'linear', 'extrap');
                ymodi = min(1, max(0, ymodi));
                CDFdist = abs(cdf - ymodi);
                avDist = mean(CDFdist); % average over data points
                costComponents.(varLabel) = mean(avDist); % average over trajectory selections
            end
            
            % Vector (size) data
            trophicLevels = unique(Data.size.trophicLevel);
            groupedByWaterOrigin = isfield(Data.size, 'waterMass');
            if groupedByWaterOrigin
                waterMasses = unique(Data.size.waterMass);
            else
                waterMasses = {[]};
            end
            % The VarsSize labels are wrong => redefine
            switch Data.size.obsInCostFunction{1}
                case 'BioVol'; VarsSize = {'BioVolDensity'};
                case 'CellConc'; VarsSize = {'cellDensity'};
            end
            allEvents = unique(Data.size.Event, 'stable');
            nVars = length(VarsSize);
            nWaterMasses = length(waterMasses);
            nTrophicLevels = length(trophicLevels);
            nESD = length(unique(Data.size.ESD));
            
            % Extract modelled output and data into arrays with sample
            % event and depth in the trailing dimension
            yobs = nan(nVars, nWaterMasses, nTrophicLevels, nESD);
            ymod = nan(nVars, nWaterMasses, nTrophicLevels, nESD);
            for i = 1:nVars
                varLabel = VarsSize{i};
                for w = 1:nWaterMasses
                    wm = waterMasses{w};
                    ind0 = ismember(Data.size.Event, allEvents(strcmp(Data.size.waterMass, wm)));
                    Events = unique(Data.size.Event(ind0));
                    nEvents = length(Events);
                    for it = 1:nTrophicLevels
                        trophicLevel = trophicLevels{it};
                        ind1 = ind0 & strcmp(Data.size.trophicLevel, trophicLevel);
                        counter = 0;
                        for ie = 1:nEvents
                            event = Events(ie);
                            ind2 = ind1 & Data.size.Event == event;
                            Depths = unique(Data.size.Depth(ind2));
                            nDepths = length(Depths);
                            for id = 1:nDepths
                                counter = counter + 1;
                                Depth = Depths(id);
                                ind3 = ind2 & Data.size.Depth == Depth;
                                yobs_ = Data.size.(varLabel)(ind3);
                                ymod_ = modData.size.(varLabel)(ind3);
                                % Fit only to relative abundance-at-size =>
                                % normalise the data & modelled output
                                yobs_ = yobs_ ./ sum(yobs_);
                                ymod_ = ymod_ ./ sum(ymod_);
                                yobs(i,w,it,:,counter) = yobs_;
                                ymod(i,w,it,:,counter) = ymod_;
                            end
                        end
                    end
                end
            end
            % Average over sample events and depths
            yobs = mean(yobs, 5);
            ymod = mean(ymod, 5);
            % Find Hellinger distances betwen model output and data for
            % these averaged spectra, for all variables, water masses and
            % trophic levels.
            HellingerDistances = nan(nVars, nWaterMasses, nTrophicLevels);
            for i = 1:nVars
                for w = 1:nWaterMasses
                    for it = 1:nTrophicLevels
                        yobs_ = yobs(i,w,it,:);
                        ymod_ = ymod(i,w,it,:);
                        HellingerDistance = (1 - sum((yobs_ .* ymod_) .^ 0.5)) .^ 0.5;
                        HellingerDistances(i,w,it) = HellingerDistance;
                    end
                end
            end
            
            % Store in costComponents
            for i = 1:nVars
                varLabel = VarsSize{i};
                for it = 1:nTrophicLevels
                    trophicLevel = trophicLevels{it};
                    if groupedByWaterOrigin
                        for w = 1:nWaterMasses
                            waterMass = waterMasses{w};
                            label = [varLabel '_' waterMass '_' trophicLevel];
                            costComponents.(label) = HellingerDistances(i,w,it);
                        end
                    else
                        label = [varLabel '_' trophicLevel];
                        costComponents.(label) = HellingerDistances(i,1,it);
                    end
                end
            end
            
            % Calculate total cost -- a scalar value
            fields = fieldnames(costComponents);
            costComponents_ = struct2array(costComponents);
            costNutrient = costComponents_(ismember(fields, Vars));
            costSize = costComponents_(~ismember(fields, Vars));
            cost_ = [mean(costNutrient), mean(costSize)];
            cost = mean(cost_);

    end
    
end


% switch selectFunction
%     case 'LSS'
%         % Least sum of squares
%         
%         % Scalar data
%         for i = 1:length(Vars)
%             varLabel = Vars{i};            
%             yobs = Data.scalar.scaled_Value(strcmp(Data.scalar.Variable, varLabel));
%             ymod = modData.scalar.scaled_Value(strcmp(modData.scalar.Variable, varLabel),:);            
%             squaredError = (yobs - ymod) .^ 2;
%             L.(varLabel) = sum(squaredError(:)) / numel(squaredError);
%         end
%         
%         % Size spectra data
%         for i = 1:length(VarsSize)
%             varLabel = VarsSize{i};
%             ind0 = strcmp(Data.size.dataBinned.Variable, varLabel);
%             % autotrophs
%             ind = ind0 & strcmp(Data.size.dataBinned.trophicLevel, 'autotroph');
%             yobs = Data.size.dataBinned.scaled_Value(ind);
%             ymod = modData.size.scaled_Value(ind,:);
%             squaredError = (yobs - ymod) .^ 2;
%             L.([varLabel '_autotroph']) = sum(squaredError(:)) / numel(squaredError);
%             % heterotrophs
%             ind = ind0 & strcmp(Data.size.dataBinned.trophicLevel, 'heterotroph');
%             yobs = Data.size.dataBinned.scaled_Value(ind);
%             ymod = modData.size.scaled_Value(ind,:);
%             squaredError = (yobs - ymod) .^ 2;
%             L.([varLabel '_heterotroph']) = sum(squaredError(:)) / numel(squaredError);
%             L.(varLabel) = L.([varLabel '_autotroph']) + L.([varLabel '_heterotroph']);
%         end
%         
%         costComponents = L;
%         cost = 0;
%         for i = 1:length(Vars)
%             cost = cost + costComponents.(Vars{i});
%         end
% %         cost = cost / length(Vars); % group all scalar variables together, with weighting equal to the size spectra
%         for i = 1:length(VarsSize)
%             cost = cost + costComponents.(VarsSize{i});
%         end
%         
%         
%     case 'RMS'
%         % Root mean square
%         
%         % Scalar data
%         for i = 1:length(Vars)
%             varLabel = Vars{i};
%             yobs = Data.scalar.scaled_Value(strcmp(Data.scalar.Variable, varLabel));
%             ymod = modData.scalar.scaled_Value(strcmp(modData.scalar.Variable, varLabel),:);
%             absError = sqrt((yobs - ymod) .^ 2);
%             L.(varLabel) = sum(absError(:)) / numel(absError);
%         end
%         
%         % Size spectra data
%         for i = 1:length(VarsSize)
%             varLabel = VarsSize{i};
%             ind0 = strcmp(Data.size.dataBinned.Variable, varLabel);
%             % autotrophs
%             ind = ind0 & strcmp(Data.size.dataBinned.trophicLevel, 'autotroph');            
%             yobs = Data.size.dataBinned.scaled_Value(ind);
%             ymod = modData.size.scaled_Value(ind,:);
%             absError = sqrt((yobs - ymod) .^ 2);
%             L.([varLabel '_autotroph']) = sum(absError(:)) / numel(absError);
%             % heterotrophs
%             ind = ind0 & strcmp(Data.size.dataBinned.trophicLevel, 'heterotroph');            
%             yobs = Data.size.dataBinned.scaled_Value(ind);
%             ymod = modData.size.scaled_Value(ind,:);
%             absError = sqrt((yobs - ymod) .^ 2);
%             L.([varLabel '_heterotroph']) = sum(absError(:)) / numel(absError);
%             L.(varLabel) = L.([varLabel '_autotroph']) + L.([varLabel '_heterotroph']);
%         end
%         
%         costComponents = L;
%         cost = 0;
%         for i = 1:length(Vars)
%             cost = cost + costComponents.(Vars{i});
%         end
%         %         cost = cost / length(Vars); % group all scalar variables together, with weighting equal to the size spectra
%         for i = 1:length(VarsSize)
%             cost = cost + costComponents.(VarsSize{i});
%         end
%         
%     case 'syntheticLikelihood_ScalarNormal_SizeSpectraLogNormalDirichlet' % fit every data point...
%         % Model misfit to data described using a 'synthetic likelihood',
%         % as described by Wood (2010), Nature Letters, 466. doi:10.1038/nature09319
%         % Standardised scalar data are approximately normally distributed.
%         % A variety of trajectory combinations are used to generate model
%         % outputs, which are transformed identically to the data. Model
%         % outputs define normal distributions used as likelihood terms --
%         % what is likelihood of observing data given model expectations?
%         % Running the model over multiple forcing data trajectories
%         % generates the output variability (process error which is probably
%         % underestimated).
%         
%         % Scalar data
%         % Each standardised data point is assigned an independent normal 
%         % distribution parameterised using model outputs over multiple 
%         % trajectories. The likelihood is the product of probabilities of 
%         % observing all data points given these model-estimated 
%         % distributions.
% 
%         log2pi = log(2*pi);
%         for i = 1:length(Vars)
%             varLabel = Vars{i};
%             yobs = Data.scalar.scaled_Value(strcmp(Data.scalar.Variable, varLabel));
%             ymod = modData.scalar.scaled_Value(strcmp(modData.scalar.Variable, varLabel),:);
%             mu = mean(ymod, 2); % expectation of each data point
%             sig2 = var(ymod, 0, 2); % variance of model output
%             sig2(sig2 == 0) = min(sig2(sig2 > 0)); % include for robustness (we only see zero variability when using single trajectories for any sampling event)
%             n = length(yobs); % sample size
% %             L = prod(1 ./ ((2*pi*sig2) .^ 0.5) .* exp(-0.5 ./ sig2 .* (yobs - mu) .^ 2)) .^ (1/n);
%             negLogLik = 0.5 .* (log2pi + 1/n .* sum(log(sig2) + (yobs - mu) .^ 2 ./ sig2));
%             L.(varLabel) = negLogLik;
%         end
% 
%         % Size data        
%         % Dirichlet distribution on relative abundance info in size
%         % spectra. Log-normal distribution on total abundance in size data        
%         for i = 1:length(VarsSize)
%             varLabel = VarsSize{i};
%             ind0 = strcmp(Data.size.dataBinned.Variable, varLabel);
%             
%             % autotrophs
%             ind = ind0 & strcmp(Data.size.dataBinned.trophicLevel, 'autotroph');
%             yobs = Data.size.dataBinned.Value(ind);
%             ymod = modData.size.Value(ind,:);
%             % Derive Dirichlet distribution parameters from model output.
%             yobsTot = sum(yobs);
%             ymodTot = sum(ymod);
%             pobs = yobs ./ yobsTot; % observed relative abundance -- simplex
%             pmod = ymod ./ ymodTot;
%             alpha = fitDirichlet(pmod); % estimate concentration parameter
%             % Dirichlet likelihood for simplex
% %            L = gamma(alpha0) ./ prod(gamma(alpha)) .* prod(p_obs .^ (alpha-1));
%             negLogLik = sum(gammaln(alpha)) - gammaln(sum(alpha)) - sum((alpha - 1) .* log(pobs));            
%             
%             % Lognormal likelihood for total
%             yobsTot_log = log(yobsTot);
%             ymodTot_log = log(ymodTot);
%             mu = mean(ymodTot_log);
%             sig2 = var(ymodTot_log);
%             negLogLik2 = 0.5 .* (log2pi + sum(log(sig2) + 1 ./ sig2 .* (yobsTot_log - mu) .^ 2));
%             L.([varLabel '_Rel_autotroph']) = negLogLik;
%             L.([varLabel '_Tot_autotroph']) = negLogLik2;
% 
%             % heterotrophs
%             ind = ind0 & strcmp(Data.size.dataBinned.trophicLevel, 'heterotroph');
%             yobs = Data.size.dataBinned.Value(ind);
%             ymod = modData.size.Value(ind,:);
%             
%             if FixedParams.nZP_size ~= 1
%                 
%                 % Derive Dirichlet distribution parameters from model output.
%                 yobsTot = sum(yobs);
%                 ymodTot = sum(ymod);
%                 pobs = yobs ./ yobsTot; % observed relative abundance -- simplex
%                 pmod = ymod ./ ymodTot;
%                 alpha = fitDirichlet(pmod); % estimate concentration parameter
%                 % Dirichlet likelihood for simplex
%                 %            L = gamma(alpha0) ./ prod(gamma(alpha)) .* prod(p_obs .^ (alpha-1));
%                 negLogLik = sum(gammaln(alpha)) - gammaln(sum(alpha)) - sum((alpha - 1) .* log(pobs));
%                 
%                 % Lognormal likelihood for total
%                 yobsTot_log = log(yobsTot);
%                 ymodTot_log = log(ymodTot);
%                 mu = mean(ymodTot_log);
%                 sig2 = var(ymodTot_log);
%                 negLogLik2 = 0.5 .* (log2pi + sum(log(sig2) + 1 ./ sig2 .* (yobsTot_log - mu) .^ 2));
%                 L.([varLabel '_Rel_heterotroph']) = negLogLik;
%                 L.([varLabel '_Tot_heterotroph']) = negLogLik2;
%                 
%                 L.(varLabel) = L.([varLabel '_Rel_autotroph']) + L.([varLabel '_Tot_autotroph']) + ... 
%                     L.([varLabel '_Rel_heterotroph']) + L.([varLabel '_Tot_heterotroph']);
% 
%             else
%                 
%                 yobsTot = sum(yobs);
%                 ymodTot = ymod(1,:);
% 
%                 % Lognormal likelihood for total
%                 yobsTot_log = log(yobsTot);
%                 ymodTot_log = log(ymodTot);
%                 mu = mean(ymodTot_log);
%                 sig2 = var(ymodTot_log);
%                 negLogLik = 0.5 .* (log2pi + sum(log(sig2) + 1 ./ sig2 .* (yobsTot_log - mu) .^ 2));
%                 L.([varLabel '_Tot_heterotroph']) = negLogLik;
%                 
%                 L.(varLabel) = L.([varLabel '_Rel_autotroph']) + L.([varLabel '_Tot_autotroph']) + ...
%                     L.([varLabel '_Tot_heterotroph']);
% 
%             end
%             
%         end
%         
%         costComponents = L;
%         cost = 0;
%         for i = 1:length(Vars)
%             cost = cost + costComponents.(Vars{i});
%         end
% %         cost = cost / length(Vars); % group all scalar variables together, with weighting equal to the size spectra
%         for i = 1:length(VarsSize)
%             cost = cost + costComponents.(VarsSize{i});
%         end
% %         cost = sum(struct2array(costComponents));
%     
%     
%     case 'syntheticLikelihood_ScalarNormal_SizeSpectraLogNormal_logisticNormal'
%         
%         % Scalar data
%         
%         log2pi = log(2*pi);
%         for i = 1:length(Vars)
%             varLabel = Vars{i};
%             yobs = Data.scalar.scaled_Value(strcmp(Data.scalar.Variable, varLabel));
%             ymod = modData.scalar.scaled_Value(strcmp(modData.scalar.Variable, varLabel),:);
%             n = size(ymod, 1); % sample size
%             m = size(ymod, 2); % number of replicates
%             mu = mean(ymod); % modelled expectations and 
%             sig = std(ymod); % standard deviations for selected trajectories
%             sig2 = sig .^ 2;
%             negLogLik = 0.5 .* (log2pi + log(sig2) + (yobs - mu) .^ 2 ./ sig2);
%             negLogLik = sum(negLogLik(:)) / n / m;
%             L.(varLabel) = negLogLik;
%         end
%         
%         % Size data
%         
%         for i = 1:length(VarsSize)
%             varLabel = VarsSize{i};
%             ind0 = strcmp(Data.size.dataBinned.Variable, varLabel);
%             
%             % autotrophs
%             ind = ind0 & strcmp(Data.size.dataBinned.trophicLevel, 'autotroph');
%             
%             yobs = Data.size.dataBinned.Value(ind);
%             ymod = modData.size.Value_allEvents(ind,:,:); % modelled values for each separate sampling event
%             ymodMean = modData.size.Value(ind,:); % modelled values averaged over sampling events
%             % dimension of ymod is [size, sampling event, replicate (trajectory choice), sample number (year or cruise)]
%             nsize = size(ymod, 1); % number of sizes
% %             ne = size(ymod, 2); % number of sampling events
%             m = size(ymod, 3); % number of replicates (trajectory choices)
% %             n = size(ymod, 4); % number of samples/years/cruises etc
%             
%             % decompose spectra into total and relative abundance (a single number and a simplex)
%             yobsTot = sum(yobs); % totals
%             ymodTot = sum(ymod);
%             ymodMeanTot = sum(ymodMean);
%             pobs = yobs ./ yobsTot; % simplices
%             pmod = ymod ./ ymodTot;
%             pmodMean = ymodMean ./ ymodMeanTot;
%             
%             % Use a logistic-normal distribution to calculate likelihood of
%             % observed relative abundance (simplex) given the modelled
%             % equivalents.
%             % Methods detailed in: Francis, R.I.C.C. (2014) Fish.Res.151:70-84. doi:10.1016/j.fishres.2013.12.015
%             
%             % Assuming observed simplex, pobs, has a logistic-normal
%             % distribution with expectation, pmod => X is multi-variate normal
%             % if pobs = exp(X) / sum(exp(X)), where X has expectation, log(pmod),
%             % and covariance, C.
%             % The pobs-X transform is not one-to-one, so transform pobs
%             % by reducing the size dimension by one.
%             Yobs = log(pobs(1:end-1) ./ pobs(end));
%             % Now assume that Yobs is multi-variate normal with
%             % expectation, mu, and covariance, V.
% %             mu = log(pmod(1:end-1,:,:) ./ pmod(end,:,:)); % expectations for [nsize-1] multi-variate normal Yobs (all sampling events)
%             Mu = log(pmodMean(1:end-1,:) ./ pmodMean(end,:)); % expectations for [nsize-1] multi-variate normal Yobs (averaged sampling events)
% 
%             log_pmod = log(pmod); % expectations (for each separate sampling event) of [nsize] multivariate-normal distribution for X
%             
%             K = [eye(nsize-1) -ones(nsize-1,1)];
%             
%             % Use modelled output for each separate sampling event to
%             % estimate (co)variance parameter for the event-averaged data
%             
%             C = nan(nsize, nsize, m); % Covaraince of multi-variate normal distribution
%             V = nan(nsize-1, nsize-1, m); % Transformed covariance -- useful form for likelihood
%             
%             covarianceTypes = {'var', 'CVtridiag', 'CVfull'}; % variance, tri-diagonal covariance matrix, full covariance matrix
%             covarianceType = covarianceTypes{1};
%             
%             switch covarianceType
%                 case 'var'
%                     for j = 1:m
%                         C(:,:,j) = diag(diag(cov(log_pmod(:,:,j)')));
%                         V(:,:,j) = K * C(:,:,j) * K';
%                     end
%                 case 'CVfull'
%                     for j = 1:m
%                         C(:,:,j) = cov(log_pmod(:,:,j)');
%                         V(:,:,j) = K * C(:,:,j) * K';
%                     end
%                 case 'CVtridiag'
%                     for j = 1:m
%                         C(:,:,j) = cov(log_pmod(:,:,j)');
%                         C(:,:,j) = diag(diag(C(:,:,j))) + ... 
%                             diag(diag(C(:,:,j), -1), -1) + ... 
%                             diag(diag(C(:,:,j), 1), 1);
%                         V(:,:,j) = K * C(:,:,j) * K';
%                     end
%             end
%             
%             w = Yobs - Mu; % error vectors
%             
%             for j = 1:m
%                 negLogLik(j) = log(det(V(:,:,j))) + (w(:,j)' / V(:,:,j) * w(:,j));
%             end
% %             negLogLik = 0.5 .* (negLogLik + (nsize - 1) .* log2pi) + sum(log(pobs));
%             negLogLik = 0.5 .* (negLogLik + (nsize - 1) .* log2pi);
%             negLogLik_rel = sum(negLogLik) / m;
% %             negLogLik_rel = sum(negLogLik) / m / (nsize-1);
%             
%             % Lognormal likelihood for total abundance
%             yobsTot_log = log(yobsTot);
%             ymodTot_log = log(ymodTot);
%             mu = squeeze(mean(ymodTot_log));
%             sig = squeeze(std(ymodTot_log));
%             sig2 = sig .^ 2;
%             
% %             Lik = prod(1 ./ ((2*pi*sig2) .^ 0.5) .* exp(-0.5 .* (yobsTot_log - mu) .^ 2 ./ sig2)) .^ (1/m)            
%             negLogLik_tot = 0.5 .* (m .* log2pi + sum(log(sig2) + (yobsTot_log - mu) .^ 2 ./ sig2)) ./ m;
%             
%             L.([varLabel '_Rel_autotroph']) = negLogLik_rel;
%             L.([varLabel '_Tot_autotroph']) = negLogLik_tot;
%             L.([varLabel '_autotroph']) = negLogLik_rel + negLogLik_tot;
% 
%             
%             % heterotrophs
%             ind = ind0 & strcmp(Data.size.dataBinned.trophicLevel, 'heterotroph');
%             
%             yobs = Data.size.dataBinned.Value(ind);
%             ymod = modData.size.Value_allEvents(ind,:,:); % modelled values for each separate sampling event
% %             ymodMean = modData.size.Value(ind,:); % modelled values averaged over sampling events
%             
%             singleSizeClass = FixedParams.nZP_size == 1;
%             
%             if singleSizeClass
%                 ymod = ymod(1,:,:);
% %                 ymodMean = ymodMean(1,:);
%                 yobs = sum(yobs); % sum observations over size classes
%                 
%                 % dimension of ymod is [size, sampling event, replicate (trajectory choice), sample number (year or cruise)]
%                 m = size(ymod, 3); % number of replicates (trajectory choices)
%                 
%                 % Lognormal likelihood for total abundance
%                 yobsTot_log = log(yobs);
%                 ymodTot_log = log(ymod);
%                 mu = squeeze(mean(ymodTot_log));
%                 sig2 = squeeze(var(ymodTot_log));
%                 
%                 negLogLik_tot = 0.5 .* (m .* log2pi + sum(log(sig2) + (yobsTot_log - mu) .^ 2 ./ sig2)) ./ m;
% 
%                 L.([varLabel '_heterotroph']) = negLogLik_tot;
%                 
%             else
%                 
%                 warning('still need to set up the multiple size class model!')
%                 
%                 
%             end
%             
%             L.(varLabel) = L.([varLabel '_autotroph']) + L.([varLabel '_heterotroph']);
%             
%         end
%         
%         costComponents = L;
%         cost = 0;
%         for i = 1:length(Vars)
%             cost = cost + costComponents.(Vars{i});
%         end
% %         cost = cost / length(Vars); % group all scalar variables together, with weighting equal to the size spectra
%         for i = 1:length(VarsSize)
%             cost = cost + costComponents.(VarsSize{i});
%         end
% %         cost = sum(struct2array(costComponents));
% 
% 
% 
%     case 'syntheticLikelihood_ScalarNormalShape_SizeSpectraLogNormalDirichlet'
%         % Model misfit to data described using a 'synthetic likelihood',
%         % as described by Wood (2010), Nature Letters, 466. doi:10.1038/nature09319
%         % Standardised scalar data are approximately normally distributed.
%         % A variety of trajectory combinations are used to generate model
%         % outputs, which are transformed identically to the data. Model
%         % outputs define normal distributions used as likelihood terms --
%         % what is likelihood of observing data given model expectations?
%         % Running the model over multiple forcing data trajectories
%         % generates the output variability (process error which is probably
%         % underestimated).
%         
%         % Scalar data
%         % Each standardised data point is assigned an independent normal 
%         % distribution parameterised using model outputs over multiple 
%         % trajectories. The likelihood is the product of probabilities of 
%         % observing all data points given these model-estimated 
%         % distributions.
%         log2pi = log(2*pi);
%         for i = 1:length(Vars)
%             varLabel = Vars{i};
%             yobs = Data.scalar.scaled_Value(strcmp(Data.scalar.Variable, varLabel));
%             ymod = modData.scalar.scaled_Value(strcmp(modData.scalar.Variable, varLabel),:);
% %             n = length(yobs);
%             n = size(ymod, 2);
%             mu_obs = mean(yobs);
%             v_obs = var(yobs);
%             Mu = mean(ymod);
%             V = var(ymod);
%             % likelihood of mu_obs
%             mu = mean(Mu);
%             v = var(Mu);
%             negLogLik_mu = log2pi + log(v) + (mu_obs - mu) .^ 2 ./ v;
%             % likelihood of v_obs
%             mu = mean(V);
%             v = var(V);
%             negLogLik_v = log2pi + log(v) + (v_obs - mu) .^ 2 ./ v;
%             % combine
%             negLogLik = 0.5 * n * (negLogLik_mu + negLogLik_v);
%             L.(varLabel) = negLogLik;
%         end
%         
% 
%         % Size data        
%         % Dirichlet distribution on relative abundance info in size
%         % spectra. Log-normal distribution on total abundance in size data        
%         for i = 1:length(VarsSize)
%             varLabel = VarsSize{i};
%             yobs = Data.size.dataBinned.Value(strcmp(Data.size.dataBinned.Variable, varLabel));
%             ymod = modData.size.Value(strcmp(modData.size.Variable, varLabel),:);
%             % Derive Dirichlet distribution parameters from model output.
%             yobsTot = sum(yobs);
%             ymodTot = sum(ymod);
%             pobs = yobs ./ yobsTot; % observed relative abundance -- simplex
%             pmod = ymod ./ ymodTot;
%             alpha = fitDirichlet(pmod); % estimate concentration parameter
%             % Dirichlet likelihood for simplex
% %            L = gamma(alpha0) ./ prod(gamma(alpha)) .* prod(p_obs .^ (alpha-1));
%             negLogLik = sum(gammaln(alpha)) - gammaln(sum(alpha)) - sum((alpha - 1) .* log(pobs));            
%             
%             % Lognormal likelihood for total
%             yobsTot_log = log(yobsTot);
%             ymodTot_log = log(ymodTot);
%             mu = mean(ymodTot_log);
%             sig = std(ymodTot_log);
%             sig2 = sig .^ 2;
%             negLogLik2 = 0.5 .* (log2pi + sum(log(sig2) + 1 ./ sig2 .* (yobsTot_log - mu) .^ 2));
%             L.([varLabel '_Rel']) = negLogLik;
%             L.([varLabel '_Tot']) = negLogLik2;
%             L.(varLabel) = negLogLik + negLogLik2;
%         end
%         
%         costComponents = L;
%         cost = 0;
%         for i = 1:length(Vars)
%             cost = cost + costComponents.(Vars{i});
%         end
% %         cost = cost / length(Vars); % group all scalar variables together, with weighting equal to the size spectra
%         for i = 1:length(VarsSize)
%             cost = cost + costComponents.(VarsSize{i});
%         end
% %         cost = sum(struct2array(costComponents));
% 
% 
%     case 'N_LN-Dir_groupWaterOrigin'
%         % Synthetic likelihoods where variability parameters are estimated
%         % using variability in model outputs from the various trajectories.
%         % Normal distributions for scalar data points.
%         % Size spectra vectors decomposed into totals and simplexes --
%         % lognormal distributions for the totals, Dirichlet distributions
%         % for the simplexes.
%         % Size spectra data for autotrophs and heterotrophs, and for water 
%         % of Arctic and of Atlantic origin -- 4 separate data vectors.
%         
%         % Scalar data
%         log2pi = log(2*pi);
%         for i = 1:length(Vars)
%             varLabel = Vars{i};
%             yobs = Data.scalar.scaled_Value(strcmp(Data.scalar.Variable, varLabel));
%             ymod = modData.scalar.scaled_Value(strcmp(modData.scalar.Variable, varLabel),:);
%             mu = mean(ymod, 2); % expectation of each data point
%             sig2 = var(ymod, 0, 2); % variance of model output
%             sig2(sig2 == 0) = min(sig2(sig2 > 0)); % include for robustness (we only see zero variability when using single trajectories for any sampling event)
%             n = length(yobs); % sample size
%             %             L = prod(1 ./ ((2*pi*sig2) .^ 0.5) .* exp(-0.5 ./ sig2 .* (yobs - mu) .^ 2)) .^ (1/n);
%             negLogLik = 0.5 .* (log2pi + 1/n .* sum(log(sig2) + (yobs - mu) .^ 2 ./ sig2));
%             L.(varLabel) = negLogLik;
%         end
%         
%         % Vector (size) data
%         for i = 1:length(VarsSize)
%             varLabel = VarsSize{i};
%             ind0 = strcmp(Data.sizeFull.dataBinned.groupedByOrigin.Variable, varLabel);
%             waterMasses = unique(Data.sizeFull.dataBinned.groupedByOrigin.waterMass);
%             
%             for w = 1:length(waterMasses)
%                 wm = waterMasses{w};
%                 ind1 = ind0 & strcmp(Data.sizeFull.dataBinned.groupedByOrigin.waterMass, wm);
%                 
%                 % autotrophs
%                 ind = ind1 & strcmp(Data.sizeFull.dataBinned.groupedByOrigin.trophicLevel, 'autotroph');
%                 yobs = Data.sizeFull.dataBinned.groupedByOrigin.Value(ind);
%                 ymod = modData.sizeFull.(['Value_' wm])(ind,:);
%                 % Derive Dirichlet distribution parameters from model output.
%                 yobsTot = sum(yobs);
%                 ymodTot = sum(ymod);
%                 pobs = yobs ./ yobsTot; % observed relative abundance -- simplex
%                 pmod = ymod ./ ymodTot;
%                 alpha = fitDirichlet(pmod); % estimate concentration parameter
%                 % Dirichlet likelihood for simplex
%                 negLogLik = sum(gammaln(alpha)) - gammaln(sum(alpha)) - sum((alpha - 1) .* log(pobs));
% 
%                 % Lognormal likelihood for total
%                 yobsTot_log = log(yobsTot);
%                 ymodTot_log = log(ymodTot);
%                 mu = mean(ymodTot_log);
%                 sig2 = var(ymodTot_log);
%                 negLogLik2 = 0.5 .* (log2pi + sum(log(sig2) + 1 ./ sig2 .* (yobsTot_log - mu) .^ 2));
%                 
%                 L.([varLabel '_' wm '_autotroph_Rel']) = negLogLik;
%                 L.([varLabel '_' wm '_autotroph_Tot']) = negLogLik2;
%                 
%                 % heterotrophs
%                 ind = ind1 & strcmp(Data.sizeFull.dataBinned.groupedByOrigin.trophicLevel, 'heterotroph');
%                 yobs = Data.sizeFull.dataBinned.groupedByOrigin.Value(ind);
%                 ymod = modData.sizeFull.(['Value_' wm])(ind,:);
%                 
%                 % Derive Dirichlet distribution parameters from model output.
%                 yobsTot = sum(yobs);
%                 ymodTot = sum(ymod);
%                 pobs = yobs ./ yobsTot; % observed relative abundance -- simplex
%                 pmod = ymod ./ ymodTot;
%                 alpha = fitDirichlet(pmod); % estimate concentration parameter
%                 % Dirichlet likelihood for simplex
%                 negLogLik = sum(gammaln(alpha)) - gammaln(sum(alpha)) - sum((alpha - 1) .* log(pobs));
%                 
%                 % Lognormal likelihood for total
%                 yobsTot_log = log(yobsTot);
%                 ymodTot_log = log(ymodTot);
%                 mu = mean(ymodTot_log);
%                 sig2 = var(ymodTot_log);
%                 negLogLik2 = 0.5 .* (log2pi + sum(log(sig2) + 1 ./ sig2 .* (yobsTot_log - mu) .^ 2));
%                 
%                 L.([varLabel '_' wm '_heterotroph_Rel']) = negLogLik;
%                 L.([varLabel '_' wm '_heterotroph_Tot']) = negLogLik2;
%             end
%         end
%         
%         costComponents = L;
%         
%         % Apply any weightings to data types here, before summing
%         % costComponents to find total cost...
%         cost = sum(struct2array(costComponents));
% 
%     
%     case 'Hellinger_groupWaterOrigin'
%         % Calulate Hellinger distances between observations and their 
%         % equivalent model outputs.
%         % For consistency/comparability use this metric for all data types.
%         
%         % Scalar data
%         for i = 1:length(Vars)
%             varLabel = Vars{i};
%             yobs = Data.scalar.scaled_Value(strcmp(Data.scalar.Variable, varLabel));
%             ymod = modData.scalar.scaled_Value(strcmp(modData.scalar.Variable, varLabel),:);
%             n = size(ymod, 1);
%             m = size(ymod, 2);
%             % As the standardised data are approximately normally distributed,
%             % we could assume that identically transformed model outputs
%             % are also normally distributed, then find their means and
%             % variances to calculate Hellinger distances analytically.
%             % However, assuming normality of transformed model output is
%             % inappropriate because it could be any shape...
%             % Thus, derive cdfs and pdfs directly without assuming that
%             % transformed model outputs follow normal distributions.
%             % To compare observed and modelled distributions, the pdfs will
%             % need to be evaluated across identical domains.
%             cdf = (1:n)' ./ n;
%             yobsc = sort(yobs);                        
%             ymodc = sort(ymod);            
%             % remove any duplicate values (more likely in model output than
%             % in data), retaining the largest probability values
%             keep = ~[diff(yobsc) == 0; false];
%             cdfobs = cdf(keep);
%             yobsc = yobsc(keep);
%             % store cdfmod as cell-array because vector lengths may differ
%             % after removing duplicates
%             keep = ~[diff(ymodc) == 0; false(1, m)];
%             cdfmod = cell(1, m);
%             ymodc_ = cell(1, m);
%             for ij = 1:m
%                 cdfmod{:,ij} = cdf(keep(:,ij));
%                 ymodc_{:,ij} = ymodc(keep(:,ij),ij);
%             end
%             % define regular grid across measurement space -- define number
%             % of grid nodes using the number of observations
%             vrange = [min([yobsc(:); ymodc(:)]), max([yobsc(:); ymodc(:)])];
%             grid = nan(n, 1);
%             grid(2:end) = linspace(vrange(1), vrange(2), n-1)';
%             sp = diff(grid(2:3));
%             grid(1) = grid(2) - sp;
%             % interpolate cdfs and derive pdfs
%             cdfobsi = interp1(yobsc, cdfobs, grid, 'linear', 'extrap');
%             cdfobsi(cdfobsi < min(cdfobs)) = 0;
%             cdfobsi(cdfobsi > 1) = 1;
%             pdfobsi = diff(cdfobsi);            
%             cdfmodi = nan(n, m);
%             pdfmodi = nan(n-1, m);
%             for ij = 1:m
%                 cdfmodi(:,ij) = interp1(ymodc_{ij}, cdfmod{ij}, grid, 'linear', 'extrap');
%                 cdfmodi(cdfmodi(:,ij) < min(cdfmod{ij}), ij) = 0;
%                 cdfmodi(cdfmodi(:,ij) > 1, ij) = 1;
%                 pdfmodi(:,ij) = diff(cdfmodi(:,ij));
%             end
%                         
%             hellingerDistance = (1 - sum((pdfobsi .* pdfmodi) .^ 0.5)) .^ 0.5;
%             
%             costComponents.(varLabel) = mean(hellingerDistance); % average over trajectory selections
%             
%         end
%         
%         % Vector (size) data
%         % Hellinger distance is appropriate for relative abundance-at-size,
%         % but I'm not sure that it can easily be used for absolute 
%         % abundance at size...
%         % The pmf (probability mass function) for relative abundance is
%         % simply a vector equalling the relative abundances.
%         % However, a pmf for absolute abundance at size would need to be
%         % derived by assuming the data follow some distribution (e.g.
%         % multivariate-normal). This complicates using Hellinger distance
%         % for absolute abundance-at-size because the necessary(?)
%         % distributional assumptions require estimation of variability
%         % parameters for BOTH the observed AND modelled values...
%         % We could calculate measured abundance-at-size variability in
%         % Arctic and Atlantic waters using data from all sampling events,
%         % and variability in modelled values could be calculated across
%         % trajectories... then, by assuming that abundance-at-size is
%         % multivariate normal, we could derive probabilities for absolute
%         % values.        
%         for i = 1:length(VarsSize)
%             varLabel = VarsSize{i};
%             ind0 = strcmp(Data.sizeFull.dataBinned.groupedByOrigin.Variable, varLabel);
%             waterMasses = unique(Data.sizeFull.dataBinned.groupedByOrigin.waterMass);
%             
%             for w = 1:length(waterMasses)
%                 wm = waterMasses{w};
%                 ind1 = ind0 & strcmp(Data.sizeFull.dataBinned.groupedByOrigin.waterMass, wm);
%                 
%                 % autotrophs
%                 ind = ind1 & strcmp(Data.sizeFull.dataBinned.groupedByOrigin.trophicLevel, 'autotroph');
%                 yobs = Data.sizeFull.dataBinned.groupedByOrigin.Value(ind);
%                 ymod = modData.sizeFull.(['Value_' wm])(ind,:);
% %                 n = size(ymod, 1);
%                 m = size(ymod, 2);
%                 cdfobs = cumsum(yobs) ./ sum(yobs);
%                 pdfobs = diff([0; cdfobs]);
%                 cdfmod = cumsum(ymod) ./ sum(ymod);
%                 pdfmod = diff([zeros(1, m); cdfmod]);
%                 hellingerDistance = (1 - sum((pdfobs .* pdfmod) .^ 0.5)) .^ 0.5;
%                 costComponents.([varLabel '_' wm '_autotroph']) = mean(hellingerDistance); % average over trajectory selections
% 
%                 % heterotrophs
%                 ind = ind1 & strcmp(Data.sizeFull.dataBinned.groupedByOrigin.trophicLevel, 'heterotroph');
%                 yobs = Data.sizeFull.dataBinned.groupedByOrigin.Value(ind);
%                 ymod = modData.sizeFull.(['Value_' wm])(ind,:);                
% %                 n = size(ymod, 1);
%                 m = size(ymod, 2);
%                 cdfobs = cumsum(yobs) ./ sum(yobs);
%                 pdfobs = diff([0; cdfobs]);
%                 cdfmod = cumsum(ymod) ./ sum(ymod);
%                 pdfmod = diff([zeros(1, m); cdfmod]);
%                 hellingerDistance = (1 - sum((pdfobs .* pdfmod) .^ 0.5)) .^ 0.5;
%                 costComponents.([varLabel '_' wm '_heterotroph']) = mean(hellingerDistance);
% 
%             end
%         end
%         
%         % Take averages across data-types to assign equal weightings to the
%         % scalar data and size data.
%         cost = zeros(1,2);
%         fields = fieldnames(costComponents);
%         for i = 1:length(fields)
%             if ismember(fields{i}, Vars)
%                 % scalar data
%                cost(1) = cost(1) + costComponents.(fields{i});
%             else
%                 % size data
%                 cost(2) = cost(2) + costComponents.(fields{i});
%             end
%         end
%         cost(1) = cost(1) ./ length(Vars);
%         cost(2) = cost(2) ./ (length(fields) - length(Vars));
%         cost = mean(cost);
%         
% %         cost = sum(struct2array(costComponents));
% 
% %         disp(costComponents)
% %         disp(cost)
% 
%     case 'Hellinger2_groupWaterOrigin'
%         % Calculate non-parameteric Hellinger distances for the scalar
%         % (nutrient) data and for size spectra vectors of relative
%         % abundance.
%         % Use a different, but comparable, metric for the total abundance
%         % data contained in the size spectra. This is necessary because
%         % these data do not form distributions so Hellinger distance is not
%         % applicable. Like the Hellinger distance, the metric takes values
%         % in the in the interval [0,1], with 0 represetning a perfect fit 
%         % and 1 being an extremely poor fit.
%         % This method makes zero assumptions about variability -- unlike
%         % the other methods which use variability across particle
%         % trajectories to create distributions. The Hellinger distances are
%         % calculated numerically instead of using formulas based on
%         % assumptions of Gaussian distributions. The scalar data are
%         % continuous and approximately Gaussian whereas the size spectra
%         % data are discrete, so Hellinger distance calculations are
%         % different for each...
%         
%         % Scalar data
%         for i = 1:length(Vars)
%             varLabel = Vars{i};
%             yobs = Data.scalar.scaled_Value(strcmp(Data.scalar.Variable, varLabel));
%             ymod = modData.scalar.scaled_Value(strcmp(modData.scalar.Variable, varLabel),:);
%             n = size(ymod, 1);
%             m = size(ymod, 2);
%             % As the standardised data are approximately normally distributed,
%             % we could assume that identically transformed model outputs
%             % are also normally distributed, then find their means and
%             % variances to calculate Hellinger distances analytically.
%             % However, assuming normality of transformed model output is
%             % inappropriate because it could be any shape...
%             % Thus, derive cdfs and pdfs directly without assuming that
%             % transformed model outputs follow normal distributions.
%             % To compare observed and modelled distributions, the pdfs will
%             % need to be evaluated across identical domains.
%             cdf = (1:n)' ./ n;
%             yobsc = sort(yobs); % sort observations
%             ymodc = sort(ymod);
%             % remove any duplicate values (more likely in model output than
%             % in data), retaining the largest probability values
%             keep = ~[diff(yobsc) == 0; false];
%             cdfobs = cdf(keep);
%             yobsc = yobsc(keep);
%             % store cdfmod as cell-array because vector lengths may differ
%             % after removing duplicates
%             keep = ~[diff(ymodc) == 0; false(1, m)];
%             cdfmod = cell(1, m);
%             ymodc_ = cell(1, m);
%             for ij = 1:m
%                 cdfmod{:,ij} = cdf(keep(:,ij));
%                 ymodc_{:,ij} = ymodc(keep(:,ij),ij);
%             end
%             % define regular grid across measurement space -- define number
%             % of grid nodes, ng, relative to number of observations, n.
%             vrange = [min([yobsc(:); ymodc(:)]), max([yobsc(:); ymodc(:)])];
%             nr = 1; % choice of grid resolution influences cost values... is there a better way to do this, or a principled way to choose grid resolution?
%             ng = ceil(nr * n);
%             grid = nan(ng, 1);
%             grid(2:end) = linspace(vrange(1), vrange(2), ng-1)';
%             sp = diff(grid(2:3));
%             grid(1) = grid(2) - sp;
%             % interpolate cdfs and derive pdfs
%             cdfobsi = interp1(yobsc, cdfobs, grid, 'linear', 'extrap');
%             cdfobsi(cdfobsi < min(cdfobs)) = 0;
%             cdfobsi(cdfobsi > 1) = 1;
%             pdfobsi = diff(cdfobsi);
%             cdfmodi = nan(ng, m);
%             pdfmodi = nan(ng-1, m);
%             for ij = 1:m
%                 cdfmodi(:,ij) = interp1(ymodc_{ij}, cdfmod{ij}, grid, 'linear', 'extrap');
%                 cdfmodi(cdfmodi(:,ij) < min(cdfmod{ij}), ij) = 0;
%                 cdfmodi(cdfmodi(:,ij) > 1, ij) = 1;
%                 pdfmodi(:,ij) = diff(cdfmodi(:,ij));
%             end
%             hellingerDistance = (1 - sum((pdfobsi .* pdfmodi) .^ 0.5)) .^ 0.5;
%             costComponents.(varLabel) = mean(hellingerDistance); % average over trajectory selections
%         end
%         
%         % Vector (size) data
%         groupedByWaterOrigin = isfield(Data.size.dataBinned, 'waterMass');
%         a = log(3) / log(2); % steepness of cost metric for total abundance
%         for i = 1:length(VarsSize)
%             varLabel = VarsSize{i};
%             ind0 = strcmp(Data.size.dataBinned.Variable, varLabel);
%             if groupedByWaterOrigin
%                 waterMasses = unique(Data.size.dataBinned.waterMass);
%             else
%                 waterMasses = {[]};
%             end
%             for w = 1:length(waterMasses)
%                 wm = waterMasses{w};
%                 ind1 = ind0;
%                 if groupedByWaterOrigin
%                     ind1 = ind0 & strcmp(Data.size.dataBinned.waterMass, wm);
%                 end
%                 if groupedByWaterOrigin
%                     label_autoRel = [varLabel '_' wm '_autotroph_Rel'];
%                     label_autoTot = [varLabel '_' wm '_autotroph_Tot'];
%                     label_heteroRel = [varLabel '_' wm '_heterotroph_Rel'];
%                     label_heteroTot = [varLabel '_' wm '_heterotroph_Tot'];
%                 else
%                     label_autoRel = [varLabel '_autotroph_Rel'];
%                     label_autoTot = [varLabel '_autotroph_Tot'];
%                     label_heteroRel = [varLabel '_heterotroph_Rel'];
%                     label_heteroTot = [varLabel '_heterotroph_Tot'];
%                 end
%                 
%                 % autotrophs
%                 ind = ind1 & strcmp(Data.size.dataBinned.trophicLevel, 'autotroph');
%                 yobs = Data.size.dataBinned.Value(ind);
%                 ymod = modData.size.Value(ind,:);
%                 nsize = size(ymod, 1);
%                 ymod = reshape(ymod(~isnan(ymod)), nsize, []); % remove NaNs ymod contains when fitting to Arctic AND Atlantic data
%                 nsample = size(ymod, 2);
%                 yobsTot = sum(yobs);
%                 ymodTot = sum(ymod);
%                 % Cost metric for relative abundance (vector)
%                 cdfobs = cumsum(yobs) ./ yobsTot;
%                 pdfobs = diff([0; cdfobs]);
%                 cdfmod = cumsum(ymod) ./ ymodTot;
%                 pdfmod = diff([zeros(1, nsample); cdfmod]);
%                 hellingerDistance = (1 - sum((pdfobs .* pdfmod) .^ 0.5)) .^ 0.5;
%                 costComponents.(label_autoRel) = mean(hellingerDistance); % average over trajectory selections
%                 % Cost metric for total abundance (scalar)
%                 u = abs(log(ymodTot / yobsTot));
%                 u = exp(-a .* u);
%                 z = (1 - u) ./ (1 + u);
%                 costComponents.(label_autoTot) = mean(z); % average over trajectory selections
%                 
%                 % heterotrophs
%                 ind = ind1 & strcmp(Data.size.dataBinned.trophicLevel, 'heterotroph');
%                 yobs = Data.size.dataBinned.Value(ind);
%                 ymod = modData.size.Value(ind,:);
%                 nsize = size(ymod, 1);
%                 ymod = reshape(ymod(~isnan(ymod)), nsize, []); % remove NaNs ymod contains when fitting to Arctic AND Atlantic data
%                 nsample = size(ymod, 2);
%                 yobsTot = sum(yobs);
%                 ymodTot = sum(ymod);
%                 % Cost metric for relative abundance (vector)
%                 cdfobs = cumsum(yobs) ./ yobsTot;
%                 pdfobs = diff([0; cdfobs]);
%                 cdfmod = cumsum(ymod) ./ ymodTot;
%                 pdfmod = diff([zeros(1, nsample); cdfmod]);
%                 hellingerDistance = (1 - sum((pdfobs .* pdfmod) .^ 0.5)) .^ 0.5;
%                 costComponents.(label_heteroRel) = mean(hellingerDistance);
%                 % Cost metric for total abundance (scalar)
%                 u = abs(log(ymodTot / yobsTot));
%                 u = exp(-a .* u);
%                 z = (1 - u) ./ (1 + u);
%                 costComponents.(label_heteroTot) = mean(z); % average over trajectory selections
%             end
%         end
%         
%         % Within each size data group, weight relative abundance-at-size
%         % relative to total abundance.
%         weight_relVsTot = 3; % weighting factor of relative vs total abundance
%         costSize = zeros(2,length(waterMasses)); % store weighted costs for all size data groups
%         for i = 1:length(waterMasses)
%             wm = waterMasses{i};
%             if groupedByWaterOrigin
%                 cta = costComponents.([varLabel '_' wm '_autotroph_Tot']);
%                 cra = costComponents.([varLabel '_' wm '_autotroph_Rel']);
%                 cth = costComponents.([varLabel '_' wm '_heterotroph_Tot']);
%                 crh = costComponents.([varLabel '_' wm '_heterotroph_Rel']);
%             else
%                 cta = costComponents.([varLabel '_autotroph_Tot']);
%                 cra = costComponents.([varLabel '_autotroph_Rel']);
%                 cth = costComponents.([varLabel '_heterotroph_Tot']);
%                 crh = costComponents.([varLabel '_heterotroph_Rel']);
%             end
%             costSize(1,i) = mean(2 .* [1, weight_relVsTot] ./ (weight_relVsTot+1) .* [cta, cra]);
%             costSize(2,i) = mean(2 .* [1, weight_relVsTot] ./ (weight_relVsTot+1) .* [cth, crh]);
%         end
%         
%         costNutrient = zeros(1,length(Vars));
%         for i = 1:length(Vars)
%             costNutrient(i) = costComponents.(Vars{i});
%         end
%         
%         % Average the cost across data-types (nutrient & size)
%         cost = [mean(costNutrient), mean(costSize(:))];
%         % Assign size vs nutrient weighting
%         weight_sizeVsNutrient = 1; % weighting factor of size vs nutrient data
%         weights = 2 .* [1, weight_sizeVsNutrient] ./ (weight_sizeVsNutrient+1);
%         cost = weights .* cost;
% 
%         cost = mean(cost); % finally, average cost over nutrient and size data components
%         
%     case 'Hellinger_MVN_groupWaterOrigin'
%         % Same cost function as 'Hellinger_groupWaterOrigin' except that
%         % information on absolute abundances-at-size is included in the
%         % cost function by assuming that abundance-at-size is normally
%         % distributed and using generic variability parameters (from the
%         % Gaussian approximation to the multinomial distribution) to derive
%         % the probability density functions used as arguments in the
%         % Hellinger distance metric.
%         % In short: this fits to abundance-at-size rather than merely
%         % fitting to relative abundance-at-size.
%         
%         % Scalar data
%         for i = 1:length(Vars)
%             varLabel = Vars{i};
%             yobs = Data.scalar.scaled_Value(strcmp(Data.scalar.Variable, varLabel));
%             ymod = modData.scalar.scaled_Value(strcmp(modData.scalar.Variable, varLabel),:);
%             n = size(ymod, 1);
%             m = size(ymod, 2);
%             % As the standardised data are approximately normally distributed,
%             % we could assume that identically transformed model outputs
%             % are also normally distributed, then find their means and
%             % variances to calculate Hellinger distances analytically.
%             % However, assuming normality of transformed model output is
%             % inappropriate because it could be any shape...
%             % Thus, derive cdfs and pdfs directly without assuming that
%             % transformed model outputs follow normal distributions.
%             % To compare observed and modelled distributions, the pdfs will
%             % need to be evaluated across identical domains.
%             cdf = (1:n)' ./ n;
%             yobsc = sort(yobs);
%             ymodc = sort(ymod);
%             % remove any duplicate values (more likely in model output than
%             % in data), retaining the largest probability values
%             keep = ~[diff(yobsc) == 0; false];
%             cdfobs = cdf(keep);
%             yobsc = yobsc(keep);
%             % store cdfmod as cell-array because vector lengths may differ
%             % after removing duplicates
%             keep = ~[diff(ymodc) == 0; false(1, m)];
%             cdfmod = cell(1, m);
%             ymodc_ = cell(1, m);
%             for ij = 1:m
%                 cdfmod{:,ij} = cdf(keep(:,ij));
%                 ymodc_{:,ij} = ymodc(keep(:,ij),ij);
%             end
%             % define regular grid across measurement space -- define number
%             % of grid nodes using the number of observations
%             vrange = [min([yobsc(:); ymodc(:)]), max([yobsc(:); ymodc(:)])];
%             grid = nan(n, 1);
%             grid(2:end) = linspace(vrange(1), vrange(2), n-1)';
%             sp = diff(grid(2:3));
%             grid(1) = grid(2) - sp;
%             % interpolate cdfs and derive pdfs
%             cdfobsi = interp1(yobsc, cdfobs, grid, 'linear', 'extrap');
%             cdfobsi(cdfobsi < min(cdfobs)) = 0;
%             cdfobsi(cdfobsi > 1) = 1;
%             pdfobsi = diff(cdfobsi);
%             cdfmodi = nan(n, m);
%             pdfmodi = nan(n-1, m);
%             for ij = 1:m
%                 cdfmodi(:,ij) = interp1(ymodc_{ij}, cdfmod{ij}, grid, 'linear', 'extrap');
%                 cdfmodi(cdfmodi(:,ij) < min(cdfmod{ij}), ij) = 0;
%                 cdfmodi(cdfmodi(:,ij) > 1, ij) = 1;
%                 pdfmodi(:,ij) = diff(cdfmodi(:,ij));
%             end
%             
%             hellingerDistance = (1 - sum((pdfobsi .* pdfmodi) .^ 0.5)) .^ 0.5;
%             
%             costComponents.(varLabel) = mean(hellingerDistance); % average over trajectory selections
%             
%         end
%         
%         % Vector (size) data
%         for i = 1:length(VarsSize)
%             varLabel = VarsSize{i};
%             ind0 = strcmp(Data.sizeFull.dataBinned.groupedByOrigin.Variable, varLabel);
%             waterMasses = unique(Data.sizeFull.dataBinned.groupedByOrigin.waterMass);
%             
%             for w = 1:length(waterMasses)
%                 wm = waterMasses{w};
%                 ind1 = ind0 & strcmp(Data.sizeFull.dataBinned.groupedByOrigin.waterMass, wm);
%                 
%                 % autotrophs
%                 ind = ind1 & strcmp(Data.sizeFull.dataBinned.groupedByOrigin.trophicLevel, 'autotroph');
%                 yobs = Data.sizeFull.dataBinned.groupedByOrigin.Value(ind);
%                 ymod = modData.sizeFull.(['Value_' wm])(ind,:);
%                 % Derive generic variabilities using normal approximation
%                 % to multinomial distribution...
%                 % This approximation is only valid for large multinomial
%                 % samples. This should be OK for cell-conc data, but maybe
%                 % not for the bio-volume data because the chosen units are
%                 % m^3 / m^3... we can just rescale the units here, although
%                 % it will be cleaner to alter the units before arranging
%                 % the data..                
%                 switch varLabel, case 'BioVol'
%                     % Alter bio-volume units from m^3/m^3 to (mm)^3/m^3
%                     yobs = 1e9 .* yobs;
%                     ymod = 1e9 .* ymod;
%                 end
%                 Nobs = sum(yobs); % totals
%                 Nmod = sum(ymod);                
%                 pobs = yobs ./ Nobs; % size-class probabilities
%                 pmod = ymod ./ Nmod;                
%                 vobs = Nobs .* pobs .* (1 - pobs); % multinomial variabilities
%                 vmod = Nmod .* pmod .* (1 - pmod);
%                 
%                 % Use analytical formula for Hellinger distance between two
%                 % normal distributions. Note that since we are omitting
%                 % covariances, determinants of covariance matrices equal
%                 % the products of the variance vectors, and matrix
%                 % inversions are not necessary.
%                 vsum = vobs + vmod;
%                 hellingerDistance = (1 - (2 .* vobs .^ 0.5 .* vmod .^ 0.5 ./ vsum) .^ 0.5 .* ... 
%                     exp(-0.25 .* (yobs - ymod) .^ 2 ./ vsum)) .^ 0.5;
%                 costComponents.([varLabel '_' wm '_autotroph']) = mean(hellingerDistance(:)); % average over sizes and trajectory selections                
%                 % The multivariate expression (commented out below) may
%                 % seem more appropriate than separately calculating 
%                 % Hellinger distances for each size class then averaging.
%                 % However, treating size classes as independently normal is
%                 % fine if covariances are omitted.
%                 % Multivariate-normal version -- use log-scale for stability
% %                 hellingerDistance = (1 - exp(0.25 .* (sum(log(vobs)) + ... 
% %                     sum(log(vmod)) - 2 .* sum(log(0.5 .* vsum)) - ...
% %                     0.5 .* sum(2 .* (yobs - ymod) .^ 2 ./ vsum)))) .^ 0.5;
% %                 costComponents.([varLabel '_' wm '_autotroph']) = mean(hellingerDistance); % average over trajectory selections
%                 
%                 % heterotrophs
%                 ind = ind1 & strcmp(Data.sizeFull.dataBinned.groupedByOrigin.trophicLevel, 'heterotroph');
%                 yobs = Data.sizeFull.dataBinned.groupedByOrigin.Value(ind);
%                 ymod = modData.sizeFull.(['Value_' wm])(ind,:);
% 
%                 switch varLabel, case 'BioVol'
%                     % Alter bio-volume units from m^3/m^3 to (mm)^3/m^3
%                     yobs = 1e9 .* yobs;
%                     ymod = 1e9 .* ymod;
%                 end
%                 Nobs = sum(yobs); % totals
%                 Nmod = sum(ymod);                
%                 pobs = yobs ./ Nobs; % size-class probabilities
%                 pmod = ymod ./ Nmod;                
%                 vobs = Nobs .* pobs .* (1 - pobs); % multinomial variabilities
%                 vmod = Nmod .* pmod .* (1 - pmod);
%                 vsum = vobs + vmod;
%                 hellingerDistance = (1 - (2 .* vobs .^ 0.5 .* vmod .^ 0.5 ./ vsum) .^ 0.5 .* ... 
%                     exp(-0.25 .* (yobs - ymod) .^ 2 ./ vsum)) .^ 0.5;
%                 costComponents.([varLabel '_' wm '_heterotroph']) = mean(hellingerDistance(:)); % average over sizes and trajectory selections
% %                 % Multivariate version
% %                 hellingerDistance = (1 - exp(0.25 .* (sum(log(vobs)) + ... 
% %                     sum(log(vmod)) - 2 .* sum(log(0.5 .* vsum)) - ...
% %                     0.5 .* sum(2 .* (yobs - ymod) .^ 2 ./ vsum)))) .^ 0.5;
% %                 costComponents.([varLabel '_' wm '_heterotroph']) = mean(hellingerDistance); % average over trajectory selections
%             end
%         end
%         
%         % Take averages across data-types to assign equal weightings to the
%         % scalar data and size data.
%         cost = zeros(1,2);
%         fields = fieldnames(costComponents);
%         for i = 1:length(fields)
%             if ismember(fields{i}, Vars)
%                 % scalar data
%                 cost(1) = cost(1) + costComponents.(fields{i});
%             else
%                 % size data
%                 cost(2) = cost(2) + costComponents.(fields{i});
%             end
%         end
%         cost(1) = cost(1) ./ length(Vars);
%         cost(2) = cost(2) ./ (length(fields) - length(Vars));
%         
% %        disp(costComponents)
% %        disp(cost)
% 
%         cost = mean(cost);
% 
% %        disp(cost)
%         
% end
% 

%%

if ~exist('cost', 'var')
    cost = nan;
    warning('Could not evaluate cost... Check that the name-value pair (selectFunction,costFunctionType) is properly specified and corresponds to a viable option within costFunction.m')
end
if ~exist('costComponents', 'var')
    costComponents = nan;
end
