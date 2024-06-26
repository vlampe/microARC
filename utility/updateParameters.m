function Params = updateParameters(Params,FixedParams,varargin)

if ~isempty(varargin)
    if nargin == 3 && isnumeric(varargin{1}) % if new param values are passed as a single vector (with labels stored in FixedParams)
        newvals = varargin{1};
        parnames = FixedParams.tunePars;
        npars = length(parnames);
        Vol = FixedParams.sizeAll;
        for i = 1:npars
            pn = parnames{i};
            Params.(pn) = newvals(i); % update scalars
            % update size-dependent parameters - constructed from scalars
            % suffixed with _a and _b
            if ~isempty(regexp(pn,'_a', 'once')) || ~isempty(regexp(pn,'_b', 'once'))
                Name = pn(1:end-2);
                param_a = Params.([Name '_a']);
                param_b = Params.([Name '_b']);
                Params.(Name) = Params.([Name '_func'])(param_a, param_b, Vol);
            end
            if strcmp(pn, 'beta1') || strcmp(pn, 'beta2') || strcmp(pn, 'beta3')
                Name = pn(1:end-1);
                param1 = Params.([Name '1']);
                param2 = Params.([Name '2']);
                param3 = Params.([Name '3']);
                Params.(Name) = Params.([Name '_func'])(param1, param2, param3, Vol);
            end
            if strcmp(pn, 'wDOM1') || strcmp(pn, 'wPOM1')
                Name = pn(1:end-1);
                param1 = Params.([Name '1']);
                Params.(Name) = Params.([Name '_func'])(param1);
            end
        end
        
        % Functions of parameters
        
        if any(strcmp('Qmin_QC_a', parnames)) || any(strcmp('Qmin_QC_b', parnames)) || ...
                any(strcmp('Qmax_delQ_a', parnames)) || any(strcmp('Qmax_delQ_b', parnames))
            Params.Qmax_QC = Params.Qmax_QC_func(Params.Qmin_QC, Params.Qmax_delQ);
            Params.delQ_QC = Params.delQ_QC_func(Params.Qmin_QC, Params.Qmax_QC);
        end
        
        if any(strcmp('Vmax_QC_a', parnames)) || any(strcmp('Vmax_QC_b', parnames)) || ...
                any(strcmp('aN_QC_a', parnames)) || any(strcmp('aN_QC_b', parnames))
            Params.kN = Params.kN_func(Params.Vmax_QC, Params.aN_QC);
        end
        
        if any(strcmp('rDOC', parnames)) || any(strcmp('rDON', parnames)) || ...
                any(strcmp('rPOC', parnames)) || any(strcmp('rPON', parnames))
            if isfield(Params, 'rDOC_func')
                Params.rDOC = Params.rDOC_func(Params.rDON);
            end
            if isfield(Params, 'rPOC_func')
                Params.rPOC = Params.rPOC_func(Params.rPON);
            end
            Params.rOM = Params.rOM_func(Params.rDOC, Params.rDON, ... 
                Params.rPOC, Params.rPON);
        end
        
        if any(strcmp('wDOM1', parnames)) || any(strcmp('wPOM1', parnames))
            Params.wk = Params.wk_func(Params.wDOM, Params.wPOM);
        end
        
        if ismember('Gmax_a', parnames) || ismember('aG', parnames)
            Params.k_G = Params.k_G_func(Params.Gmax_a, Params.aG);
        end

        if any(strcmp('sigG', parnames))
            Params.phi = Params.phi_func(Params.delta_opt, Params.sigG, FixedParams.delta);
        end
        
    else % if new params are passed as name-values pairs
        
        extractVarargin(varargin)
        
        % Scalar parameters
        parNames = Params.scalarParams;
        for i = 1:length(parNames)
            if exist(parNames{i}, 'var')
                name = parNames{i};
                value = eval(parNames{i});
                Params.(name) = value;
            end
        end
        
        % Vector parameters
        Vol = FixedParams.sizeAll(:);
        parNames = Params.vectorParams;
        suf = {'_a', '_b'};
        
        for i = 1:length(parNames)
            
            if exist(parNames{i}, 'var')
            
                name = parNames{i};
                value = eval(parNames{i});
                Params.(name) = value; % update the scalar

                % Most size-dependent terms are parameterised by power functions of
                % size using two parameters subscripted with '_a' or '_b', one
                % of which has just been updated... now recalculate the vector
                
                if any(strcmp(name(end-1:end), suf))
                    Name = name(1:length(name)-2);
                    param_a = Params.([Name '_a']);
                    param_b = Params.([Name '_b']);
                    Params.(Name) = Params.([Name '_func'])(param_a, param_b, Vol);
                end                
                if any(strcmp(name, {'beta1','beta2','beta3'}))
                    Name = name(1:length(name)-1);
                    b1 = Params.([Name '1']);
                    b2 = Params.([Name '2']);
                    b3 = Params.([Name '3']);
                    Params.(Name) = Params.([Name '_func'])(b1, b2, b3, Vol);
                end                
                if any(strcmp(name, {'wDOM1','wPOM1'}))
                    Name = name(1:length(name)-1);
                    param1 = Params.([Name '1']);
                    Params.(Name) = Params.([Name '_func'])(param1);
                end                
            end
        end
        
        
        % Functions of parameters

        if exist('Qmin_QC_a', 'var') || exist('Qmin_QC_b', 'var') || ...
                exist('Qmax_delQ_a', 'var') || exist('Qmax_delQ_b', 'var')
            Params.Qmax_QC = Params.Qmax_QC_func(Params.Qmin_QC, Params.Qmax_delQ);
            Params.delQ_QC = Params.delQ_QC_func(Params.Qmin_QC, Params.Qmax_QC);
        end
        
        if exist('Vmax_QC_a', 'var') || exist('Vmax_QC_b', 'var') || ... 
                exist('aN_QC_a', 'var') || exist('aN_QC_b', 'var')
            Params.kN = Params.kN_func(Params.Vmax_QC, Params.aN_QC);
        end
        
        if exist('rDOC', 'var') || exist('rDON', 'var') || ... 
                exist('rPOC', 'var') || exist('rPON', 'var')
            if isfield(Params, 'rDOC_func')
                Params.rDOC = Params.rDOC_func(Params.rDON);
            end
            if isfield(Params, 'rPOC_func')
                Params.rPOC = Params.rPOC_func(Params.rPON);
            end
            Params.rOM = Params.rOM_func(Params.rDOC, Params.rDON, ... 
                Params.rPOC, Params.rPON);
        end
        
        if exist('wDOM1', 'var') || exist('wPOM1', 'var')
            Params.wk = Params.wk_func(Params.wDOM, Params.wPOM);
        end
        
        if exist('Gmax_a', 'var') || exist('aG', 'var')
            Params.k_G = Params.k_G_func(Params.Gmax_a, Params.aG);
        end

        if exist('sigG', 'var')
            Params.phi = Params.phi_func(Params.delta_opt, Params.sigG, FixedParams.delta);
        end
        
    end

end

end
