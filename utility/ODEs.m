function [dvdt, out] = ODEs(t, v_in, parameterList, forc, timeStep, returnExtra)

fixedParams = parameterList.FixedParams;
params = parameterList.Params;

%% MODEL DIMENSIONS

nz = fixedParams.nz;    % number of depth layers
nPP_size = fixedParams.nPP_size;  % number of phytoplankton size classes
nPP_nut = fixedParams.nPP_nut;  % number of phytoplankton nutrient classes
nZP_size = fixedParams.nZP_size;  % number of zooplankton size classes
nZP_nut = fixedParams.nZP_nut;  % number of zooplankton nutrient classes
nOM_type = fixedParams.nOM_type;  % number of organic matter types
nOM_nut = fixedParams.nOM_nut;  % number of organic nutrient classes
phyto = fixedParams.phytoplankton;
zoo = fixedParams.zooplankton;
nsize = nPP_size + nZP_size;

%% INITIAL CONDITIONS

% Inorganic nitrogen
N = v_in(fixedParams.IN_index)';

% Plankton
B = [reshape(v_in(fixedParams.PP_index), [nPP_size nz nPP_nut]); ... % autotrophs
    cat(3, ...
    reshape(v_in(fixedParams.ZP_index), [nZP_size nz nZP_nut]), ... % heterotrophs (include extra zeros for chl-a)
    zeros(nZP_size, nz, nPP_nut - nZP_nut))];

B_C = B(:,:,fixedParams.PP_C_index); % all planktonic carbon (mmol C m^{-3})

% Organic matter
OM =reshape(v_in(fixedParams.OM_index), [nOM_type nz nOM_nut]);


%% FORCING DATA

T = forc.T(:,timeStep-1:timeStep);
K = forc.K(:,timeStep-1:timeStep);
Isurf = forc.PARsurf(:,timeStep-1:timeStep);
lat = forc.lat(:,timeStep-1:timeStep); % latitude
yd  = forc.yd(:,timeStep-1:timeStep); % day of year (year day)


% Linearly interpolate between days
T = (T(:,1) + diff(T,1,2) .* t)';
K = K(:,1) + diff(K,1,2) .* t;
Isurf = (Isurf(:,1) + diff(Isurf,1,2) .* t)';
lat = lat(1) + diff(lat) .* t;
yd = yd(1) + diff(yd) .* t;

% % Calculate light levels at depth -- within each depth layer light
% % attenuates over half the layer width plus the combined widths of all
% % shallower layers.
% att0    = fixedParams.attSW + fixedParams.attP .* sum(B(:,:,fixedParams.PP_Chl_index)); 
% att     = att0 .* fixedParams.zwidth';
% att     = 0.5 * att + [0 cumsum(att(1:nz-1))];
% out.I   = Isurf * exp(-att); % irradiance at depth layer midpoints

% Light attenuates exponentially with depth.
att0    = fixedParams.attSW + fixedParams.attP .* sum(B(:,:,fixedParams.PP_Chl_index)); 
att1 = att0 .* abs(fixedParams.zw(1:nz))';
att2 = att0 .* abs(fixedParams.zw(2:nz+1))';
out.I = Isurf ./ att0 ./ fixedParams.zwidth' .* (exp(-att1) - exp(-att2)); % irradiance averaged within depth layers


%% MODEL EQUATIONS

%~~~~~~~~~~~
% Physiology
%~~~~~~~~~~~

% nutrient quotas
out.Q = B ./ B_C; % mmol / mmol C

% Nutrient limitation
out.gammaN = max(0, min(1, (out.Q(:,:,fixedParams.PP_N_index) - params.Qmin_QC) ./ params.delQ_QC));

% Uptake regulation
out.Qstat = 1 - out.gammaN .^ params.h; %%% ATTENTION: in paper.pdf, it says 1 - gammaN^(1/h). But then, curve looks unsuitable 

% For testing: Nutrient uptake regulation after Schartau et al 2007, Eq. B7: 
% out.Qstat = 1- exp(-1000 .*(abs(out.Q(:,:,fixedParams.PP_N_index) - 0.2) - (out.Q(:,:,fixedParams.PP_N_index) - 0.2)).^2);

% Temperature dependence
out.gammaT = exp(params.A .* (T - params.Tref));

% Background mortality
% out.mortality = params.m .* B;
out.mortality = (params.m + params.m2 .* B) .* B; % mmol /m^3 / d
if isfield(params, 'K_m') && any(params.K_m > 0)
    % hyperbolic term reduces mortality at very low abundance
    out.mortality = (B ./ (params.K_m + B)) .* out.mortality;
end
    


%~~~~~~~~~~~
% Autotrophy
%~~~~~~~~~~~

out.V = zeros(nsize, nz, nPP_nut); % all uptake rates

% Nutrient uptake
out.V(:,:,fixedParams.PP_N_index) = ... 
    MichaelisMenton(params.Vmax_QC, params.kN, N) .* out.gammaT .* out.Qstat; % d^-1

% Photosynthesis
zeroLight = out.I(1) == 0;

if ~zeroLight
    out.psat = params.pmax .* out.gammaT .* out.gammaN; % light saturated photosynthetic rate (N-dependent)
    
    %aP_Q_I = (params.aP .* out.I) .* out.Q(:,:,fixedParams.PP_Chl_index); % photon capture rate (Chl-dependent)
    %out.pc = out.psat .* (1 - exp(-aP_Q_I ./ out.psat)); % photosynthetic (carbon production) rate (1 / day)  
        
    a_Chl = params.aP .* out.Q(:,:,fixedParams.PP_Chl_index); % Chl-dependent initial slope of photon capture rate (day * m^2 / μEinstein)
    
    [out.pc, out.I_lim, out.dl] = light_lim(a_Chl, out.psat, Isurf, att0, lat, yd, fixedParams.zwidth, 'dist', 'flat'); % pc = photosynthetic (carbon production) rate (1 / day), I_lim = light limitation (depth- and, optionally, time averaged), and dl = daylength
  
    aP_Q_I  = a_Chl .* out.I; 
    out.rho = params.theta .* min(1, out.pc ./ aP_Q_I, 'includenan'); % proportion of new nitrogen prodcution diverted to chlorophyll (mg Chl / mmol N)
    
    % The min function in rho should only be required to correct numerical
    % inaccuracies when I -> 0.
    out.V(:,:,fixedParams.PP_C_index) = max(0, ... 
        out.pc - params.xi .* out.V(:,:,fixedParams.PP_N_index)); % photosynthesis minus metabolic cost: carbon production rate (d^-1)
    out.V(:,:,fixedParams.PP_Chl_index) = out.rho .* out.V(:,:,fixedParams.PP_N_index); % chlorophyll production rate (d^-1) 
end

out.V(isnan(out.V)) = 0;

out.uptake = B_C .* out.V; % NPP mmol / m^3 / d

out.N_uptake_losses = sum(out.uptake(:,:,fixedParams.PP_N_index)); % mmol N / m^3 / d


%~~~~~~~~~~~~~
% Heterotrophy
%~~~~~~~~~~~~~

phi_BC = params.phi .* reshape(B_C, [1 size(B_C)]);
F = max(sum(phi_BC, 2),2e-30);  % sonst produziert exp funktion bei out.G -Inf werte die sich als NaNs durchziehen!
phi_BC2 = phi_BC .^ 2;
Phi = phi_BC2 ./ sum(phi_BC2, 2); % prey preference

out.G = (reshape(out.gammaT, [1 size(out.gammaT)]) .* ...
    MichaelisMenton(params.Gmax, params.k_G, F) .* (1-exp(params.Lambda .* F))) .* Phi;  % grazing rate (1 / day)

out.predation_losses_all = reshape(out.Q, [1, nsize, nz, nPP_nut]) .* reshape(B_C(zoo,:), [nZP_size, 1, nz]) .* out.G;

out.lambda = zeros(nZP_size, 1, nz, nPP_nut);
out.lambda(:,:,:,fixedParams.ZP_C_index) = out.gammaN(zoo,:);
out.lambda(:,:,:,fixedParams.ZP_N_index) = out.Qstat(zoo,:);
out.lambda = params.lambda_max .* out.lambda;

out.predation_gains_all = out.lambda .* out.predation_losses_all;

mess = out.predation_losses_all(:,:,:,~fixedParams.PP_Chl_index) - ... 
    out.predation_gains_all(:,:,:,~fixedParams.PP_Chl_index);

if nZP_size > 1
    out.predation_losses = sum(out.predation_losses_all);  % sum over predators
end
out.predation_losses = reshape(out.predation_losses, [nsize, nz, nPP_nut]);
out.predation_gains = [zeros(nPP_size, nz, nPP_nut);  
    reshape(sum(out.predation_gains_all, 2), [nZP_size, nz, nPP_nut])];  % sum over prey

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Sources and sinks of organic matter
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Messy feeding
out.OM_mess = zeros(nOM_type, nz, nOM_nut);
beta_mess = reshape(params.beta, [1, nsize]) .* mess;
out.OM_mess(fixedParams.DOM_index,:,:) = sum(beta_mess, [1, 2]);
out.OM_mess(fixedParams.POM_index,:,:) = sum(mess - beta_mess, [1, 2]);

% Mortality
out.OM_mort = zeros(nOM_type, nz, nOM_nut);
beta_m_B = params.beta .* out.mortality(:,:,~fixedParams.PP_Chl_index);
out.OM_mort(fixedParams.DOM_index,:,:) = sum(beta_m_B);
out.OM_mort(fixedParams.POM_index,:,:) = sum(out.mortality(:,:,~fixedParams.PP_Chl_index) - beta_m_B);

% Remineralisation
out.OM_remin = params.rOM .* OM;

SOM = out.OM_mort + out.OM_mess - out.OM_remin; % (mmol / m^3 / day)
  
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Sources of inorganic nutrients
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Remineralisation
SN = sum(out.OM_remin(:,:,fixedParams.OM_N_index)); % (mmol N / m^3 / day)

%~~~~~~~~~~
% Diffusion
%~~~~~~~~~~

B_C_t = B_C';
OM_ = reshape(permute(OM, [2, 1, 3]), [nz, nOM_type * nOM_nut]);

v_diffuse = diffusion_1D([N(:), B_C_t, OM_], K, fixedParams.zwidth, fixedParams.delz);

N_diffuse = v_diffuse(:,1);
B_diffuse = out.Q .* v_diffuse(:,2:nsize+1)';

OM_diffuse = permute(reshape( ...
    v_diffuse(:,nsize+2:end), ...
    [nz, nOM_type, nOM_nut]), [2, 1, 3]);

out.B_diffuse = B_diffuse; 
out.OM_diffuse = OM_diffuse;

%~~~~~~~~
% Sinking
%~~~~~~~~

v_sink = sinking([B_C_t, OM_], [params.wp, params.wk], fixedParams.zwidth);

B_sink = out.Q .* v_sink(:,1:nsize)';
OM_sink = permute(reshape(v_sink(:,nsize+1:end), ... 
    [nz, nOM_type, nOM_nut]), [2, 1, 3]);

out.B_sink = B_sink;
out.OM_sink = OM_sink;

%~~~~~
% ODEs
%~~~~~

% Inorganic nutrients
dNdt = N_diffuse - out.N_uptake_losses(:) + SN(:);

% Plankton
dBdt = B_sink + B_diffuse + out.uptake - out.predation_losses + out.predation_gains - out.mortality;
dPPdt = dBdt(phyto,:,:);
dZPdt = dBdt(zoo,:,~fixedParams.PP_Chl_index);

% Organic matter
dOMdt = OM_diffuse + OM_sink + SOM; 

dvdt = [dNdt; dPPdt(:); dZPdt(:); dOMdt(:)];

%% AUXILIARY OUTPUTS

if (islogical(returnExtra) && returnExtra) || ... 
        (~islogical(returnExtra) && ~any(strcmp(returnExtra, 'none')))
    
    out.cellDensity = B_C ./ params.Q_C;  % (mmol C m^{-3}) ./ (mmol C cell^{-1}) = cell * m^{-3}
    out.biovolume = fixedParams.sizeAll .* out.cellDensity; % µm^3 cell^{-1} * cell m^{-3} = µm^3 m^{-3}
%     out.biovolume = 1e-18 * fixedParams.sizeAll .* out.cellDensity;

    % Extra output variables retained by default when return = true or 'all'.
    keepVars = {'I', 'Q', 'V', 'G', 'lambda', 'cellDensity', 'biovolume', 'I_lim'};
 %  keepVars = {'B_sink', 'B_diffuse', 'uptake', 'predation_losses', 'predation_gains', 'mortality'  }; 
    keepVars = [keepVars, {'OM_sink', 'B_sink', 'B_diffuse', 'OM_diffuse', 'uptake', 'predation_losses', 'predation_gains', 'mortality', 'OM_remin', 'OM_mort', 'OM_mess', 'dl', 'pc', 'B_diff_exp'} ]; 
 %keepVars = [keepVars, {'OM_sink', 'B_sink', 'B_diffuse', 'uptake', 'predation_losses', 'predation_gains', 'mortality', 'OM_remin', 'OM_mort', 'OM_mess'} ]; 
    % Any term can be included in keepVars, but it's useful to be sparing
    % with memory by only returning a few terms then deriving more extra
    % output outside this ODEs.m function.
    
    if ~islogical(returnExtra) && ~any(strcmp(returnExtra, 'all'))
        % if extra output variables have been specified explicitly...
        keepVars = returnExtra;
    end
    
    fields = fieldnames(out);
    out = rmfield(out, fields(~ismember(fields, keepVars)));
    
else
    out = struct();
end

    
end


%% Functions

function v = diffusion_1D(u,K,w,delz)
% Rate of change of u due to diffusion, assuming zero flux boundary conditions.
% Inputs: u = concentrations, size(u)=[nz nvar]
%         K = diffusivities, size(K)=[nz-1 1]
%         w = depth layer widths, size(w)=[nz 1]
%         delz = distance between depth layer centers, size(delz)=[nz-1 1]
z = zeros(1,size(u,2));
v = diff([z; (K ./ delz) .* diff(u); z]) ./ w;
end

function v = sinking(u,s,w)
% Rate of change of u due to sinking.
% Inputs: u = concentrations, size(u)=[nz nvar]
%         s = sinking speed, size(s)=[1 nvar]
%         w = depth layer widths, size(w)=[nz 1]
v = ((-s) ./ w) .* diff([zeros(1,size(u,2)); u]);
end

function v = MichaelisMenton(m,k,u)
% Uptake rate of u, given maximum m and half saturation k
u(u<0) = 0; % include for robustness... there shouldn't be any negatives
v = m .* u ./ (u + k);
end

