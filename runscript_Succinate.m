%==========================================================================
% Generate strain designs using EMILiO                               
%============================================================ Laurence Yang
%==========================================================================
% Sept 10, 2010

clear;
%addpath /home/shyam/Matlab/Emilio % EMILiO package location

organism = 'ecoli';
growthmin = 0.1;
aerobic = 1;

%==========================================================================
% EMILiO parameters
noConc = 0;         % Prevent forcing of reversible fluxes in one direction
nCuts = 3;          % Maximum number of alternate optima enumerated
KpY=100;   % Objective weight for production rate relative to KKT violation

EMmaxIter = 5;      % EMmaxIter > 1: after finding one strain design, 
                    % identify the most important mod, blacklist it, and 
                    % run EMILiO again, in order to find totally different
                    % strain designs
LPpruneFrac = 0.9;
MILPpruneFrac = 0.99;   % Avoiding putting 1 to avoid numerical issues with MILP solver
epsProd = 1e-3;
epsProdILP = 1e-2;     % Make epsProdILP>epsProd to help growth-coupling
ICEM=[];
maxFailedIter = 5;
KKTtol = 1e-5;
KKTviol = 1e-7;

%==========================================================================
% Load model
%-----------
[model,blacklistEM,R,pvind,growth,glcind,o2ind]=loadmodels(organism);
[Sm,Sn]=size(model.S);
%-----------------------------
% Define the target metabolite 
% *********
% *CAUTION*
% *********
% Some models denote compartmentalization as
% "EX_ac(e)" while others denote it as "EX_ac_e"
%-----------------------------

prodind = strmatch('EX_succ(e)',model.rxns);

%****************************
% Sync LY format with COBRA format:
% model.lb = model.vl;
% model.ub = model.vu;

%****************************
% Remove reactions:
%model = changeRxnBounds(model, {'OHACT1','OHACT2','OHACT3','OHACT4','PDHm'},0,'b');

%****************************
% Allow acetate supplemented medium:
% model = changeRxnBounds(model, {'EX_ac(e)'}, -1., 'l');

%****************************
% Set changes to model in LY format:
% model.vl = model.lb;
% model.vu = model.ub;

%****************************
% Exclude trivial solutions!!
% ALCD2irm = strmatch('ALCD2irm',model.rxns);
% ALCD2ir = strmatch('ALCD2ir',model.rxns,'exact');
% ALCD2x = strmatch('ALCD2x',model.rxns);
% CSNATr = strmatch('CSNATr',model.rxns);
% CSNATm = strmatch('CSNATm',model.rxns);
% PRPPS = strmatch('PRPPS',model.rxns);
% ASADi = strmatch('ASADi',model.rxns);
% trivialmods = [PRPPS];
trivialmods = [];
blacklistEM = union(blacklistEM, trivialmods);
%****************************


%% ======================================================================== 
% Run EMILiO
%==========================================================================
model.vu(prodind) = 100;

%==========================================================================
% Report conditions before starting EMILiO
if aerobic
    aerostr='aerobic';
else
    aerostr='anaerobic';
end
fprintf('Designing strains of %s to make %s under %s conditions\n',...
    organism,model.rxns{prodind},aerostr);

runEMILiO;

if not(isempty(allsets2))
    % Note that after removing insignificant mods, some strains might become
    % equivalent. Therefore, reassess the set of unique strains
    [strains,strainsetind]=uniquestrains(model,allsets2);
end

% Summarize the strain design strategies
strainSummary = struct;
for i=1:length(strains)
    strainName = sprintf('strain_%d',i);
    strainSummary.(strainName) = printEMILIOmods(model, strains{1,i}, 0);
end

% Plot the strain designs
plotDesigns(model,strains,growth,prodind,1,growthmin);

% Draw Phenotypic phase planes of designs
figure()
cc=hsv(12);
hold on
models = struct;
models.wt = changeRxnBounds(model,{'EX_o2(e)'},-10,'l');
productionEnvelope(models.wt,{},'black', 'TAR_malcoa','biomass_wild');

for s=1:length(strains)
    modelName = sprintf('m%d',s);
    models.(modelName) = models.wt;
    models.(modelName).vl = strains{1,s}.vl;
    models.(modelName).vu = strains{1,s}.vu;
    models.(modelName).vl(strains{1,s}.KO)=0;
    models.(modelName).vu(strains{1,s}.KO)=0;
    models.(modelName).lb = models.(modelName).vl;
    models.(modelName).ub = models.(modelName).vu;
    productionEnvelope(models.(modelName),{},cc(s,:), 'TAR_malcoa','biomass_wild');
end

hold off
lineLegend = strread(num2str(1:length(strains)),'%s');
lineLegend = ['wt';lineLegend];
legend(lineLegend)
title('Malonyl-CoA EMILIO production envelope')
xlabel('Growth rate [1/hr]')
ylabel('Malonyl-coA accumulation [mmol/gDW/hr]')
