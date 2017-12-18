%function [model,blacklistEM,R,pvind,growth,glcind,o2ind] = loadmodels
function [model,blacklistEM,R,pvind,growth,glcind,o2ind] = loadmodels(system)
% function [model,blacklistEM,R,pvind,growth,glcind,o2ind] = loadmodels(system)
% Laurence Yang, Jan 7, 2010 --> July 29, 2010
% Load model
switch lower(system)
    case 'iaf1260'
        
        load ('D:\My Documents\MATLAB\models\iAF1260\Ecoli_iAF1260.mat');             
        load ('D:\My Documents\MATLAB\models\iAF1260\R.mat'); 
        load ('D:\My Documents\MATLAB\models\iAF1260\pvind.mat');
        
        %model = iMMr;
        model.vl = model.lb;
        model.vu = model.ub;
        
          
        [Sm,Sn]=size(model.S);
        b=sparse(Sm,1,0);
        %growth=find(model.c);
        growth = find(strncmp('Ec_biomass_iAF1260_core_59p81M',model.rxns,10));
        
        glcind=strncmp('EX_glc',model.rxns,6);
        glycind=strncmp('EX_glyc',model.rxns,7);
        o2ind = strncmp('EX_o2',model.rxns,5);
        model.glcind = find(glcind);
        model.o2ind=find(o2ind);
        model.growth=growth;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Reduce bounds to +/- 100
        toobigvu=abs(model.vu)>100; % for the unlikely case that vu<-100
        model.vu(toobigvu)=sign(model.vu(toobigvu))*100;
        toobigvl=abs(model.vl)>100;
        model.vl(toobigvl)=sign(model.vl(toobigvl))*100;
        
        model.vl(glcind)=-20;
        model.vl(o2ind)=-20;
        model.vl(glycind)=0;    % Don't allow uptake of glycerol
        removerxn = [];
        model.vl(removerxn)=0; model.vu(removerxn)=0;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        blacklist.KO=[];
        blacklist.Reg=[]; 
        exchangeind=[];
        for i=1:Sn
            % If reaction name ends with pp, don't regulate this
            namelength = length(model.rxns{i});
            ppind = strfind(model.rxns{i},'pp');
            if not(isempty(ppind))
                if ppind==namelength-1;
                    exchangeind=[exchangeind; i];
                end
            end
            % If reaction name ends with tex, don't regulate this
            texind = strfind(model.rxns{i},'tex');
            if not(isempty(texind))
                if texind==namelength-2;
                    exchangeind=[exchangeind; i];
                end
            end            
        end
        EXind= strmatch('EX_',model.rxns);
        exchangeind=union(exchangeind,EXind);
        blacklist.Reg=union(blacklist.Reg,exchangeind);
        blacklist.KO=union(blacklist.KO,exchangeind);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Remove non-gene-associated reactions from target list
        NGAR=find(sum(model.rxnGeneMat,2)==0);
        blacklist.Reg=union(blacklist.Reg,NGAR);
        blacklist.KO=union(blacklist.KO,NGAR);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Remove essential reactions from target KOs
        essrxns=[];
        
        blacklist.KO=union(blacklist.KO,essrxns);
        
        % Exclude user-defined reactions
        adipaterxns=[strmatch('2OXOADPRD',model.rxns);
            strmatch('2HADPDH',model.rxns);
            strmatch('H2EDORD',model.rxns);
            strmatch('ADt',model.rxns);
            strmatch('EX_ad(e)',model.rxns)];
        blacklist.KO=union(blacklist.KO,adipaterxns);
        blacklist.Reg=union(blacklist.Reg,adipaterxns);
        
        % Exclude from targets based on subsystem:
        subsysind=[];
        % Cell envelope biosynthesis
        subsysind=strmatch('Cell Envelope Biosynthesis',model.subSystems);
        % Glycerophospholipid metabolism
        subsysind=union(subsysind,...
            strmatch('Glycerophospholipid Metabolism',model.subSystems));
        % Inorganic ion transport and metabolism
        subsysind=union(subsysind,...
            strmatch('Inorganic Ion Transport and Metabolism',model.subSystems));
        % Lipopolysaccharide biosynthesis and recycling
        subsysind=union(subsysind,...
            strmatch('Lipopolysaccharide Biosynthesis Recycling',model.subSystems));
        % Membrane lipid metabolism
        subsysind=union(subsysind,...
            strmatch('Membrane Lipid Metabolism',model.subSystems));
        % Inner membrane transport
        subsysind=union(subsysind,...
            strmatch('Transport Inner Membrane',model.subSystems));
        % Sphingolipid Metabolism
        subsysind=union(subsysind,...
            strmatch('Sphingolipid Metabolism',model.subSystems));
        % Outer membrane transport
        subsysind=union(subsysind,...
            strmatch('Transport Outer Membrane',model.subSystems));
        % Outer membrane porin transport
        subsysind=union(subsysind,...
            strmatch('Transport Outer Membrane Porin',model.subSystems));
        % Murein biosynthesis
        subsysind=union(subsysind,...
            strmatch('Murein Biosynthesis',model.subSystems));
        % Murein recycling
        subsysind=union(subsysind,...
            strmatch('Murein Recycling',model.subSystems));
        % Fatty Acid Biosynthesis
        subsysind=union(subsysind,...
            strmatch('Fatty Acid Biosynthesis',model.subSystems));
        % tRNA charging
        blacklist.KO=union(blacklist.KO,subsysind);
        blacklist.Reg=union(blacklist.Reg,subsysind);
        % Remove E_o2 from blacklist to allow choice of
        % aerobic/anaerobic/microaerobic
        
        % Provide initial modifications for OptReg:LS
        IC.KO=[]; IC.Down=[]; IC.UF=[]; IC.UR=[];
        IC.KO=strmatch('ACOAHim',model.rxns,'exact');
        %IC.KO=strmatch('GPDDA2',model.rxns,'exact');    % Makes ethanol from g3pe
        %IC.KO=strmatch('ALCD',model.rxns);  % Alcohol dehydrogenases
        %(ethanol & glycerol) 
 
end
blacklistEM = blacklist.Reg;    % EMILiO takes array, not structure







