function [model,blacklistEM,R,pvind,growth,glcind,o2ind] = loadmodels(system)
%% function [model,blacklistEM,R,pvind,growth,glcind,o2ind] = loadmodels(system)
%Inputs 
%system - Name of model/organism already present in loadmodels
%Compatiable Models 
%iMM904 - S. ereviciae
%iAF1260 - E. coli - Not added yet
%Outputs
%model - Loaded model of organism
%blacklistEM - List of all blacklisted reaction (Reactions that cannot be
%egineered)
%R - Reduced row echelon form of the S matrix. It is recommended to
%calculate R in advance and provide it as an input to EMiLio to reduce the
%time taken to run EMiLio. Calculation of the reduced row echelon form for
%S matrices of genome-scale models takes a long time
%pvind - Index of fluxes that are pivot variables used to form R. 
%growth - Biomass reaction index in model.rxns
%glcind - GLucose uptake reaction index in model.rxns
%O2ind - O2 uptake reaction index in model.rxns
%Laurence Yang, Jan 7, 2010 --> July 29, 2010
%Shyam S, May 2013 - Added Documentation
% Load model
switch lower(system)
    case 'central'
        switch lower(reference)
            case 'artificial'                
                load centralIshiiReduced
                model=modelred{1};  %1: glucose aerobic
                pntAB=strmatch('PNT2A',model.rxns,'exact');
                model.vl(pntAB)=0; % PntAB irreversible
                vmin = model.vl;
                vmax = model.vu;
                %[vmin,vmax]=fvacplexint(model.S,model.vl,model.vu);
                % Temporary reference flux
                % Or, a vref that allows for best objective
                vmid = (vmax+vmin)/2;
                vrefmin = vmin+tempFac*(vmid-vmin);
                vrefmax = vmid+(1-tempFac)*(vmax-vmid);
            case 'experimental'
                load centralIshiiModel
                model=modelR;
                vmin = model.vl; vmax = model.vu;
                load ishiifluxesorderedRev
                vexpind=fluxindexschillingishii;
                vexpL=ishiifluxes(:,end); vexpU=ishiifluxes(:,end);
                [vrefmin,vrefmax]=fvaVref(model,vexpind,vexpL,vexpU); %Ishii data
        end
        [Sm,Sn]=size(model.S);
        b=sparse(Sm,1,0);
        [R,pvind]=rref(full([model.S, b])); % RREF is fast for central mdl
       
        % Free up dG bounds for thermokinetic constraints
        model.dGl = -400*ones(Sn,1); 
        model.dGu = 400*ones(Sn,1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        
        growth=find(model.c)
        %prodind = strmatch('SUCCxtO',model.rxns)
        blacklist.KO = [];
        blacklist.Reg = strmatch('SUCCUPR',model.rxns);
        %
        % Also prevent exchange fluxes from being up/down regulated
        exchangeind = [];
        for i=1:Sn
            if not(isempty(strfind(model.rxns{i},'xtO')))
                exchangeind=[exchangeind; i];
            end
        end
        blacklist.Reg = union(blacklist.Reg, exchangeind);
        blacklist.KO = union(blacklist.KO, exchangeind);
        findessential;
        blacklist.KO=union(blacklist.KO,essrxns);
        %}
        % Indicate proton traslocating reactions for EMILiOT/TMFA
        Hin = strmatch('HIN',model.mets,'exact');
        Hext = strmatch('HEXT',model.mets,'exact');
        protrans = find(model.S(Hin,:).*model.S(Hext,:));
        model.protrans = protrans;
        % But recall ion transport across membrane/delpH already included
        % in dGo. Need to fix pH. This pH should be consistent with assumed
        % pH when calculating dGo, or correct the dGo for different pH
        % according to Eqs. 6-10, (Henry et al., 2007).
        IC.KO=[]; IC.Down=[]; IC.UF=[]; IC.UR=[];
        IC.KO=strmatch('ADHE2Ra',model.rxns);
    
    case 'imm904'
        % S. cerevisiae without adipate pathway
        load iMM904                
        load RREFiMM904 R pvind
        
        % Model with adipate pathway added
        %load iMM904ad
        %load RREFiMM904ad R pvind
        
        [Sm,Sn]=size(model.S);
        b=sparse(Sm,1,0);
        %growth=find(model.c);
        growth = strmatch('biomass',model.rxns);
        
        glcind=strmatch('EX_glc',model.rxns);
        glycind=strmatch('EX_glyc',model.rxns);
        o2ind = strmatch('EX_o2',model.rxns);
        model.glcind = glcind;
        model.o2ind=o2ind;
        model.growth=growth;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Reduce bounds to +/- 100
        toobigvu=abs(model.vu)>100; % for the unlikely case that vu<-100
        model.vu(toobigvu)=sign(model.vu(toobigvu))*100;
        toobigvl=abs(model.vl)>100;
        model.vl(toobigvl)=sign(model.vl(toobigvl))*100;
        
        model.vl(glcind)=-10;
        model.vl(o2ind)=-10;
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

    case 'imm904rev_anaer'
        % S. cerevisiae without adipate pathway
        load /home/chrisg/research/phyto/x/models/iMM904rev/iMM904rev_anaer             
        load /home/chrisg/research/phyto/x/models/iMM904rev/iMM904rev_anaer_RREF R pvind
        
        model = iMMr_a;
        model.vl = model.lb;
        model.vu = model.ub;
        
        % Model with adipate pathway added
        %load iMM904ad
        %load RREFiMM904ad R pvind
        
        [Sm,Sn]=size(model.S);
        b=sparse(Sm,1,0);
        %growth=find(model.c);
        growth = strmatch('biomass_wild',model.rxns);
        
        glcind=strmatch('EX_glc',model.rxns);
        glycind=strmatch('EX_glyc',model.rxns);
        o2ind = strmatch('EX_o2',model.rxns);
        model.glcind = glcind;
        model.o2ind=o2ind;
        model.growth=growth;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Reduce bounds to +/- 100
        toobigvu=abs(model.vu)>100; % for the unlikely case that vu<-100
        model.vu(toobigvu)=sign(model.vu(toobigvu))*100;
        toobigvl=abs(model.vl)>100;
        model.vl(toobigvl)=sign(model.vl(toobigvl))*100;
        
        model.vl(glcind)=-10;
        model.vl(o2ind)=-10;
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
        
    case 'imm904rev'
        % S. cerevisiae without adipate pathway
        load /home/chrisg/research/phyto/x/models/iMM904rev/iMM904rev             
        load /home/chrisg/research/phyto/x/models/iMM904rev/iMM904rev_RREF R pvind
        
        model = iMMr;
        model.vl = model.lb;
        model.vu = model.ub;
        
        % Model with adipate pathway added
        %load iMM904ad
        %load RREFiMM904ad R pvind
        
        [Sm,Sn]=size(model.S);
        b=sparse(Sm,1,0);
        %growth=find(model.c);
        growth = strmatch('biomass_wild',model.rxns);
        
        glcind=strmatch('EX_glc',model.rxns);
        glycind=strmatch('EX_glyc',model.rxns);
        o2ind = strmatch('EX_o2',model.rxns);
        model.glcind = glcind;
        model.o2ind=o2ind;
        model.growth=growth;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Reduce bounds to +/- 100
        toobigvu=abs(model.vu)>100; % for the unlikely case that vu<-100
        model.vu(toobigvu)=sign(model.vu(toobigvu))*100;
        toobigvl=abs(model.vl)>100;
        model.vl(toobigvl)=sign(model.vl(toobigvl))*100;
        
        model.vl(glcind)=-10;
        model.vl(o2ind)=-10;
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







