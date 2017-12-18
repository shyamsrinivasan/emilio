function [sets,plist,vl,vu,activevl,activevu,v,KKTviol,exitflag]=EMILiO(model,growth,growthmin,prodind,blacklistEM,R,pvind,noConc,KpY,ICEM,LPpruneFrac,MILPpruneFrac,epsProd,epsProdILP,ncuts,slackTol,vpminfrac,foldername)
%function [sets,plist,vl,vu,activevl,activevu,v,KKTviol,exitflag]=EMILiO(model,growth,growthmin,prodind,blacklistEM,R,pvind,noConc,KpY,ICEM,LPpruneFrac,MILPpruneFrac,epsProd,epsProdILP,ncuts,slackTol,vpminfrac,foldername)

if nargin<13
    epsProd=1e-3;
end
if nargin<14
    epsProdILP=epsProd;
end
if nargin<15
    ncuts=1;
end
if nargin<16
    slackTol=1e-5;
end
if nargin<17
    vpminfrac=LPpruneFrac;
end
if nargin<18
    foldername=[];
end

nmodsMax = 100;  % If this number is too high, then Phase 3 of EMILiO might take a long time.

[Sm,Sn]=size(model.S);
sets=[]; plist=[]; vl=[]; vu=[]; activevl=[]; activevu=[]; v=[]; KKTviol=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First, determine acceptable vProd
cprod = sparse(1,prodind,1,1,Sn);
b=sparse(Sm,1,0);
vl1=model.vl; vu1=model.vu;
vl1(growth)=growthmin;
vLP=cplexlp(-cprod(:),[],[],model.S,b,vl1,vu1);
vPtheor = vLP(prodind);
prodmin = LPpruneFrac * vPtheor;    % Allow strategies producing at least prodmin

%%%%%%%%%%%%%%
% Start EMILiO
starttime =tic;
exitflag = 0;
[v,vl,vu,activevl,activevu,KKTviol,exitflag]=doGDILP;

switch exitflag
        case 0
            fprintf('SLP couldn''t find a feasible solution meeting production requirements.\n');
            fprintf('Try lowering production requirements\n');
        case 1                   %Successful ILP
            plist=doLPpruning;   %Returns all alternate plists 
            sets=doMILPpruning;  %Generates alternate optima for all plists
        case 2
            fprintf('SLP solution inconsistent with FBA. Try increasing epsProd.\n');
        case -1
            fprintf('Convex relaxation suggests that the solution cannot be found\n');
            fprintf('Terminating SLP\n');
end

elapsed = toc(starttime);
fprintf('Total strain design(s) completed in %g seconds\n',elapsed);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The three stages of EMILiO

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [v,vl,vu,activevl,activevu,KKTviol,exitflag] = doGDILP
        % Start GDILP
        gdilpstart=tic;
        [v,vl,vu,activevl,activevu,KKTviol,exitSLP]=GDILPLS( ...
            model.S,model.vl,model.vu,growth,growthmin,prodind,blacklistEM,R,pvind,noConc,KpY,ICEM,epsProdILP,slackTol,vpminfrac);
        if exitSLP == 1
            prodmax = v(prodind);
            gdilptime=toc(gdilpstart);
            fprintf('Finished GDILP in %g seconds\n',gdilptime);
            
            % Confirm that GDILP solution works
            vl( abs(vl)<1e-9) = 0;  % Avoid numerical issues
            vu( abs(vu)<1e-9) = 0;
            
            % Find maximum growth rate
            c=sparse(1,growth,1,1,Sn);
            vLP = cplexlp(-c(:),[],[],model.S,b,vl,vu);
            growthmax = vLP(growth);
            vl2=vl; vu2=vu;
            vl2(growth)=growthmax;
            vu2(growth)=growthmax;
            
            % Min/max vprod at max growth
            [vLP,fval] = cplexlp(cprod(:),[],[],model.S,b,vl2,vu2);
            vProdMin = fval;
            [vLP,fval] = cplexlp(-cprod(:),[],[],model.S,b,vl2,vu2);
            vProdMax = -fval;
            fprintf('Production range: %g -- %g\n',[vProdMin,vProdMax]);
            if vProdMin < prodmin
                exitflag = 2;   % Flag to suggest changing epsProd
            else
                exitflag = 1;
            end
        elseif exitSLP == 0
            exitflag = -1;   % Meaning, terminate
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [plist,nLists] = doLPpruning        
        LPstart=tic;
        plist = prunedesign(model,vl,vu,growth,growthmin,prodind,prodmin);
        LPtime=toc(LPstart);                
        nLists = length(plist);
        fprintf('Identified %g active subsets using recursive LP in %g seconds\n',[nLists LPtime]);
        badplist = logical(zeros(nLists,1));
        for i=1:nLists  % Re-format plist
            % Double-check that modifications are present            
            badplist(i) = isempty(plist{i}.activevl) && isempty(plist{i}.activevu);
            plist{i}.vl=model.vl; plist{i}.vu=model.vu;
            plist{i}.vl(plist{i}.activevl)=vl(plist{i}.activevl);
            plist{i}.vu(plist{i}.activevu)=vu(plist{i}.activevu);
        end
        plist(badplist) = [];   % Remove bad plists
        nLists = length(plist);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function sets = doMILPpruning
        % MILP-based parsimonious and alternate designs        
        sets={};
        milpstart=tic;
        nLists = length(plist);
        k=0;
        for i=1:nLists
            % Determine production using all bounds
            cboth=sparse([1 1],[growth prodind],[1 -epsProd],1,Sn);
            b=sparse(Sm,1,0);
            [vLP,fval,output]=cplexlp(-cboth(:),[],[],model.S,b,vl,vu);
            vprodmin=MILPpruneFrac*vLP(prodind);
            
            plist2=plist{i};
            plist2.vld=plist2.vl(plist2.activevl);
            plist2.vud=plist2.vu(plist2.activevu);
            
            % Check that the pruned set is not too large
            nmods = length(plist2.vld)+length(plist2.vud);
            if nmods > nmodsMax
                fprintf('Pruned modification set %g contains too many (%g) mods.\n',[i,nmods]);                
                fprintf('Not aborting, but this may take a while.\n');
                %fprintf('Pruned modification set %g contains too many (%g) mods. Aborting.\n',[i,nmods]);
                %break;
            end
            
            % Alternate/Minimal strain designs using MILP+Int.Cuts
            altsets = genAltStrainsMILPCPXINT(...
                model,plist2,growth,growthmin,prodind,ncuts,vprodmin,...
                epsProdILP,foldername);
            naltsets=length(altsets);
            for j=1:naltsets
                k=k+1;
                sets{k}=altsets{j};
                fprintf('Finished MILP-based design pruning for %g design(s)\n',k);
            end            
        end
        milptime=toc(milpstart);
        fprintf('Finished MILP-based pruning in %g seconds\n',milptime);                
        
        % Sort strains by number of modifications
        nsets = length(sets);
        nMod = zeros(1,nsets);
        for i=1:nsets
            nMod(i)=length(sets{i}.activevl)+length(sets{i}.activevu);
        end
        [sortMod,sortInd]=sort(nMod);
        sets=sets(sortInd);                
        
        % Format strain design as activevl, activevu, KO
        KOtol = 1e-5;
        for i=1:nsets
            sets{i}.KO=[];
            KOs=(abs(vl(sets{i}.activevl))<KOtol) & (abs(model.vu(sets{i}.activevl))<KOtol);
            sets{i}.KO=sets{i}.activevl(KOs);
            sets{i}.activevl(KOs)=[];
            KOs=(abs(model.vl(sets{i}.activevu))<KOtol) & (abs(vu(sets{i}.activevu))<KOtol);
            sets{i}.KO=union(sets{i}.KO,sets{i}.activevu(KOs));
            sets{i}.activevu(KOs)=[];
            sets{i}.vl=vl;
            sets{i}.vu=vu;
        end
    end
end