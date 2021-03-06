function [vbest,vlbest,vubest,activevl,activevu,slack,exitSLP] = GDILPLS(S,vl,vu,growth,growthmin,prodind,fixed,R,pvind,noConc,Kp,IC,epsProd,slackTol,vpminfrac)
% function [vbest,vlbest,vubest,activevl,activevu,slack] = GDILPLS(S,vl,vu,growth,growthmin,prodind,fixed,R,pvind,noConc,Kp,IC,epsProd,slackTol,vpminfrac)
% Inout arguments:
% S, vl,vu - obtained from the in-silico model of organism
% growth 
% growthmin - User defined parameter for minimum required growth rate 
% prodind - Target metabolite's index in model.rxns
% blacklistEM - 
% R, pvind - 
% 
% Laurence Yang Dec. 10, 2009 -> August 3, 2010
%
% Added convex relaxation

if nargin<7
    fixed=[];
end
if nargin<10
    noConc=0;
end
stepTol = 0;     % Step below stepTol indicates convergence of algorithm
wTol = 1e-3;        % Duals greater than wTol are considered active constraints
Kkkt=1; % Objective weight to minimize KKT violation
if nargin<11
    Kp=1000;            % Objective weight to maximize production rate
end
if nargin<12
    IC=[];
end
if nargin<13
    epsProd= 1e-3;         % Constant to minimize vprod in primal objective, relative to max growth
end
if nargin<14
    slackTol = 1e-6;    % Complementarity constraint violation tolerance    
end
if nargin<15    
    vpminfrac=0.99;
end
vbest=[];       % initial definitions for function output arguments
vlbest=[];
vubest=[];
activevl=[];
activevu=[];
exitSLP = 1;    % exitflag starts at okay. If solution not found, changes to 0.
Kw=0;           % Objective weight to minimize active constraints
maxIter = 30;
sameTol=1e-3;   % Modified bounds must differ from original bounds by at least sameTol to be considered modified
objTol = 1e-6;
epsGrowth = 1;
resLS = 1000;    % Resolution for line search

maxflux=max([abs(vl(:)); abs(vu(:))]);
KKTviolMax = 4*maxflux*maxflux;   % Bound on s, the auxiliary variables ==> KKT violation

[Sm,Sn]=size(S);   %Length of S - Sm metabolites (rows);Sn fluxes (columns)

% First, assess maximum possible production rate = maxProd
% maxProd is assessed by solving an LP with constrainsts that subject
% minimum growth to 'growthmin'
% 
cprod = sparse(1,prodind,-1,1,Sn);
vl2=vl; vu2=vu;
vu2(prodind)=maxflux;
vl2(growth)=growthmin;
v=cplexlp(cprod(:),[],[],S,sparse(Sm,1,0),vl2,vu2);
maxProd = abs(cprod*v);
fprintf('Maximum possible production rate for growth rate of %g h-1: %g\n',[growthmin,maxProd]);

vprodmin = abs(vpminfrac*maxProd);

m = 2*Sn;
if nargin<8 %reduced row echelon form of S and a sparse column matrix b
    [R,pvind]=rref(full([S sparse(Sm,1,0)])); %length(pvind) gives rank of R
end
[TX,rX,pvind,freeind]=projRREF(R,pvind);
nfree = length(freeind); npv=length(pvind);

vflength = nfree; %
vulength = Sn;
vllength = Sn;
mulength = 2*Sn;
wlength = 2*Sn;
slength = 2*Sn;

vfstart = 0;
vustart = vfstart+vflength;
vlstart = vustart+vulength;
mustart = vlstart+vllength;
wstart = mustart+mulength;
sstart = wstart+wlength;

nX = nfree+Sn+Sn+2*Sn+2*Sn;   %[dvf; dvu; dvl; dmu; dw]
nY = vflength+vulength+vllength+mulength+wlength+slength;
   % [dvf; dvu; dvl; dmu; dw; s]
N = nY;
Nones = ones(1,nX);
pvf = sparse(1:vflength,vfstart+(1:vflength),1,vflength,N);
pvfX = sparse(1:vflength,vfstart+(1:vflength),1,vflength,nX); 
pvu = sparse(1:vulength,vustart+(1:vulength),1,vulength,N);
pvl = sparse(1:vllength,vlstart+(1:vllength),1,vllength,N);
pvuX = sparse(1:vulength,vustart+(1:vulength),1,vulength,nX);
pvlX = sparse(1:vllength,vlstart+(1:vllength),1,vllength,nX);
pmu = sparse(1:mulength,mustart+(1:mulength),1,mulength,N);
pmuX = sparse(1:mulength,mustart+(1:mulength),1,mulength,nX);
pw = sparse(1:wlength,wstart+(1:wlength),1,wlength,N);
pwX = sparse(1:wlength,wstart+(1:wlength),1,wlength,nX);
ps = sparse(1:slength,sstart+(1:slength),1,slength,N);

pd = sparse(1:nX,1:nX,1,nX,N);     % d=pd*X

%%%%%%%%%%%%%%%
% Maximin   % VERY IMPORTANT: Prevent alternate optimal production rate
% that is actually small. Hence, we are actually doing maximin
cboth=sparse([1 1],[growth prodind],[epsGrowth -epsProd],1,Sn);
c=cboth;
%%%%%%%%%%%%%%%

A = [-TX; TX];  % v = rX-TX*vf

E=sparse(1:2*Sn,nfree+2*Sn+(1:2*Sn),1,2*Sn,nX);     % dmu = E*d
F=sparse(1:2*Sn,nfree+4*Sn+(1:2*Sn),1,2*Sn,nX);     % dw = F*d

vul=vl; vll=vl;
vuu=vu; vlu=vu;
% Irreversible reactions stay irreversible
Irrev = vl==0;
vll(Irrev) = 0;
Irrev = vu==0;  %Irreversible in reverse direction
vuu(Irrev) = 0;

% If noConc, disallow forcing reversible reactions.
% This assumes we provide "true" vmin,vmax using fva
if noConc
    rev = (vu>0) & (vl<0);
    vul(rev)=0;
    vlu(rev)=0;
end

% Prevent active constraint of production flux
vul(prodind)=maxflux; vuu(prodind)=maxflux;
vll(prodind)=0; vlu(prodind)=0;

% Other fixed bounds. Note, vld>=vl & vud<=vu.
% Hence, fluxes like ATPM won't get manipulated
if not(isempty(fixed))
    vul(fixed)=vu(fixed);
    vuu(fixed)=vu(fixed);
    vll(fixed)=vl(fixed);
    vlu(fixed)=vl(fixed);
end

% Set minimum growth rate
vul(growth)=vu(growth); vuu(growth)=vu(growth);
vll(growth)=growthmin; vlu(growth)=growthmin;

xl=[vl(freeind); vul; vll; sparse(4*Sn,1,0)];
xu=[vu(freeind); vuu; vlu; 2*maxflux*ones(4*Sn,1)];

%%%%%%%%%%%%%%%%%%%%
% Initial solution %
% Use McCormick convex relaxations of bilinear constraints and find an
% initial point that is, hopefully close to optimum
f = Kkkt*ones(1,2*Sn)*ps + Kp*cprod*-TX*pvf;
[X,Y,exitflag]=McCormick(f,xl,xu);
if exitflag==1    
    fprintf('Found initial solution using McCormick relaxations\n');    
    fprintf('Relaxed optimal estimated KKT = %g\n',sum(abs(ps*Y)));
    fprintf('Relaxed optimal actual KKT = %g\n',sum(abs((E*X).*(F*X))));
    fprintf('Relaxed optimal Vprod = %g\n',cprod*(rX-TX*pvfX*X));
else
    X=(xl+xu)/2;
    fprintf('Failed to find initial solution\n');    
end
if not(isempty(IC))
    X(1:nfree)=(IC.vl(freeind)+IC.vu(freeind)) / 2;
end
%%%%%%%%%%%%%%%%%%%%

vprod = abs(cprod*(rX-TX*pvfX*X));
slack = 1000;
optimal=0;
objs=[];
set(gcf,'Color','w');
tic
while not(optimal)
    step = 100; Xbest = X; objBest = Inf; iter=0; %  Reinitialize
    f = Kkkt*ones(1,2*Sn)*ps + Kp*cprod*-TX*pvf + Kw*ones(1,2*Sn)*pw;  % Outer objective function updated
    while iter < maxIter && abs(step) > stepTol
        iter = iter+1;
                        
        % Re-construct constraint matrices at the new solution
        EX = E*X;
        FX = F*X;
        b = [pvuX*X-rX; rX-pvlX*X];
        
        Aineq = [ (FX(:,Nones).*E + EX(:,Nones).*F)*pd-ps];
        bineq = [-EX.*FX];
        
        AeqY = [A*pvf+pmu-[pvu;-pvl]; A'*pw];
        beqY = [b-A*pvfX*X-pmuX*X; (c*-TX)'-A'*pwX*X];
                       
        % Re-construct bounds on delX and s at the new solution
        yl = [xl-X; zeros(slength,1)];
        yu = [xu-X; KKTviolMax*ones(slength,1)];
        
        % Determine if we can satisfy KKT s.t. target production rate
        % Convex relaxation of bilinear constraints
        
        
        % Solve LP for SLP
        [Y,fval,exitflag,details] = cplexlp(f',Aineq,bineq,AeqY,beqY,yl,yu);
        if exitflag ~= 1
            fprintf('Stopping ILP due to the following CPLEX status:\n');
            fprintf([details.message,'\n']);
            break;
        else
            d = Y(1:nX);
            %X = X + d;
            % Line search necessary for convergence to local optima
            % Note: this depends on the objective function, f
            exfx=sum((E*X).*(F*X));
            edfd=sum((E*d).*(F*d));
            exfd=sum((E*X).*(F*d));
            fxed=sum((F*X).*(E*d));
            delvprod = cprod*-TX*pvfX*d;    % rX gets canceled out. Change in vprod
            sumdelw = sum(pwX*d);
            a0=Kkkt*edfd;
            b0=Kkkt*(exfd+fxed)-Kp*delvprod;
            c0=Kkkt*exfx-Kp*cprod*(rX-TX*pvfX*X);
            %%%%%%%%%%%%%%%
            % Line search %         
            lambdas=linspace(0,1,resLS);
            Zs=a0*(lambdas.^2) + b0*lambdas + c0;
            minZ=min(Zs);
            lambdamin=lambdas(Zs==minZ);                        
            step = lambdamin(1);            
            %%%%%%%%%%%%%%%%%%%%%%%%%
            % If step size is 0, then SLP has converged.
            % However, if KKT is still violated, then we can try to escape
            % this local optimum and continue the SLP until next
            % convergence. For this, we should perturb the solution/objective a bit.
            %if step < 1e-8 && slack > slackTol
            if step < 1e-8 && not(optimal)
                % 1. Simplest heuristic without much theoretical basis
                %step = 1;
                %fprintf('Converged to local optimum but KKT still violated. Continuing SLP\n');
                
                % 2. Use McCormick relaxation in deviation variables
                % relative to current point to determine relaxed KKT
                
                % 3. Use McCormick with bounds tightened to be near the
                % current solution and determine relaxed KKT
                nrel = 5;
                Xms = sparse(nX,nrel);
                Objs = sparse(1,nrel);
                krel=0;
                for pfrac=linspace(0.1,0.9,nrel)
                    krel=krel+1;
                    %pfrac = 0.2;    % Fraction of feasible range to consider                    
                    xl2 = X-pfrac*(X-xl);
                    xu2 = X+pfrac*(xu-X);
                    %f2 = Kkkt*ones(1,2*Sn)*ps; %+ Kp*cprod*-TX*pvf;
                    f2 = Kkkt*ones(1,2*Sn)*ps + Kp*cprod*-TX*pvf;   % cprod already has -1
                    [X2,Y2,exitflag]=McCormick(f2,xl2,xu2);
                    if exitflag==1
                        KKT=sum(abs((E*X2).*(F*X2)));
                        KKTlin=sum(abs(ps*Y2));
                        vprod2 = abs(cprod*(rX-TX*pvfX*X2));
                        Xms(:,krel)=X2;
                        Obj = Kkkt*KKT - Kp*vprod2;
                        Objs(krel) = Obj;
                        fprintf('McCormick relaxed optimal KKT near current solution is %g\n',KKT);
                        fprintf('McCormick relaxed optimal linearized KKT near current solution is %g\n',KKTlin);
                        fprintf('McCormick relaxed optimal Vprod near current solution is %g\n',vprod2);
                    else
                        fprintf('Relaxed problem is infeasible\n');
                    end
                end
                bestObj = min(Objs);
                bestrel = (Objs == bestObj);
                
                % See if relaxed objective is better than current value
                ObjNow = Kkkt*slack - Kp*vprod;
                if abs(ObjNow - bestObj) < 1e-6
                    fprintf('Objective cannot be improved further because no relaxed objective better than current SLP solution exists\n');                    
                    fprintf('Try resolving SLP from different initial solution\n');
                    fprintf('Stopping SLP\n');
                    exitSLP = 0;
                    return;
                else
                    X=Xms(:,bestrel);
                    X=full(X);
                    step = 1;
                    fprintf('Better relaxed objective found. Continuing SLP\n');
                end                
            else
                % Update to the new point
                X = X+step*d;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%
        % Print some diagnostics
        vprod = abs(cprod*(rX-TX*pvfX*X));
        slack = sum(abs((E*X).*(F*X)));
        Pviol = sum( [pvuX*X+1e-6;-(pvlX*X-1e-6)] < [-TX;TX]*pvfX*X );
        Dviols = abs((c*-TX)' - A'*pwX*X);
        Dviol = sum(Dviols);
        KKTdual = max(abs(pvfX*X).*Dviols);    % Worst eta*x violation due to dual equality constraint violation
        boundViol = sum( (X<xl-1e-6)+(X>xu+1e-6) );
        fprintf('Production rate: %g (%g%% of max)\n', [vprod 100*vprod/maxProd]);
        fprintf('Complementarity violation: %g\n',slack);
        fprintf('Worst dual complementary violation: %g\n',KKTdual);
        fprintf('Primal constraint violation: %g\n',Pviol);
        fprintf('Dual constraint violation: %g\n',Dviol);
        fprintf('Bound violation: %g\n',boundViol);
        fprintf('Norm of displacement vector: %g\n',norm(d));
        fprintf('Step size: %g\n',step);
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Evaluate objective at new solution
        %obj=slack+Pviol+abs(Dviol)-vprod;        
        %obj = Kkkt*slack - Kp*vprod;
        obj = evalObj(X,E,F);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % See the objective & step size
        objs=[objs obj];
        subplot(2,1,1); 
        cla;
        hold on;
        plot(lambdas,Zs,'g','LineWidth',2);
        scatter(lambdamin,minZ,64,'md','filled');
        line([step step],[minZ max(Zs)]);
        %myarrow([step step],[minZ max(Zs)]);
        %myarrow([step step],[max(Zs) minZ]);
        xlabel('Step size'); ylabel('Objective');
        subplot(2,1,2);
        plot(1:iter,objs,'g','LineWidth',2);
        xlabel('Iteration'); ylabel('Objective');
        drawnow;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Retain solution if it is an improvement
        if obj < objBest
            Xbest=X;
            objBest = obj;
        end
        optimal = (slack<=slackTol) & (vprod>=vprodmin);
    end
    if slack>slackTol
        fprintf('KKT not satisfied\n');
    end
    if vprod < vprodmin
        fprintf('Minimum production requirement not met\n');
    end            
end
toc
vbest = rX-TX*pvfX*Xbest;
vlbest = pvlX*Xbest;
vubest = pvuX*Xbest;
activebounds = pwX*Xbest;
activevu = activebounds(1:Sn) > wTol;       %Sn
activevl = activebounds(Sn+(1:Sn)) > wTol;  %Sn
% Only change active bounds
vlbest(not(activevl))=vl(not(activevl));
vubest(not(activevu))=vu(not(activevu));

% In fact, we can clean results a bit more
% Some bounds might be active, but were already active. e.g.,
% irreversibility constraints or substrate/oxygen uptake constraints
redundvl=abs(vlbest(activevl)-vl(activevl))<sameTol;
redundvu=abs(vubest(activevu)-vu(activevu))<sameTol;
activevl(redundvl)=0;
activevu(redundvu)=0;

vlbest(not(activevl))=vl(not(activevl));    % Re-clean active bounds
vubest(not(activevu))=vu(not(activevu));
activevl=find(activevl);
activevu=find(activevu);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function Z=evalObj(X,E,F)
        exfx=sum((E*X).*(F*X));
        vprod = abs(cprod*(rX-TX*pvfX*X));
        Z = Kkkt*exfx - Kp*vprod;
    end
    
    function [X,Y,exitflag]=McCormick(f,xl,xu)
        EXL = E*xl; FXL = F*xl;
        EXU = E*xu; FXU = F*xu;
        % McCormick relaxations
        Aineq = [(EXL(:,Nones).*F+FXL(:,Nones).*E)*pd-ps;
            (EXU(:,Nones).*F + FXU(:,Nones).*E)*pd-ps;
            ps-(EXU(:,Nones).*F+FXL(:,Nones).*E)*pd;
            ps-(EXL(:,Nones).*F+FXU(:,Nones).*E)*pd];
        bineq = [ EXL.*FXL;
            EXU.*FXU;
            -EXU.*FXL;
            -EXL.*FXU];
        AeqY = [A*pvf+pmu-[pvu;-pvl]; % Primal constraints
            A'*pw]; % Dual constraints
        beqY = [sparse(2*Sn,1,0);
            (c*-TX)'];
        yl = [xl; sparse(slength,1,0)];
        %yl = [xl; -KKTviolMax*ones(slength,1)];
        yu = [xu; KKTviolMax*ones(slength,1)];        
        [Y,fval,exitflag,details] = cplexlp(f(:),Aineq,bineq,AeqY,beqY,yl,yu);        
        if exitflag==1
            X = Y(1:nX);
            fprintf('Found solution using McCormick relaxations\n');
        else
            X=[]; Y=[];
            fprintf('Failed to find solution with McCormick relaxations\n');
        end
    end
end