function [mod,straindef,modKO,modactivevl,modactivevu,modvl,modvu] = uniquemods(model,strains)
% function [mod,straindef,modKO,modactivevl,modactivevu,modvl,modvu] =
% uniquemods(model,strains)
% Generate unique mod list from all strain designs
% 
% Laurence Yang. Feb 12, 2010
nsets = length(strains);


modKO=[]; modactivevl=[]; modactivevu=[]; modvl=[]; modvu=[];


for i=1:nsets        
    %mod.KO=unique([mod.KO; strains{i}.KO(:)]);
    %mod.activevl = unique([mod.activevl; strains{i}.activevl(:)]);
    %mod.activevu = unique([mod.activevu; strains{i}.activevu(:)]);
    modKO = [modKO; strains{i}.KO(:)];
    modactivevl = [modactivevl; strains{i}.activevl(:)];
    modactivevu = [modactivevu; strains{i}.activevu(:)];
    modvl = [modvl; strains{i}.vl(strains{i}.activevl)];
    modvu = [modvu; strains{i}.vu(strains{i}.activevu)];
end

modKO = unique(modKO);
[modactivevl,vlind] = unique(modactivevl);
[modactivevu,vuind] = unique(modactivevu);
nKO = length(modKO);
nvl = length(modactivevl);
nvu = length(modactivevu);
nmod = nKO+nvl+nvu;

for i=1:nKO
    mod(i).rxn=modKO(i); mod(i).KO=modKO(i);
    mod(i).activevl=[]; mod(i).activevu=[];
    mod(i).vl=[]; mod(i).vu=[];
end

for i=1:nvl
    mod(nKO+i).rxn=modactivevl(i); mod(nKO+i).KO=[];
    mod(nKO+i).activevl=modactivevl(i); mod(nKO+i).activevu=[];
    mod(nKO+i).vl=modvl(vlind(i)); mod(nKO+i).vu=[];
end

for i=1:nvu
    mod(nKO+nvl+i).rxn=modactivevu(i); mod(nKO+nvl+i).KO=[];
    mod(nKO+nvl+i).activevl=[]; mod(nKO+nvl+i).activevu=modactivevu(i);
    mod(nKO+nvl+i).vl=[]; mod(nKO+nvl+i).vu=modvu(vuind(i));
end

straindef = cell(nsets,1);
for i=1:nsets
    nKO1 = length(strains{i}.KO);
    nvl1 = length(strains{i}.activevl);
    nvu1 = length(strains{i}.activevu);
    straindef{i}=zeros(1,nKO1+nvl1+nvu1);
    k=0;
    for j=1:nKO1
        k=k+1;
        straindef{i}(k)=find(modKO==strains{i}.KO(j));
    end
    for j=1:nvl1
        k=k+1;
        straindef{i}(k)=nKO+find(modactivevl==strains{i}.activevl(j));
    end
    for j=1:nvu1
        k=k+1;
        straindef{i}(k)=nKO+nvl+find(modactivevu==strains{i}.activevu(j));
    end
end