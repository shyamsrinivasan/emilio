function [strainset,strainsetind]=uniquestrains(model,strains)
% function [strainset,strainsetind]=uniquestrainset(model,strains)
%
% Identify the unique strain sets
% Laurence Yang. Aug 5, 2010

nstrains = length(strains);
strainsetind = 1:nstrains;

[mods,straindef]=uniquemods(model,strains);

% If 2 identical strains are found, then they may still differ in terms of
% their quantitative vl and vu values. So, keep the better strain.
for i=1:nstrains-1
    for j= i+1:nstrains
        if isempty(setxor(straindef{i},straindef{j}))   % If both strains have identical activevl & activevu            
            strainsetind(j)=strainsetind(i);
            % Replace the i strain with the better one found so far
            if strains{j}.vprod > strains{i}.vprod                
                strains(i) = strains(j);
            end
        end        
    end
end
strainset=strains(unique(strainsetind));