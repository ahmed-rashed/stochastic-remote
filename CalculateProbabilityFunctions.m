function [Prob,ProbDensity,ProbDist]=CalculateProbabilityFunctions(samples_row,equiSpacedintervals_vec)

N=length(samples_row);
delta=equiSpacedintervals_vec(2)-equiSpacedintervals_vec(1);

num=zeros(1,length(equiSpacedintervals_vec)-1);
for ii=1:length(equiSpacedintervals_vec)-1
    num(ii)=sum(samples_row>equiSpacedintervals_vec(ii) & samples_row<=equiSpacedintervals_vec(ii+1));
end

relativeFrequency=num/N;
Prob=relativeFrequency;
ProbDensity=Prob/delta;
ProbDist=cumsum(ProbDensity)*delta;
