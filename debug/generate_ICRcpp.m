function icr = generate_ICRcpp()
%Loads the files written by the ICR::IndependentContactRegions class (if the DEBUG_ICR flag is
%set in the debug.h file).

regions=textread('icr.txt','%s','delimiter','\n');
n_C=size(regions,1);

for i=1:n_C
    ind=strread(regions{i});
    icr(i).ind=ind;
    icr(i).N=numel(ind);
end    




