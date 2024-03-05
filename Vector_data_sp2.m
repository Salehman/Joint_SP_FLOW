function [data]=Vector_data_sp2(OB,states,NT)
n=1;
for h=1:length(OB)
    for j=1:NT
    data(n)=states{j}(OB(h));
    n=n+1;
    end
end