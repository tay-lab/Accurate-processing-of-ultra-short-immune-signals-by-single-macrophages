function [Ikr,Ms]=genReceptors(mm)
rng default
Ms=gamrnd(4.32,1035,[mm 1]);
Ikr=rand(mm,1);
end