function [opt,mode] = configureLaplaceCapacitanceInterpolation(opt)
%CONFIGURELAPLACECAPACITANCEINTERPOLATION Compatibility capacitance wrapper.

[opt,mode] = configureLaplaceCmapInterpolation(opt,'capacitance');
end
