function compatible = isLaplaceCmapSignatureCompatible( ...
        actual,expected,problem)
%ISLAPLACECMAPSIGNATURECOMPATIBLE Compare current and legacy signatures.

try
    actual = normalizeLaplaceCmapInterpolationSignature(actual,problem);
    expected = normalizeLaplaceCmapInterpolationSignature(expected,problem);
    compatible = isequaln(actual,expected);
catch
    compatible = false;
end
end
