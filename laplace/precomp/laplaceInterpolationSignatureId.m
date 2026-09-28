function identifier = laplaceInterpolationSignatureId(signature)
%LAPLACEINTERPOLATIONSIGNATUREID Stable short identifier for a struct.

bytes = uint8(jsonencode(signature));
identifier = [fnv32(bytes) fnv32(fliplr(bytes))];
end

function value = fnv32(bytes)
hash = uint32(2166136261);
prime = uint64(16777619);
modulus = uint64(2)^32;
for k = 1:numel(bytes)
    hash = bitxor(hash,uint32(bytes(k)));
    hash = uint32(mod(uint64(hash)*prime,modulus));
end
value = lower(dec2hex(hash,8));
end
