function block = flushSparseBuilderBlock(block)
%FLUSHSPARSEBUILDERBLOCK Flush buffered triplets into the sparse matrix.

n = block.next - 1;
if n > 0
    block.S = block.S + sparse(block.rows(1:n),block.cols(1:n), ...
        block.vals(1:n),size(block.S,1),size(block.S,2));
    block.next = 1;
end
end
