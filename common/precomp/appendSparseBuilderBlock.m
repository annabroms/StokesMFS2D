function block = appendSparseBuilderBlock(block,row_idx,col_idx,A)
%APPENDSPARSEBUILDERBLOCK Append one dense block to a sparse accumulator.

row_idx = row_idx(:);
col_idx = col_idx(:);
nr = numel(row_idx);
nc = numel(col_idx);
n = nr*nc;
if n == 0
    return
end

if block.next+n-1 > numel(block.vals)
    block = flushSparseBuilderBlock(block);
end

if n > numel(block.vals)
    block.S = block.S + sparse(repmat(row_idx,nc,1), ...
        repelem(col_idx,nr),A(:),size(block.S,1),size(block.S,2));
    return
end

loc = block.next:block.next+n-1;
block.rows(loc) = repmat(row_idx,nc,1);
block.cols(loc) = repelem(col_idx,nr);
block.vals(loc) = A(:);
block.next = block.next + n;
end
