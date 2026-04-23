function block = initSparseBuilder(capacity,n_rows,n_cols)
%INITSPARSEBUILDER Allocate a buffered sparse triplet accumulator.

capacity = max(0,round(capacity));
block = struct();
block.rows = zeros(capacity,1);
block.cols = zeros(capacity,1);
block.vals = zeros(capacity,1);
block.next = 1;
block.S = sparse(n_rows,n_cols);
end
