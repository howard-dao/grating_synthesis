function [ max_elem, i_row, i_col ] = f_get_max_matrix_elem( A )
% returns maximum value and indices of largest element of matrix A

[ max_elem, indx ]    = max(A(:));
[ i_row, i_col ]      = ind2sub(size(A),indx);

end

