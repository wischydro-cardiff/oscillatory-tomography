function result = toepmat_vector_math(toepmat_row,operation,varargin)

%toepmat_vector_math: Function for efficiently performing computations
%utilizing distance or covariance matrices that have a Toeplitz or
%block-Toeplitz structure (such as those produced when evaluating
%covariance matrices for a regular, equispaced grid. The function
%performs 4 different types of computations:
%  1) Toeplitz matrix - vector product
%  2) Toeplitz matrix inverse - vector product (BETA: may not work in all
%  cases)
%  3) Eigenvalues of toeplitz-embedded circulant matrix (primarily for
%  testing / debugging purposes)
%  4) Generation of unconditional realizations (2 per run) with given
%  covariance structure
%
%[result] = toepmat_vector_math(toepmat_row, operation, {v_input,
%toep_level, disc_vec, inv_reg})
%
%Where:
%   -result is the output result, dependent on operation
%      1) Multiplication: result is an (n x 1) output result of
%      matrix-vector multiplication, where n is the dimension of toepmat_row
%      2) Division: result is an (n x 1) output result of matrix
%      inverse-vector multiplication.
%      3) Eigenvalues: result is a (2(i - 1) x 2(j-1) x 2(k-1)) vector where
%      i,j,k are the dimensions of the grid the toeplitz matrix is
%      representing
%      4) Realizations: result is a (n x 2) output containing 2 conditional
%      realizations of a random field with the given toeplitz covariance
%      matrix
%    -toepmat_row (1 x n vector) is the first row/column of a toeplitz or
%    block-toeplitz matrix
%    -operation (string, 1 char) is the operation to be performed (*, \, e,
%    or r for the various actions)
%    -v_input (n x 1) is the vector to be multiplied / divided by (if
%    multiplication or division is the operation)
%    -toeplevel (scalar) defines what level toeplitz structure is given. 1
%    = toeplitz matrix, 2 = toeplitz block toeplitz, 3 = toeplitz block
%    toeplitz block toeplitz. DEFAULT = 1
%    -disc_vec (d x 1), where d is the dimension of the problem (same a
%    toeplevel) defines the number of elements along each dimension.
%    -inv_reg (scalar) (scalar) is a constant used to stabilize matrix
%    inversion for cases where a positive-definite matrix has small
%    negative eigenvalues.
%
% Code by Michael Cardiff
% 8/2009, Last Updated: 12/2014 

vecsize = numel(toepmat_row);

toep_level = 1;
disc_vec = [];
inv_reg = 0;
v_input = [];

[v_input, toep_level, disc_vec, inv_reg] = process_extra_args(varargin,v_input,toep_level,disc_vec,inv_reg);

if isempty(disc_vec)
    disc_vec = round(vecsize.^(1/toep_level)).*ones(toep_level,1);
end

if (strcmp(operation,'*') || strcmp(operation,'\')) && isempty(v_input)
    error('Vector to be multiplied / divided was not supplied. Supply v_input after operation');
end

if toep_level == 1
    circmat_row = [toepmat_row toepmat_row((end-1):-1:2)];
    if strcmp(operation,'*')
        v_padded = [v_input; v_input((end-1):-1:2)]; 
        fft_product = fft(circmat_row').*fft(v_padded);
        result_padded = ifft(fft_product);
        result = result_padded(1:disc_vec(1),:);
    elseif strcmp(operation,'\')
        v_padded = [v_input; v_input((end-1):-1:2)]; 
        fft_product = fft(v_padded)./(fft(circmat_row')+inv_reg);
        result_padded = ifft(fft_product);
        result = result_padded(1:disc_vec(1),:);
    elseif strcmp(operation,'e')
        result = fft(circmat_row');
    elseif strcmp(operation,'r')
        n = size(circmat_row,2);
        e1 = randn(n,1);
        e2 = randn(n,1);
        ec = e1 + 1i*e2;
        fft_product = fft(circmat_row'.*n).^.5.*ec; 
        result_padded = ifft(fft_product);
        result = result_padded(1:disc_vec(1),:);
        result = [real(result) imag(result)];
    end
elseif toep_level == 2
    circmat_rowmat = reshape(toepmat_row,disc_vec(1),disc_vec(2));
    circmat_rowmat = cat(2,circmat_rowmat, circmat_rowmat(:,(end-1):-1:2));
    circmat_rowmat = cat(1,circmat_rowmat, circmat_rowmat((end-1):-1:2,:));
    if strcmp(operation,'*')
        v_padded = reshape(v_input,disc_vec(1),disc_vec(2));
        v_padded = cat(2,v_padded, zeros(disc_vec(1),disc_vec(2)-2));
        v_padded = cat(1,v_padded, zeros(disc_vec(1)-2,2*disc_vec(2)-2));
        fft_product = fft2(circmat_rowmat).*fft2(v_padded);
        result_padded = ifft2(fft_product);
        result = result_padded(1:disc_vec(1),1:disc_vec(2));
        result = reshape(result,numel(result),1);
    elseif strcmp(operation,'\')
        v_padded = reshape(v_input,disc_vec(1),disc_vec(2));
        v_padded = cat(2,v_padded, zeros(disc_vec(1),disc_vec(2)-2));
        v_padded = cat(1,v_padded, zeros(disc_vec(1)-2,2*disc_vec(2)-2));
        fft_product = fft2(v_padded)./(fft2(circmat_rowmat)+inv_reg);
        result_padded = ifft2(fft_product);
        result = result_padded(1:disc_vec(1),1:disc_vec(2));
        result = reshape(result,numel(result),1);
    elseif strcmp(operation,'e')
        result = fft2(circmat_rowmat);
        result = reshape(result,numel(result),1);
    elseif strcmp(operation,'r')
        n = numel(circmat_rowmat);
        cdims = size(circmat_rowmat);
        e1 = randn(cdims);
        e2 = randn(cdims);
        ec = e1 + 1i*e2;
        fft_product = fft2(circmat_rowmat.*n).^.5.*ec; 
        result_padded = ifft2(fft_product);
        result = result_padded(1:disc_vec(1),1:disc_vec(2));
        result = reshape(result,numel(result),1);
        result = [real(result) imag(result)];
    end
elseif toep_level == 3
    circmat_rowmat = reshape(toepmat_row,disc_vec(1),disc_vec(2),disc_vec(3));
    circmat_rowmat = cat(2,circmat_rowmat, circmat_rowmat(:,(end-1):-1:2,:));
    circmat_rowmat = cat(1,circmat_rowmat, circmat_rowmat((end-1):-1:2,:,:));
    circmat_rowmat = cat(3,circmat_rowmat, circmat_rowmat(:,:,(end-1):-1:2));
    if strcmp(operation,'*')
        v_padded = reshape(v_input,disc_vec(1),disc_vec(2),disc_vec(3));
        v_padded = cat(2,v_padded, zeros(disc_vec(1),disc_vec(2)-2,disc_vec(3)));
        v_padded = cat(1,v_padded, zeros(disc_vec(1)-2,2*disc_vec(2)-2,disc_vec(3)));
        v_padded = cat(3,v_padded, zeros(2*disc_vec(1)-2,2*disc_vec(2)-2,disc_vec(3)-2));
        fft_product = fftn(circmat_rowmat).*fftn(v_padded);
        result_padded = ifftn(fft_product);
        result = result_padded(1:disc_vec(1),1:disc_vec(2),1:disc_vec(3));
        result = reshape(result,numel(result),1);
    elseif strcmp(operation,'\')
        v_padded = reshape(v_input,disc_vec(1),disc_vec(2),disc_vec(3));
        v_padded = cat(2,v_padded, zeros(disc_vec(1),disc_vec(2)-2,disc_vec(3)));
        v_padded = cat(1,v_padded, zeros(disc_vec(1)-2,2*disc_vec(2)-2,disc_vec(3)));
        v_padded = cat(3,v_padded, zeros(2*disc_vec(1)-2,2*disc_vec(2)-2,disc_vec(3)-2));
        fft_product = fftn(v_padded)./(fftn(circmat_rowmat) + inv_reg);
        result_padded = ifftn(fft_product);
        result = result_padded(1:disc_vec(1),1:disc_vec(2),1:disc_vec(3));
        result = reshape(result,numel(result),1);
    elseif strcmp(operation,'e')
        result = fftn(circmat_rowmat);
        result = reshape(result,numel(result),1);
    elseif strcmp(operation,'r')
        n = numel(circmat_rowmat);
        cdims = size(circmat_rowmat);
        e1 = randn(cdims);
        e2 = randn(cdims);
        ec = e1 + 1i*e2;
        %Appears that a multiplication correction factor is needed here, unlike presented in Dietrich, Nowak (who do division).
        fft_product = fftn(circmat_rowmat.*n).^.5.*ec; 
        result_padded = ifftn(fft_product);
        result = result_padded(1:disc_vec(1),1:disc_vec(2),1:disc_vec(3));
        result = reshape(result,numel(result),1);
        result = [real(result) imag(result)];
    end
end
    