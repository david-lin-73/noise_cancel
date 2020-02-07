function [ProcMatrix, Mean, CovMatrix] = preprocessing(Matrix)
    [MatrixCent, Mean] = noises.ica.center_rows(Matrix);
    [ProcMatrix, CovMatrix] = noises.ica.whitening(MatrixCent);
end