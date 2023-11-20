function result = multi_process_inverse(matrix)
    % Create a parallel pool
    

    % Initialize the result matrix
    result = zeros(size(matrix));

    % Use parfor for parallel computation
    parfor i = 1:size(matrix, 3)
        result(:, :, i) = inv(matrix(:, :, i));
    end

    % Close the parallel pool
    %parpool('close');
end
