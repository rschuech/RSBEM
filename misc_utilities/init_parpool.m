function [] = init_parpool(max_threads, requested_threads)
% makes sure a parpool is open with at least requested_threads workers, if
% achievable


pool_h = gcp('nocreate');

if isempty(gcp('nocreate')) %no pool started already
    parpool( min(max_threads, requested_threads) );  %no point in starting more workers than there are parfor iterations
elseif pool_h.NumWorkers < min(max_threads, requested_threads) %we want more workers than we have already
    delete(gcp);
    parpool( min(max_threads, requested_threads) );  %no point in starting more workers than there are parfor iterations
end

