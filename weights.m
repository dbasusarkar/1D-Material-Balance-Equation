%% Submitted by Debajyoti Basu Sarkar (DBS)
%  For MAT 695 (Fall 2017) offered by Dr. Feng Tian, Hampton University
%-------------------------------------------------------------------------%

function w = weights(x)
    w = harmmean([x(1:end-1), x(2:end)],2);
    w = w';
end
