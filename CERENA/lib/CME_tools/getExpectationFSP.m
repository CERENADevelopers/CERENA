%getExpectationFSP.m
%Calculates the Expectation of a FSP-model
%
%Usage:
%[y]=getExpectationFSP(P,index)
%
%Input:
%P ... Probability matrix from simulationFSP.
%index ... Index-Mapping from simulateFSP.
%species ... Cell array containing names of the species regarded
%
%Output:
%y ... Expectation values.

function [y]=getExpectationFSP(P,index)
%% 1. Check inputs
if(nargin<2)
    error('This routine requires two inputs')
end
if max(size(P)) ~= size(index,1)
    error('Dimensions of P and index do not agree!');
end
    for j=1:length(index(1,:))
        [P_red,index_red]=getMarginalization(P,index,[j,j]);
        y(:,j)=P_red'*index_red;
    end
end