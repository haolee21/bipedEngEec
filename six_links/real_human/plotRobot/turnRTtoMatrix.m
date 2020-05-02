function transMat = turnRTtoMatrix (rtResult)
transMat = [rtResult.n,rtResult.o,rtResult.a,rtResult.t];
transMat = [transMat;0 0 0 1];
end