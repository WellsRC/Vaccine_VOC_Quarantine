function [pA,epsI]=VariantParameters(VOC)

pAv=[0.351 0.817 0.9637;
    0.351 0.8105 0.9624;
    0.2204 0.7092 0.9493;
    0.275 0.3388  0.7629];
epsIv=[0.74 0.9549;
       0.737 0.9544;
       0.637 0.937;
       0.138 0.705];
pA=pAv(VOC+1,:);
epsI=epsIv(VOC+1,:);
end