function [pA,epsI]=VariantParameters(VOC)

pAv=[0.351 0.9611 0.9789;
    0.351 0.9546 0.9754;
    0.2204 0.9064 0.9493;
    0.275 0.275  0.7629];
epsIv=[0.92 0.98;
       0.917 0.9792;
       0.76 0.94;
       0.38 0.82];
pA=pAv(VOC+1,:);
epsI=epsIv(VOC+1,:);
end