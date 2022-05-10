function [pA,epsI]=VariantParameters(VOC)

pAv=[0.351 0.817 0.9681;
    0.351 0.8105 0.9670;
    0.2204 0.7092 0.9493;
    0.275 0.3388  0.7629];
epsIv=[0.74 0.9585;
       0.737 0.9580;
       0.637 0.942;
       0.138 0.721];
pA=pAv(VOC+1,:);
epsI=epsIv(VOC+1,:);
end