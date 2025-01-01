function Lr = getLrPair(B1,B2,K)
%Returns Lr matrix for a pair of circles. 
   % Lr = blkdiag(B1*K',B2*K'); %Not the right thing! 

    Lr1 = B1*K'; % xx xy; yx yy
    Lr2 = B2*K';

    xx1 = Lr1(1:end/2,1:end/2);
    xx2 = Lr2(1:end/2,1:end/2);
    yy1 = Lr1(end/2+1:end,end/2+1:end);
    yy2 = Lr2(end/2+1:end,end/2+1:end);
    yx1 = Lr1(end/2+1:end,1:end/2);
    yx2 = Lr2(end/2+1:end,1:end/2);
    xy1 = Lr1(1:end/2,end/2+1:end);
    xy2 = Lr2(1:end/2,end/2+1:end);
   
    Lrxx = [xx1 zeros(size(xx1)); zeros(size(xx1)) xx2];
    Lryy = [yy1 zeros(size(xx1)); zeros(size(yy1)) yy2];
    Lrxy = [xy1 zeros(size(xx1)); zeros(size(xx1)) xy2];
    Lryx = [yx1 zeros(size(xx1)); zeros(size(yy1)) yx2];

    Lr = [Lrxx Lrxy; Lryx Lryy];
    %Is this consistent with what I do in 3D, where the
    %ordering of the unknowns is different?

   % M1 = length([q(i)+rout_base; fine_1]);
   % M2 = length([q(p2)+rout_base; fine_2]);
    %Lr = blkdiag(B1*K'/M1,B2*K'/M2); Does it help to rescale?
    
end