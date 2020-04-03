function ph = hu_invariants(A)

r = 1:size(A,1);
c = 1:size(A,2);
m10 = sum(sum(r.'.*A));
m01 = sum(sum(c.*A));
m00 = sum(sum(A));

xmean = m10/m00;
ymean = m01/m00;
mu00 = m00;

n11 = mu(1,1,xmean,ymean,A)/(mu00^2);
n12 = mu(1,2,xmean,ymean,A)/(mu00^2.5);
n21 = mu(2,1,xmean,ymean,A)/(mu00^2.5);
n20 = mu(2,0,xmean,ymean,A)/(mu00^2);
n02 = mu(0,2,xmean,ymean,A)/(mu00^2);
n30 = mu(3,0,xmean,ymean,A)/(mu00^2.5);
n03 = mu(0,3,xmean,ymean,A)/(mu00^2.5);

ph(1) = n20 + n02;
ph(2) = (n20 - n02)^2 + 4*(n11^2);
ph(3) = (n30 - 3*n12)^2 + (3*n21 - n03)^2;
ph(4) = (n30 + n12)^2 + (n21 + n03)^2;
ph(5) = (n30 - 3*n12)*(n30 + n12)*((n30 + n12)^2 - 3*(n21 + n03)^2) + ...
        (3*n21 - n03)*(n21 + n03)*(3*(n30 + n12)^2 - (n21 + n03)^2);
ph(6) = (n20 - n02)*((n30 + n12)^2 - (n21 + n03)^2) + ...
        4*n11*(n30 + n12)*(n21 + n03);
ph(7) = (3*n21 - n03)*(n30 + n12)*((n30 + n12)^2 - 3*(n21 + n03)^2) + ...
        (3*n12 - n03)*(n21 + n03)*(3*(n30 + n12)^2 - (n21 + n03)^2);

ph = sign(ph).*log(abs(ph));
end

function mupq = mu(p,q,xmean,ymean,A)
    r = 1:size(A,1);
    c = 1:size(A,2);
    f1 = (r - xmean).^p;
    f2 = (c - ymean).^q;
    aux = f1.'.*A;
    aux = f2.*aux;
    mupq = sum(sum(aux));
end

