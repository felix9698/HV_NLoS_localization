function A = matA_caliter(psi, phi, Qprev, wprev, tdoa, d1, vprev, c)

    P = length(vprev);
    A = [];

    for indp = 2:P
        asin = (-vprev(1)+d1).*triop(3,psi(1)+Qprev,phi(1)+wprev)+(vprev(indp)-d1-c.*tdoa(indp)).*triop(3,psi(indp)+Qprev,phi(indp)+wprev);
        A = [A; [asin,0]];        
    end

    for indp = 2:P
        acos = (-vprev(1)+d1).*triop(4,psi(1)+Qprev,phi(1)+wprev)+(vprev(indp)-d1-c.*tdoa(indp)).*triop(4,psi(indp)+Qprev,phi(indp)+wprev);
        A = [A; [acos,0]];
    end    

    for indp = 2:P
        aelev = (vprev(1)-d1).*sin(psi(1)+Qprev)+(-vprev(indp)+d1+c.*tdoa(indp)).*sin(psi(indp)+Qprev);
        A = [A; [aelev,0]];
    end

    for indp = 2:P
        asin = (vprev(1)-d1).*triop(2,psi(1)+Qprev,phi(1)+wprev)+(-vprev(indp)+d1+c.*tdoa(indp)).*triop(2,psi(indp)+Qprev,phi(indp)+wprev);
        A = [A; [0,asin]];
    end

    for indp = 2:P
        acos = (-vprev(1)+d1).*triop(1,psi(1)+Qprev,phi(1)+wprev)+(vprev(indp)-d1-c.*tdoa(indp)).*triop(1,psi(indp)+Qprev,phi(indp)+wprev);
        A = [A; [0,acos]];
    end
    
end
