function B = matB_caliter(alpha, theta, psi, phi, Qprev, wprev, tdoa, d1, vprev, c)

    P = length(vprev);
    B = [];
    
    for indp = 2:P
        bsin = (vprev(1)-d1).*triop(1,psi(1)+Qprev,phi(1)+wprev)+(-vprev(indp)+d1+c.*tdoa(indp)).*triop(1,psi(indp)+Qprev,phi(indp)+wprev)+vprev(1).*triop(1,alpha(1),theta(1))-vprev(indp).*triop(1,alpha(indp),theta(indp));
        B = [B; bsin];
    end

    for indp = 2:P
        bcos = (vprev(1)-d1).*triop(2,psi(1)+Qprev,phi(1)+wprev)+(-vprev(indp)+d1+c.*tdoa(indp)).*triop(2,psi(indp)+Qprev,phi(indp)+wprev)+vprev(1).*triop(2,alpha(1),theta(1))-vprev(indp).*triop(2,alpha(indp),theta(indp));
        B = [B; bcos];
    end    

    for indp = 2:P
        belev = (vprev(1)-d1).*cos(psi(1)+Qprev)+(-vprev(indp)+d1+c.*tdoa(indp)).*cos(psi(indp)+Qprev)+vprev(1).*cos(alpha(1))-vprev(indp).*cos(alpha(indp));
        B = [B; belev];
    end

    for indp = 2:P
        bsin = (vprev(1)-d1).*triop(1,psi(1)+Qprev,phi(1)+wprev)+(-vprev(indp)+d1+c.*tdoa(indp)).*triop(1,psi(indp)+Qprev,phi(indp)+wprev)+vprev(1).*triop(1,alpha(1),theta(1))-vprev(indp).*triop(1,alpha(indp),theta(indp));
        B = [B; bsin];
    end

    for indp = 2:P
        bcos = (vprev(1)-d1).*triop(2,psi(1)+Qprev,phi(1)+wprev)+(-vprev(indp)+d1+c.*tdoa(indp)).*triop(2,psi(indp)+Qprev,phi(indp)+wprev)+vprev(1).*triop(2,alpha(1),theta(1))-vprev(indp).*triop(2,alpha(indp),theta(indp));
        B = [B; bcos];
    end
    
end
