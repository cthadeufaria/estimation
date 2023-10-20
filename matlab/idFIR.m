

function [Y, Phi, sys] = idFIR(data, nb)
    
    Y = data.y((nb + 1):end);
    Phi = zeros(length(data.u) - nb, nb + 1);

    for i = 1 : length(data.u)
        Phi(i, :) = u(i + nb : -1 : i)';
    end

    Theta = Phi \ Y;
    % Theta = (transpose(Phi) * Phi) \ (transpose(Phi) * Y);

    sys = idpoly(1,Theta');

end