function check_ev(ev, tol)
    for i = 1 : size(ev)
        x = [abs(ev(i)), log(ev(i))/1i/pi];
        if x(1) > 1-tol && x(1) < 1 + tol
            disp(x)
        end
    end
end
