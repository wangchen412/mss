function check_ev(ev)
    for i = 1 : size(ev)
        x = [abs(ev(i)), log(ev(i))/1i/pi/2];
        if x(1) > 0.8 && x(1) < 1.2
            disp(x)
        end
    end
end
