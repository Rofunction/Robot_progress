%{
    tp is customized, eg. noise_t=noise(t,1,0.1,10);
%}
function noise_t=noise(t,tp,peak,fat)
noise_t=peak*exp(-fat*(t-tp).^2);
end