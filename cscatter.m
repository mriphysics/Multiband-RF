function varargout = cscatter(r,cl,lw,ls)
% Function that scatter plots on the complex plane.
% Useful for work in DSP and SLR beta-polynomials
% e.g. to plot for a filter: cscatter(roots(b))
% 23/7/2015 Samy Abo Seada

% 13/01/2016 - allow color input

% 03/03/2016 - allow linewidth and line-spec of the unit circle

nout = nargout;

if nargin == 1
    cl = 'b'; %By default use blue as color.
    lw = 2;
    ls = '-.r';
end
if nargin == 2 
    lw = 2;
    ls = '-.r';
end
if nargin == 3
    ls = '-.r';
end

% scatter(real(r),imag(r),'color',cl,'linewidth',lw);
sh = scatter(real(r),imag(r),'MarkerEdgeColor',cl,'linewidth',lw);
hold on;
grid on;

unit_circle=exp(1i*(0:0.001:2*pi));
plot(real(unit_circle),imag(unit_circle),ls)

% If output requested, return scatter handle
if nout == 1
    varargout{1} = sh;
end
hold off;
end