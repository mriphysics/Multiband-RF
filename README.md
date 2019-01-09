# Multiband RF toolbox
Design the following types of RF pulses:

Singleband pulses with the options:
* Linear-phase, minimum-phase, quadratic-phase
* Excitation, refocusing and matched-excitation refocusing
* Using Least-squares or CVX optimization

Multiband pulses:
* Phase-optimized (Wong ISMRM 2012)
* Time-shifted (Auerbach MRM 2013)
* Root-flipped (Sharma MRM 2015)
* VERSE-Multiband (Our work in Abo Seada MRM 2018)
* VERSE-Multiband with GIRF correction (Abo Seada ISMRM 2017)
* PINS (Norris MRM 2011)
* Multi-PINS (Eicher MRM 2014)


This code relies on the following repos:
* AM_multiband  (https://github.com/mriphysics/AM_multiband)
* verse-mb (https://github.com/mriphysics/verse-mb)
Both have already been incorporated in this repo, for convenience. At the time of writing they are exact copies.

Also relies on
*  CVX from (http://cvxr.com/cvx/). 
*  Pauly's RF tools from (http://rsl.stanford.edu/research/software.html).


## Multiband examples

### A usual multiband refocusing pulse (all in-phase, from precalculated RF).

<img src="Figures/MBexample.png" alt="MBexample" width="100%">

### Quadratic phase VERSE multiband pulse with GIRF correction


### PINS refocusing pulse

## Singleband examples

### A usual singleband refocusing pulse (linear-phase, designed with CVX)
The most common parameters you would want to change are near the top.
```
Nt = 256; %<-- nr of time-points
tb = 6;   %<-- Time bandwidth product

d1 = 0.01; %<-- passband ripple [%]
d2 = 0.01; %<-- stopband ripple [%]

% Set Flip angle, design mode, pulse type and phase type.
flip = pi; % Flip-angle in radians
mode = 'cvx'; % Set to cvx or ls
type = 'ref'; % Set to exc, ref or me (excitation, refocusing or matched-excitation refocusing respectively).
phase = 'linear'; % Set to linear, minimum, maximum or quadratic.
plot_fa = 0; %<-- set to 1 to plot profile in flip-angle representation
quiet = 1;  %<-- set to 1 to reduce command window output
slthick = 2*1e-3; %<-- Slice-thickness in m
```

These go into the singleband_rf function, and the output should look something like this
<img src="Figures/SBexample.png" alt="SBexample" width="100%">

### A more exotic, quadratic phase, matched-excitation refocusing pulse
To do somehting like this, set pulse `type` to 'me' (short for matched-excitation) and set `phase` to 'quadratic'.

```python
type = 'me'; % Set to exc, ref or me (excitation, refocusing or matched-excitation refocusing respectively).
phase = 'quadratic'; % Set to linear, minimum, maximum or quadratic.
```

When designing matched-excitation refocusing pulse, and a third output is requested from the `singleband_rf` function,
an matching excitation pulse is designed.
```
[rf,tb,rf_me_exc] = singleband_rf(Nt,tb,flip,mode,type,phase,d1,d2,quiet);
```

The easiest way to check the slice-profile of a matched-excitation pulse is the following CK representation
``mxy_display = @(a180,b180)(sqrt(2)*abs(b180).^4.*sqrt(1-abs(b180).^4/2))``

However, in the Bloch-simulation section there is commented where you can then simulate both excitation pulse and refocusing pulse separately.
The result should be the same (as far as I checked..).

```
 dt90 = max(abs(rf_me_exc))/(gamma_mT*b1max);
 T90 = (length(rf_me_exc))*dt90;        

 BW90 = 2*tb/T90; %<-- double BW corrects sb ripples. Why??
 Gsel90 = 2*pi*BW90/(gamma_mT*slthick)*ones(length(rf_me_exc),1);
 G90 =[0*Gsel90 0*Gsel90 Gsel90];
 % Simulate excitation pulse
 [~,~,~,~,aexc,bexc] = blochsim_CK(rf_me_exc(:)/gamma_mT/dt90,G90,pos,ones([Nz 1]),zeros([Nz 1]),'dt',dt90);
 % Simulate refocusing pulse
 [~,~,~,~,aref,bref] = blochsim_CK(rf(:)/gamma_mT/dt,G,pos,ones([Nz 1]),zeros([Nz 1]),'dt',dt);
 
 % This is the product of an CK excitation profile 2*alpha*beta and refocusing profile beta^2.
 mxy = 2*conj(aexc(:,end)).*bexc(:,end).*bref(:,end).^2;
```
The result should look like this:
<img src="/Figures/me_quad_SBexample.PNG" alt="me_quad_SBexample" width="100%">