# Multiband RF pulse design for realistic gradient performance

## Content
This repository contains code to produce multiband VERSE (MBv) and VERSE multiband RF and gradient pulses
, associated with our published work [Abo Seada et al. MRM 2018](https://onlinelibrary.wiley.com/doi/abs/10.1002/mrm.27411) 
titled "**Multiband RF pulse design for realistic gradient performance**" (DOI: 10.1002/mrm.27411) available 
 at https://onlinelibrary.wiley.com/doi/abs/10.1002/mrm.27411 

The VERSE code used is an modified version version of https://github.com/mriphysics/reVERSE-GIRF by Shaihan Malik, which in turn
is an implementation of Time-optimal VERSE pulse design [(Lee et al. MRM 2009)](http://doi.org/10.1002/mrm.21950). This has been implemented by modifying code released by Miki Lustig for designing time-optimal gradients, publically available on the  [authors' website](http://www.eecs.berkeley.edu/~mlustig/Software.html). 

Acknowledgements to the time-optimal gradient framework should be attributed to [Lustig IEEE-TMI 2008](https://www.ncbi.nlm.nih.gov/pubmed/18541493) whilst acknowledgement to the Time-optimal VERSE should be attributed to Lee et al 2009.

## Purpose
Multiband RF pulses are an essential building block in multiband/simultaneous multi-slice imaging sequences. However conventional RF 
pulse design strategies lead to long RF pulse durations, which leads to long echo-times and thus lower SNR. One way of significantly reducing
pulse durations is to use time-variable selection gradients, like those used in the ISMRM pulse design challenge in 2016 [(Grissom et al. MRM 2016)](https://onlinelibrary.wiley.com/doi/abs/10.1002/mrm.26512). 

However, time-optimal design methods such as VERSE assume RF and gradient systems have similar output bandwidth, however in practice gradient systems have much lower temporal bandwidth than RF systems.
This is down to a number of factors such as eddy currents, amplifier bandwidth and mechanical vibrations.
The output bandwidth of a gradient system can be well characterized by a Gradient Impulse Response Function (GIRF) as demonstrated by [Vannesjo et al. MRM 2013](https://onlinelibrary.wiley.com/doi/abs/10.1002/mrm.24263).

When ignored, this can lead to slice distortions and image artefacts.
This is well-appreciated by looking at a Short-Time Fourier Transform of a VERSE gradient, when VERSE is
applied to a multiband RF pulse.
As can be seen in the top-left example, such a gradient pulse with fast-oscillating ripples lead to high-frequency 
components, which are visible in the STFT in the bottom-left. This demanded gradient waveform can be beyond the
FWHM of a typical GIRF (shown in blue dashed lines).
<img src="bin/Github_Figures/fig1.png" alt="fig1" width="100%">

In the right column, we show an equivalent result when we apply VERSE to a single-band RF waveform, which does not contain such high-frequency components
and is thus "kinder" to the gradient system, and more likely to be reproduced with fidelity.

We investigated the slice-profile artefacts associated with two methods of combining multiband RF pulse design and VERSE.
* The first is by directly applying VERSE to a multiband RF pulse, and is referred to as **MBv**.
* The second method is applying VERSE to a singleband RF pulse, which is
subsequently modulated into a multiband pulse, and is referred to as **vMB**.

The latter method will retain the smooth gradient waveform shown in the right-column.

For these two approaches, we investigated the slice-profile distortions one might expect based on two different GIRF models, representative of
two different MRI systems.

## Results
Temporal profiles of RF and Gradient pulses are shown below, as well as the simulated gradient distortion and
the slice-profile effects due to the GIRF.

<img src="bin/Github_Figures/fig2.png" alt="fig2" width="100%">

For the **MBv** method (Multiband pulse followed by VERSE) the errors in reproduced gradient are shown to excite ghost slices outside the desired imaging field-of-view.
This does not happen for the **vMB** method, which is a clear benefit of this approach.

<img src="bin/Github_Figures/fig3.png" alt="fig3" width="100%">

## Contact

Thank you for taking an interest in our work, and if you have any questions or suggestion please reach out to 
contact us.

**Shaihan Malik** - Email:
shaihan.malik@kcl.ac.uk

**Samy Abo Seada** - Email:
samy.abo_seada@kcl.ac.uk

King's College London, 2018.

Initial upload: 02/02/2018

Update after publication: 08/10/2018
