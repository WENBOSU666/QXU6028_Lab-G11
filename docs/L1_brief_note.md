L1 Brief Note: Carrier Counting and Mobility Correction

We compared two carrier-counting definitions using the same linear density of states, g(E) = cE, and finite-temperature Fermi–Dirac distribution at T=300 K. The band-edge count integrates from E=0, whereas the Fermi-referenced count integrates from E=EF ; both use an upper limit of 2.5 eV. Numerical integration used 20,001 points and the trapezoidal rule. The supplied band-edge counts were used for the sweep data, and both integrals were evaluated at the exact operating point specified in the G11 parameter file.
At EF = 0.1555011006 eV and sheet resistivity ρ= 285.8450549Ω/m , the results were:
Quantity	Result
Band-edge carrier density, nband	1.82345× 1012cm−2
Fermi-referenced carrier density, nF	4.61273×1011cm−2
Correction factor, C=nband/nF	3.95308
Band-count mobility	1.19747×104cm2V−1S−1
Corrected mobility	4.73371×104cm2V−1S−1

The correction factor increases monotonically across the supplied positive EF sweep, reaching C=2 at 0.072762 eV and C=4 at 0.157347 eV. With increasing positive doping, the band-edge integral includes an expanding occupied interval below EF , which the Fermi-referenced definition excludes. Consequently, the difference between the two counts grows.
At fixed resistivity, μ=1/(enρ), so the corrected mobility is larger than the band-count estimate by the same factor, 3.95308. This is a correction to the extracted mobility, not evidence of a change in the material’s scattering mechanism. Under the single-carrier current-related Hall sensitivity model, SI=1/(en), the corrected sensitivity is likewise 3.95308 times the band-count prediction—an increase of approximately 295.3%. Sensor calibration should therefore use the corrected count within this model; this correction alone does not establish improved noise performance.
Validation gave C(0)=1.00000, C(0.10eV)=2.59885, and C(0.20eV)=5.10437. Numerical and closed-form Fermi counts agreed to within 1.2 × 10-6 relative difference. Doubling the integration grid changed the operating-point counts by less than 2.6×10-7 relatively. These establish numerical convergence, not experimental uncertainty; uncertainties in the supplied parameters were not provided.