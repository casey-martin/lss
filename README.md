# Lactococcus Normalization
We found substantial Lactococcus contamination that we suspect is coming from the purified diets.
Because there is a large increase in Lactococcus relative abundance during the time of antibiotic administration,
we suspect that this signal is inert (genetic material) rather than biologically active and growing.

In addition to contaminant filtering and total sum scaling, we are also attempting to normalize relative abundance
data to the Lactococcus relative abundance under the assumption that the rate of contamination is constant throughout
the timecourse. In doing so, we hoped create a pseudo measure of absolute abundance.

To run the bootleg core diversity analysis, execute the following:
`sh core_diversity_analysis.sh`
