# Description to explore data + compare to Kara's
loadd(all_plume_speciation_df, cache = drake::new_cache("scripts/.drake"))

# kara's data
test_filt <- all_plume_speciation_df
             
# http://wiki.seas.harvard.edu/geos-chem/index.php/Particulate_matter_in_GEOS-Chem#PM2.5_definition
# to calculate PM2.5:
# SO4 (sulfate) + NIT (nitrate) + NH4 (ammonium) + BC (black carbon) + OC (organic carbon) + Dust + Sea salt + SOA (secondary organic aerosol)
# mapping to IMPROVE variable names:
# SO4f for SO4
# NO3f for NIT
# ammNO3f + ammSO4f for the sum for NH4. 
# But then we only want the associated mass from NH4, so need to drop the mass from SO4 and NO3 (otw double counted). 
# If my math is right, then it should be 18/(18+62)*ammNO3f + 18/(18+96)*ammSO4f
# ECf for BC
# OMCf for OC+SOA
# We don’t care about seasalt here
# Dust is tricky. We would previously calculate it as proportional to iron, but this might overlap with metal species. 
  # There is a column called SOILf, but I don’t see any dust paper using it.  
  # iron /X (X being the percentage of iron in typical dust particles, e.g. 3.5%)
# http://vista.cira.colostate.edu/Improve/reconstructed-fine-mass/
# USE SOIL COLUMN AS DUST -> FINAL ANSWER
  # breakdown of what goes into the soil calculation: Soil (88348) = 2.2*Aluminum + 2.49*Silicon + 1.63*Calcium + 2.43* Iron + 1.94* Titanium
# Metal is not the main focus of most PM2.5 CTM work. 
  # I would vote to have a category as the sum of all metal,
  # and then track some individual metal components of interest (from the health perspective).


# differentiate our analysis into two parts: metal vs non-metal. 
# For the non-metal part, we could take advantage of the reconstructed formula from IMPROVE 
# [RCFM] = ammonium sulfate + ammonium nitrate + organic mass + elemental carbon + fine soil + sea salt
# Therefore, we would (by design) have a clean and mutually exclusive grouping. 
# This could be used for the analysis marshall laid out above: 
# pool monitors within a region and then run monitor-day regressions of species conc on smoke PM. 
# To make this work, smoke PM would need to use the reconstructed fine matter (RCFM). 

# Metal part would just be tracking the metal components of interest (using either the regression/non-regression method). 
# The metal analysis might need to be separate anyway since their conc is much lower than the other big groups 
