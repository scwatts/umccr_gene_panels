# Important

## hmftools panel data

* genes with ambiguous (oncogene and tsgene) or absent roles handled by:
  * deferring to Hartwig configuration where available
  * fallback to tsgene in all other cases
* impact of this is minimal since gene role is only used in driver likelihood calculation
  * not used by curators
* best solution would be to have hmftools be able to support dual role
* otherwise consult/resolve with curators if driver likelihood calculation is to be used
