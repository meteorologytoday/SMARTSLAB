α   = 3e-4     # K^-1    http://www.kayelaby.npl.co.uk/general_physics/2_7/2_7_9.html
β   = 1e-3     # Simple estimation
c_p = 3985.0   # J / kg / K
ρ   = 1027.0   # kg / m^3
g   = 9.8      # m / s^2
h_ML_min = 30.0
h_ML_max = 1000.0

"""
    printConstants()

# Description
This function prints all constants used in this module.
"""
function printConstants()
    @printf("α   = %8.2f. Logrithmic expansion of density ρ as function of temperature T.\n", α)
    @printf("β   = %8.2f. Logrithmic expansion of density ρ as function of salinity S.\n", β)
    @printf("c_p = %8.2f J/kg/K. Specific heat of seawater.\n", c_p)
    @printf("ρ   = %8.2f kg/m^3. Mass density of seawater.\n", ρ)

end


