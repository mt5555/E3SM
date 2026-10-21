#ifndef EAMXX_VIRTUAL_TEMPERATURE_HPP
#define EAMXX_VIRTUAL_TEMPERATURE_HPP

#include "share/diagnostics/abstract_diagnostic.hpp"

namespace scream
{

/*
 * This diagnostic will produce the virtual temperature, T*R*/Rd, with R* computed
 * from the dry mass (pseudo_density_dry) as in the dycore:
 *   R*·pseudo_density = Rd·pseudo_density_dry + Rv·qv·pseudo_density
 * With condensate water loading, pseudo_density_dry also excludes the condensates.
 */

class VirtualTemperature : public AbstractDiagnostic
{
public:
  // Constructors
  VirtualTemperature (const ekat::Comm& comm, const ekat::ParameterList& params,
                      const std::shared_ptr<const AbstractGrid>& grid);

  // The name of the diagnostic CLASS (not the computed field)
  std::string name () const { return "VirtualTemperature"; }

protected:
#ifdef KOKKOS_ENABLE_CUDA
public:
#endif
  void compute_impl ();
};

} //namespace scream

#endif // EAMXX_VIRTUAL_TEMPERATURE_HPP
