#include "virtual_temperature.hpp"
#include "share/physics/eamxx_common_physics_functions.hpp"

namespace scream
{

VirtualTemperature::
VirtualTemperature (const ekat::Comm& comm, const ekat::ParameterList& params,
                    const std::shared_ptr<const AbstractGrid>& grid)
 : AbstractDiagnostic(comm,params,grid)
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  const auto& grid_name = m_grid->name();

  auto scalar3d = m_grid->get_3d_scalar_layout(LEV);

  m_field_in_names.push_back("T_mid");
  m_field_in_names.push_back("qv");
  m_field_in_names.push_back("pseudo_density");
  m_field_in_names.push_back("pseudo_density_dry");

  FieldIdentifier fid (name(), scalar3d, K, grid_name);
  m_diagnostic_output = Field(fid);
  m_diagnostic_output.allocate_view();
}

void VirtualTemperature::compute_impl()
{
  using KT      = KokkosTypes<DefaultDevice>;
  using PF      = PhysicsFunctions<DefaultDevice>;
  using MDRange = Kokkos::MDRangePolicy<typename KT::ExeSpace,Kokkos::Rank<2>>;

  const auto& virtualT = m_diagnostic_output.get_view<Real**>();
  const auto& T_mid    = m_fields_in.at("T_mid").get_view<const Real**>();
  const auto& qv_mid   = m_fields_in.at("qv").get_view<const Real**>();
  const auto& dp       = m_fields_in.at("pseudo_density").get_view<const Real**>();
  const auto& dp_dry   = m_fields_in.at("pseudo_density_dry").get_view<const Real**>();

  const int ncols = m_grid->get_num_local_dofs();
  const int nlevs = m_grid->get_num_vertical_levels();

  auto lambda = KOKKOS_LAMBDA(const int icol, const int ilev) {
    // R* from the dry mass, as in the dycore. With water loading the dry mass
    // also excludes the condensates.
    virtualT(icol,ilev) = PF::calculate_virtual_temperature(T_mid(icol,ilev),qv_mid(icol,ilev),
                                                            dp(icol,ilev),dp_dry(icol,ilev));
  };
  MDRange policy ({0,0},{ncols,nlevs});
  Kokkos::parallel_for("VirtualTemperature",policy,lambda);
}

} //namespace scream
