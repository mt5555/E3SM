#include "catch2/catch.hpp"

#include "share/grid/point_grid.hpp"
#include "share/diagnostics/register_diagnostics.hpp"

#include "share/physics/physics_constants.hpp"

#include "share/core/eamxx_setup_random_test.hpp"
#include "share/physics/eamxx_common_physics_functions.hpp"
#include "share/field/field_utils.hpp"

#include <ekat_pack.hpp>
#include <ekat_team_policy_utils.hpp>
#include <ekat_view_utils.hpp>

#include <iomanip>

namespace scream {


//-----------------------------------------------------------------------------------------------//
template<typename DeviceT>
void run(std::mt19937_64& engine)
{
  using PF         = scream::PhysicsFunctions<DeviceT>;
  using PC         = scream::physics::Constants<Real>;
  using KT         = ekat::KokkosTypes<DeviceT>;
  using ExecSpace  = typename KT::ExeSpace;
  using MemberType = typename KT::MemberType;
  using TPF        = ekat::TeamPolicyFactory<ExecSpace>;

  const int packsize = SCREAM_PACK_SIZE;

  // A world comm
  ekat::Comm comm(MPI_COMM_WORLD);

  // Create a grids manager - single column for these tests
  const int ncols = 1;
  const int nlevs = packsize*2 + 1; // Number of levels to use for tests, make sure the last pack can also have some empty slots (packsize>1).
  auto grid = create_point_grid("physics",ncols,nlevs,comm);

  // Kokkos Policy
  auto policy = TPF::get_default_team_policy(ncols, nlevs);

  // Construct random input data
  using RPDF = std::uniform_real_distribution<Real>;
  RPDF pdf_qv(1e-6,1e-3),
       pdf_temp(200.0,400.0),
       pdf_dp(1.0,100.0);

  // A time stamp
  util::TimeStamp t0 ({2022,1,1},{0,0,0});

  // Construct the Diagnostic
  ekat::ParameterList params;
  register_diagnostics();
  auto& diag_factory = DiagnosticFactory::instance();
  auto diag = diag_factory.create("VirtualTemperature",comm,params, grid);
  // Set the required fields for the diagnostic.
  std::map<std::string,Field> input_fields;
  {
    using namespace ShortFieldTagsNames;
    using namespace ekat::units;
    auto scalar3d = grid->get_3d_scalar_layout(LEV);
    for (const auto& fname : diag->get_input_fields_names()) {
      FieldIdentifier fid(fname, scalar3d, Pa, grid->name());
      Field f(fid);
      f.allocate_view();
      f.get_header().get_tracking().update_time_stamp(t0);
      diag->set_input_field(f);
      input_fields.emplace(fname, f);
    }
  }

  // Initialize the diagnostic
  diag->initialize();

  // Run tests, without and with condensate in the dry mass
  for (const bool condensate : {false, true}) {
    // Construct random data to use for test
    // Get views of input data and set to random values
    const auto& T_mid_f = input_fields["T_mid"];
    const auto& T_mid_v = T_mid_f.get_view<Real**>();
    const auto& qv_mid_f = input_fields["qv"];
    const auto& qv_mid_v = qv_mid_f.get_view<Real**>();
    const auto& dp_v     = input_fields["pseudo_density"].get_view<Real**>();
    const auto& dp_dry_v = input_fields["pseudo_density_dry"].get_view<Real**>();
    for (int icol=0;icol<ncols;icol++) {
      const auto& T_sub = ekat::subview(T_mid_v,icol);
      const auto& qv_sub = ekat::subview(qv_mid_v,icol);
      const auto& dp_sub = ekat::subview(dp_v,icol);
      ekat::genRandArray(T_sub,  engine, pdf_temp);
      ekat::genRandArray(qv_sub, engine, pdf_qv);
      ekat::genRandArray(dp_sub, engine, pdf_dp);
    }
    // dp_dry = dp*(1-qv-qc), with qc=0 if there is no condensate
    {
      const auto qc = condensate ? Real(pdf_qv(engine)) : Real(0);
      Kokkos::parallel_for("", policy, KOKKOS_LAMBDA(const MemberType& team) {
        const int icol = team.league_rank();
        Kokkos::parallel_for(Kokkos::TeamVectorRange(team,nlevs), [&] (const Int& ilev) {
          dp_dry_v(icol,ilev) = dp_v(icol,ilev)*(1 - qv_mid_v(icol,ilev) - qc);
        });
      });
      Kokkos::fence();
    }

    // Run diagnostic and compare with manual calculation
    diag->compute(t0);
    const auto& diag_out = diag->get();
    Field virtualT_f = diag_out.clone();
    virtualT_f.deep_copy(0);
    const auto& virtualT_v = virtualT_f.get_view<Real**>();
    // Reference written out, Tv = T*R*/Rd = T*(q_dry + qv*Rv/Rd), rather than
    // calling the PhysicsFunctions helper the diagnostic uses
    const Real rv_rd = 1/PC::ep_2.value;
    Kokkos::parallel_for("", policy, KOKKOS_LAMBDA(const MemberType& team) {
      const int icol = team.league_rank();
      Kokkos::parallel_for(Kokkos::TeamVectorRange(team,nlevs), [&] (const Int& ilev) {
        const Real q_dry = dp_dry_v(icol,ilev)/dp_v(icol,ilev);
        virtualT_v(icol,ilev) = T_mid_v(icol,ilev)*(q_dry + rv_rd*qv_mid_v(icol,ilev));
      });
      team.team_barrier();
    });
    Kokkos::fence();
    REQUIRE(views_are_equal(diag_out,virtualT_f));

    // With condensate, the result differs from the qv-only virtual temperature
    if (condensate) {
      Field virtualT_qv_f = diag_out.clone();
      const auto& virtualT_qv_v = virtualT_qv_f.get_view<Real**>();
      Kokkos::parallel_for("", policy, KOKKOS_LAMBDA(const MemberType& team) {
        const int icol = team.league_rank();
        Kokkos::parallel_for(Kokkos::TeamVectorRange(team,nlevs), [&] (const Int& ilev) {
          virtualT_qv_v(icol,ilev) = PF::calculate_virtual_temperature(T_mid_v(icol,ilev),qv_mid_v(icol,ilev));
        });
      });
      Kokkos::fence();
      REQUIRE(not views_are_equal(diag_out,virtualT_qv_f));
    }
  }
 
} // run()

TEST_CASE("virtual_temperature_test", "virtual_temperature_test]"){
  // Run tests for both Real and Pack, and for (potentially) different pack sizes
  using scream::Real;
  using Device = scream::DefaultDevice;

  constexpr int num_runs = 5;

  auto engine = scream::setup_random_test();

  printf(" -> Number of randomized runs: %d\n\n", num_runs);

  printf(" -> Testing Pack<Real,%d> scalar type...",SCREAM_PACK_SIZE);
  for (int irun=0; irun<num_runs; ++irun) {
    run<Device>(engine);
  }
  printf("ok!\n");

  printf("\n");

} // TEST_CASE

} // namespace
