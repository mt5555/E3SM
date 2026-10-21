#include <catch2/catch.hpp>

#include "ElementOps.hpp"

#include "utilities/TestUtils.hpp"
#include "utilities/SubviewUtils.hpp"
#include "utilities/SyncUtils.hpp"
#include "utilities/ViewUtils.hpp"
#include "Types.hpp"

#include <iostream>
#include <random>

using namespace Homme;

extern "C" {
void init_f90 (const Real* hyai_ptr, const Real& ps0);
void compute_r_star_f90(const int& num_elems,
                        const bool& moist,
                        const Real*& Q,
                        Real*& R_star);
void compute_r_star_wl_f90(const int& num_elems,
                           const bool& moist,
                           const bool& water_loading,
                           const Real*& q1, const Real*& q2, const Real*& q3,
                           Real*& R_star);
} // extern "C"

// ============= ELEMENT OPS ================ //

TEST_CASE("elem_ops", "elem_ops") {

  constexpr int num_elems = 2;
  auto policy = Homme::get_default_team_policy<ExecSpace>(num_elems);

  std::random_device rd;
  const unsigned int catchRngSeed = Catch::rngSeed();
  const unsigned int seed = catchRngSeed==0 ? rd() : catchRngSeed;
  std::cout << "seed: " << seed << (catchRngSeed==0 ? " (catch rng seed was 0)\n" : "\n");
  using rngAlg = std::mt19937_64;
  rngAlg engine(seed);
  std::uniform_real_distribution<Real> pdf(0.01, 1.0);

  HybridVCoord hvcoord;
  hvcoord.random_init(seed);

  decltype(hvcoord.hybrid_ai)::host_mirror_type hyai = Kokkos::create_mirror_view(hvcoord.hybrid_ai);
  Kokkos::deep_copy(hyai,hvcoord.hybrid_ai);
  const Real* hyai_ptr = hyai.data();
  init_f90(hyai_ptr,hvcoord.ps0);

  // F90 and CXX views (and host mirrors of cxx views)
  HostViewManaged<Real*[NUM_PHYSICAL_LEV][NP][NP]>  dp_f90("",num_elems);
  HostViewManaged<Real*[NUM_PHYSICAL_LEV][NP][NP]>  theta_ref_f90("",num_elems);

  ExecViewManaged<Scalar*[NP][NP][NUM_LEV]>   dp_cxx("",num_elems);
  ExecViewManaged<Scalar*[NP][NP][NUM_LEV]>   theta_ref_cxx("",num_elems);
  ExecViewManaged<Scalar*[NP][NP][NUM_LEV]>   p_cxx("",num_elems);
  ExecViewManaged<Scalar*[NP][NP][NUM_LEV_P]> p_i_cxx("",num_elems);

  auto h_theta_ref = Kokkos::create_mirror_view(theta_ref_cxx);

  ElementOps elem_ops;
  elem_ops.init(hvcoord);

  SECTION("r_star") {
    for (bool moist : {false, true}) {

      // Use some existing buffers
      auto Q_cxx = theta_ref_cxx;
      auto Q_f90 = theta_ref_f90;
      auto R_f90 = dp_f90;
      genRandArray(Q_cxx,engine,pdf);
      sync_to_host(Q_cxx,Q_f90);

      // Run cxx version
      Kokkos::parallel_for(policy,KOKKOS_LAMBDA(const TeamMember& team) {
        KernelVariables kv(team);
        Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team,NP*NP),
                             [&](const int idx) {
          const int igp = idx / NP;
          const int jgp = idx % NP;

          auto Q = Homme::subview(Q_cxx,kv.ie,igp,jgp);
          auto R = Homme::subview(dp_cxx,kv.ie,igp,jgp);

          elem_ops.get_R_star(kv,moist,Q,R);
        });
      });

      // Run f90 version
      const Real* Q_ptr = Q_f90.data();
      Real* R_ptr = R_f90.data();
      compute_r_star_f90(num_elems,moist,Q_ptr,R_ptr);

      // Compare answers
      auto R_cxx = Kokkos::create_mirror_view(dp_cxx);
      Kokkos::deep_copy(R_cxx,dp_cxx);

      for (int ie=0; ie<num_elems; ++ie) {
        for (int igp=0; igp<NP; ++igp) {
          for (int jgp=0; jgp<NP; ++jgp) {
            for (int k=0; k<NUM_PHYSICAL_LEV; ++k) {
              const int ilev = k / VECTOR_SIZE;
              const int ivec = k % VECTOR_SIZE;
              REQUIRE (R_cxx(ie,igp,jgp,ilev)[ivec] == R_f90(ie,k,igp,jgp));
            }
          }
        }
      }
    }
  }

  SECTION("r_star_water_loading") {
    // Three tracers: vapor (0) and two condensates, listed out of order (2,1)
    std::uniform_real_distribution<Real> q_pdf(0.0, 0.02);
    ExecViewManaged<Scalar*[NP][NP][NUM_LEV]> q0_cxx("",num_elems), q1_cxx("",num_elems), q2_cxx("",num_elems);
    HostViewManaged<Real*[NUM_PHYSICAL_LEV][NP][NP]> q0_f90("",num_elems), q1_f90("",num_elems), q2_f90("",num_elems);
    HostViewManaged<Real*[NUM_PHYSICAL_LEV][NP][NP]> R_f90("",num_elems);
    ExecViewManaged<Scalar*[NP][NP][NUM_LEV]> R_cxx_d("",num_elems);

    WaterLoading wl;
    wl.num_condensates = 2;
    wl.condensate_idx[0] = 2;
    wl.condensate_idx[1] = 1;

    for (bool moist : {false, true}) {
      for (bool water_loading : {false, true}) {
        wl.enabled = water_loading;

        genRandArray(q0_cxx,engine,q_pdf);
        genRandArray(q1_cxx,engine,q_pdf);
        genRandArray(q2_cxx,engine,q_pdf);
        sync_to_host(q0_cxx,q0_f90);
        sync_to_host(q1_cxx,q1_f90);
        sync_to_host(q2_cxx,q2_f90);

        // Run cxx version
        Kokkos::parallel_for(policy,KOKKOS_LAMBDA(const TeamMember& team) {
          KernelVariables kv(team);
          Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team,NP*NP),
                               [&](const int idx) {
            const int igp = idx / NP;
            const int jgp = idx % NP;

            auto R  = Homme::subview(R_cxx_d,kv.ie,igp,jgp);
            auto q = [&](const int iq, const int ilev) -> Scalar {
              return iq==0 ? q0_cxx(kv.ie,igp,jgp,ilev) :
                     (iq==1 ? q1_cxx(kv.ie,igp,jgp,ilev) : q2_cxx(kv.ie,igp,jgp,ilev));
            };

            auto qv = [&](const int ilev) -> Scalar { return q(0,ilev); };
            elem_ops.get_q_dry(kv,moist,wl,q,R);
            elem_ops.get_R_star(kv,moist,wl,R,qv,R);
          });
        });

        // Run f90 version
        const Real* q0_ptr = q0_f90.data();
        const Real* q1_ptr = q1_f90.data();
        const Real* q2_ptr = q2_f90.data();
        Real* R_ptr = R_f90.data();
        compute_r_star_wl_f90(num_elems,moist,water_loading,q0_ptr,q1_ptr,q2_ptr,R_ptr);

        // Compare answers
        auto R_cxx = Kokkos::create_mirror_view(R_cxx_d);
        Kokkos::deep_copy(R_cxx,R_cxx_d);

        for (int ie=0; ie<num_elems; ++ie) {
          for (int igp=0; igp<NP; ++igp) {
            for (int jgp=0; jgp<NP; ++jgp) {
              for (int k=0; k<NUM_PHYSICAL_LEV; ++k) {
                const int ilev = k / VECTOR_SIZE;
                const int ivec = k % VECTOR_SIZE;
                REQUIRE (R_cxx(ie,igp,jgp,ilev)[ivec] == R_f90(ie,k,igp,jgp));
              }
            }
          }
        }
      }
    }
  }

}
