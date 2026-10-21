#ifndef HOMMEXX_ELEMENT_OPS_HPP
#define HOMMEXX_ELEMENT_OPS_HPP

#include "Types.hpp"
#include "KernelVariables.hpp"
#include "HybridVCoord.hpp"
#include "ColumnOps.hpp"
#include "EquationOfState.hpp"
#include "PhysicalConstants.hpp"
#include "SimulationParams.hpp"

#include "utilities/BfbUtils.hpp"

namespace Homme {

class ElementOps {
public:
  ElementOps () = default;

  ~ElementOps () = default;

  void init (const HybridVCoord& hvcoord) {
    m_hvcoord = hvcoord;
  }

  // Dry mass of a layer: the total mass dp minus vapor and condensates.
  // With qcond=0 this is bitwise identical to dp*(1-qv).
  template<typename ScalarT>
  KOKKOS_INLINE_FUNCTION static
  ScalarT compute_dp_dry (const ScalarT& dp, const ScalarT& qv, const ScalarT& qcond) {
    return dp*(1.0 - qv - qcond);
  }

  // Dry mass of each level. dp is total mass (dry air + vapor + condensates), and
  // q(iq,ilev) returns the mixing ratio (per unit dp) of tracer iq, with vapor in
  // tracer 0. The condensates are the tracers listed in wl (none if disabled).
  template<typename QProvider>
  KOKKOS_INLINE_FUNCTION
  void get_dry_mass (const KernelVariables& kv,
                     const WaterLoading& wl,
                     const ExecViewUnmanaged<const Scalar[NUM_LEV]>& dp,
                     const QProvider& q,
                     const ExecViewUnmanaged<Scalar[NUM_LEV]>& dp_dry) const {
    Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team,NUM_LEV),
                         [&](const int ilev) {
      Scalar qcond(0.0);
      if (wl.enabled) {
        for (int ic=0; ic<wl.num_condensates; ++ic) {
          qcond += q(wl.condensate_idx[ic],ilev);
        }
      }
      dp_dry(ilev) = compute_dp_dry(dp(ilev),q(0,ilev),qcond);
    });
  }

  // R* involves only dry air and vapor: R*·dp = Rd·dp_dry + Rv·qv·dp,
  // with dp_dry from get_dry_mass. R may alias dp_dry.
  template<typename InputProvider>
  KOKKOS_INLINE_FUNCTION
  void get_R_star (const KernelVariables& kv,
                   const ExecViewUnmanaged<const Scalar[NUM_LEV]>& dp_dry,
                   const ExecViewUnmanaged<const Scalar[NUM_LEV]>& dp,
                   const InputProvider& qv,
                   const ExecViewUnmanaged<Scalar[NUM_LEV]>& R) const {
    using namespace PhysicalConstants;
    Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team,NUM_LEV),
                         [&](const int ilev) {
      R(ilev) = (Rgas*dp_dry(ilev) + Rwater_vapor*qv(ilev)*dp(ilev)) / dp(ilev);
    });
  }

  // R* for a column. With water loading, dry mass removes the condensates from dp
  // (get_dry_mass), and R* follows from the dry mass. q(iq,ilev) is as in
  // get_dry_mass; it is not called if use_moisture is false.
  template<typename QProvider>
  KOKKOS_INLINE_FUNCTION
  void compute_R_star (const KernelVariables& kv,
                       const bool use_moisture,
                       const WaterLoading& wl,
                       const ExecViewUnmanaged<const Scalar[NUM_LEV]>& dp,
                       const QProvider& q,
                       const ExecViewUnmanaged<Scalar[NUM_LEV]>& R) const {
    auto qv = [&](const int ilev)->Scalar { return q(0,ilev); };
    if (use_moisture && wl.enabled) {
      get_dry_mass(kv,wl,dp,q,R);
      get_R_star(kv,R,dp,qv,R);
    } else {
      // BFB-only: legacy qv-only form, kept so water_loading=false matches baselines; remove in follow-up PR
      get_R_star(kv,use_moisture,qv,R);
    }
  }

  template<typename InputProvider>
  KOKKOS_INLINE_FUNCTION
  void get_R_star (const KernelVariables& kv,
                   const bool use_moisture,
                   const InputProvider& Q,
                   const ExecViewUnmanaged<Scalar[NUM_LEV]>& R) const {
    using namespace PhysicalConstants;
    if (use_moisture) {
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team,NUM_LEV),
                           [&](const int ilev) {
        R(ilev) = (Rgas + (Rwater_vapor-Rgas)*Q(ilev));
      });
    } else {
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team,NUM_LEV),
                           [&](const int ilev) {
        R(ilev) = Rgas;
      });
    }
  }

  KOKKOS_INLINE_FUNCTION
  void compute_hydrostatic_p (const KernelVariables& kv,
                              const ExecViewUnmanaged<const Scalar[NUM_LEV  ]>& dp,
                              const ExecViewUnmanaged<      Scalar[NUM_LEV_P]>& p_i,
                              const ExecViewUnmanaged<      Scalar[NUM_LEV  ]>& pi) const
  {
    p_i(0)[0] = m_hvcoord.hybrid_ai0*m_hvcoord.ps0;
    ColumnOps::column_scan_mid_to_int<true>(kv,dp,p_i);
    ColumnOps::compute_midpoint_values(kv,p_i,pi);
  }

  KOKKOS_FUNCTION
  void get_temperature (const KernelVariables& kv,
                        const EquationOfState& eos,
                        const bool use_moisture,
                        const ExecViewUnmanaged<const Scalar[NUM_LEV]>& dp,
                        const ExecViewUnmanaged<const Scalar[NUM_LEV]>& exner,
                        const ExecViewUnmanaged<const Scalar[NUM_LEV]>& vtheta_dp,
                        const ExecViewUnmanaged<const Scalar[NUM_LEV]>& qv,
                        const ExecViewUnmanaged<Scalar[NUM_LEV]>& Rstar,
                        const ExecViewUnmanaged<Scalar[NUM_LEV]>& T) const {
    using namespace PhysicalConstants;
    get_R_star(kv, use_moisture, qv, Rstar);
    Kokkos::parallel_for(
      Kokkos::ThreadVectorRange(kv.team, NUM_LEV),
      [&] (const int k) { T(k) = Rgas * vtheta_dp(k) * exner(k) / (Rstar(k) * dp(k)); });
  }

  // As above, but with R* from compute_R_star, so that water loading is honored.
  // q(iq,ilev) is as in get_dry_mass.
  template<typename QProvider>
  KOKKOS_INLINE_FUNCTION
  void get_temperature (const KernelVariables& kv,
                        const EquationOfState& eos,
                        const bool use_moisture,
                        const WaterLoading& wl,
                        const ExecViewUnmanaged<const Scalar[NUM_LEV]>& dp,
                        const ExecViewUnmanaged<const Scalar[NUM_LEV]>& exner,
                        const ExecViewUnmanaged<const Scalar[NUM_LEV]>& vtheta_dp,
                        const QProvider& q,
                        const ExecViewUnmanaged<Scalar[NUM_LEV]>& Rstar,
                        const ExecViewUnmanaged<Scalar[NUM_LEV]>& T) const {
    using namespace PhysicalConstants;
    compute_R_star(kv, use_moisture, wl, dp, q, Rstar);
    Kokkos::parallel_for(
      Kokkos::ThreadVectorRange(kv.team, NUM_LEV),
      [&] (const int k) { T(k) = Rgas * vtheta_dp(k) * exner(k) / (Rstar(k) * dp(k)); });
  }

private:

  HybridVCoord    m_hvcoord;
};

} // namespace Homme

#endif // HOMMEXX_ELEMENT_OPS_HPP
