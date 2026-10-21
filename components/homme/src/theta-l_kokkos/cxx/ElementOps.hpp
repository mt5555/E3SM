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

  // Dry-air mass fraction of a layer, q_dry = dp_dry/dp. Like qv and qcond, it is
  // per unit total mass dp (dry air + vapor + condensates).
  template<typename ScalarT>
  KOKKOS_INLINE_FUNCTION static
  ScalarT compute_q_dry (const ScalarT& qv, const ScalarT& qcond) {
    return 1.0 - qv - qcond;
  }

  // Dry-air mass fraction q_dry of each level. q(iq,ilev) returns the mixing ratio
  // (per unit total mass dp) of tracer iq, with vapor in tracer 0. The condensates
  // are the tracers listed in wl (none if disabled). If use_moisture is false, the
  // air is all dry (q is not called, since qsize might be 0).
  template<typename QProvider>
  KOKKOS_INLINE_FUNCTION
  void get_q_dry (const KernelVariables& kv,
                  const bool use_moisture,
                  const WaterLoading& wl,
                  const QProvider& q,
                  const ExecViewUnmanaged<Scalar[NUM_LEV]>& q_dry) const {
    if (not use_moisture) {
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team,NUM_LEV),
                           [&](const int ilev) {
        q_dry(ilev) = 1.0;
      });
      return;
    }
    Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team,NUM_LEV),
                         [&](const int ilev) {
      Scalar qcond(0.0);
      if (wl.enabled) {
        for (int ic=0; ic<wl.num_condensates; ++ic) {
          qcond += q(wl.condensate_idx[ic],ilev);
        }
      }
      q_dry(ilev) = compute_q_dry(q(0,ilev),qcond);
    });
  }

  // R* involves only dry air and vapor: R* = Rd*q_dry + Rv*qv, with q_dry from
  // get_q_dry and qv both per unit total mass. With water loading, q_dry also
  // excludes the condensates. R may alias q_dry. qv is not called if
  // use_moisture is false.
  template<typename InputProvider>
  KOKKOS_INLINE_FUNCTION
  void get_R_star (const KernelVariables& kv,
                   const bool use_moisture,
                   const WaterLoading& wl,
                   const ExecViewUnmanaged<const Scalar[NUM_LEV]>& q_dry,
                   const InputProvider& qv,
                   const ExecViewUnmanaged<Scalar[NUM_LEV]>& R) const {
    using namespace PhysicalConstants;
    if (use_moisture && wl.enabled) {
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team,NUM_LEV),
                           [&](const int ilev) {
        R(ilev) = Rgas*q_dry(ilev) + Rwater_vapor*qv(ilev);
      });
    } else {
      // BFB-only: legacy qv-only form, kept so water_loading=false matches baselines; remove in follow-up PR
      get_R_star(kv,use_moisture,qv,R);
    }
  }

  // Legacy qv-only R*, R* = Rd + (Rv-Rd)*qv.
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

  // Temperature from vtheta_dp, with R* from the dry-air mass fraction (get_q_dry,
  // then get_R_star). q(iq,ilev) is as in get_q_dry. Rstar is scratch.
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
    // Rstar holds q_dry, then R* (computed in place)
    get_q_dry(kv, use_moisture, wl, q, Rstar);
    const auto qv = [&](const int ilev)->Scalar { return q(0,ilev); };
    get_R_star(kv, use_moisture, wl, Rstar, qv, Rstar);
    Kokkos::parallel_for(
      Kokkos::ThreadVectorRange(kv.team, NUM_LEV),
      [&] (const int k) { T(k) = Rgas * vtheta_dp(k) * exner(k) / (Rstar(k) * dp(k)); });
  }

private:

  HybridVCoord    m_hvcoord;
};

} // namespace Homme

#endif // HOMMEXX_ELEMENT_OPS_HPP
