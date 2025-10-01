#ifndef FLOWFUNCTIONS_H
#define FLOWFUNCTIONS_H

#include <vector>
#include <array>
#include <utility>
#include <cmath>
#include "TComplex.h"
//#include "TMath.h"

namespace flow_functions
{

  // Constants related to Q-vector components
  constexpr int max_harmonic = 10;
  constexpr int max_power    = 10;

  // Calculate the Q-vectors
  TComplex get_flow_vector(const std::vector<double>& phi_angles, int harmonic);
  TComplex get_weighted_flow_vector(const std::vector<std::pair<double,double>>& phi_weight, int harmonic);
  std::array<TComplex, max_harmonic> get_flow_vectors(const std::vector<double>& phi_angles);
  std::array<TComplex, max_harmonic> get_weighted_flow_vectors(const std::vector<std::pair<double,double>>& phi_weight);
  std::array<std::array<TComplex, max_harmonic>, max_power> get_power_weighted_flow_vectors(const std::vector<std::pair<double,double>>& phi_weight);

  // Some basic correlation functions
  double calccosevent(const std::array<TComplex, max_harmonic>& allQ, int harmonic);
  double calccsinevent(const std::array<TComplex, max_harmonic>& allQ, int harmonic);
  double calc2event(const std::array<TComplex, max_harmonic>& allQ, int harmonic);
  double calcSPevent(const std::array<TComplex, max_harmonic>& allQA, const std::array<TComplex, max_harmonic>& allQB, int harmonic);
  double calccossum2event(const std::array<TComplex, max_harmonic>& allQ, int harmonic);
  double calcsinsum2event(const std::array<TComplex, max_harmonic>& allQ, int harmonic);
  double calccos3event(const std::array<TComplex, max_harmonic>& allQ, int harmonic);
  double calcsin3event(const std::array<TComplex, max_harmonic>& allQ, int harmonic);
  double calc4event(const std::array<TComplex, max_harmonic>& allQ, int harmonic);
  double calc6_event_jamie(TComplex& qn, TComplex& q2n, TComplex& q3n, float M);
  double calc6event(const std::array<TComplex, max_harmonic>& allQ, int harmonic);

} // namespace flow_functions

#endif // FLOWFUNCTIONS_H
