#include "flow_functions.h"

namespace flow_functions
{

TComplex get_flow_vector(const std::vector<double>& phi_angles, const int harmonic)
{
  TComplex Q(0.0,0.0);
  for ( auto it = phi_angles.begin(); it != phi_angles.end(); ++it )
    {
      double phi = *it;
      TComplex u(cos(harmonic*phi),sin(harmonic*phi));
      Q += u;
    }
  return Q;
}

TComplex get_weighted_flow_vector(const std::vector<std::pair<double,double>>& phi_weight, const int harmonic)
{
  TComplex Q(0.0,0.0);
  for ( auto it = phi_weight.begin(); it != phi_weight.end(); ++it )
    {
      double phi = it->first;
      double wgt = it->second;
      TComplex u(cos(harmonic*phi),sin(harmonic*phi));
      Q += wgt*u;
    }
  return Q;
}

// --- wondering if I want to change this...
std::array<TComplex,max_harmonic> get_flow_vectors(const std::vector<double>& phi_angles)
{
  std::array<TComplex,max_harmonic> allQ{};
  for ( int i = 0; i < max_harmonic; ++i )
    {
      allQ[i] = get_flow_vector(phi_angles,i);
    }
  return allQ;
}

// --- because I think I want to do this differently
std::array<std::array<TComplex,max_harmonic>,max_power> get_weighted_flow_vectors(const std::vector<std::pair<double,double>>& phi_weight)
{
  std::vector<std::pair<double,double>> new_phi_weight = phi_weight;
  std::array<std::array<TComplex,max_harmonic>,max_power> allQ{{}};
  for ( int i = 0; i < max_harmonic; ++i )
    {
      for ( int j = 0; j < max_power; ++j )
        {
          // weight needs to be raised to power j
          auto it_old = phi_weight.begin();
          auto it_new = new_phi_weight.begin();
          for ( ; it_old != phi_weight.end() && it_new != new_phi_weight.end(); ++it_old, ++it_new )
            {
              it_new->second = pow(it_old->second,j);
            }
          // use the new weight to get the weighted flow vector
          allQ[i][j] = get_weighted_flow_vector(new_phi_weight,i);
        }
    }
  return allQ;
}

// <cos(nphi)>
double calccosevent(const std::array<TComplex,max_harmonic>& allQ, int harmonic)
{
  double M = allQ[0].Re();
  if ( M < 1 ) return -9999;
  return allQ[harmonic].Re()/M;
}

// <sin(nphi)>
double calccsinevent(const std::array<TComplex,max_harmonic>& allQ, int harmonic)
{
  double M = allQ[0].Re();
  if ( M < 1 ) return -9999;
  return allQ[harmonic].Im()/M;
}

// <cos(n(phi1-phi2))>
double calc2event(const std::array<TComplex,max_harmonic>& allQ, int harmonic)
{
  double M = allQ[0].Re();
  if ( M < 2 ) return -9999;
  // ---
  TComplex Q = allQ[harmonic];
  // ---
  double numerator = Q.Rho2()-M;
  double denominator = M*(M-1);
  // ---
  return numerator/denominator;
}

// <cos(n(phi1a-phi2b))>
double calcSPevent(const std::array<TComplex,max_harmonic>& allQA, const std::array<TComplex,max_harmonic>& allQB, int harmonic)
{
  double MA = allQA[0].Re();
  double MB = allQA[0].Re();
  if ( MA < 1 || MB < 1 ) return -9999;
  // ---
  TComplex QA = allQA[harmonic];
  TComplex QBstar = TComplex::Conjugate(allQB[harmonic]);
  TComplex tc_numerator = QA*QBstar;
  // ---
  double numerator = tc_numerator.Re();
  double denominator = MA*MB;
  // ---
  return numerator/denominator;
}

// <cos(n(phi1+phi2))>
double calccossum2event(const std::array<TComplex,max_harmonic>& allQ, int harmonic)
{
  double M = allQ[0].Re();
  if ( M < 2 ) return -9999;
  // ---
  TComplex Q = allQ[harmonic];
  TComplex Q2 = allQ[2*harmonic];
  TComplex result = Q*Q - Q2;
  // ---
  double numerator = result.Re();
  double denominator = M*(M-1);
  // ---
  return numerator/denominator;
}

// <sin(n(phi1+phi2))>
double calcsinsum2event(const std::array<TComplex,max_harmonic>& allQ, int harmonic)
{
  double M = allQ[0].Re();
  if ( M < 2 ) return -9999;
  // ---
  TComplex Q = allQ[harmonic];
  TComplex Q2 = allQ[2*harmonic];
  TComplex result = Q*Q - Q2;
  // ---
  double numerator = result.Im();
  double denominator = M*(M-1);
  // ---
  return numerator/denominator;
}

// <cos(n(phi1-phi2-phi3))>
double calccos3event(const std::array<TComplex,max_harmonic>& allQ, int harmonic)
{
  double M = allQ[0].Re();
  if ( M < 3 ) return -9999;
  // ---
  TComplex Q = allQ[harmonic];
  TComplex Q2 = allQ[2*harmonic];
  TComplex Qstar = TComplex::Conjugate(allQ[harmonic]);
  TComplex Q2star = TComplex::Conjugate(allQ[2*harmonic]);
  TComplex result = Q*Qstar*Qstar - Q*Q2star;
  // ---
  double numerator = result.Re() - 2*(M-1)*Qstar.Re();
  double denominator = M*(M-1)*(M-2);
  // ---
  return numerator/denominator;
}

// <sin(n(phi1-phi2-phi3))>
double calcsin3event(const std::array<TComplex,max_harmonic>& allQ, int harmonic)
{
  double M = allQ[0].Re();
  if ( M < 3 ) return -9999;
  // ---
  TComplex Q = allQ[harmonic];
  TComplex Q2 = allQ[2*harmonic];
  TComplex Qstar = TComplex::Conjugate(allQ[harmonic]);
  TComplex Q2star = TComplex::Conjugate(allQ[2*harmonic]);
  TComplex result = Q*Qstar*Qstar - Q*Q2star;
  // ---
  double numerator = result.Im() - 2*(M-1)*Qstar.Im();
  double denominator = M*(M-1)*(M-2);
  // ---
  return numerator/denominator;
}

// <cos(n(phi1+phi2-phi3-phi4))>
double calc4event(const std::array<TComplex,max_harmonic>& allQ, int harmonic)
{
  double M = allQ[0].Re();
  if ( M < 4 ) return -9999;
  // ---
  TComplex Q = allQ[harmonic];
  TComplex Q2 = allQ[2*harmonic];
  TComplex Qstar = TComplex::Conjugate(allQ[harmonic]);
  TComplex Q2star = TComplex::Conjugate(allQ[2*harmonic]);
  TComplex tc_three = 2*Q2*Qstar*Qstar;
  // ---
  double one   = pow(Q.Rho2(),2);
  double two   = Q2.Rho2();
  double three = tc_three.Re();
  double four  = 2*(2*(M-2)*Q.Rho2());
  double five  = 2*(M*(M-3));
  // ---
  double numerator = one + two - three - four + five;
  double denominator = M*(M-1)*(M-2)*(M-3);
  // ---
  return numerator/denominator;
}

// <cos(n(phi1+phi2+phi3-phi4-phi5-phi6))>
double calc6_event_jamie(TComplex& qn, TComplex& q2n, TComplex& q3n, float M)
{

  if ( M < 6 ) return -9999;

  // TComplex qn, q2n, q3n;
  // qn = TComplex(Q2x,Q2y);
  // q2n = TComplex(Q4x,Q4y);
  // q3n = TComplex(Q6x,Q6y);

  TComplex temp1;

  // first term
  // |Qn|^6 + 9*|Q2n|^2|Qn|^2 - 6 x Re[Q2n x Qn x Qn* x Qn* x Qn*] / (Mx(M-1)x(M-2)x(M-3)x(M-4)x(M-5)
  double term1a = TMath::Power((qn*TComplex::Conjugate(qn)),3);
  double term1b = 9.0 * q2n*TComplex::Conjugate(q2n) * qn*TComplex::Conjugate(qn);
  temp1 = q2n * qn * TComplex::Conjugate(qn) * TComplex::Conjugate(qn) * TComplex::Conjugate(qn);
  double term1c = -6.0 * temp1.Re();
  double term1 = (term1a+term1b+term1c)/(M*(M-1)*(M-2)*(M-3)*(M-4)*(M-5));

  // second term
  // 4 * [Re[Q3nQn*Qn*Qn*] - 3 Re[Q3nQ2n*Qn*]] / (M(M-1)(M-2)(M-3)(M-4)(M-5)
  temp1 = q3n * TComplex::Conjugate(qn) * TComplex::Conjugate(qn) * TComplex::Conjugate(qn);
  double term2a = temp1.Re();
  temp1 = q3n * TComplex::Conjugate(q2n) * TComplex::Conjugate(qn);
  double term2b = -3.0 * temp1.Re();
  double term2 = 4.0 * (term2a+term2b)/(M*(M-1)*(M-2)*(M-3)*(M-4)*(M-5));

  // third term
  // +2 * (9*(M-4)*Re[Q2nQn*qn*] + 2 |Q3n|^2) / ((M(M-1)(M-2)(M-3)(M-4)(M-5))
  temp1 = q2n*TComplex::Conjugate(qn)*TComplex::Conjugate(qn);
  double term3a = 9.0*(M-4)*temp1.Re();
  double term3b = 2.0*q3n*TComplex::Conjugate(q3n);
  double term3 = 2.0 * (term3a + term3b) / (M*(M-1)*(M-2)*(M-3)*(M-4)*(M-5));

  // fourth term
  //double term4 = -9.0 * (TMath::Power(qn*TComplex::Conjugate(qn),2)+q2n*TComplex::Conjugate(q2n)) / (M*(M-1)*(M-2)*(M-3)*(M-5));
  double term4 = -9.0 * (TMath::Power(qn*TComplex::Conjugate(qn),2)+q2n*TComplex::Conjugate(q2n)) ;
  term4 /= (M*(M-1)*(M-2)*(M-3)*(M-5));

  // fifth term
  //double term5 = 18.0 * qn*TComplex::Conjugate(qn) / (M*(M-1)*(M-3)*(M-4));
  double term5 = 18.0 * qn*TComplex::Conjugate(qn) ;
  term5 /=  (M*(M-1)*(M-3)*(M-4));

  // sixth term
  double term6 = -6.0/((M-1)*(M-2)*(M-3));

  // cos(n(phi1+phi2+phi3-phi4-phi5-phi6))
  double six = term1 + term2 + term3 + term4 + term5 + term6;

  return six;

}

// <cos(n(phi1+phi2+phi3-phi4-phi5-phi6))>
double calc6event(const std::array<TComplex,max_harmonic>& allQ, int harmonic)
{

  double M = allQ[0].Re();
  if ( M < 6 ) return -9999;
  // ---
  TComplex Q = allQ[harmonic];
  TComplex Q2 = allQ[2*harmonic];
  TComplex Q3 = allQ[3*harmonic];
  TComplex Qstar = TComplex::Conjugate(allQ[harmonic]);
  TComplex Q2star = TComplex::Conjugate(allQ[2*harmonic]);
  TComplex Q3star = TComplex::Conjugate(allQ[3*harmonic]);

  TComplex temp1;

  // first term
  // |Qn|^6 + 9*|Q2n|^2|Qn|^2 - 6 x Re[Q2n x Qn x Qn* x Qn* x Qn*] / (Mx(M-1)x(M-2)x(M-3)x(M-4)x(M-5)
  double term1a = TMath::Power((Q*Qstar),3);
  double term1b = 9.0 * Q2*Q2star * Q*Qstar;
  temp1 = Q2 * Q * Qstar * Qstar * Qstar;
  double term1c = -6.0 * temp1.Re();
  double term1 = (term1a+term1b+term1c)/(M*(M-1)*(M-2)*(M-3)*(M-4)*(M-5));

  // second term
  // 4 * [Re[Q3nQn*Qn*Qn*] - 3 Re[Q3nQ2n*Qn*]] / (M(M-1)(M-2)(M-3)(M-4)(M-5)
  temp1 = Q3 * Qstar * Qstar * Qstar;
  double term2a = temp1.Re();
  temp1 = Q3 * Q2star * Qstar;
  double term2b = -3.0 * temp1.Re();
  double term2 = 4.0 * (term2a+term2b)/(M*(M-1)*(M-2)*(M-3)*(M-4)*(M-5));

  // third term
  // +2 * (9*(M-4)*Re[Q2nQn*qn*] + 2 |Q3n|^2) / ((M(M-1)(M-2)(M-3)(M-4)(M-5))
  temp1 = Q2*Qstar*Qstar;
  double term3a = 9.0*(M-4)*temp1.Re();
  double term3b = 2.0*Q3*Q3star;
  double term3 = 2.0 * (term3a + term3b) / (M*(M-1)*(M-2)*(M-3)*(M-4)*(M-5));

  // fourth term
  //double term4 = -9.0 * (TMath::Power(Q*Qstar,2)+Q2*Q2star) / (M*(M-1)*(M-2)*(M-3)*(M-5));
  double term4 = -9.0 * (TMath::Power(Q*Qstar,2)+Q2*Q2star) ;
  term4 /= (M*(M-1)*(M-2)*(M-3)*(M-5));

  // fifth term
  //double term5 = 18.0 * Q*Qstar / (M*(M-1)*(M-3)*(M-4));
  double term5 = 18.0 * Q*Qstar ;
  term5 /=  (M*(M-1)*(M-3)*(M-4));

  // sixth term
  double term6 = -6.0/((M-1)*(M-2)*(M-3));

  // cos(n(phi1+phi2+phi3-phi4-phi5-phi6))
  double six = term1 + term2 + term3 + term4 + term5 + term6;

  return six;

} // double calc6event

} // namespace
