#ifndef JetFlavour_H
#define JetFlavour_H

#include <vector>
#include <cmath>

#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
#include "DataFormats/Math/interface/LorentzVector.h"

  /**
   * Class containing the main MC information on a jet, as extracted from the HepMC event.
   * Some information, see below, depends on the chosen definition, 'algorithmic' or 'physical'
   */
namespace BTagMCTools {

class JetFlavour {
public:

  JetFlavour(void);

  /**
   * Returns the flavour of the jet. This changes according of the definition, 'algorithmic' or 'physical'
   */
  int flavour(void) const { return m_flavour; }

  /**
   * The momentum of the underlying parton. This parton changes according of the
   * definition, 'algorithmic' or 'physical'
   */
  const math::XYZTLorentzVector & underlyingParton4Vec (void) const { return m_underlyingParton4Vec; }

  int  mainFlavour                  (void) const { return m_mainFlavour; }
  int  heaviestFlavour              (void) const { return m_heaviestFlavour; }
  int  minimumDeltaRFlavour         (void) const { return m_minimumDeltaRFlavour; }

  int  mainOrigFlavour              (void) const { return m_mainOrigFlavour; }
  int  originFlavour                (void) const { return m_originFlavour; }
  int  initialFlavour               (void) const { return m_initialFlavour; }

  bool initialPartonSplitsToC       (void) const { return m_initialPartonSplitsToC; }
  bool initialPartonSplitsToB       (void) const { return m_initialPartonSplitsToB; }

  double pMainParton                (void) const { return m_pMainParton; }
  double pClosestParton             (void) const { return m_pClosestParton; }
  double pHeaviestParton            (void) const { return m_pHeaviestParton; }
  double deltaRMainParton           (void) const { return m_deltaRMainParton; }
  double deltaRClosestParton        (void) const { return m_deltaRClosestParton; }
  const math::XYZTLorentzVector & vec4MainParton     (void) const { return m_Vec4MainParton; }
  const math::XYZTLorentzVector & vec4ClosestParton  (void) const { return m_Vec4ClosestParton; }
  const math::XYZTLorentzVector & vec4HeaviestParton (void) const { return m_Vec4HeaviestParton; }

  const math::XYZTLorentzVector & vec4SummedPartons  (void) const { return m_Vec4SummedPartons; }
  const math::XYZTLorentzVector & vec4OriginParton   (void) const { return m_Vec4OriginParton; }

  // if based on heavy hadrons
  bool hasBottomHadronInCone        (void) const { return ( m_mainFlavour == 5 ); }
  bool hasCharmHadronInCone         (void) const { return ( m_mainFlavour == 4 ); }
  bool hasStrangeHadronInCone       (void) const { return ( m_mainFlavour == 3 ); }
  bool hasHeavyHadronInCone         (void) const { return ( m_mainFlavour == 5 || m_mainFlavour == 4 ); }

  bool hasDown                      (void) const { return m_hasDown; }
  bool hasUp                        (void) const { return m_hasUp; }
  bool hasStrange                   (void) const { return m_hasStrange; }
  bool hasCharm                     (void) const { return m_hasCharm; }
  bool hasBottom                    (void) const { return m_hasBottom; }
  bool hasGluon                     (void) const { return m_hasGluon; }

  int  nDown                        (void) const { return m_nDown; }
  int  nUp                          (void) const { return m_nUp; }
  int  nStrange                     (void) const { return m_nStrange; }
  int  nCharm                       (void) const { return m_nCharm; }
  int  nBottom                      (void) const { return m_nBottom; }
  int  nGluon                       (void) const { return m_nGluon; }

  int  numberOfSources              (void) const { return m_numberOfSources; }
  const std::vector<int> & flavourSources (void) const { return m_flavourSources; }

//set methods:

  void flavour(const int a) {m_flavour = a;}
  void underlyingParton4Vec(const math::XYZTLorentzVector & a) { m_underlyingParton4Vec = a;}

  void mainFlavour              (const int a)    { m_mainFlavour = a; }
  void heaviestFlavour          (const int a)    { m_heaviestFlavour = a; }
  void minimumDeltaRFlavour     (const int a)    { m_minimumDeltaRFlavour = a; }

  void mainOrigFlavour          (const int a)    { m_mainOrigFlavour = a; }
  void originFlavour            (const int a)    { m_originFlavour = a; }
  void initialFlavour           (const int a)    { m_initialFlavour = a; }

  void initialPartonSplitsToC   (const bool a)   { m_initialPartonSplitsToC = a; }
  void initialPartonSplitsToB   (const bool a)   { m_initialPartonSplitsToB = a; }

  void pMainParton              (const double a) { m_pMainParton = a; }
  void pClosestParton           (const double a) { m_pClosestParton = a; }
  void pHeaviestParton          (const double a) { m_pHeaviestParton = a; }
  void deltaRMainParton         (const double a) { m_deltaRMainParton = a; }
  void deltaRClosestParton      (const double a) { m_deltaRClosestParton = a; }
  void vec4MainParton           (const math::XYZTLorentzVector & a) { m_Vec4MainParton = a; }
  void vec4ClosestParton        (const math::XYZTLorentzVector & a) { m_Vec4ClosestParton = a; }
  void vec4HeaviestParton       (const math::XYZTLorentzVector & a) { m_Vec4HeaviestParton = a; }

  void vec4SummedPartons        (const math::XYZTLorentzVector & a) { m_Vec4SummedPartons = a; }
  void vec4OriginParton         (const math::XYZTLorentzVector & a) { m_Vec4OriginParton = a; }

  // if based on heavy hadrons

  void hasDown          (const bool a)  { m_hasDown = a; }
  void hasUp            (const bool a)  { m_hasUp = a; }
  void hasStrange       (const bool a)  { m_hasStrange = a; }
  void hasCharm         (const bool a)  { m_hasCharm = a; }
  void hasBottom        (const bool a)  { m_hasBottom = a; }
  void hasGluon         (const bool a)  { m_hasGluon = a; }

  void nDown            (const int a)   { m_nDown = a; }
  void nUp              (const int a)   { m_nUp = a; }
  void nStrange         (const int a)   { m_nStrange = a; }
  void nCharm           (const int a)   { m_nCharm = a; }
  void nBottom          (const int a)   { m_nBottom = a; }
  void nGluon           (const int a)   { m_nGluon = a; }

  void numberOfSources  (const int a)   { m_numberOfSources = a; }
  void flavourSources   (const std::vector<int> & a) { m_flavourSources = a; }

private:

  int    m_flavour;             // THE flavour, depends on whether the physics of the algo def is chosen
  math::XYZTLorentzVector m_underlyingParton4Vec;
  int    m_mainFlavour;         // main is the hardest
  int    m_heaviestFlavour;     // heaviest does not include gluons
  int    m_minimumDeltaRFlavour;

  int    m_originFlavour;       // the one of the status==3 mother
  int    m_mainOrigFlavour;     // mother flavour of main parton

  double m_pMainParton;         // mom. of main parton
  double m_pClosestParton;      // mom. of closest parton
  double m_pHeaviestParton;     // mom. of heaviest parton
  double m_deltaRMainParton;    // deltaR of main parton
  double m_deltaRClosestParton; // deltaR of closest parton

  math::XYZTLorentzVector m_Vec4MainParton;
  math::XYZTLorentzVector m_Vec4ClosestParton;
  math::XYZTLorentzVector m_Vec4HeaviestParton;

  math::XYZTLorentzVector m_Vec4SummedPartons;
  math::XYZTLorentzVector m_Vec4OriginParton;

  bool m_heavyHadronBased;
  bool m_partonBased;

  bool m_hasDown;
  bool m_hasUp;
  bool m_hasStrange;
  bool m_hasCharm;
  bool m_hasBottom;
  bool m_hasGluon;

  int m_nDown;
  int m_nUp;
  int m_nStrange;
  int m_nCharm;
  int m_nBottom;
  int m_nGluon;

  int m_numberOfSources;

  std::vector<int> m_flavourSources;

  int m_initialFlavour;        // when matching to initial partons only

  bool m_initialPartonSplitsToC;
  bool m_initialPartonSplitsToB;

};
}
#endif
