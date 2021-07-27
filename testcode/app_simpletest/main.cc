#include "NCPhysicsModel.hh"
#include <iostream>
#include "NCrystal/internal/NCMath.hh"//for NC::linspace

int main()
{
  //Very simple test which instantiates our model and calculates a few cross
  //sections and samples a few scattering events:
  
  auto pm = NCP::PhysicsModel(66.40901846,  1.50279313, 14.64271592,  4.14778786, 2.13191984, 0.016, 5.551);

  for ( auto wl : NC::linspace(0.01, 8.0, 20) ) {
    std::cout << "cross section @ " << wl << " Aa is "
              << pm.calcCrossSection( NC::wl2ekin(wl) ) <<" barn" << std::endl;
  }

  auto rng = NC::getRNG();

  for ( auto wl : NC::linspace(9, 11, 2) ) {
    for (unsigned i = 0; i < 10; ++i) {
      auto outcome = pm.sampleScatteringEvent( *rng, NC::wl2ekin(wl) );
      if (std::isnan(outcome.mu)) {
        std::cout << "scattering @ " << wl << " Aa gives neutron with wl = "
                << NC::ekin2wl(outcome.ekin_final) <<" Aa and a scattering   "
                << std::acos(outcome.mu)*NC::kToDeg<<" degrees. The " 
                << std::endl;
      }
    }
  }



  return 0;
}
