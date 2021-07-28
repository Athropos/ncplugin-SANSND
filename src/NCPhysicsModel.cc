#include "NCPhysicsModel.hh"

//Include various utilities from NCrystal's internal header files:
#include "NCrystal/internal/NCString.hh"
#include "NCrystal/internal/NCRandUtils.hh"

#include <boost/math/special_functions/gamma.hpp>
#include <math.h>
#include <iostream>

bool NCP::PhysicsModel::isApplicable( const NC::Info& info )
{
  //Accept if input is NCMAT data with @CUSTOM_SANSND section:
  return info.countCustomSections(pluginNameUpperCase()) > 0;
}

NCP::PhysicsModel NCP::PhysicsModel::createFromInfo( const NC::Info& info )
{
  //Parse the content of our custom section. In case of syntax errors, it 
  //raises BadInput exceptions
  //Get the relevant custom section data (and verify that there are not multiple
  //such sections in the input data):
  if ( info.countCustomSections( pluginNameUpperCase() ) != 1 )
    NCRYSTAL_THROW2(BadInput,"Multiple @CUSTOM_"<<pluginNameUpperCase()<<" sections are not allowed");
  auto data = info.getCustomSection( pluginNameUpperCase() );

  // data is here a vector of lines, and each line is a vector of words. In our
  // case, we want to accept sections of the form (units are barn and angstrom as
  // is usual in NCrystal):
  //
  // @CUSTOM_SANSND
  //    v                   # plugin version
  //    A s rg m p          # piecewise GP and power law parameters
  //    q1 sigma_0          # boundary parameter and absolute cross section 
  //              

  //Verify we have exactly three lines and 1-5-2 numbers:
  if ( data.size() != 3 || data.at(0).size()!=1 || data.at(1).size()!=5 || data.at(2).size()!=2 )
    NCRYSTAL_THROW2(BadInput,"Data in the @CUSTOM_"<<pluginNameUpperCase()
                    <<" section should be three lines with 1-5-2 numbers");

  //Parse and validate values:
  double supp_version = 2.0;
  double version, A, s, rg, m, p, q1, sigma0 ;
  if ( ! NC::safe_str2dbl( data.at(0).at(0), version )
       || ! NC::safe_str2dbl( data.at(1).at(0), A )
       || ! NC::safe_str2dbl( data.at(1).at(1), s )
       || ! NC::safe_str2dbl( data.at(1).at(2), rg )
       || ! NC::safe_str2dbl( data.at(1).at(3), m )
       || ! NC::safe_str2dbl( data.at(1).at(4), p )
       || ! NC::safe_str2dbl( data.at(2).at(0), q1 )
       || ! NC::safe_str2dbl( data.at(2).at(1), sigma0 )
       || !(p>-2) || !(m>-2) || !(q1>0) || !(sigma0>0)) //condition for a safe integral
    NCRYSTAL_THROW2( BadInput,"Invalid values specified in the @CUSTOM_"<<pluginNameUpperCase()
                     <<" section (see the plugin readme for more info)" );
  //special warning for wrong version
  if ( ! (version==supp_version) ) 
    NCRYSTAL_THROW2( BadInput,"Invalid version specified for the "<<pluginNameUpperCase()
                     <<" plugin. Only the version "<<supp_version<<" is supported." );

  //Parsing done! Create and return our model:
  return PhysicsModel(A, s, rg, m, p, q1, sigma0);
}

NCP::PhysicsModel::PhysicsModel( double A, double s,double rg,double m,double p,double q1,double sigma0 )
  : m_A(A),
    m_s(s),
    m_rg(rg),
    m_m(m),
    m_p(p),
    m_q1(q1),
    m_sigma0(sigma0)
{
  //Important note to developers who are using the infrastructure in the
  //testcode/ subdirectory: If you change the number or types of the arguments
  //for the constructor here, you should make sure to perform a corresponding
  //change in three files in the testcode/ directory: _cbindings.py,
  //__init__.py, and NCForPython.cc - that way you can still instantiate your
  //model directly from your python test code).

  nc_assert( m_s < 3 );
  nc_assert( m_s < m_m );
  nc_assert( m_m > 2 );
  nc_assert( m_p > 2 );
  nc_assert( m_q1 > 0.0 );
  nc_assert( m_sigma0 > 0.0 );
}

double NCP::PhysicsModel::calcCrossSection( double neutron_ekin ) const
{

  double q2 = 1.0/m_rg*sqrt((m_m-m_s)*(3-m_s)/2);
  double B = std::pow(q2,m_m-m_s)*exp((-q2*q2*m_rg*m_rg)/(3-m_s));
  //double C = std::pow(m_q1,m_p-m_s)*exp((-m_q1*m_q1*m_rg*m_rg)/(3-m_s));
  
  double k =  NC::k2Pi/ NC::ekin2wl(neutron_ekin); //wavevector
  //definite power guinier integral from 0 to q1
  double r = m_rg*m_rg/(3-m_s);
  double defint_guinier = std::pow(r,m_s/2-1)/2*boost::math::tgamma_lower(1-m_s/2,r*q2*q2);
  double defint_porod = B*(std::pow(2*k,2-m_m)/(2-m_m) - std::pow(q2,2-m_m)/(2-m_m));
  double total_sigma = (m_sigma0/(2*k*k))*m_A*(defint_guinier+defint_porod);
  return total_sigma;
}

double NCP::PhysicsModel::sampleScatteringVector( NC::RNG& rng, double neutron_ekin ) const 
{
  double rand = rng.generate();
  double Q;
  double k =  NC::k2Pi/NC::ekin2wl(neutron_ekin); //wavevector
  //sample a random scattering vector Q from the inverse CDF (see plugin readme)
  double ratio_sigma = (m_sigma0/(2*k*k))/calcCrossSection(neutron_ekin); //cross section over total cross section ratio
  double CDF_Q0 = (m_A*std::pow(m_q1, m_p+2)/(m_p+2))*ratio_sigma;
  if(rand < CDF_Q0){
    Q = std::pow(((m_p+2)*rand/m_A)/ratio_sigma, 1/(m_p+2));
  } else {
    Q = std::pow((rand/ratio_sigma - (m_A/(m_p+2))*std::pow(m_q1,m_p+2) + (m_A/(m_p+2))*std::pow(m_q1,m_p+2))*(m_p+2)/m_A, 1/(m_p+2));
  }

  return Q;
}
NCP::PhysicsModel::ScatEvent NCP::PhysicsModel::sampleScatteringEvent( NC::RNG& rng, double neutron_ekin ) const
{
  ScatEvent result;

  //section for energy limits
  /*if ( ! (neutron_ekin > 1) ) {
    result.ekin_final = neutron_ekin;
    result.mu = 1.0;
    return result;
  }*/

  //Implement our actual model here:
  result.ekin_final = neutron_ekin;//Elastic
  double Q = sampleScatteringVector(rng, neutron_ekin);
  double ksquared = NC::k4PiSq*NC::ekin2wlsqinv(neutron_ekin);
  result.mu = 1-0.5*(Q*Q/ksquared);

  return result;
}


