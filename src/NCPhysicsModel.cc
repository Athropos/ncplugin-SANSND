#include "NCPhysicsModel.hh"

//Include various utilities from NCrystal's internal header files:
#include "NCrystal/internal/NCString.hh"
#include "NCrystal/internal/NCRandUtils.hh"

#include <boost/math/special_functions/gamma.hpp>
#include <sys/stat.h>
#include <algorithm>
#include <fstream>
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
  //    filename            # read file mode
  //   
  //    or
  //
  //    v                   # plugin version
  //    A s rg m p          # piecewise GP and power law parameters
  //    sigma_0          # boundary parameter and absolute cross section 
  //              

  bool file_mode=false;
  //Verify we have exactly 1-2 (file_mode on) or 1-5-1 input :
  if ( data.size() == 2 ) {
    if ( data.at(0).size()!=1 || data.at(1).size()!=1 ) {
    NCRYSTAL_THROW2(BadInput,"Bad input for file mode in the @CUSTOM_"<<pluginNameUpperCase()
                    <<" section (see the plugin readme for more info).");
    } else { file_mode=true; }   
  } else if ( data.size() == 3 ) {
    if ( data.at(0).size()!=1 || data.at(1).size()!=5 ||  data.at(2).size()!=1  ) 
    NCRYSTAL_THROW2(BadInput,"Bad input for fit mode in the @CUSTOM_"<<pluginNameUpperCase()
                    <<" section (see the plugin readme for more info).")   
  } else {
    NCRYSTAL_THROW2(BadInput,"Bad input in the @CUSTOM_"<<pluginNameUpperCase()
                    <<" section (see the plugin readme for more info).");    
  }
  //Parse and validate values:
  double supp_version = 2.0;
  double version, A, s, rg, m, p, sigma0;
  if ( ! NC::safe_str2dbl( data.at(0).at(0), version ))
    NCRYSTAL_THROW2( BadInput,"Invalid version input in the @CUSTOM_"<<pluginNameUpperCase()
                    <<" section (see the plugin readme for more info).");  
  //special warning for wrong version
  if ( ! (version==supp_version) ) 
    NCRYSTAL_THROW2( BadInput,"Invalid version specified for the "<<pluginNameUpperCase()
                      <<" plugin. Only the version "<<supp_version<<" is supported." );
  
  if (file_mode){
    std::string filename= data.at(1).at(0);
    std::string root_rel = "../data/";
    std::string rel_path = root_rel+filename;
    //check existence
    struct stat buffer;   
    if (!(stat (rel_path.c_str(), &buffer) == 0))
      NCRYSTAL_THROW2( BadInput,"The filename specified for the "<<pluginNameUpperCase()
                      <<" plugin is invalid or the file could not be found in the data/ directory. " );
    //Checks done! Create and return our model:
    return PhysicsModel(filename);
  } else {
      if ( ! NC::safe_str2dbl( data.at(2).at(0), A )
      || ! NC::safe_str2dbl( data.at(2).at(1), s )
      || ! NC::safe_str2dbl( data.at(2).at(2), rg )
      || ! NC::safe_str2dbl( data.at(2).at(3), m )
      || ! NC::safe_str2dbl( data.at(2).at(4), p )
      || ! NC::safe_str2dbl( data.at(3).at(1), sigma0 )
      || !(sigma0>0)) 
      NCRYSTAL_THROW2( BadInput,"Invalid values specified in the @CUSTOM_"<<pluginNameUpperCase()
                      <<" section (see the plugin readme for more info)" );

    //Parsing done! Create and return our model:
    return PhysicsModel(A, s, rg, m, p, sigma0);
  }
  
}

NCP::PhysicsModel::PhysicsModel(std::string filename)
  : m_filename(filename),
    m_q(),
    m_I_q(),
    m_CDF(),
    m_f_integral(),
    m_sigma0(5.551),
    extr_porod_par(),
    m_xs_max_k(),
    m_xs_max()
{   
  //Parse the input file and create the vector with the q, Iq info
  std::string root_rel = "data/";
  std::string rel_path = root_rel+m_filename;
  std::ifstream input_file(rel_path);
  if(input_file) {
    double temp_x, temp_y;
    std::string line;
    while (std::getline(input_file, line)){
          std::istringstream iss(line);
          iss >> temp_x >> temp_y;
          m_q.push_back(temp_x);
          m_I_q.push_back(temp_y);
    }
  } else {
    NCRYSTAL_THROW2( BadInput,rel_path << ": Invalid data file for the "<<pluginNameUpperCase()
                    <<" plugin" );   
  }
  //DEBUG ONLY
  /*for(int i=0;i<=m_q.size()-1;i++){
    std::cout<< m_q[i] << std::endl;
  }*/

  //Fits porod a*x^b for high-q, after given q_max, with the last 5 points
  int i;
  double m=0, x_bar=0, t=0, a=0, b=0, num=0, den=0;
  for(i=m_q.size()-6;i<=m_q.size()-1;i++){
    m += (m_I_q[i]*m_I_q[i]);
    x_bar += (m_I_q[i]*m_I_q[i])*log(m_q[i]);
  }
  x_bar /= m;
  for(i=m_q.size()-6;i<=m_q.size()-1;i++){
    t += (m_I_q[i]*m_I_q[i])*log(m_q[i] - x_bar)*log(m_q[i] - x_bar);
  }
  for(i=m_q.size()-6;i<=m_q.size()-1;i++){
      b += -(m_I_q[i]*m_I_q[i])*log(m_I_q[i])*log(m_q[i] - x_bar);
  }
  b /= t;
  extr_porod_par.push_back(b);
  for(i=m_q.size()-6;i<=m_q.size()-1;i++){
      num += std::pow(m_q[i],b)*m_I_q[i];
      den += std::pow(m_q[i],2*b);
  }
  a = num/den;
  extr_porod_par.push_back(a);
  std::cout << "Extrapolation after q_max with Porod a*x^b with coefficients:\n a=" << a << " and b=" << b << std::endl;

  //Assign integral function
  m_f_integral = [this](double k) {
    double sum=0;
    double xs_corr = 0;
    // upper_bound is the index of first m_q not smaller than k 
    int upper_bound =  std::distance(m_q.begin(),std::lower_bound(m_q.begin(),m_q.end(), k));
    int lower_bound;
    //check also if k is equal to any q bin boundary : 
    int is_equal = std::count(m_q.begin(),m_q.end(), k);
    //Find the interval that contains k assuming q is made of increasing, not duplicated q. 
    if (upper_bound == 0){
      // k smaller or equal to then the smallest q, upper_bound = 0
      // STILL TO FIND A SOLUTION
    } else if (upper_bound > m_q.size()-1){
      // k bigger then the biggest q, upper_bound = size of q + 1
      double a=extr_porod_par.at(1);
      double b=extr_porod_par.at(0);
      lower_bound = m_q.size()-1;
      //Evaluate the integral qI(q) from q_min to q_max
      for (size_t i = 0; i < lower_bound; i++)
        {
          // sum from q_min and not 0, small mistake in total xs (?). 
          sum += (m_q.at(i+1)*m_I_q.at(i+1) + m_q.at(i)*m_I_q.at(i))*(m_q.at(i+1) - m_q.at(i));
        }
      // Add missing xs integrating qI(q) from q_max to k with Porod extrapolation
      xs_corr = a*(std::pow(k,2+b)/(2+b) - std::pow(m_q.at(m_q.size()-1),2+b)/(2+b));
    } else {
      // found interval
      if (is_equal==1){
        // k equal to one of the edges
        //Evaluate the integral qI(q) from q_min to exactly k
        for (size_t i = 0; i < upper_bound; i++)
        {
          // sum from q_min and not 0, small mistake in total xs (?). 
          sum += (m_q.at(i+1)*m_I_q.at(i+1) + m_q.at(i)*m_I_q.at(i))*(m_q.at(i+1) - m_q.at(i));
        }
      } else {
        // k is in the middle
        lower_bound=upper_bound-1;
        // find linear interpolation ratio for k on the q-axis
        double ratio = (k-m_q.at(lower_bound))/(m_q.at(upper_bound)-m_q.at(lower_bound));
        // use the same ratio on the y-axis to calculate the interpolated value of I(k)
        double I_q_k=m_I_q.at(upper_bound)*ratio+m_I_q.at(lower_bound)*(1-ratio);
        //Evaluate the integral qI(q) from q_min to k
        for (size_t i = 0; i < lower_bound; i++)
        {
          // sum from q_min and not 0, small mistake in total xs (?). 
          sum += (m_q.at(i+1)*m_I_q.at(i+1) + m_q.at(i)*m_I_q.at(i))*(m_q.at(i+1) - m_q.at(i));
        }
        // evaluate the correction to the xs as integral of qI(q) from lower_bound to k 
        xs_corr = ((k)*I_q_k + m_q.at(lower_bound)*m_I_q.at(lower_bound))*(k - m_q.at(lower_bound));
      }
    }
    //try to add a correction for low q constant based on the fit (temporary for debugging)
    //the correction is the integrated guinier in the range 0-q_min
    double r = 14.64*14.64/(3-1.5);
    sum += 66*std::pow(r,1.5/2-1)/2*(boost::math::tgamma_lower(1-1.5/2,r*m_q.at(0)*m_q.at(0)));
    
    return (m_sigma0/2)*(sum+xs_corr);
    
  };
  //Evaluate CDF as xs(q)/xs(q_max) at the same q as m_q
  m_xs_max = m_f_integral(m_q.back());
  for (size_t i = 0; i < m_q.size()-1; i++){
    //Here we make the approximation of negligible changes in the CDF for q>q_max, CDF(q_max) = 1
    m_CDF.push_back(m_f_integral(m_q.at(i))/m_xs_max);
  }
  std::cout << "Initialization SANSND done\n";
};

NCP::PhysicsModel::PhysicsModel( double A, double s,double rg,double m,double p,double sigma0 )
  : m_A(A),
    m_s(s),
    m_rg(rg),
    m_m(m),
    m_p(p),
    m_sigma0(sigma0),
    m_q2(),
    m_B(),
    m_r(),
    m_q(),
    m_CDF(),
    m_xs_max_k(),
    m_xs_max(),
    m_def_int_guinier(),
    m_f_integral()
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
  nc_assert( m_sigma0 > 0.0 );

  //Parameters
  m_q2 = 1.0/m_rg*sqrt((m_m-m_s)*(3-m_s)/2);
  m_B = std::pow(m_q2,m_m-m_s)*exp((-m_q2*m_q2*m_rg*m_rg)/(3-m_s));
  m_r = m_rg*m_rg/(3-m_s);
  m_def_int_guinier = std::pow(m_r,m_s/2-1)/2*boost::math::tgamma_lower(1-m_s/2,m_r*m_q2*m_q2);
  //Assign integral function
  m_f_integral = [this](double k) {
    double def_int_porod = m_B*(std::pow(k,2-m_m)/(2-m_m) - std::pow(m_q2,2-m_m)/(2-m_m));
    return (m_sigma0/2)*m_A*(m_def_int_guinier+def_int_porod);
  };
  //TO DO!
  //generate a vector of q to sample numerically the CDF (log or linear? use NCMath.hh )
  //Evaluate CDF as xs(q)/xs(q_max) at the same q as m_q
}

double NCP::PhysicsModel::calcCrossSection( double neutron_ekin ) const
{
  
  double k =  NC::k2Pi/ NC::ekin2wl(neutron_ekin); //wavevector
  //Calculate total sigma using the suitable m_f_integral
  m_xs_max_k = m_f_integral(2*k);
  double total_sigma = 1.0/(k*k)*m_xs_max_k;
  return total_sigma;
}

double NCP::PhysicsModel::sampleScatteringVector( NC::RNG& rng, double neutron_ekin ) const 
{
  double rand = rng.generate();
  double Q;
  double k =  NC::k2Pi/NC::ekin2wl(neutron_ekin); //wavevector
  //If 2*k < q_max then adjust rand to give always CDF^-1(1) = 2*k
  if (2*k<m_q.back()) { rand *= m_xs_max_k/m_xs_max; }

  //sample a random scattering vector Q
  //Find the CDF interval that contains rand
  int lower_bound_CDF;
  // upper_bound_CDF is the first m_CDF index with a CDF not smaller than rand
  int upper_bound_CDF = std::distance(m_CDF.begin(), 
                                  std::lower_bound (m_CDF.begin(), 
                                                    m_CDF.end(), rand) ) - 1 ;
  if (upper_bound_CDF==0){
    // rand smaller then the smallest CDF bin, then no info: Q can be any vecotor between 0 and q_min
    Q = m_q.at(0)*rng.generate();
  } else if (upper_bound_CDF == m_q.size()-1){
    // rand = 1 and 2*k>q_max, so no info: Q can be any vector between q_max and 2*k
    Q = (2*k-m_q.back())*rng.generate() + m_q.back();
  } else {
    // found interval
    lower_bound_CDF=upper_bound_CDF-1;
    // find fractional location of rand relative to lower and higher
    double ratio = (rand - m_CDF.at(lower_bound_CDF))/(m_CDF.at(upper_bound_CDF)-m_CDF.at(lower_bound_CDF));
    // use the same ratio on the q-axis to calculate the log-interpolated value of q to scatter
    Q =std::pow(m_q.at(upper_bound_CDF),ratio)*std::pow(m_q.at(lower_bound_CDF),1-ratio);
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


