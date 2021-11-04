#include "NCPhysicsModel.hh"

// Include various utilities from NCrystal's internal header files:
#include "NCrystal/internal/NCString.hh"
#include "NCrystal/internal/NCRandUtils.hh"

#include <algorithm>
#include <sys/stat.h>
#include <fstream>
#include <iostream>

bool NCP::PhysicsModel::isApplicable(const NC::Info &info)
{
  // Accept if input is NCMAT data with @CUSTOM_SANSND section:
  return info.countCustomSections(pluginNameUpperCase()) > 0;
}

NCP::PhysicsModel NCP::PhysicsModel::createFromInfo(const NC::Info &info)
{
  // Parse the content of our custom section. In case of syntax errors, it
  // raises BadInput exceptions
  // Get the relevant custom section data (and verify that there are not multiple
  // such sections in the input data):
  if (info.countCustomSections(pluginNameUpperCase()) != 1)
    NCRYSTAL_THROW2(BadInput, "Multiple @CUSTOM_" << pluginNameUpperCase() << " sections are not allowed");
  auto data = info.getCustomSection(pluginNameUpperCase());

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
  //    model               # which model
  //    parametrs           # paramters for the selected model
  //

  bool file_mode = false;
  // Verify we have exactly 1-2 (file_mode on) or 1-1 input (fit mode) :
  if (data.size() == 2)
  {
    if (data.at(0).size() != 1 || data.at(1).size() != 1)
    {
      NCRYSTAL_THROW2(BadInput, "Bad input for file mode in the @CUSTOM_" << pluginNameUpperCase()
                                                                          << " section (see the plugin readme for more info).");
    }
    else
    {
      file_mode = true;
    }
  }
  else if (data.size() == 3)
  {
    if (data.at(0).size() != 1 || data.at(1).size() != 1)
      NCRYSTAL_THROW2(BadInput, "Bad input for fit mode in the @CUSTOM_" << pluginNameUpperCase()
                                                                         << " section (see the plugin readme for more info).")
  }
  else
  {
    NCRYSTAL_THROW2(BadInput, "Bad input in the @CUSTOM_" << pluginNameUpperCase()
                                                          << " section (see the plugin readme for more info).");
  }
  // Parse and validate values:
  double supp_version = 2.0;
  double version;
  // NC::VectD param;
  double p0, p1, p2, p3, p4;
  if (!NC::safe_str2dbl(data.at(0).at(0), version))
    NCRYSTAL_THROW2(BadInput, "Invalid version input in the @CUSTOM_" << pluginNameUpperCase()
                                                                      << " section (see the plugin readme for more info).");
  // special warning for wrong version
  if (!(version == supp_version))
    NCRYSTAL_THROW2(BadInput, "Invalid version specified for the " << pluginNameUpperCase()
                                                                   << " plugin. Only the version " << supp_version << " is supported.");

  if (file_mode)
  {
    std::string filename = data.at(1).at(0);
    std::string root_rel = "data/";
    std::string rel_path = root_rel + filename;
    std::cout << "Input file: " << rel_path << std::endl;
    // check existence
    struct stat buffer;
    if (!(stat(rel_path.c_str(), &buffer) == 0))
      NCRYSTAL_THROW2(BadInput, "The filename specified for the " << pluginNameUpperCase()
                                                                  << " plugin is invalid or the file could not be found in the data/ directory. ");
    // Checks done! Create and return our model:
    return PhysicsModel(filename);
  }
  else
  {
    std::string model = data.at(1).at(0);
    if (model == "GP")
    {
      if (!NC::safe_str2dbl(data.at(2).at(0), p0) || !NC::safe_str2dbl(data.at(2).at(1), p1) || !NC::safe_str2dbl(data.at(2).at(2), p2) || !NC::safe_str2dbl(data.at(2).at(3), p3) || !NC::safe_str2dbl(data.at(2).at(4), p4))
      {
        NCRYSTAL_THROW2(BadInput, "Invalid values specified for " << model << " model in the @CUSTOM_" << pluginNameUpperCase()
                                                                  << " section (see the plugin readme for more info)");
      }
      else
      {
        // CHECK THE INPUT PARAM

        // param.insert(param.end(), {A,s,rg,m});
        // Parsing done! Create and return our model:
        return PhysicsModel(model, p0, p1, p2, p3, p4);
      }
    }
    else if (model == "PPF")
    {
      if (!NC::safe_str2dbl(data.at(2).at(0), p0) || !NC::safe_str2dbl(data.at(2).at(1), p1) || !NC::safe_str2dbl(data.at(2).at(2), p2) || !NC::safe_str2dbl(data.at(2).at(3), p3) || !NC::safe_str2dbl(data.at(2).at(4), p4)
          //|| ! NC::safe_str2dbl( data.at(2).at(5), sigma0 )
      )
      {
        NCRYSTAL_THROW2(BadInput, "Invalid values specified for " << model << " model in the @CUSTOM_" << pluginNameUpperCase()
                                                                  << " section (see the plugin readme for more info)");
      }
      else
      {
        // CHECK THE INPUT PARAM
        // nc_assert(Q0>0);
        // nc_assert(sigma0>0);
        // param.insert(param.end(), {A1,b1,A2,b2,Q0,sigma0});
        // Parsing done! Create and return our model:
        return PhysicsModel(model, p0, p1, p2, p3, p4);
      }
    }
    else if (model == "NP-FBA")
    {
      if (NC::safe_str2dbl(data.at(2).at(0), p0)){
        return PhysicsModel(model , p0 );
      } else {
        std::string filename = data.at(2).at(0);
        std::string root_rel = "data/";
        std::string rel_path = root_rel + filename;
        std::cout << "Input file: " << rel_path << std::endl;
        // check existence
        struct stat buffer;
        if (!(stat(rel_path.c_str(), &buffer) == 0))
        {
          NCRYSTAL_THROW2(BadInput, "The filename specified for the " << pluginNameUpperCase()
                                                                      << " plugin is invalid or the file could not be found in the data/ directory. ");
        }
        else
        {
          // CHECK THE INPUT PARAM
          // nc_assert(R>0);
          // param.insert(param.end(), {R});}
          return PhysicsModel(model , filename );
        }    
      }  
    }
    else if (model == "simple")
    {
      if (NC::safe_str2dbl(data.at(2).at(0), p0)){
        // Parsing done! Create and return our model:
        return PhysicsModel(model, p0);
      }
    }
    else
    {
      NCRYSTAL_THROW2(BadInput, "Invalid model specified in the @CUSTOM_" << pluginNameUpperCase()
                                                                          << " section (see the plugin readme for more info)");
    }
  }
}

NCP::PhysicsModel::PhysicsModel(std::string model, double p0, double p1, double p2, double p3, double p4)
    : m_model(model),
      m_param({p0, p1, p2, p3, p4}),
      m_helper(([this]() -> NC::IofQHelper
                { 
      
      //Generate vector of data q and IofQ
      NC::VectD q = NC::logspace(-6,1,100000);
      NC::VectD IofQ = q;
      if (m_model == "GP") {
        double A=m_param.value().at(0);
        double s=m_param.value().at(1);
        double rg=m_param.value().at(2);
        double m=m_param.value().at(3);
        double p=m_param.value().at(4);
        //Q1 is when IofQ stops being evaluated as power law and Guinier starts (maybe parameter?)
        double Q1=0.016;
        auto it_q1 = std::lower_bound(IofQ.begin(),IofQ.end(), Q1);
        //Approximation valid as long as we have high sampling. Otherwise interpolation needed
        if (it_q1==IofQ.end()) 
          NCRYSTAL_THROW2( BadInput,"Invalid parameters, Q1 bigger then 10 AA-1 in the @CUSTOM_"<<pluginNameUpperCase()
                            <<" section (see the plugin readme for more info)" );
        if (it_q1==IofQ.begin()) 
          NCRYSTAL_THROW2( BadInput,"Invalid parameters, Q1 smaller then 1e-5 AA-1 in the @CUSTOM_"<<pluginNameUpperCase()
                            <<" section (see the plugin readme for more info)" );
        double C = A*std::pow(Q1,p-s)*std::exp((-Q1*Q1*rg*rg)/(3-s));
        std::for_each(IofQ.begin(),it_q1,
                      [C,p](double &x) { x = C*std::pow(x,-p);}
                      ); 
        std::advance(it_q1,1);
        //evaluate Q2 where IofQ stops being evaluated as Guinier and Porod starts
        double Q2 = 1.0/rg*std::sqrt((m-s)*(3-s)/2);
        //IofQ still filled with q here
        auto it_q2 = std::lower_bound(it_q1,IofQ.end(), Q2);
        //Approximation valid as long as we have high sampling. Otherwise interpolation needed
        if (it_q2==IofQ.end()) 
          NCRYSTAL_THROW2( BadInput,"Invalid parameters, Q2 bigger then 10 AA-1 in the @CUSTOM_"<<pluginNameUpperCase()
                            <<" section (see the plugin readme for more info)" );
        if (it_q2==it_q1) 
          NCRYSTAL_THROW2( BadInput,"Invalid parameters, Q1==Q2 in the @CUSTOM_"<<pluginNameUpperCase()
                            <<" section (see the plugin readme for more info)" );
        double B = A*std::pow(Q2,m-s)*std::exp((-Q2*Q2*rg*rg)/(3-s));
        std::for_each(it_q1,it_q2,
                      [A,s,rg](double &x) { x = A*std::pow(x,-s)*std::exp((-x*x*rg*rg)/(3-s));}
                      ); 
        std::advance(it_q2,1);
        std::for_each(it_q2,IofQ.end(),
                      [B,m](double &x) { x = B*std::pow(x,-m);}
                      );
      } else if (m_model == "PPF") {
        double A1=m_param.value().at(0);
        double b1=m_param.value().at(1);
        double A2=m_param.value().at(2);
        double b2=m_param.value().at(3);
        double Q0=m_param.value().at(4);

        auto it_q0 = std::lower_bound(IofQ.begin(),IofQ.end(), Q0);
        if (it_q0==IofQ.end()) 
          NCRYSTAL_THROW2( BadInput,"Invalid parameters, Q0 bigger then 10 AA-1 in the @CUSTOM_"<<pluginNameUpperCase()
                            <<" section (see the plugin readme for more info)" );
        std::for_each(IofQ.begin(),it_q0,
                      [A1,b1](double &x) { x = A1*std::pow(x,-b1);}
                      ); 
        std::advance(it_q0,1);
        std::for_each(it_q0,IofQ.end(),
                      [A2,b2](double &x) { x = A2*std::pow(x,-b2);}
                      );    
      } 
      //Initialize the helper           
      NC::IofQHelper helper(q,IofQ);
      return helper; })())
{
  // Important note to developers who are using the infrastructure in the
  // testcode/ subdirectory: If you change the number or types of the arguments
  // for the constructor here, you should make sure to perform a corresponding
  // change in three files in the testcode/ directory: _cbindings.py,
  //__init__.py, and NCForPython.cc - that way you can still instantiate your
  // model directly from your python test code).
}

NCP::PhysicsModel::PhysicsModel(std::string model, double p)
    : m_model(model),
      m_p(p)
{     
};

NCP::PhysicsModel::PhysicsModel(std::string model, std::string filename)
    : m_model(model),
      m_R(),
      m_freq()
{
  // read R distribution information
  NC::VectD R;
  NC::VectD freq;
  std::string root_rel = "data/";
  std::string rel_path = root_rel + filename;
  std::ifstream input_file(rel_path);
  if (input_file)
  {
    double temp_x, temp_y;
    std::string line;
    while (std::getline(input_file, line))
    {
      std::istringstream iss(line);
      iss >> temp_x >> temp_y;
      R.push_back(temp_x);
      freq.push_back(temp_y);
    }
    m_R = R;
    m_freq = freq;
  }
  else
  {
    NCRYSTAL_THROW2(BadInput, rel_path << ": Invalid data file for the " << pluginNameUpperCase()
                                       << " plugin");
  }
};

NCP::PhysicsModel::PhysicsModel(std::string filename)
    : m_model("file"),
      m_helper(([filename]() -> NC::IofQHelper
                {
      
      //Parse the input file and create the vector with the q, Iq info
      NC::VectD q;
      NC::VectD IofQ;
      std::string root_rel = "data/";
      std::string rel_path = root_rel+filename;
      std::ifstream input_file(rel_path);
      if(input_file) {
        double temp_x, temp_y;
        std::string line;
        while (std::getline(input_file, line)){
              std::istringstream iss(line);
              iss >> temp_x >> temp_y;
              q.push_back(temp_x);
              IofQ.push_back(temp_y);
        }
      } else {
        NCRYSTAL_THROW2( BadInput,rel_path << ": Invalid data file for the "<<pluginNameUpperCase()
                        <<" plugin" );   
      }
      /*for(auto i :q){
        std::cout<<i<<std::endl;;
      }*/ 
      NC::IofQHelper helper(q,IofQ);
      return helper; })())
      {

      };

double NCP::PhysicsModel::calcCrossSection(double neutron_ekin) const
{
  NC::NeutronEnergy ekin(neutron_ekin);
  double k = NC::k2Pi / NC::ekin2wl(neutron_ekin); // wavevector
  double SANS_xs;
  if (m_model == "simple")
  {
    if (m_p.has_value()) SANS_xs = m_p.value() * 25.3e-3 / neutron_ekin;
  }
  else if (m_model == "NP-FBA")
  {
    // h_bar = 1.054571817e-34 ;//[J*s] <- Reduced planck constant
    // m = 1.674927498e-27;// [Kg] <- neutron mass
    // V = 290.0e-9*1.60218e-19;// [J]
    // 2*pi*(2mV/h_bar^2)^2 =   1.2307e+33 [1/m^4]
    double physical_constant = 1.2307e+33;
    SANS_xs = 0;
    double R, freq, _2kr, I;
    if (m_R.has_value() && m_freq.has_value()){
      for (size_t i = 0; i < m_R.value().size(); ++i)
      {
        R = m_R.value().at(i) * 10; // converto to AA
        freq = m_freq.value().at(i);
        _2kr = 2 * k * R;
        I = 0.25 * (1 - 1 / (_2kr * _2kr) + sin(2 * _2kr) / (_2kr * _2kr * _2kr) - sin(_2kr) * sin(_2kr) / (_2kr * _2kr * _2kr * _2kr));
        SANS_xs += freq * std::pow(R * 1e-10, 6) / (k * k * R * R) * I * 1e+28; //[barn]
      }
      SANS_xs *= physical_constant;
    } else if (m_p.has_value()){
        R = m_p.value() * 10; // converto to AA
        _2kr = 2 * k * R;
        I = 0.25 * (1 - 1 / (_2kr * _2kr) + sin(2 * _2kr) / (_2kr * _2kr * _2kr) - sin(_2kr) * sin(_2kr) / (_2kr * _2kr * _2kr * _2kr));
        SANS_xs = physical_constant*std::pow(R * 1e-10, 6) / (k * k * R * R) * I * 1e+28; //[barn]
    }
  }
  else
  {
    if (m_helper.has_value())
    {
      SANS_xs = 1 / (2 * k * k) * m_helper.value().calcQIofQIntegral(ekin);
    }
    else
    {
      NCRYSTAL_THROW2(BadInput, "Attempt to use IofQHelper when not defined in the " << pluginNameUpperCase()
                                                                                     << " plugin");
    }
  }
  return SANS_xs;
}

double NCP::PhysicsModel::sampleScatteringVector(NC::RNG &rng, double neutron_ekin) const
{
  NC::NeutronEnergy ekin(neutron_ekin);
  double Q;
  if (m_model == "simple" || m_model == "NP-FBA")
  {
    double rand = rng.generate();
    double k = NC::k2Pi / NC::ekin2wl(neutron_ekin); // wavevector
    // sample a random scattering vector Q from the inverse PPF CDF (see plugin readme)
    double A1 = 132.869;
    double b1 = 1.33605;
    double A2 = 0.0519763;
    double b2 = 3.97314;
    double Q0 = 0.0510821;
    double ratio_sigma = (1 / (2 * k * k)) / calcCrossSection(neutron_ekin); // cross section over total cross section ratio
    double CDF_Q0 = (A1 * std::pow(Q0, -b1 + 2) / (-b1 + 2)) * ratio_sigma;
    if (rand < CDF_Q0)
    {
      Q = std::pow(((-b1 + 2) * rand / A1) / ratio_sigma, 1 / (-b1 + 2));
    }
    else
    {
      Q = std::pow((rand / ratio_sigma - (A1 / (-b1 + 2)) * std::pow(Q0, -b1 + 2) + (A2 / (-b2 + 2)) * std::pow(Q0, -b2 + 2)) * (-b2 + 2) / A2, 1 / (-b2 + 2));
    }
  }
  else if (m_helper.has_value())
  {
    Q = m_helper.value().sampleQValue(rng, ekin);
  }
  else
  {
    NCRYSTAL_THROW2(BadInput, "Attempt to use IofQHelper when not defined in the " << pluginNameUpperCase()
                                                                                   << " plugin");
  }
  return Q;
}
NCP::PhysicsModel::ScatEvent NCP::PhysicsModel::sampleScatteringEvent(NC::RNG &rng, double neutron_ekin) const
{
  ScatEvent result;

  // section for energy limits
  /*if ( ! (neutron_ekin > 1) ) {
    result.ekin_final = neutron_ekin;
    result.mu = 1.0;
    return result;
  }*/

  // Implement our actual model here:
  result.ekin_final = neutron_ekin; // Elastic
  double Q = sampleScatteringVector(rng, neutron_ekin);
  double ksquared = NC::k4PiSq * NC::ekin2wlsqinv(neutron_ekin);
  result.mu = 1 - 0.5 * (Q * Q / ksquared);

  return result;
}
