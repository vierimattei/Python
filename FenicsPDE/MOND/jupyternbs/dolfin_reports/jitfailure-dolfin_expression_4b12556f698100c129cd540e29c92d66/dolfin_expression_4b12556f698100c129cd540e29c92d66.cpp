
// Based on https://gcc.gnu.org/wiki/Visibility
#if defined _WIN32 || defined __CYGWIN__
    #ifdef __GNUC__
        #define DLL_EXPORT __attribute__ ((dllexport))
    #else
        #define DLL_EXPORT __declspec(dllexport)
    #endif
#else
    #define DLL_EXPORT __attribute__ ((visibility ("default")))
#endif

#include <dolfin/function/Expression.h>
#include <dolfin/math/basic.h>
#include <Eigen/Dense>


// cmath functions
using std::cos;
using std::sin;
using std::tan;
using std::acos;
using std::asin;
using std::atan;
using std::atan2;
using std::cosh;
using std::sinh;
using std::tanh;
using std::exp;
using std::frexp;
using std::ldexp;
using std::log;
using std::log10;
using std::modf;
using std::pow;
using std::sqrt;
using std::ceil;
using std::fabs;
using std::floor;
using std::fmod;
using std::max;
using std::min;

const double pi = DOLFIN_PI;


namespace dolfin
{
  class dolfin_expression_4b12556f698100c129cd540e29c92d66 : public Expression
  {
     public:
       double a0;
double ms;
double mgb;
double G;
double ly;
double kp;
double radius_tot;
double volume_out;
double stand_dev;
double p;
double source_number;
double source_mass;


       dolfin_expression_4b12556f698100c129cd540e29c92d66()
       {
            
       }

       void eval(Eigen::Ref<Eigen::VectorXd> values, Eigen::Ref<const Eigen::VectorXd> x) const override
       {
          values[0] = (((pow((x[0]-position_x[0]),2)+pow((x[1]-position_y[0]),2)+pow((x[2]-position_z[0]),2)) <=pow(radius_tot, 2)) ||((pow((x[0]-position_x[1]),2)+pow((x[1]-position_y[1]),2)+pow((x[2]-position_z[1]),2)) <=pow(radius_tot, 2)) ||((pow((x[0]-position_x[2]),2)+pow((x[1]-position_y[2]),2)+pow((x[2]-position_z[2]),2)) <=pow(radius_tot, 2)) ||((pow((x[0]-position_x[3]),2)+pow((x[1]-position_y[3]),2)+pow((x[2]-position_z[3]),2)) <=pow(radius_tot, 2)) ||((pow((x[0]-position_x[4]),2)+pow((x[1]-position_y[4]),2)+pow((x[2]-position_z[4]),2)) <=pow(radius_tot, 2)))? a0*4*pi*G*source_mass/volume_out : 0;

       }

       void set_property(std::string name, double _value) override
       {
          if (name == "a0") { a0 = _value; return; }          if (name == "ms") { ms = _value; return; }          if (name == "mgb") { mgb = _value; return; }          if (name == "G") { G = _value; return; }          if (name == "ly") { ly = _value; return; }          if (name == "kp") { kp = _value; return; }          if (name == "radius_tot") { radius_tot = _value; return; }          if (name == "volume_out") { volume_out = _value; return; }          if (name == "stand_dev") { stand_dev = _value; return; }          if (name == "p") { p = _value; return; }          if (name == "source_number") { source_number = _value; return; }          if (name == "source_mass") { source_mass = _value; return; }
       throw std::runtime_error("No such property");
       }

       double get_property(std::string name) const override
       {
          if (name == "a0") return a0;          if (name == "ms") return ms;          if (name == "mgb") return mgb;          if (name == "G") return G;          if (name == "ly") return ly;          if (name == "kp") return kp;          if (name == "radius_tot") return radius_tot;          if (name == "volume_out") return volume_out;          if (name == "stand_dev") return stand_dev;          if (name == "p") return p;          if (name == "source_number") return source_number;          if (name == "source_mass") return source_mass;
       throw std::runtime_error("No such property");
       return 0.0;
       }

       void set_generic_function(std::string name, std::shared_ptr<dolfin::GenericFunction> _value) override
       {

       throw std::runtime_error("No such property");
       }

       std::shared_ptr<dolfin::GenericFunction> get_generic_function(std::string name) const override
       {

       throw std::runtime_error("No such property");
       }

  };
}

extern "C" DLL_EXPORT dolfin::Expression * create_dolfin_expression_4b12556f698100c129cd540e29c92d66()
{
  return new dolfin::dolfin_expression_4b12556f698100c129cd540e29c92d66;
}

