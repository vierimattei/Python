
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
  class dolfin_expression_875b5d46befb4f814cfcf7a25e11518d : public Expression
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
double x_close;
double y_close;
double z_close;
double radius_close;
double eps;
double marmar;


       dolfin_expression_875b5d46befb4f814cfcf7a25e11518d()
       {
            
       }

       void eval(Eigen::Ref<Eigen::VectorXd> values, Eigen::Ref<const Eigen::VectorXd> x) const override
       {
          values[0] = a0*4*pi*G*mgb/pow((marmar*sqrt(2*pi)),3)*pow(-(1/2)*((pow(x[0],2) + pow(x[1],2)+ pow(x[2],2))/pow(marmar,2)));

       }

       void set_property(std::string name, double _value) override
       {
          if (name == "a0") { a0 = _value; return; }          if (name == "ms") { ms = _value; return; }          if (name == "mgb") { mgb = _value; return; }          if (name == "G") { G = _value; return; }          if (name == "ly") { ly = _value; return; }          if (name == "kp") { kp = _value; return; }          if (name == "radius_tot") { radius_tot = _value; return; }          if (name == "volume_out") { volume_out = _value; return; }          if (name == "x_close") { x_close = _value; return; }          if (name == "y_close") { y_close = _value; return; }          if (name == "z_close") { z_close = _value; return; }          if (name == "radius_close") { radius_close = _value; return; }          if (name == "eps") { eps = _value; return; }          if (name == "marmar") { marmar = _value; return; }
       throw std::runtime_error("No such property");
       }

       double get_property(std::string name) const override
       {
          if (name == "a0") return a0;          if (name == "ms") return ms;          if (name == "mgb") return mgb;          if (name == "G") return G;          if (name == "ly") return ly;          if (name == "kp") return kp;          if (name == "radius_tot") return radius_tot;          if (name == "volume_out") return volume_out;          if (name == "x_close") return x_close;          if (name == "y_close") return y_close;          if (name == "z_close") return z_close;          if (name == "radius_close") return radius_close;          if (name == "eps") return eps;          if (name == "marmar") return marmar;
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

extern "C" DLL_EXPORT dolfin::Expression * create_dolfin_expression_875b5d46befb4f814cfcf7a25e11518d()
{
  return new dolfin::dolfin_expression_875b5d46befb4f814cfcf7a25e11518d;
}

