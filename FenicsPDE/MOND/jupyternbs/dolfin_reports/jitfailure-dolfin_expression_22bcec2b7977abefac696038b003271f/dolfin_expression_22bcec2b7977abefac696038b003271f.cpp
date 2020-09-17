
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
  class dolfin_expression_22bcec2b7977abefac696038b003271f : public Expression
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
double position_x0;
double position_y0;
double position_z0;
double position_x1;
double position_y1;
double position_z1;
double source_number;
double source_mass;


       dolfin_expression_22bcec2b7977abefac696038b003271f()
       {
            
       }

       void eval(Eigen::Ref<Eigen::VectorXd> values, Eigen::Ref<const Eigen::VectorXd> x) const override
       {
          values[0] = (((pow((x[0]-position_x_0),2)+pow((x[1]-position_y_0),2)+pow((x[2]-position_z_0),2)) <=pow(radius_tot, 2)) ||((pow((x[0]-position_x_1),2)+pow((x[1]-position_y_1),2)+pow((x[2]-position_z_1),2)) <=pow(radius_tot, 2)) ||((pow((x[0]-position_x_2),2)+pow((x[1]-position_y_2),2)+pow((x[2]-position_z_2),2)) <=pow(radius_tot, 2)) ||((pow((x[0]-position_x_3),2)+pow((x[1]-position_y_3),2)+pow((x[2]-position_z_3),2)) <=pow(radius_tot, 2)) ||((pow((x[0]-position_x_4),2)+pow((x[1]-position_y_4),2)+pow((x[2]-position_z_4),2)) <=pow(radius_tot, 2)))? a0*4*pi*G*source_mass/volume_out : 0;

       }

       void set_property(std::string name, double _value) override
       {
          if (name == "a0") { a0 = _value; return; }          if (name == "ms") { ms = _value; return; }          if (name == "mgb") { mgb = _value; return; }          if (name == "G") { G = _value; return; }          if (name == "ly") { ly = _value; return; }          if (name == "kp") { kp = _value; return; }          if (name == "radius_tot") { radius_tot = _value; return; }          if (name == "volume_out") { volume_out = _value; return; }          if (name == "stand_dev") { stand_dev = _value; return; }          if (name == "p") { p = _value; return; }          if (name == "position_x0") { position_x0 = _value; return; }          if (name == "position_y0") { position_y0 = _value; return; }          if (name == "position_z0") { position_z0 = _value; return; }          if (name == "position_x1") { position_x1 = _value; return; }          if (name == "position_y1") { position_y1 = _value; return; }          if (name == "position_z1") { position_z1 = _value; return; }          if (name == "source_number") { source_number = _value; return; }          if (name == "source_mass") { source_mass = _value; return; }
       throw std::runtime_error("No such property");
       }

       double get_property(std::string name) const override
       {
          if (name == "a0") return a0;          if (name == "ms") return ms;          if (name == "mgb") return mgb;          if (name == "G") return G;          if (name == "ly") return ly;          if (name == "kp") return kp;          if (name == "radius_tot") return radius_tot;          if (name == "volume_out") return volume_out;          if (name == "stand_dev") return stand_dev;          if (name == "p") return p;          if (name == "position_x0") return position_x0;          if (name == "position_y0") return position_y0;          if (name == "position_z0") return position_z0;          if (name == "position_x1") return position_x1;          if (name == "position_y1") return position_y1;          if (name == "position_z1") return position_z1;          if (name == "source_number") return source_number;          if (name == "source_mass") return source_mass;
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

extern "C" DLL_EXPORT dolfin::Expression * create_dolfin_expression_22bcec2b7977abefac696038b003271f()
{
  return new dolfin::dolfin_expression_22bcec2b7977abefac696038b003271f;
}

