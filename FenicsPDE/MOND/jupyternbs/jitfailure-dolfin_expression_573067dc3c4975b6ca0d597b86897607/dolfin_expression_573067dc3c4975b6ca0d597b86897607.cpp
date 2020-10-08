
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
  class dolfin_expression_573067dc3c4975b6ca0d597b86897607 : public Expression
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
double position_x_0;
double position_y_0;
double position_z_0;
double position_x_1;
double position_y_1;
double position_z_1;
double position_x_2;
double position_y_2;
double position_z_2;
double position_x_3;
double position_y_3;
double position_z_3;
double position_x_4;
double position_y_4;
double position_z_4;
double position_x_5;
double position_y_5;
double position_z_5;
double position_x_6;
double position_y_6;
double position_z_6;
double position_x_7;
double position_y_7;
double position_z_7;
double position_x_8;
double position_y_8;
double position_z_8;
double position_x_9;
double position_y_9;
double position_z_9;
double position_x_10;
double position_y_10;
double position_z_10;
double position_x_11;
double position_y_11;
double position_z_11;
double position_x_12;
double position_y_12;
double position_z_12;
double position_x_13;
double position_y_13;
double position_z_13;
double position_x_14;
double position_y_14;
double position_z_14;
double position_x_15;
double position_y_15;
double position_z_15;
double position_x_16;
double position_y_16;
double position_z_16;
double position_x_17;
double position_y_17;
double position_z_17;
double position_x_18;
double position_y_18;
double position_z_18;
double position_x_19;
double position_y_19;
double position_z_19;
double position_x_20;
double position_y_20;
double position_z_20;
double position_x_21;
double position_y_21;
double position_z_21;
double position_x_22;
double position_y_22;
double position_z_22;
double position_x_23;
double position_y_23;
double position_z_23;
double position_x_24;
double position_y_24;
double position_z_24;
double position_x_25;
double position_y_25;
double position_z_25;
double position_x_26;
double position_y_26;
double position_z_26;
double position_x_27;
double position_y_27;
double position_z_27;
double position_x_28;
double position_y_28;
double position_z_28;
double position_x_29;
double position_y_29;
double position_z_29;
double position_x_30;
double position_y_30;
double position_z_30;
double position_x_31;
double position_y_31;
double position_z_31;
double position_x_32;
double position_y_32;
double position_z_32;
double position_x_33;
double position_y_33;
double position_z_33;
double position_x_34;
double position_y_34;
double position_z_34;
double position_x_35;
double position_y_35;
double position_z_35;
double position_x_36;
double position_y_36;
double position_z_36;
double position_x_37;
double position_y_37;
double position_z_37;
double position_x_38;
double position_y_38;
double position_z_38;
double position_x_39;
double position_y_39;
double position_z_39;
double position_x_40;
double position_y_40;
double position_z_40;
double position_x_41;
double position_y_41;
double position_z_41;
double position_x_42;
double position_y_42;
double position_z_42;
double position_x_43;
double position_y_43;
double position_z_43;
double position_x_44;
double position_y_44;
double position_z_44;
double position_x_45;
double position_y_45;
double position_z_45;
double position_x_46;
double position_y_46;
double position_z_46;
double position_x_47;
double position_y_47;
double position_z_47;
double position_x_48;
double position_y_48;
double position_z_48;
double position_x_49;
double position_y_49;
double position_z_49;


       dolfin_expression_573067dc3c4975b6ca0d597b86897607()
       {
            
       }

       void eval(Eigen::Ref<Eigen::VectorXd> values, Eigen::Ref<const Eigen::VectorXd> x) const override
       {
          values[0] = (((pow((x[0]-position_x_0),2)+pow((x[1]-position_y_0),2)+pow((x[2]-position_z_0),2)) <=pow(radius_tot, 2)) ||((pow((x[0]-position_x_1),2)+pow((x[1]-position_y_1),2)+pow((x[2]-position_z_1),2)) <=pow(radius_tot, 2)) ||((pow((x[0]-position_x_2),2)+pow((x[1]-position_y_2),2)+pow((x[2]-position_z_2),2)) <=pow(radius_tot, 2)) ||((pow((x[0]-position_x_3),2)+pow((x[1]-position_y_3),2)+pow((x[2]-position_z_3),2)) <=pow(radius_tot, 2)) ||((pow((x[0]-position_x_4),2)+pow((x[1]-position_y_4),2)+pow((x[2]-position_z_4),2)) <=pow(radius_tot, 2)) ||((pow((x[0]-position_x_5),2)+pow((x[1]-position_y_5),2)+pow((x[2]-position_z_5),2)) <=pow(radius_tot, 2)) ||((pow((x[0]-position_x_6),2)+pow((x[1]-position_y_6),2)+pow((x[2]-position_z_6),2)) <=pow(radius_tot, 2)) ||((pow((x[0]-position_x_7),2)+pow((x[1]-position_y_7),2)+pow((x[2]-position_z_7),2)) <=pow(radius_tot, 2)) ||((pow((x[0]-position_x_8),2)+pow((x[1]-position_y_8),2)+pow((x[2]-position_z_8),2)) <=pow(radius_tot, 2)) ||((pow((x[0]-position_x_9),2)+pow((x[1]-position_y_9),2)+pow((x[2]-position_z_9),2)) <=pow(radius_tot, 2)) ||((pow((x[0]-position_x_10),2)+pow((x[1]-position_y_10),2)+pow((x[2]-position_z_10),2)) <=pow(radius_tot, 2)) ||((pow((x[0]-position_x_11),2)+pow((x[1]-position_y_11),2)+pow((x[2]-position_z_11),2)) <=pow(radius_tot, 2)) ||((pow((x[0]-position_x_12),2)+pow((x[1]-position_y_12),2)+pow((x[2]-position_z_12),2)) <=pow(radius_tot, 2)) ||((pow((x[0]-position_x_13),2)+pow((x[1]-position_y_13),2)+pow((x[2]-position_z_13),2)) <=pow(radius_tot, 2)) ||((pow((x[0]-position_x_14),2)+pow((x[1]-position_y_14),2)+pow((x[2]-position_z_14),2)) <=pow(radius_tot, 2)) ||((pow((x[0]-position_x_15),2)+pow((x[1]-position_y_15),2)+pow((x[2]-position_z_15),2)) <=pow(radius_tot, 2)) ||((pow((x[0]-position_x_16),2)+pow((x[1]-position_y_16),2)+pow((x[2]-position_z_16),2)) <=pow(radius_tot, 2)) ||((pow((x[0]-position_x_17),2)+pow((x[1]-position_y_17),2)+pow((x[2]-position_z_17),2)) <=pow(radius_tot, 2)) ||((pow((x[0]-position_x_18),2)+pow((x[1]-position_y_18),2)+pow((x[2]-position_z_18),2)) <=pow(radius_tot, 2)) ||((pow((x[0]-position_x_19),2)+pow((x[1]-position_y_19),2)+pow((x[2]-position_z_19),2)) <=pow(radius_tot, 2)) ||((pow((x[0]-position_x_20),2)+pow((x[1]-position_y_20),2)+pow((x[2]-position_z_20),2)) <=pow(radius_tot, 2)) ||((pow((x[0]-position_x_21),2)+pow((x[1]-position_y_21),2)+pow((x[2]-position_z_21),2)) <=pow(radius_tot, 2)) ||((pow((x[0]-position_x_22),2)+pow((x[1]-position_y_22),2)+pow((x[2]-position_z_22),2)) <=pow(radius_tot, 2)) ||((pow((x[0]-position_x_23),2)+pow((x[1]-position_y_23),2)+pow((x[2]-position_z_23),2)) <=pow(radius_tot, 2)) ||((pow((x[0]-position_x_24),2)+pow((x[1]-position_y_24),2)+pow((x[2]-position_z_24),2)) <=pow(radius_tot, 2)) ||((pow((x[0]-position_x_25),2)+pow((x[1]-position_y_25),2)+pow((x[2]-position_z_25),2)) <=pow(radius_tot, 2)) ||((pow((x[0]-position_x_26),2)+pow((x[1]-position_y_26),2)+pow((x[2]-position_z_26),2)) <=pow(radius_tot, 2)) ||((pow((x[0]-position_x_27),2)+pow((x[1]-position_y_27),2)+pow((x[2]-position_z_27),2)) <=pow(radius_tot, 2)) ||((pow((x[0]-position_x_28),2)+pow((x[1]-position_y_28),2)+pow((x[2]-position_z_28),2)) <=pow(radius_tot, 2)) ||((pow((x[0]-position_x_29),2)+pow((x[1]-position_y_29),2)+pow((x[2]-position_z_29),2)) <=pow(radius_tot, 2)) ||((pow((x[0]-position_x_30),2)+pow((x[1]-position_y_30),2)+pow((x[2]-position_z_30),2)) <=pow(radius_tot, 2)) ||((pow((x[0]-position_x_31),2)+pow((x[1]-position_y_31),2)+pow((x[2]-position_z_31),2)) <=pow(radius_tot, 2)) ||((pow((x[0]-position_x_32),2)+pow((x[1]-position_y_32),2)+pow((x[2]-position_z_32),2)) <=pow(radius_tot, 2)) ||((pow((x[0]-position_x_33),2)+pow((x[1]-position_y_33),2)+pow((x[2]-position_z_33),2)) <=pow(radius_tot, 2)) ||((pow((x[0]-position_x_34),2)+pow((x[1]-position_y_34),2)+pow((x[2]-position_z_34),2)) <=pow(radius_tot, 2)) ||((pow((x[0]-position_x_35),2)+pow((x[1]-position_y_35),2)+pow((x[2]-position_z_35),2)) <=pow(radius_tot, 2)) ||((pow((x[0]-position_x_36),2)+pow((x[1]-position_y_36),2)+pow((x[2]-position_z_36),2)) <=pow(radius_tot, 2)) ||((pow((x[0]-position_x_37),2)+pow((x[1]-position_y_37),2)+pow((x[2]-position_z_37),2)) <=pow(radius_tot, 2)) ||((pow((x[0]-position_x_38),2)+pow((x[1]-position_y_38),2)+pow((x[2]-position_z_38),2)) <=pow(radius_tot, 2)) ||((pow((x[0]-position_x_39),2)+pow((x[1]-position_y_39),2)+pow((x[2]-position_z_39),2)) <=pow(radius_tot, 2)) ||((pow((x[0]-position_x_40),2)+pow((x[1]-position_y_40),2)+pow((x[2]-position_z_40),2)) <=pow(radius_tot, 2)) ||((pow((x[0]-position_x_41),2)+pow((x[1]-position_y_41),2)+pow((x[2]-position_z_41),2)) <=pow(radius_tot, 2)) ||((pow((x[0]-position_x_42),2)+pow((x[1]-position_y_42),2)+pow((x[2]-position_z_42),2)) <=pow(radius_tot, 2)) ||((pow((x[0]-position_x_43),2)+pow((x[1]-position_y_43),2)+pow((x[2]-position_z_43),2)) <=pow(radius_tot, 2)) ||((pow((x[0]-position_x_44),2)+pow((x[1]-position_y_44),2)+pow((x[2]-position_z_44),2)) <=pow(radius_tot, 2)) ||((pow((x[0]-position_x_45),2)+pow((x[1]-position_y_45),2)+pow((x[2]-position_z_45),2)) <=pow(radius_tot, 2)) ||((pow((x[0]-position_x_46),2)+pow((x[1]-position_y_46),2)+pow((x[2]-position_z_46),2)) <=pow(radius_tot, 2)) ||((pow((x[0]-position_x_47),2)+pow((x[1]-position_y_47),2)+pow((x[2]-position_z_47),2)) <=pow(radius_tot, 2)) ||((pow((x[0]-position_x_48),2)+pow((x[1]-position_y_48),2)+pow((x[2]-position_z_48),2)) <=pow(radius_tot, 2)) ||((pow((x[0]-position_x_49),2)+pow((x[1]-position_y_49),2)+pow((x[2]-position_z_49),2)) <=pow(radius_tot, 2)))? a*4*pi*G*source_mass/volume_out : 0;

       }

       void set_property(std::string name, double _value) override
       {
          if (name == "a0") { a0 = _value; return; }          if (name == "ms") { ms = _value; return; }          if (name == "mgb") { mgb = _value; return; }          if (name == "G") { G = _value; return; }          if (name == "ly") { ly = _value; return; }          if (name == "kp") { kp = _value; return; }          if (name == "radius_tot") { radius_tot = _value; return; }          if (name == "volume_out") { volume_out = _value; return; }          if (name == "stand_dev") { stand_dev = _value; return; }          if (name == "p") { p = _value; return; }          if (name == "source_number") { source_number = _value; return; }          if (name == "source_mass") { source_mass = _value; return; }          if (name == "position_x_0") { position_x_0 = _value; return; }          if (name == "position_y_0") { position_y_0 = _value; return; }          if (name == "position_z_0") { position_z_0 = _value; return; }          if (name == "position_x_1") { position_x_1 = _value; return; }          if (name == "position_y_1") { position_y_1 = _value; return; }          if (name == "position_z_1") { position_z_1 = _value; return; }          if (name == "position_x_2") { position_x_2 = _value; return; }          if (name == "position_y_2") { position_y_2 = _value; return; }          if (name == "position_z_2") { position_z_2 = _value; return; }          if (name == "position_x_3") { position_x_3 = _value; return; }          if (name == "position_y_3") { position_y_3 = _value; return; }          if (name == "position_z_3") { position_z_3 = _value; return; }          if (name == "position_x_4") { position_x_4 = _value; return; }          if (name == "position_y_4") { position_y_4 = _value; return; }          if (name == "position_z_4") { position_z_4 = _value; return; }          if (name == "position_x_5") { position_x_5 = _value; return; }          if (name == "position_y_5") { position_y_5 = _value; return; }          if (name == "position_z_5") { position_z_5 = _value; return; }          if (name == "position_x_6") { position_x_6 = _value; return; }          if (name == "position_y_6") { position_y_6 = _value; return; }          if (name == "position_z_6") { position_z_6 = _value; return; }          if (name == "position_x_7") { position_x_7 = _value; return; }          if (name == "position_y_7") { position_y_7 = _value; return; }          if (name == "position_z_7") { position_z_7 = _value; return; }          if (name == "position_x_8") { position_x_8 = _value; return; }          if (name == "position_y_8") { position_y_8 = _value; return; }          if (name == "position_z_8") { position_z_8 = _value; return; }          if (name == "position_x_9") { position_x_9 = _value; return; }          if (name == "position_y_9") { position_y_9 = _value; return; }          if (name == "position_z_9") { position_z_9 = _value; return; }          if (name == "position_x_10") { position_x_10 = _value; return; }          if (name == "position_y_10") { position_y_10 = _value; return; }          if (name == "position_z_10") { position_z_10 = _value; return; }          if (name == "position_x_11") { position_x_11 = _value; return; }          if (name == "position_y_11") { position_y_11 = _value; return; }          if (name == "position_z_11") { position_z_11 = _value; return; }          if (name == "position_x_12") { position_x_12 = _value; return; }          if (name == "position_y_12") { position_y_12 = _value; return; }          if (name == "position_z_12") { position_z_12 = _value; return; }          if (name == "position_x_13") { position_x_13 = _value; return; }          if (name == "position_y_13") { position_y_13 = _value; return; }          if (name == "position_z_13") { position_z_13 = _value; return; }          if (name == "position_x_14") { position_x_14 = _value; return; }          if (name == "position_y_14") { position_y_14 = _value; return; }          if (name == "position_z_14") { position_z_14 = _value; return; }          if (name == "position_x_15") { position_x_15 = _value; return; }          if (name == "position_y_15") { position_y_15 = _value; return; }          if (name == "position_z_15") { position_z_15 = _value; return; }          if (name == "position_x_16") { position_x_16 = _value; return; }          if (name == "position_y_16") { position_y_16 = _value; return; }          if (name == "position_z_16") { position_z_16 = _value; return; }          if (name == "position_x_17") { position_x_17 = _value; return; }          if (name == "position_y_17") { position_y_17 = _value; return; }          if (name == "position_z_17") { position_z_17 = _value; return; }          if (name == "position_x_18") { position_x_18 = _value; return; }          if (name == "position_y_18") { position_y_18 = _value; return; }          if (name == "position_z_18") { position_z_18 = _value; return; }          if (name == "position_x_19") { position_x_19 = _value; return; }          if (name == "position_y_19") { position_y_19 = _value; return; }          if (name == "position_z_19") { position_z_19 = _value; return; }          if (name == "position_x_20") { position_x_20 = _value; return; }          if (name == "position_y_20") { position_y_20 = _value; return; }          if (name == "position_z_20") { position_z_20 = _value; return; }          if (name == "position_x_21") { position_x_21 = _value; return; }          if (name == "position_y_21") { position_y_21 = _value; return; }          if (name == "position_z_21") { position_z_21 = _value; return; }          if (name == "position_x_22") { position_x_22 = _value; return; }          if (name == "position_y_22") { position_y_22 = _value; return; }          if (name == "position_z_22") { position_z_22 = _value; return; }          if (name == "position_x_23") { position_x_23 = _value; return; }          if (name == "position_y_23") { position_y_23 = _value; return; }          if (name == "position_z_23") { position_z_23 = _value; return; }          if (name == "position_x_24") { position_x_24 = _value; return; }          if (name == "position_y_24") { position_y_24 = _value; return; }          if (name == "position_z_24") { position_z_24 = _value; return; }          if (name == "position_x_25") { position_x_25 = _value; return; }          if (name == "position_y_25") { position_y_25 = _value; return; }          if (name == "position_z_25") { position_z_25 = _value; return; }          if (name == "position_x_26") { position_x_26 = _value; return; }          if (name == "position_y_26") { position_y_26 = _value; return; }          if (name == "position_z_26") { position_z_26 = _value; return; }          if (name == "position_x_27") { position_x_27 = _value; return; }          if (name == "position_y_27") { position_y_27 = _value; return; }          if (name == "position_z_27") { position_z_27 = _value; return; }          if (name == "position_x_28") { position_x_28 = _value; return; }          if (name == "position_y_28") { position_y_28 = _value; return; }          if (name == "position_z_28") { position_z_28 = _value; return; }          if (name == "position_x_29") { position_x_29 = _value; return; }          if (name == "position_y_29") { position_y_29 = _value; return; }          if (name == "position_z_29") { position_z_29 = _value; return; }          if (name == "position_x_30") { position_x_30 = _value; return; }          if (name == "position_y_30") { position_y_30 = _value; return; }          if (name == "position_z_30") { position_z_30 = _value; return; }          if (name == "position_x_31") { position_x_31 = _value; return; }          if (name == "position_y_31") { position_y_31 = _value; return; }          if (name == "position_z_31") { position_z_31 = _value; return; }          if (name == "position_x_32") { position_x_32 = _value; return; }          if (name == "position_y_32") { position_y_32 = _value; return; }          if (name == "position_z_32") { position_z_32 = _value; return; }          if (name == "position_x_33") { position_x_33 = _value; return; }          if (name == "position_y_33") { position_y_33 = _value; return; }          if (name == "position_z_33") { position_z_33 = _value; return; }          if (name == "position_x_34") { position_x_34 = _value; return; }          if (name == "position_y_34") { position_y_34 = _value; return; }          if (name == "position_z_34") { position_z_34 = _value; return; }          if (name == "position_x_35") { position_x_35 = _value; return; }          if (name == "position_y_35") { position_y_35 = _value; return; }          if (name == "position_z_35") { position_z_35 = _value; return; }          if (name == "position_x_36") { position_x_36 = _value; return; }          if (name == "position_y_36") { position_y_36 = _value; return; }          if (name == "position_z_36") { position_z_36 = _value; return; }          if (name == "position_x_37") { position_x_37 = _value; return; }          if (name == "position_y_37") { position_y_37 = _value; return; }          if (name == "position_z_37") { position_z_37 = _value; return; }          if (name == "position_x_38") { position_x_38 = _value; return; }          if (name == "position_y_38") { position_y_38 = _value; return; }          if (name == "position_z_38") { position_z_38 = _value; return; }          if (name == "position_x_39") { position_x_39 = _value; return; }          if (name == "position_y_39") { position_y_39 = _value; return; }          if (name == "position_z_39") { position_z_39 = _value; return; }          if (name == "position_x_40") { position_x_40 = _value; return; }          if (name == "position_y_40") { position_y_40 = _value; return; }          if (name == "position_z_40") { position_z_40 = _value; return; }          if (name == "position_x_41") { position_x_41 = _value; return; }          if (name == "position_y_41") { position_y_41 = _value; return; }          if (name == "position_z_41") { position_z_41 = _value; return; }          if (name == "position_x_42") { position_x_42 = _value; return; }          if (name == "position_y_42") { position_y_42 = _value; return; }          if (name == "position_z_42") { position_z_42 = _value; return; }          if (name == "position_x_43") { position_x_43 = _value; return; }          if (name == "position_y_43") { position_y_43 = _value; return; }          if (name == "position_z_43") { position_z_43 = _value; return; }          if (name == "position_x_44") { position_x_44 = _value; return; }          if (name == "position_y_44") { position_y_44 = _value; return; }          if (name == "position_z_44") { position_z_44 = _value; return; }          if (name == "position_x_45") { position_x_45 = _value; return; }          if (name == "position_y_45") { position_y_45 = _value; return; }          if (name == "position_z_45") { position_z_45 = _value; return; }          if (name == "position_x_46") { position_x_46 = _value; return; }          if (name == "position_y_46") { position_y_46 = _value; return; }          if (name == "position_z_46") { position_z_46 = _value; return; }          if (name == "position_x_47") { position_x_47 = _value; return; }          if (name == "position_y_47") { position_y_47 = _value; return; }          if (name == "position_z_47") { position_z_47 = _value; return; }          if (name == "position_x_48") { position_x_48 = _value; return; }          if (name == "position_y_48") { position_y_48 = _value; return; }          if (name == "position_z_48") { position_z_48 = _value; return; }          if (name == "position_x_49") { position_x_49 = _value; return; }          if (name == "position_y_49") { position_y_49 = _value; return; }          if (name == "position_z_49") { position_z_49 = _value; return; }
       throw std::runtime_error("No such property");
       }

       double get_property(std::string name) const override
       {
          if (name == "a0") return a0;          if (name == "ms") return ms;          if (name == "mgb") return mgb;          if (name == "G") return G;          if (name == "ly") return ly;          if (name == "kp") return kp;          if (name == "radius_tot") return radius_tot;          if (name == "volume_out") return volume_out;          if (name == "stand_dev") return stand_dev;          if (name == "p") return p;          if (name == "source_number") return source_number;          if (name == "source_mass") return source_mass;          if (name == "position_x_0") return position_x_0;          if (name == "position_y_0") return position_y_0;          if (name == "position_z_0") return position_z_0;          if (name == "position_x_1") return position_x_1;          if (name == "position_y_1") return position_y_1;          if (name == "position_z_1") return position_z_1;          if (name == "position_x_2") return position_x_2;          if (name == "position_y_2") return position_y_2;          if (name == "position_z_2") return position_z_2;          if (name == "position_x_3") return position_x_3;          if (name == "position_y_3") return position_y_3;          if (name == "position_z_3") return position_z_3;          if (name == "position_x_4") return position_x_4;          if (name == "position_y_4") return position_y_4;          if (name == "position_z_4") return position_z_4;          if (name == "position_x_5") return position_x_5;          if (name == "position_y_5") return position_y_5;          if (name == "position_z_5") return position_z_5;          if (name == "position_x_6") return position_x_6;          if (name == "position_y_6") return position_y_6;          if (name == "position_z_6") return position_z_6;          if (name == "position_x_7") return position_x_7;          if (name == "position_y_7") return position_y_7;          if (name == "position_z_7") return position_z_7;          if (name == "position_x_8") return position_x_8;          if (name == "position_y_8") return position_y_8;          if (name == "position_z_8") return position_z_8;          if (name == "position_x_9") return position_x_9;          if (name == "position_y_9") return position_y_9;          if (name == "position_z_9") return position_z_9;          if (name == "position_x_10") return position_x_10;          if (name == "position_y_10") return position_y_10;          if (name == "position_z_10") return position_z_10;          if (name == "position_x_11") return position_x_11;          if (name == "position_y_11") return position_y_11;          if (name == "position_z_11") return position_z_11;          if (name == "position_x_12") return position_x_12;          if (name == "position_y_12") return position_y_12;          if (name == "position_z_12") return position_z_12;          if (name == "position_x_13") return position_x_13;          if (name == "position_y_13") return position_y_13;          if (name == "position_z_13") return position_z_13;          if (name == "position_x_14") return position_x_14;          if (name == "position_y_14") return position_y_14;          if (name == "position_z_14") return position_z_14;          if (name == "position_x_15") return position_x_15;          if (name == "position_y_15") return position_y_15;          if (name == "position_z_15") return position_z_15;          if (name == "position_x_16") return position_x_16;          if (name == "position_y_16") return position_y_16;          if (name == "position_z_16") return position_z_16;          if (name == "position_x_17") return position_x_17;          if (name == "position_y_17") return position_y_17;          if (name == "position_z_17") return position_z_17;          if (name == "position_x_18") return position_x_18;          if (name == "position_y_18") return position_y_18;          if (name == "position_z_18") return position_z_18;          if (name == "position_x_19") return position_x_19;          if (name == "position_y_19") return position_y_19;          if (name == "position_z_19") return position_z_19;          if (name == "position_x_20") return position_x_20;          if (name == "position_y_20") return position_y_20;          if (name == "position_z_20") return position_z_20;          if (name == "position_x_21") return position_x_21;          if (name == "position_y_21") return position_y_21;          if (name == "position_z_21") return position_z_21;          if (name == "position_x_22") return position_x_22;          if (name == "position_y_22") return position_y_22;          if (name == "position_z_22") return position_z_22;          if (name == "position_x_23") return position_x_23;          if (name == "position_y_23") return position_y_23;          if (name == "position_z_23") return position_z_23;          if (name == "position_x_24") return position_x_24;          if (name == "position_y_24") return position_y_24;          if (name == "position_z_24") return position_z_24;          if (name == "position_x_25") return position_x_25;          if (name == "position_y_25") return position_y_25;          if (name == "position_z_25") return position_z_25;          if (name == "position_x_26") return position_x_26;          if (name == "position_y_26") return position_y_26;          if (name == "position_z_26") return position_z_26;          if (name == "position_x_27") return position_x_27;          if (name == "position_y_27") return position_y_27;          if (name == "position_z_27") return position_z_27;          if (name == "position_x_28") return position_x_28;          if (name == "position_y_28") return position_y_28;          if (name == "position_z_28") return position_z_28;          if (name == "position_x_29") return position_x_29;          if (name == "position_y_29") return position_y_29;          if (name == "position_z_29") return position_z_29;          if (name == "position_x_30") return position_x_30;          if (name == "position_y_30") return position_y_30;          if (name == "position_z_30") return position_z_30;          if (name == "position_x_31") return position_x_31;          if (name == "position_y_31") return position_y_31;          if (name == "position_z_31") return position_z_31;          if (name == "position_x_32") return position_x_32;          if (name == "position_y_32") return position_y_32;          if (name == "position_z_32") return position_z_32;          if (name == "position_x_33") return position_x_33;          if (name == "position_y_33") return position_y_33;          if (name == "position_z_33") return position_z_33;          if (name == "position_x_34") return position_x_34;          if (name == "position_y_34") return position_y_34;          if (name == "position_z_34") return position_z_34;          if (name == "position_x_35") return position_x_35;          if (name == "position_y_35") return position_y_35;          if (name == "position_z_35") return position_z_35;          if (name == "position_x_36") return position_x_36;          if (name == "position_y_36") return position_y_36;          if (name == "position_z_36") return position_z_36;          if (name == "position_x_37") return position_x_37;          if (name == "position_y_37") return position_y_37;          if (name == "position_z_37") return position_z_37;          if (name == "position_x_38") return position_x_38;          if (name == "position_y_38") return position_y_38;          if (name == "position_z_38") return position_z_38;          if (name == "position_x_39") return position_x_39;          if (name == "position_y_39") return position_y_39;          if (name == "position_z_39") return position_z_39;          if (name == "position_x_40") return position_x_40;          if (name == "position_y_40") return position_y_40;          if (name == "position_z_40") return position_z_40;          if (name == "position_x_41") return position_x_41;          if (name == "position_y_41") return position_y_41;          if (name == "position_z_41") return position_z_41;          if (name == "position_x_42") return position_x_42;          if (name == "position_y_42") return position_y_42;          if (name == "position_z_42") return position_z_42;          if (name == "position_x_43") return position_x_43;          if (name == "position_y_43") return position_y_43;          if (name == "position_z_43") return position_z_43;          if (name == "position_x_44") return position_x_44;          if (name == "position_y_44") return position_y_44;          if (name == "position_z_44") return position_z_44;          if (name == "position_x_45") return position_x_45;          if (name == "position_y_45") return position_y_45;          if (name == "position_z_45") return position_z_45;          if (name == "position_x_46") return position_x_46;          if (name == "position_y_46") return position_y_46;          if (name == "position_z_46") return position_z_46;          if (name == "position_x_47") return position_x_47;          if (name == "position_y_47") return position_y_47;          if (name == "position_z_47") return position_z_47;          if (name == "position_x_48") return position_x_48;          if (name == "position_y_48") return position_y_48;          if (name == "position_z_48") return position_z_48;          if (name == "position_x_49") return position_x_49;          if (name == "position_y_49") return position_y_49;          if (name == "position_z_49") return position_z_49;
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

extern "C" DLL_EXPORT dolfin::Expression * create_dolfin_expression_573067dc3c4975b6ca0d597b86897607()
{
  return new dolfin::dolfin_expression_573067dc3c4975b6ca0d597b86897607;
}

