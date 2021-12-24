/*
  Copyright (C) 2016 - 2018 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  ASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/

#include <aspect/initial_temperature/potential_temperature.h>
#include <aspect/geometry_model/interface.h>
#include <aspect/gravity_model/ascii_data.h>
#include <aspect/utilities.h>
#include <deal.II/base/signaling_nan.h>


namespace aspect
{
  namespace InitialTemperature
  {

   template <int dim>
   PotentialTemperature<dim>::PotentialTemperature()
   :
   surface_boundary_id(5)
   {}
  
   template <int dim>
   void
   PotentialTemperature<dim>::initialize ()
   {
    // initialize AsciiDataInitial to read tomography file
    Utilities::AsciiDataInitial<dim>::initialize(1);
    // Find the boundary indicator that represents the surface
    surface_boundary_id = this->get_geometry_model().translate_symbolic_boundary_name_to_id("outer");

    std::set<types::boundary_id> surface_boundary_set;
    surface_boundary_set.insert(surface_boundary_id);
    // The input ascii table contains two components, the crust depth and the LAB depth
    ascii_data_lab.initialize(surface_boundary_set,1);
   }


   template <int dim>
   double
   PotentialTemperature<dim>::
   ascii_grid_vs (const Point<dim> &position) const
   {
       return  Utilities::AsciiDataInitial<dim>::get_data_component(position,0);
   }

    
    template <int dim>
    double
    PotentialTemperature<dim>::initial_temperature (const Point<dim> &position) const
    {
        const double depth = this->get_geometry_model().depth(position);
        double vs_perturbation = 0.0;

       //   Read ascii data file containing depth of lithospheric layers. 
        double LAB_depth;
        
           LAB_depth =    ascii_data_lab.get_data_component(surface_boundary_id, position, 0);
           

        
        double temperature_perturbation = 0.0;
//        double lab_temperature;
        aspect::Utilities::tk::spline s;
        s.set_points(spline_depths, reference_Vs);
         
        if (use_reference_profile)
        {
               vs_perturbation = (ascii_grid_vs(position) - s(depth))/s(depth);
        }
        else
               vs_perturbation = ascii_grid_vs(position)*0.01;

    	 // scale the density perturbation into a temperature perturbation. Velocity perturbation make sense only above 400Km depth. Do not apply perturbation below that
        if (depth > LAB_depth)
           {
               // scale the density perturbation into a temperature perturbation
              temperature_perturbation = -1./thermal_alpha * vs_to_density * vs_perturbation;
	// normalize temperature_perturbation between [-100,100]. 100 is the maximum excess temperature in the RVP (1450-1350). Where 1450 is the maximum mantle potential temperature for an ambient mantle temperature of 1350 (Rooney et al., 2012). Thus for interval [a,b], dTnormalized = (b-a)*(dT - dTmin)/(dTmax-dTmin) + a.
          std::cout << std::max(temperature_perturbation) << std::endl;

	      temperature_perturbation = 200*(temperature_perturbation - std::min(temperature_perturbation, 0.0) / (std::max(temperature_perturbation, 0.0) - std::min(temperature_perturbation))) - 100;
	      std::cout << "The minimum excess temperature if:" << std::min(temperature_perturbation, 0.0) << std::endl;
	      std::cout << "The maximum excess temperature if:" << std::max(temperature_perturbation, 0.0) << std::endl;


           }
         else
           {
              // set heterogeneity to zero down to a specified depth
               temperature_perturbation = 0.0;
           }

           if (add_temperature_perturbation == false)
              temperature_perturbation = 0.0;
////////////////////////////////////       
        double temperature;
	double LAB_temperature;
//	const gravity = this->get_gravity_model().g(position);
	const double gravity = 9.81;

		LAB_temperature = mantle_potential_temperature * std::exp(gravity * thermal_alpha * LAB_depth / specific_heat);

        if (depth < LAB_depth)
		// approximating a conductive geotherm in the lithosphere

		temperature  =  surface_temperature + (depth/LAB_depth)*(LAB_temperature - surface_temperature);
         else 
//              // the sublithospheric mantle is fully adiabatic

		temperature = mantle_potential_temperature * std::exp(gravity * thermal_alpha * depth / specific_heat) + temperature_perturbation;
//		temperature = mantle_potential_temperature * std::exp(gravity * thermal_alpha * depth / specific_heat);

        return temperature;


/////////////////////////////////////
/*
	double temperature;

        if (depth < isotherm_depth)
           temperature  =  surface_temperature + (depth/isotherm_depth)*(lab_isotherm_temperature +(isotherm_depth * temperature_gradient)  - surface_temperature);
         else
           temperature   =  lab_isotherm_temperature + (depth * temperature_gradient) + temperature_perturbation;
        return temperature;

*/
////////////////////////////////////////
//        return temperature_perturbation;
    }

    template <int dim>
    void
    PotentialTemperature<dim>::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection ("Initial temperature model");
      {
        prm.enter_subsection("Potential temperature");
        {
             prm.enter_subsection("LAB file");
             {
                 Utilities::AsciiDataBoundary<dim>::declare_parameters(prm,
                                                            "$ASPECT_SOURCE_DIR/data/initial-temperature/ascii-data/",
//                                                             "Emry.litho.aspect.input.txt");
                                                             "litho.africa.fishwick.txt");
//                                                             "litho.africa.200km.txt");
//                                                             "new.litho.africa.plus10.txt");
             }
            prm.leave_subsection();
	    
  
            prm.enter_subsection("Tomography file");
             {
 
               Utilities::AsciiDataBase<dim>::declare_parameters(prm,
                                                              "$ASPECT_SOURCE_DIR/data/initial-temperature/ascii-data/",
//                                                              "fishwick_shear_wave_velocity.txt");
                                                              "Emry_shear_wave_velocity.txt");
             }
           prm.leave_subsection();

          prm.declare_entry ("Specific heat", "1250",
                             Patterns::Double (0),
                             "The value of the specific heat of the mantle $C_p$. "
                             "Units: $J/kg/K$.");

            prm.declare_entry ("Thermal expansion coefficient", "3e-5",
          	                 Patterns::Double (0),
                             "The value of the thermal expansion coefficient $\\alpha$. "
                             "Units: $1/K$.");
//            prm.declare_entry ("Vs to density", "0.2",
            prm.declare_entry ("Vs to density", "0.15",
                      	     Patterns::Double (0),
                      	     "Vs to density scaling factor");
            prm.declare_entry ("Mantle potential temperature", "1600",
                              Patterns::Double (0),
                              "Temperature the mantle will have if decompressed to the surface without melting ");
            prm.declare_entry ("Surface temperature", "293",
                             Patterns::Double (0),
                             "Temperature at the surface or top boundary of the model. ");
//            prm.declare_entry ("LAB temperature", "1673.15",
//                              Patterns::Double (0),
//                              "Temperatue at the LAB. ");
            prm.declare_entry ("Adiabatic temperature gradient", "0.0004",
                             Patterns::Double (0),
                             "The value of the adiabatic temperature gradient. Units: $K m^{-1}$.");
            prm.declare_entry ("Data dir", "$ASPECT_SOURCE_DIR/data/initial-temperature/ascii-data/",
                             Patterns::DirectoryName (),
                             "Directory where the reference profile file is located.");
            prm.declare_entry ("Reference profile filename",
                             "AK135F_AVG.txt",
                             Patterns::FileName (),
                             "File from which the isotherm depth data is read.");
            prm.declare_entry ("Use reference profile", "true",
                             Patterns::Bool (),
                             "Option to take the thermal expansion coefficient from the "
                             "material model instead of from what is specified in this "
                             "section.");
            prm.declare_entry ("Use temperature perturbation", "true",
                             Patterns::Bool (),
                             "Option to take the thermal expansion coefficient from the "
                             "material model instead of from what is specified in this "
                             "section.");
        }
      prm.leave_subsection();
      }
      prm.leave_subsection();

    }

    template <int dim>
    void
    PotentialTemperature<dim>::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection ("Initial temperature model");
      {
        prm.enter_subsection("Potential temperature");
        {
         prm.enter_subsection("LAB file");
         {
            ascii_data_lab.initialize_simulator (this->get_simulator());
            ascii_data_lab.parse_parameters(prm);
         }
         prm.leave_subsection();
          
         prm.enter_subsection("Tomography file");
         {
               Utilities::AsciiDataBase<dim>::parse_parameters(prm);
         }
         prm.leave_subsection();

	  specific_heat				= prm.get_double ("Specific heat");
          vs_to_density                         = prm.get_double ("Vs to density");
          thermal_alpha                         = prm.get_double ("Thermal expansion coefficient");
          mantle_potential_temperature          = prm.get_double ("Mantle potential temperature");
//          LAB_temperature	                = prm.get_double ("LAB temperature");
          surface_temperature                   = prm.get_double ("Surface temperature");
          reference_profile_filename            = prm.get ("Reference profile filename");
          data_dir                             = Utilities::expand_ASPECT_SOURCE_DIR(prm.get ("Data dir"));
          add_temperature_perturbation          = prm.get_bool ("Use temperature perturbation");
          use_reference_profile                 = prm.get_bool ("Use reference profile");
          temperature_gradient                  = prm.get_double("Adiabatic temperature gradient");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
     
//      if (use_reference_profile)
//      { 
      const std::string filename = data_dir+reference_profile_filename;
//      std::cout<<data_dir<<"   "<<reference_profile_filename<<std::endl; 
      /**
 *        * Read data from disk and distribute among processes
 *               */
      std::istringstream in(Utilities::read_and_distribute_file_content(filename, this->get_mpi_communicator()));

      /**
       * Reading data lines.
       */
      double depths, ref_vs;
      while (in >> depths >> ref_vs)
        {
          spline_depths.push_back(depths);
          reference_Vs.push_back(ref_vs);
        }
  
//      }
    
       /*Read craton file */
//      const std::string filename = data_dir+plate_boundaries_file_name;

      /**
       * Read data from disk and distribute among processes
       */
//       std::istringstream in(Utilities::read_and_distribute_file_content(filename, this->get_mpi_communicator()));

      /**
       * Reading data lines
       */
//       double longitude, latitude;

//      while (in >> longitude >> latitude)
//      {
//          Point<2> point (longitude,latitude);
//          boundaries_point_lists.push_back(point);
  //    }

    
     }

  }

}


namespace aspect
{
  namespace InitialTemperature
  {
    ASPECT_REGISTER_INITIAL_TEMPERATURE_MODEL(PotentialTemperature,
                                              "potential temperature",
                                              "An initial temperature condition in which the LAB is an adiabatic "
					      "defined by the mantle potential temperature "
                                              "is designed specifically for the ellipsoidal chunk geometry model.")
  }
}
