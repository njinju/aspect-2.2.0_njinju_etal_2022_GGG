/*
  Copyright (C) 2015 - 2019 by the authors of the ASPECT code.

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


#include <aspect/material_model/edc_melting.h>
#include <aspect/adiabatic_conditions/interface.h>

#include <deal.II/base/parameter_handler.h>
#include <deal.II/numerics/fe_field_function.h>
#include <aspect/melt.h>


namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    std::vector<double>
    EdcMelting<dim>::
    compute_volume_fractions( const std::vector<double> &compositional_fields) const
    {
      std::vector<double> volume_fractions( compositional_fields.size()+1);

      // clip the compositional fields so they are between zero and one
      std::vector<double> x_comp = compositional_fields;
      for ( unsigned int i=0; i < x_comp.size(); ++i)
        x_comp[i] = std::min(std::max(x_comp[i], 0.0), 1.0);

      // sum the compositional fields for normalization purposes
      double sum_composition = 0.0;
      for ( unsigned int i=0; i < x_comp.size(); ++i)
        sum_composition += x_comp[i];

      if (sum_composition >= 1.0)
        {
          volume_fractions[0] = 0.0;  // background mantle
          for ( unsigned int i=1; i <= x_comp.size(); ++i)
            volume_fractions[i] = x_comp[i-1]/sum_composition;
        }
      else
        {
          volume_fractions[0] = 1.0 - sum_composition; // background mantle
          for ( unsigned int i=1; i <= x_comp.size(); ++i)
            volume_fractions[i] = x_comp[i-1];
        }
      return volume_fractions;
    }

    template <int dim>
    double
    EdcMelting<dim>::
    average_value ( const std::vector<double> &composition,
                    const std::vector<double> &parameter_values,
                    const enum averaging_scheme &average_type) const
    {
      double averaged_parameter = 0.0;
      const std::vector<double> volume_fractions = compute_volume_fractions(composition);

      switch (average_type)
        {
          case arithmetic:
          {
            for (unsigned int i=0; i< volume_fractions.size(); ++i)
              averaged_parameter += volume_fractions[i]*parameter_values[i];
            break;
          }
          case harmonic:
          {
            for (unsigned int i=0; i< volume_fractions.size(); ++i)
              averaged_parameter += volume_fractions[i]/(parameter_values[i]);
            averaged_parameter = 1.0/averaged_parameter;
            break;
          }
          case geometric:
          {
            for (unsigned int i=0; i < volume_fractions.size(); ++i)
              averaged_parameter += volume_fractions[i]*std::log(parameter_values[i]);
            averaged_parameter = std::exp(averaged_parameter);
            break;
          }
          case maximum_composition:
          {
            const unsigned int i = (unsigned int)(std::max_element( volume_fractions.begin(),
                                                                    volume_fractions.end() )
                                                  - volume_fractions.begin());
            averaged_parameter = parameter_values[i];
            break;

          }
          default:
          {
            AssertThrow( false, ExcNotImplemented() );
            break;
          }
        }
      return averaged_parameter;
    }

    template <int dim>
    double
    EdcMelting<dim>::
    diffusion_creep (const double &pressure,
                     const double &temperature) const
    {
      return  0.5 * std::pow(prefactor_diffusion,-1/stress_exponent_diffusion) *
              std::exp(std::max((activation_energie_diffusion + pressure*activation_volume_diffusion),0.0)/
                       (constants::gas_constant*temperature*stress_exponent_diffusion)) *
              std::pow(grain_size, grain_size_exponent_diffusion);
    }

    template <int dim>
    double
    EdcMelting<dim>::
    dislocation_creep (const double &pressure,
                       const double &temperature,
                       const SymmetricTensor<2,dim> &strain_rate) const
    {
      const double edot_ii = ( (this->get_timestep_number() == 0 && strain_rate.norm() <= std::numeric_limits<double>::min())
                               ?
                               ref_strain_rate
                               :
                               std::max(std::sqrt(std::fabs(second_invariant(deviator(strain_rate)))),
                                        min_strain_rate) );

      return 0.5 * std::pow(prefactor_dislocation,-1/stress_exponent_dislocation) *
             std::exp(std::max((activation_energie_dislocation + pressure*activation_volume_dislocation),0.0)/
                      (constants::gas_constant*temperature*stress_exponent_dislocation)) *
             std::pow(edot_ii,((1. - stress_exponent_dislocation)/stress_exponent_dislocation));
    }

    template <int dim>
    double
    EdcMelting<dim>::
    composite_viscosity (double &dislocation_creep,
                         double &diffusion_creep) const
    {

      return (diffusion_creep * dislocation_creep)/(diffusion_creep + dislocation_creep);
    }

    template <>
    void
    EdcMelting<2>::evaluate(const MaterialModel::MaterialModelInputs<2> &in,
                            MaterialModel::MaterialModelOutputs<2> &out) const
    {
      Assert (false, ExcNotImplemented());
      //return 0;
    }

    template <int dim>
    double
    EdcMelting<dim>::
    reference_viscosity () const
    {
      return eta_0;
    }

    template <int dim>
    double
    EdcMelting<dim>::
    reference_darcy_coefficient () const
    {
      // 0.01 = 1% melt
      return reference_permeability * std::pow(0.01,3.0) / eta_f;
    }

    template <int dim>
    bool
    EdcMelting<dim>::
    is_compressible () const
    {
//      return false;
      return model_is_compressible;
    }

    template <int dim>
    double
    EdcMelting<dim>::
    melt_fraction (const double temperature,
                   const double pressure) const

    {
      // anhydrous melting of peridotite after Katz, 2003
      const double T_solidus  = A1 + 273.15
                                + A2 * pressure
                                + A3 * pressure * pressure;
      const double T_lherz_liquidus = B1 + 273.15
                                      + B2 * pressure
                                      + B3 * pressure * pressure;
      const double T_liquidus = C1 + 273.15
                                + C2 * pressure
                                + C3 * pressure * pressure;

      // melt fraction for peridotite with clinopyroxene
      double melt_fraction;
      if (temperature < T_solidus || pressure > 1.3e10)
        melt_fraction = 0.0;
      else if (temperature > T_lherz_liquidus)
        melt_fraction = 1.0;
      else
        melt_fraction = std::pow((temperature - T_solidus) / (T_lherz_liquidus - T_solidus),beta);

      // melt fraction after melting of all clinopyroxene
      const double R_cpx = r1 + r2 * std::max(0.0, pressure);
      const double F_max = M_cpx / R_cpx;

      if (melt_fraction > F_max && temperature < T_liquidus)
        {
          const double T_max = std::pow(F_max,1/beta) * (T_lherz_liquidus - T_solidus) + T_solidus;
          melt_fraction = F_max + (1 - F_max) * pow((temperature - T_max) / (T_liquidus - T_max),beta);
        }
      return melt_fraction;
    }


    template <int dim>
    double
    EdcMelting<dim>::
    entropy_change (const double temperature,
                    const double pressure,
                    const double maximum_melt_fraction,
                    const NonlinearDependence::Dependence dependence) const
    {
      double entropy_gradient = 0.0;

      // calculate latent heat of melting
      // we need the change of melt fraction in dependence of pressure and temperature

      // for peridotite after Katz, 2003
      const double T_solidus        = A1 + 273.15
                                      + A2 * pressure
                                      + A3 * pressure * pressure;
      const double T_lherz_liquidus = B1 + 273.15
                                      + B2 * pressure
                                      + B3 * pressure * pressure;
      const double T_liquidus       = C1 + 273.15
                                      + C2 * pressure
                                      + C3 * pressure * pressure;
      const double dT_solidus_dp        = A2 + 2 * A3 * pressure;
      const double dT_lherz_liquidus_dp = B2 + 2 * B3 * pressure;
      const double dT_liquidus_dp       = C2 + 2 * C3 * pressure;
                                                                                                                           if (temperature > T_solidus && temperature < T_liquidus && pressure < 1.3e10)
        {
      // melt fraction when clinopyroxene is still present
      double melt_fraction_derivative_temperature
        = beta * pow((temperature - T_solidus)/(T_lherz_liquidus - T_solidus),beta-1)
          / (T_lherz_liquidus - T_solidus);
                                                                                                                           double melt_fraction_derivative_pressure
        = beta * pow((temperature - T_solidus)/(T_lherz_liquidus - T_solidus),beta-1)
          * (dT_solidus_dp * (temperature - T_lherz_liquidus)
             + dT_lherz_liquidus_dp * (T_solidus - temperature))
          / pow(T_lherz_liquidus - T_solidus,2);
                                                                                                                           // melt fraction after melting of all clinopyroxene
       const double R_cpx = r1 + r2 * std::max(0.0, pressure);
       const double F_max = M_cpx / R_cpx;
                                                                                                                            if (melt_fraction(temperature, pressure) > F_max)
         {
           const double T_max = std::pow(F_max,1.0/beta) * (T_lherz_liquidus - T_solidus) + T_solidus;
           const double dF_max_dp = - M_cpx * std::pow(r1 + r2 * pressure,-2) * r2;
           const double dT_max_dp = dT_solidus_dp
                                     + 1.0/beta * std::pow(F_max,1.0/beta - 1.0) * dF_max_dp * (T_lherz_liquidus - T_solidus)
                                     + std::pow(F_max,1.0/beta) * (dT_lherz_liquidus_dp - dT_solidus_dp);

           melt_fraction_derivative_temperature
             = (1.0 - F_max) * beta * std::pow((temperature - T_max)/(T_liquidus - T_max),beta-1)
               / (T_liquidus - T_max);
                                                                                                                                melt_fraction_derivative_pressure
             = dF_max_dp
               - dF_max_dp * std::pow((temperature - T_max)/(T_liquidus - T_max),beta)
               + (1.0 - F_max) * beta * std::pow((temperature - T_max)/(T_liquidus - T_max),beta-1)
               * (dT_max_dp * (T_max - T_liquidus) - (dT_liquidus_dp - dT_max_dp) * (temperature - T_max)) / std::pow(T_liquidus - T_max, 2);
         }
      
        double melt_fraction_derivative = 0;
        if (dependence == NonlinearDependence::temperature)
            melt_fraction_derivative = melt_fraction_derivative_temperature;
          else if (dependence == NonlinearDependence::pressure)
            melt_fraction_derivative = melt_fraction_derivative_pressure;
          else
            AssertThrow(false, ExcMessage("not implemented"));

          if (melt_fraction(temperature, pressure) >= maximum_melt_fraction)
            entropy_gradient = melt_fraction_derivative * peridotite_melting_entropy_change;
        }
      return entropy_gradient;
    }


    template <int dim>
    void
    EdcMelting<dim>::
    melt_fractions (const MaterialModel::MaterialModelInputs<dim> &in,
                    std::vector<double> &melt_fractions) const
    {

      for (unsigned int q=0; q<in.temperature.size(); ++q)
        {
          if (this->include_melt_transport())
            {
              const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");
              const unsigned int peridotite_idx = this->introspection().compositional_index_for_name("peridotite");
//              depletion = in.composition[q][peridotite_idx] - in.composition[q][porosity_idx];
            }
          melt_fractions[q] = this->melt_fraction(in.temperature[q],
                                                this->get_adiabatic_conditions().pressure(in.position[q]));

        }
    }


    template <int dim>
    void
    EdcMelting<dim>::
    evaluate(const typename Interface<dim>::MaterialModelInputs &in, typename Interface<dim>::MaterialModelOutputs &out) const
    {
      std::vector<double> old_porosity(in.position.size());

      ReactionRateOutputs<dim> *reaction_rate_out = out.template get_additional_output<ReactionRateOutputs<dim> >();

      // we want to get the porosity field from the old solution here,
      // because we need a field that is not updated in the nonlinear iterations
      if (this->include_melt_transport() && in.current_cell.state() == IteratorState::valid
          && this->get_timestep_number() > 0 && !this->get_parameters().use_operator_splitting)
        {
          // Prepare the field function
          Functions::FEFieldFunction<dim, DoFHandler<dim>, LinearAlgebra::BlockVector>
          fe_value(this->get_dof_handler(), this->get_old_solution(), this->get_mapping());

          const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");

          fe_value.set_active_cell(in.current_cell);
          fe_value.value_list(in.position,
                              old_porosity,
                              this->introspection().component_indices.compositional_fields[porosity_idx]);
        }
      else if (this->get_parameters().use_operator_splitting)
        for (unsigned int i=0; i<in.position.size(); ++i)
          {
            const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");
            old_porosity[i] = in.composition[i][porosity_idx];
///// ////// ///////////// ///////////

          if (this->include_melt_transport() && in.strain_rate.size() > 0)
            {
              const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");
              const unsigned int peridotite_idx = this->introspection().compositional_index_for_name("peridotite");
              const double old_porosity = in.composition[i][porosity_idx];
              const double maximum_melt_fraction = in.composition[i][peridotite_idx];
              // calculate the melting rate as difference between the equilibrium melt fraction
              // and the solution of the previous time step
              double porosity_change = 0.0;
//              double freezing_rate = 0.0;
              if (fractional_melting)
                {
                 // solidus is lowered by previous melting events (fractional melting)
                  const double solidus_change = (maximum_melt_fraction - old_porosity) * depletion_solidus_change;
                  const double eq_melt_fraction = melt_fraction(in.temperature[i] - solidus_change, this->get_adiabatic_conditions().pressure(in.position[i]));
                  porosity_change = eq_melt_fraction - old_porosity;
                }
              else
                {
                 // batch melting
                  porosity_change = melt_fraction(in.temperature[i], this->get_adiabatic_conditions().pressure(in.position[i]))
                                    - std::max(maximum_melt_fraction, 0.0);
                  porosity_change = std::max(porosity_change, 0.0);                 

                 // freezing of melt below the solidus
                  {
                    const double eq_melt_fraction = melt_fraction(in.temperature[i], this->get_adiabatic_conditions().pressure(in.position[i]));
 
                 // porosity reaches the equilibrium melt fraction:
                    const double porosity_change_wrt_melt_fraction = std::min(eq_melt_fraction - old_porosity - porosity_change,0.0);
                 // depletion reaches the equilibrium melt fraction:
                    const double porosity_change_wrt_depletion = std::min((eq_melt_fraction - std::max(maximum_melt_fraction, 0.0))
                                                                          * (1.0 - old_porosity) / (1.0 - maximum_melt_fraction),0.0);
                    double freezing_amount = std::max(porosity_change_wrt_melt_fraction, porosity_change_wrt_depletion);

                    if (eq_melt_fraction == 0.0)
                      freezing_amount = - old_porosity;

                    porosity_change += freezing_amount;

                    if (porosity_change < 0 )
                      porosity_change *= freezing_rate * melting_time_scale;
                  }
                }
               // remove melt that gets to the base of the lithosphere
              if (this->get_geometry_model().depth(in.position[i]) < extraction_depth)
                porosity_change = -old_porosity * (in.position[i](1) - (this->get_geometry_model().maximal_depth() - extraction_depth))/extraction_depth;

              // do not allow negative porosity  
              porosity_change = std::max(porosity_change, -old_porosity);

	    }       
          }

      for (unsigned int i=0; i<in.position.size(); ++i)
        {
/////////////////////////////
          // Calculate Density
          const Point<3> pos = in.position[i];
          const double depth = this->get_geometry_model().depth(in.position[i]);
          const std::vector<double> &composition = in.composition[i];
          std::vector<double> composition_densities (composition.size()+1);
          std::vector<double> composition_viscosities (composition.size()+1);
          std::vector<double> composition_thermal_conductivities (composition.size()+1);
          AssertThrow(this->introspection().compositional_name_exists("crust"),
                      ExcMessage("Material model Edc melting works only works if there is a "
                                 "compositional field called crust."));

          AssertThrow(this->introspection().compositional_name_exists("mantle_lithosphere"),
                      ExcMessage("Material model Edc melting works only works if there is a "
                                 "compositional field called mantle_lithosphere."));
          // The default index for the first compositional fields is 0. During the compositioning
          // the crust needs have index 1 because the index 0 is by default the background mantle.
          const unsigned int crust_idx = this->introspection().compositional_index_for_name("crust") + 1;
          const unsigned int mantle_lithosphere_idx = this->introspection().compositional_index_for_name("mantle_lithosphere") + 1;
          const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");
          const unsigned int peridotite_idx = this->introspection().compositional_index_for_name("peridotite");

          // Compositional densities.
          for (unsigned int c=0; c <= in.composition[i].size() ; ++c)
            {
              if (c == crust_idx)
               {
		 composition_densities[c] = 2700.0;
                 composition_thermal_conductivities[c] = 2.5;
//                composition_densities[c] = 2700.0 * (1 - thermal_expansivity * (in.temperature[i] - reference_T));
	       }
              else if (c == mantle_lithosphere_idx)
	       {
                 composition_densities[c] = 3400.0;
                 composition_thermal_conductivities[c] = 3.5;
//              else if (c == porosity_idx || c == peridotite_idx)
	       }
              else
                {
                 composition_thermal_conductivities[c] = 4.7;
                  // calculate density first, we need it for the reaction term
                  // temperature dependence of density is 1 - alpha * (T - T(adiabatic))
                  double temperature_dependence = 1.0;
                  if (this->include_adiabatic_heating ())
                    temperature_dependence -= (in.temperature[i] - this->get_adiabatic_conditions().temperature(in.position[i]))
                                              * thermal_expansivity;
                  else
                    temperature_dependence -= (in.temperature[i] - reference_T) * thermal_expansivity;

                  // calculate composition dependence of density
                  const double delta_rho = this->introspection().compositional_name_exists("peridotite")
                                           ?
                                           depletion_density_change * in.composition[i][this->introspection().compositional_index_for_name("peridotite")]
                                           :
                                           0.0;

                  composition_densities[c] = (reference_rho_s + delta_rho) * temperature_dependence
                                             * std::exp(compressibility * (in.pressure[i] - this->get_surface_pressure()));
                }
          out.thermal_conductivities[i] = composition_thermal_conductivities[c];
            }
//                  std::cout<<density_averaging<<std::endl;
          out.densities[i] = average_value (composition, composition_densities, density_averaging);
//          out.thermal_conductivities[i] = average_value (composition, composition_thermal_conductivities, density_averaging);

          // Calculate Viscosity
          const double temperature = in.temperature[i];
          const double pressure= in.pressure[i];

          double dislo_creep  = std::min(std::max(dislocation_creep (pressure, temperature, in.strain_rate[i]), min_visc), max_visc);
          double diff_creep = std::min(std::max(diffusion_creep (pressure, temperature), min_visc), max_visc);


//          double dislo_creep = dislocation_creep (pressure, temperature, in.strain_rate[i]);
//          double diff_creep = diffusion_creep (pressure, temperature);

          for (unsigned int c=0; c <= in.composition[i].size() ; ++c)
            {
              if (c == crust_idx)
                {
//                composition_viscosities[c] =  crust_out.viscosities[i];
                  composition_viscosities[c] = 1.0e23;
                }
//             else if (c == mantle_lithosphere_idx)
//                 composition_densities[c] = 3400.0 * (1 - thermal_alpha * (in.temperature[i] - reference_T));

              else if (c == mantle_lithosphere_idx)
                {
//                composition_viscosities[c] =  crust_out.viscosities[i];
                  composition_viscosities[c] = 1.0e23;
//                composition_viscosities[c] = std::min(std::max(dislocation_creep (pressure, temperature, in.strain_rate[i]), min_visc), max_visc);
                }

//            else if (c == porosity_idx || c == peridotite_idx)
              else
                {
//                composition_viscosities[c] = std::min(std::max(diffusion_creep (pressure, temperature), min_visc), max_visc);

                  double visc_temperature_dependence = 1.0;
                  if (this->include_adiabatic_heating ())
                    {
                      const double delta_temp = in.temperature[i]-this->get_adiabatic_conditions().temperature(in.position[i]);
                      visc_temperature_dependence = std::max(std::min(std::exp(-thermal_viscosity_exponent*delta_temp/this->get_adiabatic_conditions().temperature(in.position[i])),1e4),1e-4);
                    }
                  else
                    {
                      const double delta_temp = in.temperature[i]-reference_T;
                      visc_temperature_dependence = std::max(std::min(std::exp(-thermal_viscosity_exponent*delta_temp/reference_T),1e4),1e-4);
                    }
//          out.viscosities[i] *= visc_temperature_dependence;

                  composition_viscosities[c] = (dislo_creep * diff_creep) / (dislo_creep + diff_creep);
                }
            }
          out.viscosities[i] = average_value (composition, composition_viscosities, viscosity_averaging);


          for (unsigned int c=0; c<in.composition[i].size(); ++c)
            {
              out.reaction_terms[i][c] = 0.0;

              if (this->get_parameters().use_operator_splitting && reaction_rate_out != nullptr)
                reaction_rate_out->reaction_rates[i][c] = 0.0;
            }

          if (this->include_melt_transport())
            {
              const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");
              const double porosity = std::min(1.0, std::max(in.composition[i][porosity_idx],0.0));

              // calculate viscosity based on local melt
//              out.viscosities[i] *= exp(- alpha_phi * porosity);

              if (include_melting_and_freezing && in.strain_rate.size())
                {
                  const unsigned int peridotite_idx = this->introspection().compositional_index_for_name("peridotite");

                  // Calculate the melting rate as difference between the equilibrium melt fraction
                  // and the solution of the previous time step (or the current solution, in case
                  // operator splitting is used).
                  // The solidus is lowered by previous melting events (fractional melting).
                  const double eq_melt_fraction = melt_fraction(in.temperature[i],
                                                                this->get_adiabatic_conditions().pressure(in.position[i]));
                  double porosity_change = eq_melt_fraction - old_porosity[i];
                  // do not allow negative porosity
                  if (old_porosity[i] + porosity_change < 0)
                    porosity_change = -old_porosity[i];

                  for (unsigned int c=0; c<in.composition[i].size(); ++c)
                    {
                      if (c == peridotite_idx && this->get_timestep_number() > 1)
                        out.reaction_terms[i][c] = porosity_change
                                                   - in.composition[i][peridotite_idx] * trace(in.strain_rate[i]) * this->get_timestep();
                      else if (c == porosity_idx && this->get_timestep_number() > 1)
                        out.reaction_terms[i][c] = porosity_change
                                                   * out.densities[i] / this->get_timestep();
                      else
                        out.reaction_terms[i][c] = 0.0;

                      // fill reaction rate outputs if the model uses operator splitting
                      if (this->get_parameters().use_operator_splitting)
                        {
                          if (reaction_rate_out != nullptr)
                            {
                              if (c == peridotite_idx && this->get_timestep_number() > 0)
                                reaction_rate_out->reaction_rates[i][c] = porosity_change / melting_time_scale
                                                                          - in.composition[i][peridotite_idx] * trace(in.strain_rate[i]);
                              else if (c == porosity_idx && this->get_timestep_number() > 0)
                                reaction_rate_out->reaction_rates[i][c] = porosity_change / melting_time_scale;
                              else
                                reaction_rate_out->reaction_rates[i][c] = 0.0;
                            }
                          out.reaction_terms[i][c] = 0.0;
                        }
                    }

                  // find depletion = peridotite, which might affect shear viscosity:
                  const double depletion_visc = std::min(1.0, std::max(in.composition[i][peridotite_idx],0.0));

                  // calculate strengthening due to depletion:
                  const double depletion_strengthening = std::min(exp(alpha_depletion*depletion_visc),delta_eta_depletion_max);

                  // calculate viscosity change due to local melt and depletion:
//                out.viscosities[i] *= depletion_strengthening;
                }
            }

          out.entropy_derivative_pressure[i]    = 0.0;
          out.entropy_derivative_temperature[i] = 0.0;
          out.thermal_expansion_coefficients[i] = thermal_expansivity;
          out.specific_heat[i] = reference_specific_heat;
//          out.thermal_conductivities[i] = thermal_conductivity;
          out.compressibilities[i] = 0.0;

          double visc_temperature_dependence = 1.0;
          if (this->include_adiabatic_heating ())
            {
              const double delta_temp = in.temperature[i]-this->get_adiabatic_conditions().temperature(in.position[i]);
              visc_temperature_dependence = std::max(std::min(std::exp(-thermal_viscosity_exponent*delta_temp/this->get_adiabatic_conditions().temperature(in.position[i])),1e4),1e-4);
            }
          else if (thermal_viscosity_exponent != 0.0)
            {
              const double delta_temp = in.temperature[i]-reference_T;
              visc_temperature_dependence = std::max(std::min(std::exp(-thermal_viscosity_exponent*delta_temp/reference_T),1e4),1e-4);
            }
//          out.viscosities[i] *= visc_temperature_dependence;
        }

      // fill melt outputs if they exist
      MeltOutputs<dim> *melt_out = out.template get_additional_output<MeltOutputs<dim> >();

      if (melt_out != nullptr)
        {
          const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");

          for (unsigned int i=0; i<in.position.size(); ++i)
            {
              double porosity = std::max(in.composition[i][porosity_idx],0.0);

              melt_out->fluid_viscosities[i] = eta_f;
              melt_out->permeabilities[i] = reference_permeability * std::pow(porosity,3) * std::pow(1.0-porosity,2);
              melt_out->fluid_density_gradients[i] = Tensor<1,dim>();

              // temperature dependence of density is 1 - alpha * (T - T(adiabatic))
              double temperature_dependence = 1.0;
              if (this->include_adiabatic_heating ())
                temperature_dependence -= (in.temperature[i] - this->get_adiabatic_conditions().temperature(in.position[i]))
                                          * thermal_expansivity;
              else
                temperature_dependence -= (in.temperature[i] - reference_T) * thermal_expansivity;
              melt_out->fluid_densities[i] = reference_rho_f * temperature_dependence
                                             * std::exp(melt_compressibility * (in.pressure[i] - this->get_surface_pressure()));

              melt_out->compaction_viscosities[i] = xi_0 * exp(- alpha_phi * porosity);

              double visc_temperature_dependence = 1.0;
              if (this->include_adiabatic_heating ())
                {
                  const double delta_temp = in.temperature[i]-this->get_adiabatic_conditions().temperature(in.position[i]);
                  visc_temperature_dependence = std::max(std::min(std::exp(-thermal_bulk_viscosity_exponent*delta_temp/this->get_adiabatic_conditions().temperature(in.position[i])),1e4),1e-4);
                }
              else if (thermal_viscosity_exponent != 0.0)
                {
                  const double delta_temp = in.temperature[i]-reference_T;
                  visc_temperature_dependence = std::max(std::min(std::exp(-thermal_bulk_viscosity_exponent*delta_temp/reference_T),1e4),1e-4);
                }
              melt_out->compaction_viscosities[i] *= visc_temperature_dependence;
            }
        }
    }



    template <int dim>
    void
    EdcMelting<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Edc melting");
        {
          prm.declare_entry ("Reference strain rate","1.0e-15",Patterns::Double(0),
                             "Reference strain rate for first time step. Units: $1 / s$");
          prm.declare_entry ("Reference solid density", "3000",
                             Patterns::Double (0),
                             "Reference density of the solid $\\rho_{s,0}$. Units: $kg/m^3$.");
          prm.declare_entry ("Reference melt density", "2500",
                             Patterns::Double (0),
                             "Reference density of the melt/fluid$\\rho_{f,0}$. Units: $kg/m^3$.");
          prm.declare_entry ("Reference temperature", "293",
                             Patterns::Double (0),
                             "The reference temperature $T_0$. The reference temperature is used "
                             "in both the density and viscosity formulas. Units: $\\si{K}$.");
          prm.declare_entry ("Reference shear viscosity", "5e20",
                             Patterns::Double (0),
                             "The value of the constant viscosity $\\eta_0$ of the solid matrix. "
                             "This viscosity may be modified by both temperature and porosity "
                             "dependencies. Units: $Pa \\, s$.");
          prm.declare_entry ("Reference bulk viscosity", "1e22",
                             Patterns::Double (0),
                             "The value of the constant bulk viscosity $\\xi_0$ of the solid matrix. "
                             "This viscosity may be modified by both temperature and porosity "
                             "dependencies. Units: $Pa \\, s$.");
          prm.declare_entry ("Reference melt viscosity", "10",
                             Patterns::Double (0),
                             "The value of the constant melt viscosity $\\eta_f$. Units: $Pa \\, s$.");
          prm.declare_entry ("Exponential melt weakening factor", "27",
                             Patterns::Double (0),
                             "The porosity dependence of the viscosity. Units: dimensionless.");
          prm.declare_entry ("Thermal viscosity exponent", "0.0",
                             Patterns::Double (0),
                             "The temperature dependence of the shear viscosity. Dimensionless exponent. "
                             "See the general documentation "
                             "of this model for a formula that states the dependence of the "
                             "viscosity on this factor, which is called $\\beta$ there.");
          prm.declare_entry ("Thermal bulk viscosity exponent", "0.0",
                             Patterns::Double (0),
                             "The temperature dependence of the bulk viscosity. Dimensionless exponent. "
                             "See the general documentation "
                             "of this model for a formula that states the dependence of the "
                             "viscosity on this factor, which is called $\\beta$ there.");
          prm.declare_entry ("Use full compressibility", "true",
                             Patterns::Bool (),
                             "If the compressibility should be used everywhere in the code "
                             "(if true), changing the volume of material when the density changes, "
                             "or only in the momentum conservation and advection equations "
                             "(if false).");
          prm.declare_entry ("Thermal conductivity", "4.7",
                             Patterns::Double (0),
                             "The value of the thermal conductivity $k$. "
                             "Units: $W/m/K$.");
          prm.declare_entry ("Reference specific heat", "1250",
                             Patterns::Double (0),
                             "The value of the specific heat $C_p$. "
                             "Units: $J/kg/K$.");
          prm.declare_entry ("Thermal expansion coefficient", "2e-5",
                             Patterns::Double (0),
                             "The value of the thermal expansion coefficient $\\beta$. "
                             "Units: $1/K$.");
          prm.declare_entry ("Reference permeability", "1e-8",
                             Patterns::Double(),
                             "Reference permeability of the solid host rock."
                             "Units: $m^2$.");
          prm.declare_entry ("Depletion density change", "0.0",
                             Patterns::Double (),
                             "The density contrast between material with a depletion of 1 and a "
                             "depletion of zero. Negative values indicate lower densities of "
                             "depleted material. Depletion is indicated by the compositional "
                             "field with the name peridotite. Not used if this field does not "
                             "exist in the model. "
                             "Units: $kg/m^3$.");
          prm.declare_entry ("Surface solidus", "1300",
                             Patterns::Double (0),
                             "Solidus for a pressure of zero. "
                             "Units: $\\si{K}$.");
          prm.declare_entry ("Depletion solidus change", "200.0",
                             Patterns::Double (),
                             "The solidus temperature change for a depletion of 100\\%. For positive "
                             "values, the solidus gets increased for a positive peridotite field "
                             "(depletion) and lowered for a negative peridotite field (enrichment). "
                             "Scaling with depletion is linear. Only active when fractional melting "
                             "is used. "
                             "Units: $\\si{K}$.");
          prm.declare_entry ("Pressure solidus change", "6e-8",
                             Patterns::Double (),
                             "The linear solidus temperature change with pressure. For positive "
                             "values, the solidus gets increased for positive pressures. "
                             "Units: $1/Pa$.");
          prm.declare_entry ("Solid compressibility", "0.0",
                             Patterns::Double (0),
                             "The value of the compressibility of the solid matrix. "
                             "Units: $1/Pa$.");
          prm.declare_entry ("Melt compressibility", "0.0",
                             Patterns::Double (0),
                             "The value of the compressibility of the melt. "
                             "Units: $1/Pa$.");
          prm.declare_entry ("Melt bulk modulus derivative", "0.0",
                             Patterns::Double (0),
                             "The value of the pressure derivative of the melt bulk "
                             "modulus. "
                             "Units: None.");
          prm.declare_entry ("Include melting and freezing", "true",
                             Patterns::Bool (),
                             "Whether to include melting and freezing (according to a simplified "
                             "linear melting approximation in the model (if true), or not (if "
                             "false).");
          prm.declare_entry ("Melting time scale for operator splitting", "1e3",
                             Patterns::Double (0),
                             "In case the operator splitting scheme is used, the porosity field can not "
                             "be set to a new equilibrium melt fraction instantly, but the model has to "
                             "provide a melting time scale instead. This time scale defines how fast melting "
                             "happens, or more specifically, the parameter defines the time after which "
                             "the deviation of the porosity from the equilibrium melt fraction will be "
                             "reduced to a fraction of $1/e$. So if the melting time scale is small compared "
                             "to the time step size, the reaction will be so fast that the porosity is very "
                             "close to the equilibrium melt fraction after reactions are computed. Conversely, "
                             "if the melting time scale is large compared to the time step size, almost no "
                             "melting and freezing will occur."
                             "\n\n"
                             "Also note that the melting time scale has to be larger than or equal to the reaction "
                             "time step used in the operator splitting scheme, otherwise reactions can not be "
                             "computed. If the model does not use operator splitting, this parameter is not used. "
                             "Units: yr or s, depending on the ``Use years "
                             "in output instead of seconds'' parameter.");
          prm.declare_entry ("Exponential depletion strengthening factor", "0.0",
                             Patterns::Double (0),
                             "$\\alpha_F$: exponential dependency of viscosity on the depletion "
                             "field $F$ (called peridotite). "
                             "Dimensionless factor. With a value of 0.0 (the default) the "
                             "viscosity does not depend on the depletion. The effective viscosity increase"
                             "due to depletion is defined as $exp( \\alpha_F * F)$. "
                             "Rationale: melting dehydrates the source rock by removing most of the volatiles,"
                             "and makes it stronger. Hirth and Kohlstedt (1996) report typical values around a "
                             "factor 100 to 1000 viscosity contrast between wet and dry rocks, although some "
                             "experimental studies report a smaller (factor 10) contrast (e.g. Fei et al., 2013).");
          prm.declare_entry ("Maximum Depletion viscosity change", "1.0e3",
                             Patterns::Double (0),
                             "$\\Delta \\eta_{F,max}$: maximum depletion strengthening of viscosity. "
                             "Rationale: melting dehydrates the source rock by removing most of the volatiles,"
                             "and makes it stronger. Hirth and Kohlstedt (1996) report typical values around a "
                             "factor 100 to 1000 viscosity contrast between wet and dry rocks, although some "
                             "experimental studies report a smaller (factor 10) contrast (e.g. Fei et al., 2013).");
/////// ////// /////////

          prm.declare_entry ("Use fractional melting", "false",
                             Patterns::Bool (),
                             "If fractional melting should be used (if true), including a solidus "
                             "change based on depletion (in this case, the amount of melt that has "
                             "migrated away from its origin), and freezing of melt when it has moved "
                             "to a region with temperatures lower than the solidus; or if batch "
                             "melting should be used (if false), assuming that the melt fraction only "
                             "depends on temperature and pressure, and how much melt has already been "
                             "generated at a given point, but not considering movement of melt in "
                             "the melting parameterization."
                             "\n\n"
                             "Note that melt does not freeze unless the 'Freezing rate' parameter is set "
                             "to a value larger than 0.");
          prm.declare_entry ("Freezing rate", "0.0",
                             Patterns::Double (0),
                             "Freezing rate of melt when in subsolidus regions. "
                             "If this parameter is set to a number larger than 0.0, it specifies the "
                             "fraction of melt that will freeze per year (or per second, depending on the "
                             "``Use years in output instead of seconds'' parameter), as soon as the porosity "
                             "exceeds the equilibrium melt fraction, and the equilibrium melt fraction "
                             "falls below the depletion. In this case, melt will freeze according to the "
                             "given rate until one of those conditions is not fulfilled anymore. The "
                             "reasoning behind this is that there should not be more melt present than "
                             "the equilibrium melt fraction, as melt production decreases with increasing "
                             "depletion, but the freezing process of melt also reduces the depletion by "
                             "the same amount, and as soon as the depletion falls below the equilibrium "
                             "melt fraction, we expect that material should melt again (no matter how "
                             "much melt is present). This is quite a simplification and not a realistic "
                             "freezing parameterization, but without tracking the melt composition, there "
                             "is no way to compute freezing rates accurately. "
                             "If this parameter is set to zero, no freezing will occur. "
                             "Note that freezing can never be faster than determined by the "
                             "``Melting time scale for operator splitting''. The product of the "
                             "``Freezing rate'' and the ``Melting time scale for operator splitting'' "
                             "defines how fast freezing occurs with respect to melting (if the "
                             "product is 0.5, melting will occur twice as fast as freezing). "
                             "Units: 1/yr or 1/s, depending on the ``Use years "
                             "in output instead of seconds'' parameter.");

          prm.declare_entry ("Melt extraction depth", "100000.0",
                             Patterns::Double(0),
                             "Depth above that melt will be extracted from the model, "
                             "which is done by a negative reaction term proportional to the "
                             "porosity field. "
                             "Units: $m$.");
/////// /////// ////////

          // parameters of anhydrous melting of peridotite, Katz 2003
           prm.declare_entry ("A1", "1085.7",
                             Patterns::Double (),
                             "Constant parameter in the quadratic "
                             "function that approximates the solidus "
                             "of peridotite. "
                             "Units: $\\degree C$.");
          prm.declare_entry ("A2", "1.329e-7",
                             Patterns::Double (),
                             "Prefactor of the linear pressure term "
                             "in the quadratic function that approximates "
                             "the solidus of peridotite. "
                             "Units: $\\degree C/Pa$.");
          prm.declare_entry ("A3", "-5.1e-18",
                             Patterns::Double (),
                             "Prefactor of the quadratic pressure term "
                             "in the quadratic function that approximates "
                             "the solidus of peridotite. "
                             "Units: $\\degree C/(Pa^2)$.");
          prm.declare_entry ("B1", "1475.0",
                             Patterns::Double (),
                             "Constant parameter in the quadratic "
                             "function that approximates the lherzolite "
                             "liquidus used for calculating the fraction "
                             "of peridotite-derived melt. "
                             "Units: $\\degree C$.");
          prm.declare_entry ("B2", "8.0e-8",
                             Patterns::Double (),
                             "Prefactor of the linear pressure term "
                             "in the quadratic function that approximates "
                             "the  lherzolite liquidus used for "
                             "calculating the fraction of peridotite-"
                             "derived melt. "
                             "Units: $\\degree C/Pa$.");
          prm.declare_entry ("B3", "-3.2e-18",
                             Patterns::Double (),
                             "Prefactor of the quadratic pressure term "
                             "in the quadratic function that approximates "
                             "the  lherzolite liquidus used for "
                             "calculating the fraction of peridotite-"
                             "derived melt. "
                             "Units: $\\degree C/(Pa^2)$.");
          prm.declare_entry ("C1", "1780.0",
                             Patterns::Double (),
                             "Constant parameter in the quadratic "
                             "function that approximates the liquidus "
                             "of peridotite. "
                             "Units: $\\degree C$.");
          prm.declare_entry ("C2", "4.50e-8",
                             Patterns::Double (),
                             "Prefactor of the linear pressure term "
                             "in the quadratic function that approximates "
                             "the liquidus of peridotite. "
                             "Units: $\\degree C/Pa$.");
          prm.declare_entry ("C3", "-2.0e-18",
                             Patterns::Double (),
                             "Prefactor of the quadratic pressure term "
                             "in the quadratic function that approximates "
                             "the liquidus of peridotite. "
                             "Units: $\\degree C/(Pa^2)$.");
          prm.declare_entry ("r1", "0.5",
                             Patterns::Double (),
                             "Constant in the linear function that "
                             "approximates the clinopyroxene reaction "
                             "coefficient. "
                             "Units: non-dimensional.");
          prm.declare_entry ("r2", "8e-11",
                             Patterns::Double (),
                             "Prefactor of the linear pressure term "
                             "in the linear function that approximates "
                             "the clinopyroxene reaction coefficient. "
                             "Units: $1/Pa$.");
          prm.declare_entry ("beta", "1.5",
                             Patterns::Double (),
                             "Exponent of the melting temperature in "
                             "the melt fraction calculation. "
                             "Units: non-dimensional.");
          prm.declare_entry ("Peridotite melting entropy change", "-300",
                             Patterns::Double (),
                             "The entropy change for the phase transition "
                             "from solid to melt of peridotite. "
                             "Units: $J/(kg K)$.");
          prm.declare_entry ("Mass fraction cpx", "0.15",
                             Patterns::Double (),
                             "Mass fraction of clinopyroxene in the "
                             "peridotite to be molten. "
                             "Units: non-dimensional.");

         //rheology parameters
          prm.declare_entry ("Grain size", "1e-3", Patterns::Double(0), "Units: $m$");
          prm.declare_entry ("Minimum strain rate", "1.4e-18", Patterns::Double(0),
                             "Stabilizes strain dependent viscosity. Units: $1 / s$");
          prm.declare_entry ("Maximum viscosity", "1e34", Patterns::Double(0),
                             "Upper cutoff for effective viscosity. Units: $Pa s$");
          prm.declare_entry ("Minimum viscosity", "1e17", Patterns::Double(0),
                             "Lower cutoff for effective viscosity. Units: $Pa s$");

          //difusion creep parameters
          prm.declare_entry ("Activation energie for diffusion creep", "300e3",
                             Patterns::List(Patterns::Double(0)),
                             "Aactivation energies, $E_a$, for background mantle Units: $J / mol$");
          prm.declare_entry ("Activation volume for diffusion creep", "6e-6",
                             Patterns::List(Patterns::Double(0)),
                             "List of activation volumes, $V_a$, for background mantle Units: $m^3 / mol$");
          prm.declare_entry ("Stress exponent for diffusion creep", "1",
                             Patterns::List(Patterns::Double(0)),
                             "$n_dislocation$, for background mantle and compositional fields, "
                             "Units: None");
          prm.declare_entry ("Grain size exponent for diffusion creep", "2.5",
                             Patterns::List(Patterns::Double(0)),
                             "$m_diffusion$, for background mantle "
                             "Units: None");
          prm.declare_entry ("Prefactor for diffusion creep", "0.11e-15",
                             Patterns::List(Patterns::Double(0)),
                             "Viscosity prefactor, $A$, for background mantle,  "
                             "Units: $Pa^{-n_{diffusion}} m^{n_{diffusion}/m_{diffusion}} s^{-1}$");

          //dislocation creep parameters
          prm.declare_entry ("Activation energie for dislocation creep", "430e3",
                             Patterns::List(Patterns::Double(0)),
                             "Aactivation energies, $E_a$, for background mantle Units: $J / mol$");
          prm.declare_entry ("Activation volume for dislocation creep", "10e-6",
                             Patterns::List(Patterns::Double(0)),
                             "List of activation volumes, $V_a$, for background mantle Units: $m^3 / mol$");
          prm.declare_entry ("Stress exponent for dislocation creep", "3",
                             Patterns::List(Patterns::Double(0)),
                             "Stress exponent, $n_dislocation$, for background mantle, "
                             "Units: None");
          prm.declare_entry ("Prefactor for dislocation creep", "1.1e-16",
                             Patterns::List(Patterns::Double(0)),
                             "Viscosity prefactor, $A$, for background mantle, "
                             "Units: $Pa^{-n_{dislocation}} m^{n_{dislocation}/m_{dislocation}} s^{-1}$");
          prm.declare_entry ("Water content", "1000",
                             Patterns::List(Patterns::Double(0)),
                             "Water conten for dry olivine, $E_a$, for background mantle Units: $J / mol$");
          prm.declare_entry ("Viscosity averaging scheme", "harmonic",
                             Patterns::Selection("arithmetic|harmonic|geometric|maximum composition"),
                             "When more than one compositional field is present at a point "
                             "with different viscosities, we need to come up with an average "
                             "viscosity at that point.  Select a weighted harmonic, arithmetic, "
                             "geometric, or maximum composition.");
          prm.declare_entry ("Reference compressibility", "4e-12",
                             Patterns::Double (0),
                             "The value of the reference compressibility. "
                             "Units: $1/Pa$.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    EdcMelting<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Edc melting");
        {
          reference_rho_s                   = prm.get_double ("Reference solid density");
          reference_rho_f                   = prm.get_double ("Reference melt density");
          reference_T                       = prm.get_double ("Reference temperature");
          eta_0                             = prm.get_double ("Reference shear viscosity");
          xi_0                              = prm.get_double ("Reference bulk viscosity");
          eta_f                             = prm.get_double ("Reference melt viscosity");
          reference_permeability            = prm.get_double ("Reference permeability");
          thermal_viscosity_exponent        = prm.get_double ("Thermal viscosity exponent");
          thermal_bulk_viscosity_exponent   = prm.get_double ("Thermal bulk viscosity exponent");
          thermal_conductivity              = prm.get_double ("Thermal conductivity");
          reference_specific_heat           = prm.get_double ("Reference specific heat");
          thermal_expansivity               = prm.get_double ("Thermal expansion coefficient");
          alpha_phi                         = prm.get_double ("Exponential melt weakening factor");
          depletion_density_change          = prm.get_double ("Depletion density change");
          surface_solidus                   = prm.get_double ("Surface solidus");
          depletion_solidus_change          = prm.get_double ("Depletion solidus change");
          pressure_solidus_change           = prm.get_double ("Pressure solidus change");
          compressibility                   = prm.get_double ("Solid compressibility");
          melt_compressibility              = prm.get_double ("Melt compressibility");
          model_is_compressible             = prm.get_bool ("Use full compressibility");
          include_melting_and_freezing      = prm.get_bool ("Include melting and freezing");
          melting_time_scale                = prm.get_double ("Melting time scale for operator splitting");
          alpha_depletion                   = prm.get_double ("Exponential depletion strengthening factor");
          delta_eta_depletion_max           = prm.get_double ("Maximum Depletion viscosity change");
/////// ////////// //////
          extraction_depth           = prm.get_double ("Melt extraction depth");
          fractional_melting         = prm.get_bool ("Use fractional melting");
          freezing_rate              = prm.get_double ("Freezing rate");
/////// /////// ////////
          //rheology parameters

          A1              = prm.get_double ("A1");
          A2              = prm.get_double ("A2");
          A3              = prm.get_double ("A3");
          B1              = prm.get_double ("B1");
          B2              = prm.get_double ("B2");
          B3              = prm.get_double ("B3");
          C1              = prm.get_double ("C1");
          C2              = prm.get_double ("C2");
          C3              = prm.get_double ("C3");
          r1              = prm.get_double ("r1");
          r2              = prm.get_double ("r2");
          beta            = prm.get_double ("beta");
          peridotite_melting_entropy_change
            = prm.get_double ("Peridotite melting entropy change");
          M_cpx           = prm.get_double ("Mass fraction cpx");

          //rheology parameters
          grain_size                      = prm.get_double("Grain size");
          min_strain_rate                 = prm.get_double("Minimum strain rate");
          ref_strain_rate                 = prm.get_double("Reference strain rate");
          max_visc                        = prm.get_double ("Maximum viscosity");
          min_visc                        = prm.get_double ("Minimum viscosity");

          //diffusion creep parameters
          activation_energie_diffusion    = prm.get_double ("Activation energie for diffusion creep");
          activation_volume_diffusion     = prm.get_double ("Activation volume for diffusion creep");
          stress_exponent_diffusion       = prm.get_double ("Stress exponent for diffusion creep");
          grain_size_exponent_diffusion   = prm.get_double ("Grain size exponent for diffusion creep");
          prefactor_diffusion             = prm.get_double ("Prefactor for diffusion creep");
          C_OH                            = prm.get_double ("Water content");

          //diffusion creep parameters;
          activation_energie_dislocation  = prm.get_double ("Activation energie for dislocation creep");
          activation_volume_dislocation   = prm.get_double ("Activation volume for dislocation creep");
          stress_exponent_dislocation     = prm.get_double ("Stress exponent for dislocation creep");
          prefactor_dislocation           = prm.get_double ("Prefactor for dislocation creep");
          C_OH                            = prm.get_double ("Water content");
//          reference_compressibility  = prm.get_double ("Reference compressibility");

          // Rheological parameters
          if (prm.get ("Viscosity averaging scheme") == "harmonic")
            viscosity_averaging = harmonic;
          else if (prm.get ("Viscosity averaging scheme") == "arithmetic")
            viscosity_averaging = arithmetic;
          else if (prm.get ("Viscosity averaging scheme") == "geometric")
            viscosity_averaging = geometric;
          else if (prm.get ("Viscosity averaging scheme") == "maximum composition")
            viscosity_averaging = maximum_composition;
          else
            AssertThrow(false, ExcMessage("Not a valid viscosity averaging scheme"));

          density_averaging = arithmetic;





          if (thermal_viscosity_exponent!=0.0 && reference_T == 0.0)
            AssertThrow(false, ExcMessage("Error: Material model Edc melting with Thermal viscosity exponent can not have reference_T=0."));

          if (this->convert_output_to_years() == true)
            melting_time_scale *= year_in_seconds;

          if (this->get_parameters().use_operator_splitting)
            {
              AssertThrow(melting_time_scale >= this->get_parameters().reaction_time_step,
                          ExcMessage("The reaction time step " + Utilities::to_string(this->get_parameters().reaction_time_step)
                                     + " in the operator splitting scheme is too large to compute melting rates! "
                                     "You have to choose it in such a way that it is smaller than the 'Melting time scale for "
                                     "operator splitting' chosen in the material model, which is currently "
                                     + Utilities::to_string(melting_time_scale) + "."));
              AssertThrow(melting_time_scale > 0,
                          ExcMessage("The Melting time scale for operator splitting must be larger than 0!"));
              AssertThrow(this->introspection().compositional_name_exists("porosity"),
                          ExcMessage("Material model Edc melting with melt transport only "
                                     "works if there is a compositional field called porosity."));
            }

          if (this->include_melt_transport())
            {
              AssertThrow(this->introspection().compositional_name_exists("porosity"),
                          ExcMessage("Material model Edc melting with melt transport only "
                                     "works if there is a compositional field called porosity."));
              if (include_melting_and_freezing)
                {
                  AssertThrow(this->introspection().compositional_name_exists("peridotite"),
                              ExcMessage("Material model Edc melting only works if there is a "
                                         "compositional field called peridotite."));
                }
            }

        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    EdcMelting<dim>::create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      if (this->get_parameters().use_operator_splitting
          && out.template get_additional_output<ReactionRateOutputs<dim> >() == nullptr)
        {
          const unsigned int n_points = out.viscosities.size();
          out.additional_outputs.push_back(
            std_cxx14::make_unique<MaterialModel::ReactionRateOutputs<dim>> (n_points, this->n_compositional_fields()));
        }
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(EdcMelting,
                                   "edc melting",
                                   "A material model that implements a simple formulation of the "
                                   "material parameters required for the modelling of melt transport, "
                                   "including a source term for the porosity according to a simplified "
                                   "linear melting model similar to \\cite{schmeling2006}:\n"
                                   "$\\phi_{\\text{equilibrium}} = \\frac{T-T_{\\text{sol}}}{T_{\\text{liq}}-T_{\\text{sol}}}$\n"
                                   "with "
                                   "$T_{\\text{sol}} = T_{\\text{sol,0}} + \\Delta T_p \\, p + \\Delta T_c \\, C$ \n"
                                   "$T_{\\text{liq}} = T_{\\text{sol}}  + \\Delta T_{\\text{sol-liq}}$.")
  }
}
