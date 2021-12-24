/*
  Copyright (C) 2015 - 2018 by the authors of the ASPECT code.

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

#ifndef _aspect_material_model_edc_melting_h
#define _aspect_material_model_edc_melting_h

#include <aspect/material_model/interface.h>
#include <aspect/geometry_model/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/postprocess/melt_statistics.h>
#include <aspect/melt.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    /**
     * A material model that implements a simple formulation of the
     * material parameters required for the modelling of melt transport
     * in a global model, including a source term for the porosity according
     * a simplified linear melting model.
     *
     * The model is considered incompressible, following the definition
     * described in Interface::is_compressible.
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class EdcMelting : public MaterialModel::MeltInterface<dim>, public ::aspect::SimulatorAccess<dim>, public MaterialModel::MeltFractionModel<dim>
    {
      public:
        /**
         * Return whether the model is compressible or not.  Incompressibility
         * does not necessarily imply that the density is constant; rather, it
         * may still depend on temperature or pressure. In the current
         * context, compressibility means whether we should solve the continuity
         * equation as $\nabla \cdot (\rho \mathbf u)=0$ (compressible Stokes)
         * or as $\nabla \cdot \mathbf{u}=0$ (incompressible Stokes).
         */
        bool is_compressible () const override;

        void evaluate(const typename Interface<dim>::MaterialModelInputs &in,
                      typename Interface<dim>::MaterialModelOutputs &out) const override;

        /**
         * Compute the equilibrium melt fractions for the given input conditions.
         * @p in and @p melt_fractions need to have the same size.
         *
         * @param in Object that contains the current conditions.
         * @param melt_fractions Vector of doubles that is filled with the
         * equilibrium melt fraction for each given input conditions.
         */
        void melt_fractions (const MaterialModel::MaterialModelInputs<dim> &in,
                             std::vector<double> &melt_fractions) const override;

        /**
         * @name Reference quantities
         * @{
         */
        double reference_viscosity () const override;

        double reference_darcy_coefficient () const override;


        /**
         * @}
         */

        /**
         * @name Functions used in dealing with run-time parameters
         * @{
         */
        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        void
        parse_parameters (ParameterHandler &prm) override;
        /**
         * @}
         */

        void
        create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const override;

        enum averaging_scheme
        {
          harmonic,
          arithmetic,
          geometric,
          maximum_composition
        } viscosity_averaging, density_averaging;

        /**
        * Thermal conductivity
        */
        double compositional_delta_rho;

        std::vector<double> compute_volume_fractions( const std::vector<double> &compositional_fields) const;

        double average_value (const std::vector<double> &composition,
                              const std::vector<double> &parameter_values,
                              const averaging_scheme &average_type) const;


        double diffusion_creep (const double &pressure,
                                const double &temperature) const;

        double dislocation_creep (const double &pressure,
                                  const double &temperature,
                                  const SymmetricTensor<2,dim> &strain_rate) const;


        double composite_viscosity (double &dislocation_creep,
                                    double &diffusion_creep) const;
        /**
         * Pointer to the material model used as the base model
         */
        std::unique_ptr<MaterialModel::Interface<dim> > base_model;

      private:
        double reference_rho_s;
        double reference_rho_f;
        double reference_T;
        double eta_0;
        double xi_0;
        double eta_f;
        double thermal_viscosity_exponent;
        double thermal_bulk_viscosity_exponent;
        double thermal_expansivity;
        double reference_specific_heat;
        double thermal_conductivity;
        double reference_permeability;
        double alpha_phi;
        double depletion_density_change;
        double depletion_solidus_change;
        double pressure_solidus_change;
        double surface_solidus;
        double compressibility;
        double melt_compressibility;
        bool include_melting_and_freezing;
        double melting_time_scale;
        double alpha_depletion;
        double extraction_depth;
        double fractional_melting;
        double delta_eta_depletion_max;
        double freezing_rate;
	bool model_is_compressible;
        double reference_rho;
        double eta;
        double composition_viscosity_prefactor;
        double thermal_alpha;

        /*
         * Rheology parameters
         */
        double grain_size;
        double activation_energie_diffusion;
        double activation_volume_diffusion;
        double stress_exponent_diffusion;
        double grain_size_exponent_diffusion;
        double activation_energie_dislocation;
        double activation_volume_dislocation;
        double stress_exponent_dislocation;
        double prefactor_diffusion;
        double prefactor_dislocation;
        double min_strain_rate;
        double ref_strain_rate;
        double min_visc;
        double max_visc;
        double C_OH;

        /**
         * Parameters for anhydrous melting of peridotite after Katz, 2003
         */

        // for the solidus temperature
        double A1;   // °C
        double A2; // °C/Pa
        double A3; // °C/(Pa^2)

        // for the lherzolite liquidus temperature
        double B1;   // °C
        double B2;   // °C/Pa
        double B3; // °C/(Pa^2)

        // for the liquidus temperature
        double C1;   // °C
        double C2;  // °C/Pa
        double C3; // °C/(Pa^2)

        // for the reaction coefficient of pyroxene
        double r1;     // cpx/melt
        double r2;     // cpx/melt/GPa
        double M_cpx;  // mass fraction of pyroxene

        // melt fraction exponent
        double beta;


        // entropy change upon melting
        double peridotite_melting_entropy_change;

        virtual
        double
        melt_fraction (const double temperature,
                       const double pressure) const;

    
        /**
          * Compute the change in entropy due to melting for a given @p temperature
          * and @p pressure, and under the assumption that a fraction
          * @p maximum_melt_fraction of the material has already been molten
          * previously. The entropy change is computed with respect to temperature
          * or pressure, depending on @p dependence.
          * This is needed to calculate the latent heat of melt.
          */
        virtual
        double
        entropy_change (const double temperature,
                        const double pressure,
                        const double maximum_melt_fraction,
                        const NonlinearDependence::Dependence dependence) const;

    };

  }
}

#endif
