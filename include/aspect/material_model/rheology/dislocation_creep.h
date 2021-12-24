/*
  Copyright (C) 2019 by the authors of the ASPECT code.

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

#ifndef _aspect_material_model_rheology_dislocation_creep_h
#define _aspect_material_model_rheology_dislocation_creep_h

#include <aspect/global.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    namespace Rheology
    {
      template <int dim>
      class DislocationCreep : public ::aspect::SimulatorAccess<dim>
      {
        public:
          /**
           * Constructor.
           */
          DislocationCreep();

          /**
           * Declare the parameters this function takes through input files.
           */
          static
          void
          declare_parameters (ParameterHandler &prm);

          /**
           * Read the parameters from the parameter file.
           */
          void
          parse_parameters (ParameterHandler &prm);

          /**
           * Compute the viscosity based on the dislocation creep law.
           */
          double
          compute_viscosity (const double strain_rate,
                             const double pressure,
                             const double temperature,
                             const unsigned int composition) const;

        private:

          /**
           * List of dislocation creep prefactors A.
           */
          std::vector<double> prefactors_dislocation;

          /**
           * List of dislocation creep stress exponents n.
           */
          std::vector<double> stress_exponents_dislocation;

          /**
           * List of dislocation creep activation energies E.
           */
          std::vector<double> activation_energies_dislocation;

          /**
           * List of dislocation creep activation volumes V.
           */
          std::vector<double> activation_volumes_dislocation;

      };
    }
  }
}
#endif
