/*
  Copyright (C) 2011 - 2016 by the authors of the ASPECT code.

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
  along with ASPECT; see the file doc/COPYING.  If not see
  <http://www.gnu.org/licenses/>.
*/


#ifndef _aspect_initial_composition_lithosphere_asthenosphere_h
#define _aspect_initial_composition_lithosphere_asthenosphere_h

#include <aspect/initial_composition/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/utilities.h>
#include <aspect/geometry_model/interface.h>

namespace aspect
{
  namespace InitialComposition
  {
    using namespace dealii;

    /**
     * A class that implements initial conditions for the compositional fields
     * based on a functional description provided in the input file.The porosity field
     * is implemented by computing the equilibrium melt fraction for the given initial
     * condition and reference pressure profile. Note that this plugin only
     * works if there is a compositional field called 'porosity', and the
     * used material model implements the 'MeltFractionModel' interface.
     * All compositional fields except porosity are not changed by this plugin.
     *
     * @ingroup InitialCompositionModels
     */
    template <int dim>
    class LithosphereAsthenosphere:  public Utilities::AsciiDataBoundary<dim>, public Interface<dim>
//    class LithosphereAsthenosphere:  public Utilities::AsciiDataInitial<dim>, public Interface<dim>
//    class LithosphereAsthenosphere : public Interface<dim>,  public Utilities::AsciiDataBase<dim>//, public SimulatorAccess<dim>
//    class LithosphereAsthenosphere : public Interface<dim>,  public SimulatorAccess<dim>
    {
      public:
        /**
         * Constructor.
         */
        // avoid -Woverloaded-virtual:
//        using Utilities::AsciiDataBoundary<dim>::initialize;

        LithosphereAsthenosphere ();

        void
        initialize ();

	using Utilities::AsciiDataBoundary<dim>::initialize;
//	using Utilities::AsciiDataInitial<dim>::initialize;
        /**
         * Return the initial composition as a function of position and number
         * of compositional field.
         */
        virtual
        double initial_composition (const Point<dim> &position, const unsigned int n_comp) const;

        /**
         * Declare the parameters this class takes through input files. The
         * default implementation of this function does not describe any
         * parameters. Consequently, derived classes do not have to overload
         * this function if they do not take any runtime parameters.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         * The default implementation of this function does not read any
         * parameters. Consequently, derived classes do not have to overload
         * this function if they do not take any runtime parameters.
         */
        virtual
        void
        parse_parameters (ParameterHandler &prm);

      private:

//        double  LAB_isotherm;
        types::boundary_id surface_boundary_id;
        double  moho;

        Utilities::AsciiDataBoundary<dim> ascii_data_lab;





    };
  }
}


#endif
