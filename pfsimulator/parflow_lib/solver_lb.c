/*BHEADER**********************************************************************
*
*  Copyright (c) 1995-2024, Lawrence Livermore National Security,
*  LLC. Produced at the Lawrence Livermore National Laboratory. Written
*  by the Parflow Team (see the CONTRIBUTORS file)
*  <parflow@lists.llnl.gov> CODE-OCEC-08-103. All rights reserved.
*
*  This file is part of Parflow. For details, see
*  http://www.llnl.gov/casc/parflow
*
*  Please read the COPYRIGHT file or Our Notice and the LICENSE file
*  for the GNU Lesser General Public License.
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License (as published
*  by the Free Software Foundation) version 2.1 dated February 1999.
*
*  This program is distributed in the hope that it will be useful, but
*  WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the terms
*  and conditions of the GNU General Public License for more details.
*
*  You should have received a copy of the GNU Lesser General Public
*  License along with this program; if not, write to the Free Software
*  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
*  USA
**********************************************************************EHEADER*/

/****************************************************************************
 *
 * Top level Multi-Phase Fractional-Flow solver module:
 *              solves pressure equation,
 *              advects saturations,
 *              advects contaminants.
 *           Implicit Saturation Explicit Pressure method is used.
 *
 *-----------------------------------------------------------------------------
 *
 *****************************************************************************/

#include "parflow.h"

#include <string.h>


/*--------------------------------------------------------------------------
 * Structures
 *--------------------------------------------------------------------------*/

typedef struct {
  PFModule          *discretize_pressure;
  PFModule          *diag_scale;
  PFModule          *linear_solver;
  PFModule          *permeability_face;
  PFModule          *phase_velocity_face;
  PFModule          *total_velocity_face;
  PFModule          *advect_satur;
  PFModule          *advect_concen;
  PFModule          *set_problem_data;

  Problem           *problem;

  int sadvect_order;
  int advect_order;
  double CFL;
  int max_iterations;
  double rel_tol;                            /* relative tolerance */
  double abs_tol;                            /* absolute tolerance */
  double drop_tol;                           /* drop tolerance */
  int print_subsurf_data;                    /* print permeability/porosity? */
  int print_press;                           /* print pressures? */
  int print_satur;                           /* print saturations? */
  int print_concen;                          /* print concentrations? */
  int print_wells;                           /* print well data? */

  int write_pdi_subsurf_data;                /* print subsurf data using PDI */
  int write_pdi_press;                       /* print pressure using PDI */
  int write_pdi_satur;                       /* print saturation using PDI */
  int write_pdi_concen;                      /* print concentration using PDI */
  // int write_pdi_wells;                       /* print wells using PDI */

/* New member variables added for diffusion solver */
/* double             CFL; */
/* int                max_iterations; */
  double beta_perm;                    /* permeability compressibility */
  double beta_fluid;                   /* fluid compressibility */
  double beta_pore;                    /* pore compressibility */
  double beta_fracture;                /* permeability compressibility */
  int comp_compress_flag;
} PublicXtra;

typedef struct {
  PFModule          *discretize_pressure;
  PFModule          *diag_scale;
  PFModule          *linear_solver;
  PFModule          *permeability_face;
  PFModule          *phase_velocity_face;
  PFModule          *total_velocity_face;
  PFModule          *advect_satur;
  PFModule          *advect_concen;
  PFModule          *set_problem_data;

  PFModule          *retardation;
  PFModule          *phase_density;
  PFModule          *phase_mobility;
  PFModule          *ic_phase_satur;
  PFModule          *ic_phase_concen;
  PFModule          *bc_phase_saturation;
  PFModule          *constitutive;

  Grid              *grid;
  Grid              *x_grid;
  Grid              *y_grid;
  Grid              *z_grid;

  ProblemData       *problem_data;

  double            *temp_data;
} InstanceXtra;


/*--------------------------------------------------------------------------
 * SolverDiffusion
 *--------------------------------------------------------------------------*/

void      SolverDiffusion()
{
  PFModule      *this_module = ThisPFModule;
  PublicXtra    *public_xtra = (PublicXtra*)PFModulePublicXtra(this_module);
  InstanceXtra  *instance_xtra = (InstanceXtra*)PFModuleInstanceXtra(this_module);

  Problem      *problem = (public_xtra->problem);

  int sadvect_order = (public_xtra->sadvect_order);
  int advect_order = (public_xtra->advect_order);
  double CFL = (public_xtra->CFL);
  int max_iterations = (public_xtra->max_iterations);
  double abs_tol = (public_xtra->abs_tol);
  double drop_tol = (public_xtra->drop_tol);
  int print_subsurf_data = (public_xtra->print_subsurf_data);
  int print_press = (public_xtra->print_press);
  int print_satur = (public_xtra->print_satur);
  int print_concen = (public_xtra->print_concen);
  int print_wells = (public_xtra->print_wells);

  PFModule     *discretize_pressure = (instance_xtra->discretize_pressure);
  PFModule     *diag_scale = (instance_xtra->diag_scale);
  PFModule     *linear_solver = (instance_xtra->linear_solver);
  PFModule     *permeability_face = (instance_xtra->permeability_face);
  PFModule     *phase_velocity_face = (instance_xtra->phase_velocity_face);
  PFModule     *total_velocity_face = (instance_xtra->total_velocity_face);
  PFModule     *advect_satur = (instance_xtra->advect_satur);
  PFModule     *phase_density = (instance_xtra->phase_density);
  PFModule     *advect_concen = (instance_xtra->advect_concen);
  PFModule     *set_problem_data = (instance_xtra->set_problem_data);

  PFModule     *retardation = (instance_xtra->retardation);
  PFModule     *phase_mobility = (instance_xtra->phase_mobility);
  PFModule     *ic_phase_satur = (instance_xtra->ic_phase_satur);
  PFModule     *ic_phase_concen = (instance_xtra->ic_phase_concen);
  PFModule     *bc_phase_saturation = (instance_xtra->bc_phase_saturation);
  PFModule     *constitutive = (instance_xtra->constitutive);

  ProblemData  *problem_data = (instance_xtra->problem_data);

  Grid         *grid = (instance_xtra->grid);
  Grid         *x_grid = (instance_xtra->x_grid);
  Grid         *y_grid = (instance_xtra->y_grid);
  Grid         *z_grid = (instance_xtra->z_grid);

  Vector       *temp_mobility_x = NULL;
  Vector       *temp_mobility_y = NULL;
  Vector       *temp_mobility_z = NULL;
  Vector       *stemp = NULL;
  Vector       *ctemp = NULL;

  int start_count = ProblemStartCount(problem);
  double start_time = ProblemStartTime(problem);
  double stop_time = ProblemStopTime(problem);
  double dump_interval = ProblemDumpInterval(problem);

  GrGeomSolid  *gr_domain;

  Vector       *pressure;
  Vector       *total_mobility_x;
  Vector       *total_mobility_y;
  Vector       *total_mobility_z;
  Vector      **saturations;
  Vector      **concentrations;

  Vector      **phase_x_velocity, **phase_y_velocity, **phase_z_velocity;
  Vector       *total_x_velocity, *total_y_velocity, *total_z_velocity;
  Vector       *z_permeability;

  Vector       *solidmassfactor;

  Matrix       *A;
  Vector       *f;

  EvalStruct   *eval_struct;

  int iteration_number, number_logged, file_number, dump_index;
  int indx, phase, concen;
  int transient, recompute_pressure, still_evolving, any_file_dumped;
  int dump_files, evolve_saturations, evolve_concentrations;
  int is_multiphase;

  double t;
  double dt;
  double       *phase_dt, min_phase_dt, total_dt, print_dt, well_dt, bc_dt;
  double phase_maximum, total_maximum;
  double dtmp, *phase_densities;

  VectorUpdateCommHandle   *handle;

  char dt_info;
  char file_prefix[PATH_MAX], file_postfix[64];

  double       *time_log = NULL, *dt_log = NULL;
  int          *seq_log = NULL, *dumped_log = NULL;
  char         *recomp_log = NULL, *dt_info_log = NULL;

  /* New variables added by B. Bosl */
  Lattice      *lattice;
  Vector       *pwork;
  CharVector   *cellType;
  int a, pressure_file_number;
  char name[128];


  is_multiphase = ProblemNumPhases(problem) > 1;

  /*-------------------------------------------------------------------
   * Set up temporary vectors
   *-------------------------------------------------------------------*/
  if (is_multiphase)
  {
    temp_mobility_x = NewVectorType(grid, 1, 1, vector_cell_centered);
    temp_mobility_y = NewVectorType(grid, 1, 1, vector_cell_centered);
    temp_mobility_z = NewVectorType(grid, 1, 1, vector_cell_centered);
    stemp = NewVectorType(grid, 1, 3, vector_cell_centered);
  }
  ctemp = NewVectorType(grid, 1, 3, vector_cell_centered);


  IfLogging(1)
  {
    seq_log = talloc(int, max_iterations + 1);
    time_log = talloc(double, max_iterations + 1);
    dt_log = talloc(double, max_iterations + 1);
    dt_info_log = talloc(char, max_iterations + 1);
    dumped_log = talloc(int, max_iterations + 1);
    recomp_log = talloc(char, max_iterations + 1);
    number_logged = 0;
  }

  sprintf(file_prefix, "%s", GlobalsOutFileName);

  /* do turning bands (and other stuff maybe) */
  PFModuleInvokeType(SetProblemDataInvoke, set_problem_data, (problem_data));
  gr_domain = ProblemDataGrDomain(problem_data);

  phase_densities = talloc(double, ProblemNumPhases(problem));

  for (phase = 0; phase < ProblemNumPhases(problem) - 1; phase++)
  {
    PFModuleInvokeType(PhaseDensityInvoke, phase_density,
                       (phase, NULL, NULL, &dtmp, &phase_densities[phase],
                        CALCFCN));
  }

  if (print_subsurf_data)
  {
    sprintf(file_postfix, "perm_x");
    WritePFBinary(file_prefix, file_postfix,
                  ProblemDataPermeabilityX(problem_data));

    sprintf(file_postfix, "perm_y");
    WritePFBinary(file_prefix, file_postfix,
                  ProblemDataPermeabilityY(problem_data));

    sprintf(file_postfix, "perm_y");
    WritePFBinary(file_prefix, file_postfix,
                  ProblemDataPermeabilityY(problem_data));

    sprintf(file_postfix, "porosity");
    WritePFBinary(file_prefix, file_postfix,
                  ProblemDataPorosity(problem_data));
  }

  if (public_xtra->write_pdi_subsurf_data)
  {
    sprintf(file_postfix, "perm_x");
    WritePDI(file_prefix, file_postfix, 0,
             ProblemDataPermeabilityX(problem_data), 0, 0);

    sprintf(file_postfix, "perm_y");
    WritePDI(file_prefix, file_postfix, 0,
             ProblemDataPermeabilityY(problem_data), 0, 0);

    sprintf(file_postfix, "perm_y");
    WritePDI(file_prefix, file_postfix, 0,
             ProblemDataPermeabilityY(problem_data), 0, 0);

    sprintf(file_postfix, "porosity");
    WritePDI(file_prefix, file_postfix, 0,
             ProblemDataPorosity(problem_data), 0, 0);
  }

  if (!amps_Rank(amps_CommWorld))
  {
    PrintWellData(ProblemDataWellData(problem_data), (WELLDATA_PRINTPHYSICAL | WELLDATA_PRINTVALUES));
  }

  /* Check to see if time dependency is requested and set flags */
  if (start_count < 0)
  {
    transient = 0;
    still_evolving = 0;
    recompute_pressure = 1;
    evolve_saturations = 0;
    evolve_concentrations = 0;
  }
  else
  {
    transient = 1;
    still_evolving = 1;
    recompute_pressure = 1;
    evolve_saturations = 1;
    evolve_concentrations = 1;
  }

/***********  New for lattice solver  ***********/
  transient = 0;
  still_evolving = 0;
  recompute_pressure = 1;
  evolve_saturations = 0;
  evolve_concentrations = 0;
  is_multiphase = 0;
/************************************************/

  t = start_time;

  pressure = NewVectorType(grid, 1, 1, vector_cell_centered);
  InitVectorAll(pressure, 0.0);

  total_mobility_x = NewVectorType(grid, 1, 1, vector_cell_centered);
  InitVectorAll(total_mobility_x, 0.0);

  total_mobility_y = NewVectorType(grid, 1, 1, vector_cell_centered);
  InitVectorAll(total_mobility_y, 0.0);

  total_mobility_z = NewVectorType(grid, 1, 1, vector_cell_centered);
  InitVectorAll(total_mobility_z, 0.0);

  /*-------------------------------------------------------------------
   * Allocate and set up initial saturations
   *-------------------------------------------------------------------*/

  saturations = ctalloc(Vector *, ProblemNumPhases(problem));

  if (is_multiphase)
  {
    for (phase = 0; phase < ProblemNumPhases(problem) - 1; phase++)
    {
      saturations[phase] = NewVectorType(grid, 1, 3, vector_cell_centered);
      InitVectorAll(saturations[phase], 0.0);

      PFModuleInvokeType(ICPhaseSaturInvoke,
                         ic_phase_satur,
                         (saturations[phase], phase, problem_data));
      PFModuleInvokeType(BCPhaseSaturationInvoke, bc_phase_saturation,
                         (saturations[phase], phase, gr_domain));

      handle = InitVectorUpdate(saturations[phase], VectorUpdateGodunov);
      FinalizeVectorUpdate(handle);
    }

    saturations[ProblemNumPhases(problem) - 1] = NewVectorType(grid, 1, 3, vector_cell_centered);
    InitVectorAll(saturations[ProblemNumPhases(problem) - 1], 0.0);

    PFModuleInvokeType(SaturationConstitutiveInvoke, constitutive, (saturations));

    handle = InitVectorUpdate(saturations[ProblemNumPhases(problem) - 1],
                              VectorUpdateGodunov);
    FinalizeVectorUpdate(handle);
  }
  else
  {
    saturations[0] = NULL;
  }

  /*-------------------------------------------------------------------
   * If (transient); initialize some stuff
   *-------------------------------------------------------------------*/

  if (transient)
  {
    iteration_number = file_number = start_count;
    dump_index = 1;
    if ((t >= stop_time) || (iteration_number > max_iterations))
    {
      recompute_pressure = 0;

      transient = 0;
      still_evolving = 0;
      evolve_saturations = 0;
      evolve_concentrations = 0;

      print_press = 0;
      print_satur = 0;
      print_concen = 0;
      print_wells = 0;
    }
    else
    {
      dt = 0.0;

      dump_files = 1;

      if (is_multiphase)
        eval_struct = NewEvalStruct(problem);

      solidmassfactor = NewVectorType(grid, 1, 2, vector_cell_centered);

      /*-------------------------------------------------------------------
       * Allocate and set up initial concentrations
       *-------------------------------------------------------------------*/

      concentrations = ctalloc(Vector *,
                               ProblemNumPhases(problem) * ProblemNumContaminants(problem));

      indx = 0;
      for (phase = 0; phase < ProblemNumPhases(problem); phase++)
      {
        for (concen = 0; concen < ProblemNumContaminants(problem); concen++)
        {
          concentrations[indx] = NewVectorType(grid, 1, 3, vector_cell_centered);
          InitVectorAll(concentrations[indx], 0.0);

          PFModuleInvokeType(ICPhaseConcenInvoke, ic_phase_concen,
                             (concentrations[indx], phase,
                              concen, problem_data));

          handle = InitVectorUpdate(concentrations[indx], VectorUpdateGodunov);
          FinalizeVectorUpdate(handle);
          indx++;
        }
      }

      /*-------------------------------------------------------------------
       * Allocate phase/total velocities and total mobility
       *-------------------------------------------------------------------*/

      phase_x_velocity = ctalloc(Vector *, ProblemNumPhases(problem));
      for (phase = 0; phase < ProblemNumPhases(problem); phase++)
      {
        phase_x_velocity[phase] = NewVectorType(x_grid, 1, 1, vector_cell_centered);
        InitVectorAll(phase_x_velocity[phase], 0.0);
      }

      phase_y_velocity = ctalloc(Vector *, ProblemNumPhases(problem));
      for (phase = 0; phase < ProblemNumPhases(problem); phase++)
      {
        phase_y_velocity[phase] = NewVectorType(y_grid, 1, 1, vector_cell_centered);
        InitVectorAll(phase_y_velocity[phase], 0.0);
      }

      phase_z_velocity = ctalloc(Vector *, ProblemNumPhases(problem));
      for (phase = 0; phase < ProblemNumPhases(problem); phase++)
      {
        phase_z_velocity[phase] = NewVectorType(z_grid, 1, 2, vector_cell_centered);
        InitVectorAll(phase_z_velocity[phase], 0.0);
      }

      phase_dt = ctalloc(double, ProblemNumPhases(problem));

      if (is_multiphase)
      {
        total_x_velocity = NewVectorType(x_grid, 1, 1, vector_cell_centered);
        InitVectorAll(total_x_velocity, 0.0);

        total_y_velocity = NewVectorType(y_grid, 1, 1, vector_cell_centered);
        InitVectorAll(total_y_velocity, 0.0);

        total_z_velocity = NewVectorType(z_grid, 1, 2, vector_cell_centered);
        InitVectorAll(total_z_velocity, 0.0);

        /*-------------------------------------------------------------------
         * Allocate and set up edge permeabilities
         *-------------------------------------------------------------------*/
        z_permeability = NewVectorType(z_grid, 1, 2, vector_cell_centered);
        InitVectorAll(z_permeability, 0.0);

        PFModuleInvokeType(PermeabilityFaceInvoke, permeability_face,
                           (z_permeability,
                            ProblemDataPermeabilityZ(problem_data)));
      }

      /*****************************************************************/
      /*          Print out any of the requested initial data          */
      /*****************************************************************/

      any_file_dumped = 0;

      /*-------------------------------------------------------------------
       * Print out the initial saturations?
       *-------------------------------------------------------------------*/

      if (print_satur && dump_files && is_multiphase)
      {
        for (phase = 0; phase < ProblemNumPhases(problem); phase++)
        {
          sprintf(file_postfix, "satur.%01d.%05d", phase, file_number);
          WritePFBinary(file_prefix, file_postfix, saturations[phase]);
        }
        any_file_dumped = 1;
      }

      if (public_xtra->write_pdi_satur && dump_files && is_multiphase)
      {
        for (phase = 0; phase < ProblemNumPhases(problem); phase++)
        {
          sprintf(file_postfix, "satur.%01d.%05d", phase, file_number);
          WritePDI(file_prefix, file_postfix, file_number, saturations[phase], 0, 0);
        }
        any_file_dumped = 1;
      }

      /*-------------------------------------------------------------------
       * Print out the initial concentrations?
       *-------------------------------------------------------------------*/

      if (print_concen && dump_files)
      {
        indx = 0;
        for (phase = 0; phase < ProblemNumPhases(problem); phase++)
        {
          for (concen = 0; concen < ProblemNumContaminants(problem); concen++)
          {
            sprintf(file_postfix, "concen.%01d.%02d.%05d", phase, concen,
                    file_number);
            WritePFSBinary(file_prefix, file_postfix, concentrations[indx],
                           drop_tol);
            indx++;
          }
        }
        any_file_dumped = 1;
      }

      if (public_xtra->write_pdi_concen && dump_files)
      {
        indx = 0;
        for (phase = 0; phase < ProblemNumPhases(problem); phase++)
        {
          for (concen = 0; concen < ProblemNumContaminants(problem); concen++)
          {
            sprintf(file_postfix, "concen.%01d.%02d.%05d", phase, concen,
                    file_number);
            WritePDI(file_prefix, file_postfix, file_number,
                     concentrations[indx], 1, drop_tol);
            indx++;
          }
        }
        any_file_dumped = 1;
      }

      /*-------------------------------------------------------------------
       * Print out the initial well data?
       *-------------------------------------------------------------------*/

      if (print_wells && dump_files)
      {
        WriteWells(file_prefix,
                   problem,
                   ProblemDataWellData(problem_data),
                   t,
                   WELLDATA_WRITEHEADER);
      }

      /*-------------------------------------------------------------------
       * Log this step
       *-------------------------------------------------------------------*/

      IfLogging(1)
      {
        seq_log[number_logged] = iteration_number;
        time_log[number_logged] = t;
        dt_log[number_logged] = dt;
        dt_info_log[number_logged] = 'i';
        if (any_file_dumped)
          dumped_log[number_logged] = file_number;
        else
          dumped_log[number_logged] = -1;
        if (recompute_pressure)
          recomp_log[number_logged] = 'y';
        else
          recomp_log[number_logged] = 'n';
        number_logged++;
      }

      if (any_file_dumped)
        file_number++;
    }
  }
  else
  {
    file_number = 1;

    /*-------------------------------------------------------------------
     * Print out the initial well data?
     *-------------------------------------------------------------------*/

    if (print_wells)
    {
      WriteWells(file_prefix,
                 problem,
                 ProblemDataWellData(problem_data),
                 t,
                 WELLDATA_WRITEHEADER);
    }
  }

  /***********************************************************************/
  /*                                                                     */
  /*                Begin the main computational section                 */
  /*                                                                     */
  /***********************************************************************/

  do
  {
    if (recompute_pressure)
    {
      /*******************************************************************/
      /*                  Compute the total mobility                     */
      /*******************************************************************/
      InitVectorAll(total_mobility_x, 0.0);
      InitVectorAll(total_mobility_y, 0.0);
      InitVectorAll(total_mobility_z, 0.0);

      PFModuleInvokeType(PhaseMobilityInvoke, phase_mobility,
                         (total_mobility_x,
                          total_mobility_y,
                          total_mobility_z,
                          ProblemDataPermeabilityX(problem_data),
                          ProblemDataPermeabilityY(problem_data),
                          ProblemDataPermeabilityZ(problem_data),
                          0,
                          saturations[0],
                          ProblemPhaseViscosity(problem, 0)));

      for (phase = 1; phase < ProblemNumPhases(problem); phase++)
      {
        /* Get the mobility of this phase */
        PFModuleInvokeType(PhaseMobilityInvoke, phase_mobility,
                           (temp_mobility_x,
                            temp_mobility_y,
                            temp_mobility_z,
                            ProblemDataPermeabilityX(problem_data),
                            ProblemDataPermeabilityY(problem_data),
                            ProblemDataPermeabilityZ(problem_data),
                            phase,
                            saturations[phase],
                            ProblemPhaseViscosity(problem, phase)));

        Axpy(1.0, temp_mobility_x, total_mobility_x);
        Axpy(1.0, temp_mobility_y, total_mobility_y);
        Axpy(1.0, temp_mobility_z, total_mobility_z);
      }

      handle = InitVectorUpdate(total_mobility_x, VectorUpdateAll);
      FinalizeVectorUpdate(handle);

      handle = InitVectorUpdate(total_mobility_y, VectorUpdateAll);
      FinalizeVectorUpdate(handle);

      handle = InitVectorUpdate(total_mobility_z, VectorUpdateAll);
      FinalizeVectorUpdate(handle);

      /*=========================================================================
       *  LATTICE BOLTZMANN SECTION
       *  =========================================================================*/

      /***********************************************************************/
      /*  Allocate memory                                                    */
      /***********************************************************************/
      sprintf(file_prefix, "%s", GlobalsOutFileName);

      lattice = ctalloc(Lattice, 1);

      cellType = NewCharVector(grid, 1, 1);
      InitCharVectorAll(cellType, 1);

      /***********************************************************************/
      /*  Compute permeability, porosity and initialize lattice variables    */
      /***********************************************************************/
      PFModuleInvokeType(SetProblemDataInvoke, set_problem_data, (problem_data));
      gr_domain = ProblemDataGrDomain(problem_data);

      pwork = NewVectorType(grid, 1, 1, vector_cell_centered);
      InitVectorAll(pwork, 0.0);

      /* Assign lattice variables to the lattice structure
       * for ease in passing arguments to functions. */
      (lattice->grid) = grid;
      (lattice->pressure) = pressure;
      /* Assume a scalar perm field given with a kx value of 1.0 */
      (lattice->perm) = ProblemDataPermeabilityX(problem_data);;
      (lattice->phi) = ProblemDataPorosity(problem_data);
      (lattice->pwork) = pwork;
      (lattice->cellType) = cellType;
      (lattice->cfl) = (public_xtra->CFL);
      (lattice->beta_perm) = (public_xtra->beta_perm);
      (lattice->beta_pore) = (public_xtra->beta_pore);
      (lattice->beta_fracture) = (public_xtra->beta_fracture);
      (lattice->beta_fluid) = (public_xtra->beta_fluid);
      (lattice->viscosity) = ProblemPhaseViscosity(problem, 0);
      (lattice->comp_compress_flag) = (public_xtra->comp_compress_flag);
      (lattice->t) = 0.0;
      (lattice->start) = start_time;
      (lattice->dump) = dump_interval;
      (lattice->stop) = stop_time;
      (lattice->density) = RHO;

      InitVectorAll((lattice->phi), 0.1);

      /***********************************************************************/
      /*                  Solve for the base pressure                        */
      /***********************************************************************/

      if (SOLVE_EQUILIBRIUM)
      {
        /* Discretize and solve the pressure equation */
        PFModuleInvokeType(DiscretizePressureInvoke,
                           discretize_pressure, (&A, &f, problem_data,
                                                 t,
                                                 total_mobility_x,
                                                 total_mobility_y,
                                                 total_mobility_z,
                                                 saturations));
        PFModuleInvokeType(MatrixDiagScaleInvoke, diag_scale, (pressure, A, f, DO));
        PFModuleReNewInstanceType(LinearSolverInitInstanceXtraInvoke,
                                  linear_solver, (NULL, NULL, problem_data, A, NULL));
        PFModuleInvokeType(LinearSolverInvoke,
                           linear_solver, (pressure, f, abs_tol, 0));
        PFModuleInvokeType(MatrixDiagScaleInvoke, diag_scale, (pressure, A, f, UNDO));
      }    /* end if(SOLVE_EQUILIBRIUM) */

      /***********************************************************************/
      /*                Begin the main computational section                 */
      /***********************************************************************/

      /* Print initial permeability, porosity and pressure */
      sprintf(file_postfix, "perm_x");
      WritePFBinary(file_prefix, file_postfix,
                    ProblemDataPermeabilityX(problem_data));

      sprintf(file_postfix, "perm_y");
      WritePFBinary(file_prefix, file_postfix,
                    ProblemDataPermeabilityY(problem_data));

      sprintf(file_postfix, "perm_z");
      WritePFBinary(file_prefix, file_postfix,
                    ProblemDataPermeabilityZ(problem_data));

      sprintf(file_postfix, "porosity");
      WritePFBinary(file_prefix, file_postfix, ProblemDataPorosity(problem_data));

      /* Initialize the pressure field and initialize */
      LBInitializeBC(lattice, problem, problem_data);
      LBWells(lattice, problem, problem_data);
      PerturbSystem(lattice, problem);
      LatticeFlowInit(lattice, problem);

      pressure_file_number = 0;
      sprintf(file_postfix, "press.%05d", pressure_file_number);
      strcpy(name, file_prefix);
      strcat(name, file_postfix);
      strcat(name, ".pfb");
      amps_Printf("Writing %s at t=%f\n", name, (lattice->t));
      WritePFBinary(file_prefix, file_postfix, pressure);

      DiffuseLB(lattice, problem, max_iterations, file_prefix);

      FreeVector(pwork);
      FreeCharVector(cellType);
      for (a = 0; a < nDirections; a++)
        free((lattice->e)[a]);
      free((lattice->e));
      free(lattice->Kscale);
      for (a = 0; a < 2; a++)
        free((lattice->Ktensor)[a]);
      free(lattice->Ktensor);
      free(lattice->bforce);
      free(lattice->c);
      free(lattice);


/*=========================================================================
 *  =========================================================================
 *  END OF LATTICE BOLTZMANN SECTION
 *  =========================================================================
 *  =========================================================================*/


      if (transient)
      {
        /***********************************************************************/
        /*                      Compute the velocities                         */
        /***********************************************************************/
        for (phase = 0; phase < ProblemNumPhases(problem); phase++)
        {
          /* Compute the velocity for this phase */
          PFModuleInvokeType(PhaseVelocityFaceInvoke,
                             phase_velocity_face,
                             (phase_x_velocity[phase],
                              phase_y_velocity[phase],
                              phase_z_velocity[phase],
                              problem_data,
                              pressure,
                              saturations,
                              phase,
                              t));

          phase_maximum = MaxPhaseFieldValue(phase_x_velocity[phase],
                                             phase_y_velocity[phase],
                                             phase_z_velocity[phase],
                                             ProblemDataPorosity(problem_data));

          phase_dt[phase] = CFL / phase_maximum;
          if (phase == 0)
          {
            min_phase_dt = phase_dt[0];
          }
          else
          {
            if (phase_dt[phase] < min_phase_dt)
            {
              min_phase_dt = phase_dt[phase];
            }
          }
        }

        /* Compute the total velocity */
        if (is_multiphase)
        {
          PFModuleInvokeType(TotalVelocityFaceInvoke,
                             total_velocity_face,
                             (total_x_velocity,
                              total_y_velocity,
                              total_z_velocity,
                              problem_data,
                              total_mobility_x,
                              total_mobility_y,
                              total_mobility_z,
                              pressure,
                              saturations));

          total_maximum = MaxTotalFieldValue(problem, eval_struct,
                                             saturations[0],
                                             total_x_velocity,
                                             total_y_velocity,
                                             total_z_velocity,
                                             z_permeability,
                                             ProblemDataPorosity(problem_data));

          total_dt = CFL / total_maximum;
        }

        /***********************************************************************/
        /*                       Print the pressure                            */
        /***********************************************************************/

        /* Dump the pressure values at this time-step */
        if (print_press)
        {
          sprintf(file_postfix, "press.%05d", file_number - 1);
          WritePFBinary(file_prefix, file_postfix, pressure);

#if 0
          for (phase = 0; phase < ProblemNumPhases(problem); phase++)
          {
            sprintf(file_postfix, "phasex.%01d.%05d", phase, file_number - 1);
            WritePFBinary(file_prefix, file_postfix, phase_x_velocity[phase]);

            sprintf(file_postfix, "phasey.%01d.%05d", phase, file_number - 1);
            WritePFBinary(file_prefix, file_postfix, phase_y_velocity[phase]);

            sprintf(file_postfix, "phasez.%01d.%05d", phase, file_number - 1);
            WritePFBinary(file_prefix, file_postfix, phase_z_velocity[phase]);
          }

          if (is_multiphase)
          {
            sprintf(file_postfix, "totalx.%05d", file_number - 1);
            WritePFBinary(file_prefix, file_postfix, total_x_velocity);

            sprintf(file_postfix, "totaly.%05d", file_number - 1);
            WritePFBinary(file_prefix, file_postfix, total_y_velocity);

            sprintf(file_postfix, "totalz.%05d", file_number - 1);
            WritePFBinary(file_prefix, file_postfix, total_z_velocity);
          }
#endif
        }

        if (public_xtra->write_pdi_press)
        {
          sprintf(file_postfix, "press.%05d", file_number - 1);
          WritePDI(file_prefix, file_postfix, file_number - 1, pressure, 0, 0);
        }
      }
    }

    if (transient)
    {
      iteration_number++;

      if (iteration_number >= max_iterations)
      {
        still_evolving = 0;
      }

      /*------------------------------------------------------------------
       *
       * Get delta t's for all wells and boundary conditions.
       *
       *------------------------------------------------------------------*/

      well_dt = TimeCycleDataComputeNextTransition(problem,
                                                   t,
                                                   WellDataTimeCycleData(ProblemDataWellData(problem_data)));

      bc_dt = TimeCycleDataComputeNextTransition(problem,
                                                 t,
                                                 BCPressureDataTimeCycleData(ProblemDataBCPressureData(problem_data)));

      /*------------------------------------------------------------------
       *
       * Compute the new dt value based on the dt's from the velocity
       * field computations and any time stepping criterion imposed
       * by the user.
       *
       * Modified by RDF
       *
       *------------------------------------------------------------------*/

      {
        if (is_multiphase)
        {
          dt = pfmin(total_dt, min_phase_dt);
        }
        else
        {
          dt = phase_dt[0];
        }

        if (well_dt < 0.0)
        {
          well_dt = dt;
        }

        if (bc_dt < 0.0)
        {
          bc_dt = dt;
        }

        /*--------------------------------------------------------------
         * Determine what needs to be computed/evolved on this
         * iteration.  Indicate what determined the value of `dt'.
         *--------------------------------------------------------------*/

        if (is_multiphase)
        {
          if (dt == total_dt)
          {
            recompute_pressure = 1;
            evolve_saturations = 1;
            evolve_concentrations = 1;

            dt_info = 's';
          }
          else if (dt > well_dt)
          {
            dt = well_dt;

            recompute_pressure = 1;
            evolve_saturations = 1;
            evolve_concentrations = 1;

            dt_info = 'w';
          }
          else if (dt > bc_dt)
          {
            dt = bc_dt;

            recompute_pressure = 1;
            evolve_saturations = 1;
            evolve_concentrations = 1;

            dt_info = 'b';
          }
          else
          {
            recompute_pressure = 0;
            evolve_saturations = 0;
            evolve_concentrations = 1;

            dt_info = 'c';
          }
        }
        else
        {
          if (dt > well_dt)
          {
            dt = well_dt;

            recompute_pressure = 1;
            evolve_concentrations = 1;

            dt_info = 'w';
          }
          else if (dt > bc_dt)
          {
            dt = bc_dt;

            recompute_pressure = 1;
            evolve_concentrations = 1;

            dt_info = 'b';
          }
          else
          {
            recompute_pressure = 0;
            evolve_concentrations = 1;

            dt_info = 'c';
          }
        }

        /*--------------------------------------------------------------
         * If we are printing out results, then determine if we need
         * to print them after this time step.
         *
         * If we are dumping output at real time intervals, the value
         * of dt may be changed.  If this happens, we want to
         * compute/evolve all values.  We also set `dump_info' to `p'
         * to indicate that the dump interval decided the time step for
         * this iteration.
         *--------------------------------------------------------------*/

        if (print_press || print_satur || print_concen || print_wells ||
            public_xtra->write_pdi_press ||
            public_xtra->write_pdi_satur ||
            public_xtra->write_pdi_concen)
        {
          dump_files = 0;

          if (dump_interval > 0)
          {
            print_dt = start_time + dump_index * dump_interval - t;

            if (dt >= print_dt)
            {
              dt = print_dt;

              if (is_multiphase)
              {
                recompute_pressure = 1;
                evolve_saturations = 1;
              }
              else
              {
                recompute_pressure = 0;
              }

              evolve_concentrations = 1;

              dt_info = 'p';

              dump_files = 1;
              dump_index++;
            }
          }
          else if (dump_interval < 0)
          {
            if ((iteration_number % (-(int)dump_interval)) == 0)
            {
              dump_files = 1;
            }
          }
          else
          {
            dump_files = 0;
          }
        }

        /*--------------------------------------------------------------
         * If this is the last iteration, set appropriate variables.
         *--------------------------------------------------------------*/

        if ((t + dt) >= stop_time)
        {
          still_evolving = 0;

          dt = stop_time - t;
          dt_info = 'f';
          recompute_pressure = 0;
          evolve_saturations = 1;
          evolve_concentrations = 1;
        }

        t += dt;
      }

      any_file_dumped = 0;

      /***********************************************************************/
      /*              Solve for and print the saturations                    */
      /***********************************************************************/

      /* Solve for the saturation values at this time-step if necessary. */
      if (is_multiphase)
      {
        if (evolve_saturations)
        {
          for (phase = 0; phase < ProblemNumPhases(problem) - 1; phase++)
          {
            InitVectorAll(stemp, 0.0);
            Copy(saturations[phase], stemp);
            PFModuleInvokeType(BCPhaseSaturationInvoke, bc_phase_saturation,
                               (stemp, phase, gr_domain));


            /* Evolve to the new time */
            PFModuleInvokeType(AdvectionSaturationInvoke, advect_satur,
                               (problem_data, phase,
                                stemp, saturations[phase],
                                total_x_velocity,
                                total_y_velocity,
                                total_z_velocity,
                                z_permeability,
                                ProblemDataPorosity(problem_data),
                                ProblemPhaseViscosities(problem),
                                phase_densities,
                                ProblemGravity(problem),
                                t, dt, sadvect_order));
          }
          InitVectorAll(saturations[ProblemNumPhases(problem) - 1], 0.0);
          PFModuleInvokeType(SaturationConstitutiveInvoke, constitutive, (saturations));

          handle = InitVectorUpdate(saturations[ProblemNumPhases(problem) - 1],
                                    VectorUpdateGodunov);
          FinalizeVectorUpdate(handle);

          /* Print the saturation values at this time-step? */
          if (print_satur && dump_files)
          {
            for (phase = 0; phase < ProblemNumPhases(problem); phase++)
            {
              sprintf(file_postfix, "satur.%01d.%05d", phase, file_number);
              WritePFBinary(file_prefix, file_postfix, saturations[phase]);
            }
            any_file_dumped = 1;
          }

          if (public_xtra->write_pdi_satur && dump_files)
          {
            for (phase = 0; phase < ProblemNumPhases(problem); phase++)
            {
              sprintf(file_postfix, "satur.%01d.%05d", phase, file_number);
              WritePDI(file_prefix, file_postfix, file_number, saturations[phase], 0, 0);
            }
            any_file_dumped = 1;
          }
        }
      }

      /***********************************************************************/
      /*            Solve for and print the concentrations                   */
      /***********************************************************************/

      if (evolve_concentrations)
      {
        /* Solve for the concentration values at this time-step. */
        indx = 0;
        for (phase = 0; phase < ProblemNumPhases(problem); phase++)
        {
          for (concen = 0; concen < ProblemNumContaminants(problem); concen++)
          {
            PFModuleInvokeType(RetardationInvoke, retardation,
                               (solidmassfactor,
                                concen,
                                problem_data));
            handle = InitVectorUpdate(solidmassfactor, VectorUpdateAll2);
            FinalizeVectorUpdate(handle);

            InitVectorAll(ctemp, 0.0);
            Copy(concentrations[indx], ctemp);

            PFModuleInvokeType(AdvectionConcentrationInvoke, advect_concen,
                               (problem_data, phase, concen,
                                ctemp, concentrations[indx],
                                phase_x_velocity[phase],
                                phase_y_velocity[phase],
                                phase_z_velocity[phase],
                                solidmassfactor,
                                t, dt, advect_order));
            indx++;
          }
        }

        /* Print the concentration values at this time-step? */
        if (print_concen && dump_files)
        {
          indx = 0;
          for (phase = 0; phase < ProblemNumPhases(problem); phase++)
          {
            for (concen = 0; concen < ProblemNumContaminants(problem); concen++)
            {
              sprintf(file_postfix, "concen.%01d.%02d.%05d", phase, concen,
                      file_number);
              WritePFSBinary(file_prefix, file_postfix,
                             concentrations[indx], drop_tol);
              indx++;
            }
          }
          any_file_dumped = 1;
        }

        if (public_xtra->write_pdi_concen && dump_files)
        {
          indx = 0;
          for (phase = 0; phase < ProblemNumPhases(problem); phase++)
          {
            for (concen = 0; concen < ProblemNumContaminants(problem); concen++)
            {
              sprintf(file_postfix, "concen.%01d.%02d.%05d", phase, concen,
                      file_number);
              WritePDI(file_prefix, file_postfix, file_number,
                       concentrations[indx], 1, drop_tol);
              indx++;
            }
          }
          any_file_dumped = 1;
        }
      }

      /**********************************************************************/
      /*                      Print the Well Data                           */
      /**********************************************************************/

      if (evolve_concentrations || evolve_saturations)
      {
        /* Print the well data values at this time-step? */
        if (print_wells && dump_files)
        {
          WriteWells(file_prefix,
                     problem,
                     ProblemDataWellData(problem_data),
                     t,
                     WELLDATA_DONTWRITEHEADER);
        }
      }

      /*-------------------------------------------------------------------
       * Log this step
       *-------------------------------------------------------------------*/

      IfLogging(1)
      {
        seq_log[number_logged] = iteration_number;
        time_log[number_logged] = t;
        dt_log[number_logged] = dt;
        dt_info_log[number_logged] = dt_info;
        if (any_file_dumped)
          dumped_log[number_logged] = file_number;
        else
          dumped_log[number_logged] = -1;
        if (recompute_pressure)
          recomp_log[number_logged] = 'y';
        else
          recomp_log[number_logged] = 'n';
        number_logged++;
      }

      if (any_file_dumped)
        file_number++;
    }
    else
    {
      if (print_press)
      {
        sprintf(file_postfix, "press");
        WritePFBinary(file_prefix, file_postfix, pressure);
      }

      if (public_xtra->write_pdi_press)
      {
        sprintf(file_postfix, "press");
        WritePDI(file_prefix, file_postfix, 0, pressure, 0, 0);
      }
    }
  }
  while (still_evolving);

  if (transient)
  {
    free(phase_dt);

    if (is_multiphase)
    {
      FreeVector(z_permeability);

      FreeVector(total_z_velocity);
      FreeVector(total_y_velocity);
      FreeVector(total_x_velocity);

      FreeEvalStruct(eval_struct);
    }

    for (phase = 0; phase < ProblemNumPhases(problem); phase++)
    {
      FreeVector(phase_z_velocity[phase]);
    }
    tfree(phase_z_velocity);

    for (phase = 0; phase < ProblemNumPhases(problem); phase++)
    {
      FreeVector(phase_y_velocity[phase]);
    }
    tfree(phase_y_velocity);

    for (phase = 0; phase < ProblemNumPhases(problem); phase++)
    {
      FreeVector(phase_x_velocity[phase]);
    }
    tfree(phase_x_velocity);

    indx = 0;
    for (phase = 0; phase < ProblemNumPhases(problem); phase++)
    {
      for (concen = 0; concen < ProblemNumContaminants(problem); concen++)
      {
        FreeVector(concentrations[indx]);
        indx++;
      }
    }
    tfree(concentrations);

    FreeVector(solidmassfactor);
  }

  if (is_multiphase)
  {
    for (phase = 0; phase < ProblemNumPhases(problem); phase++)
    {
      FreeVector(saturations[phase]);
    }
  }
  tfree(saturations);
  tfree(phase_densities);

  FreeVector(total_mobility_x);
  FreeVector(total_mobility_y);
  FreeVector(total_mobility_z);
  FreeVector(pressure);

  if (!amps_Rank(amps_CommWorld))
  {
    PrintWellData(ProblemDataWellData(problem_data), (WELLDATA_PRINTSTATS));
  }

  /*-----------------------------------------------------------------------
   * Print log
   *-----------------------------------------------------------------------*/

  IfLogging(1)
  {
    FILE*  log_file;
    int k;

    log_file = OpenLogFile("SolverMulti");

    if (transient)
    {
      fprintf(log_file, "Transient Problem Solved.\n");
      fprintf(log_file, "-------------------------\n");
      fprintf(log_file, "Sequence #       Time         \\Delta t         Dumpfile #   Recompute?\n");
      fprintf(log_file, "----------   ------------   ------------ -     ----------   ----------\n");

      for (k = 0; k < number_logged; k++)
      {
        if (dumped_log[k] == -1)
          fprintf(log_file, "  %06d     %8e   %8e %1c                       %1c\n",
                  k, time_log[k], dt_log[k], dt_info_log[k], recomp_log[k]);
        else
          fprintf(log_file, "  %06d     %8e   %8e %1c       %06d          %1c\n",
                  k, time_log[k], dt_log[k], dt_info_log[k], dumped_log[k], recomp_log[k]);
      }
    }
    else
    {
      fprintf(log_file, "Non-Transient Problem Solved.\n");
      fprintf(log_file, "-----------------------------\n");
    }

    CloseLogFile(log_file);

    tfree(seq_log);
    tfree(time_log);
    tfree(dt_log);
    tfree(dt_info_log);
    tfree(dumped_log);
    tfree(recomp_log);
  }
}

/*--------------------------------------------------------------------------
 * SolverDiffusionInitInstanceXtra
 *--------------------------------------------------------------------------*/

PFModule *SolverDiffusionInitInstanceXtra()
{
  PFModule      *this_module = ThisPFModule;
  PublicXtra    *public_xtra = (PublicXtra*)PFModulePublicXtra(this_module);
  InstanceXtra  *instance_xtra;

  Problem      *problem = (public_xtra->problem);

  Grid         *grid;
  Grid         *grid2d = NULL;
  Grid         *x_grid;
  Grid         *y_grid;
  Grid         *z_grid;

  SubgridArray *new_subgrids;
  SubgridArray *all_subgrids, *new_all_subgrids;

  Subgrid      *subgrid, *new_subgrid;

  double       *temp_data, *temp_data_placeholder;
  int total_mobility_sz, pressure_sz, velocity_sz, satur_sz = 0,
      concen_sz, temp_data_size, sz;
  int is_multiphase;

  int i;


  is_multiphase = ProblemNumPhases(problem) > 1;

  if (PFModuleInstanceXtra(this_module) == NULL)
    instance_xtra = ctalloc(InstanceXtra, 1);
  else
    instance_xtra = (InstanceXtra*)PFModuleInstanceXtra(this_module);

  /*-------------------------------------------------------------------
   * Create the grids
   *-------------------------------------------------------------------*/

  /* Create the flow grid */
  grid = CreateGrid(GlobalsUserGrid);

  /* Create the x velocity grid */

  all_subgrids = GridAllSubgrids(grid);

  /***** Set up a new subgrid grown by one in the x-direction *****/

  new_all_subgrids = NewSubgridArray();
  ForSubgridI(i, all_subgrids)
  {
    subgrid = SubgridArraySubgrid(all_subgrids, i);
    new_subgrid = DuplicateSubgrid(subgrid);
    SubgridNX(new_subgrid) += 1;
    AppendSubgrid(new_subgrid, new_all_subgrids);
  }
  new_subgrids = GetGridSubgrids(new_all_subgrids);
  x_grid = NewGrid(new_subgrids, new_all_subgrids);
  CreateComputePkgs(x_grid);

  /* Create the y velocity grid */

  all_subgrids = GridAllSubgrids(grid);

  /***** Set up a new subgrid grown by one in the y-direction *****/

  new_all_subgrids = NewSubgridArray();
  ForSubgridI(i, all_subgrids)
  {
    subgrid = SubgridArraySubgrid(all_subgrids, i);
    new_subgrid = DuplicateSubgrid(subgrid);
    SubgridNY(new_subgrid) += 1;
    AppendSubgrid(new_subgrid, new_all_subgrids);
  }
  new_subgrids = GetGridSubgrids(new_all_subgrids);
  y_grid = NewGrid(new_subgrids, new_all_subgrids);
  CreateComputePkgs(y_grid);

  /* Create the z velocity grid */

  all_subgrids = GridAllSubgrids(grid);

  /***** Set up a new subgrid grown by one in the z-direction *****/

  new_all_subgrids = NewSubgridArray();
  ForSubgridI(i, all_subgrids)
  {
    subgrid = SubgridArraySubgrid(all_subgrids, i);
    new_subgrid = DuplicateSubgrid(subgrid);
    SubgridNZ(new_subgrid) += 1;
    AppendSubgrid(new_subgrid, new_all_subgrids);
  }
  new_subgrids = GetGridSubgrids(new_all_subgrids);
  z_grid = NewGrid(new_subgrids, new_all_subgrids);
  CreateComputePkgs(z_grid);

  (instance_xtra->grid) = grid;
  (instance_xtra->x_grid) = x_grid;
  (instance_xtra->y_grid) = y_grid;
  (instance_xtra->z_grid) = z_grid;

  /*-------------------------------------------------------------------
   * Create problem_data
   *-------------------------------------------------------------------*/

  (instance_xtra->problem_data) = NewProblemData(grid, grid2d);

  /*-------------------------------------------------------------------
   * Initialize module instances
   *-------------------------------------------------------------------*/

  if (PFModuleInstanceXtra(this_module) == NULL)
  {
    (instance_xtra->discretize_pressure) =
      PFModuleNewInstanceType(DiscretizePressureInitInstanceXtraInvoke,
                              (public_xtra->discretize_pressure),
                              (problem, grid, NULL));
    (instance_xtra->diag_scale) =
      PFModuleNewInstanceType(MatrixDiagScaleInitInstanceXtraInvoke,
                              (public_xtra->diag_scale),
                              (grid));
    (instance_xtra->linear_solver) =
      PFModuleNewInstanceType(LinearSolverInitInstanceXtraInvoke,
                              (public_xtra->linear_solver),
                              (problem, grid, NULL, NULL, NULL));
    (instance_xtra->phase_velocity_face) =
      PFModuleNewInstanceType(PhaseVelocityFaceInitInstanceXtraInvoke,
                              (public_xtra->phase_velocity_face),
                              (problem, grid, x_grid, y_grid, z_grid, NULL));
    (instance_xtra->advect_concen) =
      PFModuleNewInstanceType(AdvectionConcentrationInitInstanceXtraType,
                              (public_xtra->advect_concen),
                              (problem, grid, NULL));
    (instance_xtra->set_problem_data) =
      PFModuleNewInstanceType(SetProblemDataInitInstanceXtraInvoke,
                              (public_xtra->set_problem_data),
                              (problem, grid, NULL, NULL));
    (instance_xtra->retardation) =
      PFModuleNewInstanceType(RetardationInitInstanceXtraInvoke,
                              ProblemRetardation(problem), (NULL));
    (instance_xtra->phase_mobility) =
      PFModuleNewInstance(ProblemPhaseMobility(problem), ());
    (instance_xtra->ic_phase_concen) =
      PFModuleNewInstance(ProblemICPhaseConcen(problem), ());
    (instance_xtra->phase_density) =
      PFModuleNewInstance(ProblemPhaseDensity(problem), ());

    if (is_multiphase)
    {
      (instance_xtra->permeability_face) =
        PFModuleNewInstanceType(PermeabilityFaceInitInstanceXtraInvoke,
                                (public_xtra->permeability_face),
                                (z_grid));
      (instance_xtra->total_velocity_face) =
        PFModuleNewInstanceType(TotalVelocityFaceInitInstanceXtraInvoke,
                                (public_xtra->total_velocity_face),
                                (problem, grid, x_grid, y_grid, z_grid, NULL));
      (instance_xtra->advect_satur) =
        PFModuleNewInstanceType(AdvectionSaturationInitInstanceXtraInvoke,
                                (public_xtra->advect_satur),
                                (problem, grid, NULL));
      (instance_xtra->ic_phase_satur) =
        PFModuleNewInstance(ProblemICPhaseSatur(problem), ());
      (instance_xtra->bc_phase_saturation) =
        PFModuleNewInstance(ProblemBCPhaseSaturation(problem), ());
      (instance_xtra->constitutive) =
        PFModuleNewInstanceType(SaturationConstitutiveInitInstanceXtraInvoke,
                                ProblemSaturationConstitutive(problem), (grid));
    }
  }
  else
  {
    PFModuleReNewInstanceType(DiscretizePressureInitInstanceXtraInvoke,
                              (instance_xtra->discretize_pressure),
                              (problem, grid, NULL));
    PFModuleReNewInstanceType(MatrixDiagScaleInitInstanceXtraInvoke, (instance_xtra->diag_scale),
                              (grid));
    PFModuleReNewInstanceType(LinearSolverInitInstanceXtraInvoke,
                              (instance_xtra->linear_solver),
                              (problem, grid, NULL, NULL, NULL));
    PFModuleReNewInstanceType(PhaseVelocityFaceInitInstanceXtraInvoke,
                              (instance_xtra->phase_velocity_face),
                              (problem, grid, x_grid, y_grid, z_grid, NULL));
    PFModuleReNewInstanceType(AdvectionConcentrationInitInstanceXtraType,
                              (instance_xtra->advect_concen),
                              (problem, grid, NULL));
    PFModuleReNewInstanceType(SetProblemDataInitInstanceXtraInvoke,
                              (instance_xtra->set_problem_data),
                              (problem, grid, NULL, NULL));
    PFModuleReNewInstanceType(RetardationInitInstanceXtraInvoke,
                              (instance_xtra->retardation), (NULL));
    PFModuleReNewInstance((instance_xtra->phase_mobility), ());
    PFModuleReNewInstance((instance_xtra->ic_phase_concen), ());
    PFModuleReNewInstance((instance_xtra->phase_density), ());

    if (is_multiphase)
    {
      PFModuleReNewInstanceType(PermeabilityFaceInitInstanceXtraInvoke,
                                (instance_xtra->permeability_face),
                                (z_grid));
      PFModuleReNewInstanceType(TotalVelocityFaceInitInstanceXtraInvoke,
                                (instance_xtra->total_velocity_face),
                                (problem, grid, x_grid, y_grid, z_grid, NULL));
      PFModuleReNewInstanceType(AdvectionSaturationInitInstanceXtraInvoke,
                                (instance_xtra->advect_satur),
                                (problem, grid, NULL));
      PFModuleReNewInstance((instance_xtra->ic_phase_satur), ());
      PFModuleReNewInstance((instance_xtra->bc_phase_saturation), ());
      PFModuleReNewInstanceType(SaturationConstitutiveInitInstanceXtraInvoke,
                                (instance_xtra->constitutive), (grid));
    }
  }

  /*-------------------------------------------------------------------
   * Set up temporary data
   *-------------------------------------------------------------------*/

  if (is_multiphase)
  {
    /* compute size for total mobility computation */
    sz = 0;
    total_mobility_sz = sz;

    /* compute size for saturation advection */
    sz = 0;
    sz += PFModuleSizeOfTempData(instance_xtra->advect_satur);
    satur_sz = sz;
  }

  /* compute size for pressure solve */
  sz = 0;
  sz = pfmax(sz, PFModuleSizeOfTempData(instance_xtra->discretize_pressure));
  sz = pfmax(sz, PFModuleSizeOfTempData(instance_xtra->linear_solver));
  pressure_sz = sz;

  /* compute size for velocity computation */
  sz = 0;
  sz = pfmax(sz, PFModuleSizeOfTempData(instance_xtra->phase_velocity_face));
  if (is_multiphase)
  {
    sz = pfmax(sz, PFModuleSizeOfTempData(instance_xtra->total_velocity_face));
  }
  velocity_sz = sz;

  /* compute size for concentration advection */
  sz = 0;
  sz = pfmax(sz, PFModuleSizeOfTempData(instance_xtra->retardation));
  sz += pfmax(sz, PFModuleSizeOfTempData(instance_xtra->advect_concen));
  concen_sz = sz;

  /* set temp_data size to max of pressure_sz, satur_sz, and concen_sz*/
  temp_data_size = pfmax(pfmax(pressure_sz, velocity_sz), concen_sz);
  if (is_multiphase)
  {
    temp_data_size = pfmax(temp_data_size, pfmax(total_mobility_sz, satur_sz));
  }
/*     temp_data_size = total_mobility_sz + pressure_sz + velocity_sz
 *                      + satur_sz + concen_sz;  */

  /* allocate temporary data */
  temp_data = NewTempData(temp_data_size);
  (instance_xtra->temp_data) = temp_data;


  /* renew set_problem_data module */
  PFModuleReNewInstanceType(SetProblemDataInitInstanceXtraInvoke,
                            (instance_xtra->set_problem_data),
                            (NULL, NULL, NULL, temp_data));

  /* renew pressure solve modules that take temporary data */
  PFModuleReNewInstanceType(DiscretizePressureInitInstanceXtraInvoke,
                            (instance_xtra->discretize_pressure),
                            (NULL, NULL, temp_data));
/*   temp_data += PFModuleSizeOfTempData(instance_xtra -> discretize_pressure);  */
  PFModuleReNewInstanceType(LinearSolverInitInstanceXtraInvoke,
                            (instance_xtra->linear_solver),
                            (NULL, NULL, NULL, NULL, temp_data));
/*   temp_data += PFModuleSizeOfTempData(instance_xtra -> linear_solver);  */

  /* renew velocity computation modules that take temporary data */
  PFModuleReNewInstanceType(PhaseVelocityFaceInitInstanceXtraInvoke,
                            (instance_xtra->phase_velocity_face),
                            (NULL, NULL, NULL, NULL, NULL, temp_data));
  if (is_multiphase)
  {
    PFModuleReNewInstanceType(TotalVelocityFaceInitInstanceXtraInvoke,
                              (instance_xtra->total_velocity_face),
                              (NULL, NULL, NULL, NULL, NULL, temp_data));
    /* temp_data += PFModuleSizeOfTempData(instance_xtra ->
     *                                      total_velocity_face);  */

    /* renew saturation advection modules that take temporary data */
    temp_data_placeholder = temp_data;
    PFModuleReNewInstanceType(AdvectionSaturationInitInstanceXtraInvoke,
                              (instance_xtra->advect_satur),
                              (NULL, NULL, temp_data_placeholder));
    temp_data_placeholder += PFModuleSizeOfTempData(instance_xtra->advect_satur);
  }

  /* renew concentration advection modules that take temporary data */
  temp_data_placeholder = temp_data;
  PFModuleReNewInstanceType(RetardationInitInstanceXtraInvoke,
                            (instance_xtra->retardation),
                            (temp_data_placeholder));
  PFModuleReNewInstanceType(AdvectionConcentrationInitInstanceXtraType,
                            (instance_xtra->advect_concen),
                            (NULL, NULL, temp_data_placeholder));
  int size_retardation = PFModuleSizeOfTempData(instance_xtra->retardation);
  int size_advect = PFModuleSizeOfTempData(instance_xtra->advect_concen);

  temp_data_placeholder += pfmax(size_retardation,
                                 size_advect
                                 );

  temp_data += temp_data_size;

  PFModuleInstanceXtra(this_module) = instance_xtra;
  return this_module;
}

/*--------------------------------------------------------------------------
 * SolverDiffusionFreeInstanceXtra
 *--------------------------------------------------------------------------*/

void  SolverDiffusionFreeInstanceXtra()
{
  PFModule      *this_module = ThisPFModule;
  InstanceXtra  *instance_xtra = (InstanceXtra*)PFModuleInstanceXtra(this_module);

  PublicXtra    *public_xtra = (PublicXtra*)PFModulePublicXtra(this_module);
  Problem       *problem = (public_xtra->problem);
  int is_multiphase;

  is_multiphase = ProblemNumPhases(problem) > 1;

  if (instance_xtra)
  {
    FreeTempData((instance_xtra->temp_data));

    PFModuleFreeInstance((instance_xtra->ic_phase_concen));
    PFModuleFreeInstance((instance_xtra->phase_mobility));
    PFModuleFreeInstance((instance_xtra->retardation));

    PFModuleFreeInstance((instance_xtra->set_problem_data));
    PFModuleFreeInstance((instance_xtra->advect_concen));
    PFModuleFreeInstance((instance_xtra->phase_velocity_face));
    PFModuleFreeInstance((instance_xtra->linear_solver));
    PFModuleFreeInstance((instance_xtra->diag_scale));
    PFModuleFreeInstance((instance_xtra->discretize_pressure));

    PFModuleFreeInstance((instance_xtra->phase_density));

    if (is_multiphase)
    {
      PFModuleFreeInstance((instance_xtra->constitutive));
      PFModuleFreeInstance((instance_xtra->bc_phase_saturation));
      PFModuleFreeInstance((instance_xtra->ic_phase_satur));
      PFModuleFreeInstance((instance_xtra->advect_satur));
      PFModuleFreeInstance((instance_xtra->total_velocity_face));
      PFModuleFreeInstance((instance_xtra->permeability_face));
    }

    FreeProblemData((instance_xtra->problem_data));

    FreeGrid((instance_xtra->z_grid));
    FreeGrid((instance_xtra->y_grid));
    FreeGrid((instance_xtra->x_grid));
    FreeGrid((instance_xtra->grid));

    tfree(instance_xtra);
  }
}

/*--------------------------------------------------------------------------
 * SolverDiffusionNewPublicXtra
 *--------------------------------------------------------------------------*/

PFModule   *SolverDiffusionNewPublicXtra(char *name)
{
  PFModule      *this_module = ThisPFModule;
  PublicXtra    *public_xtra;

  char          *switch_name;
  int switch_value;
  NameArray switch_na;

  char key[IDB_MAX_KEY_LEN];

  NameArray diag_solver_na;
  NameArray linear_solver_na;

  switch_na = NA_NewNameArray("False True");

  public_xtra = ctalloc(PublicXtra, 1);

  diag_solver_na = NA_NewNameArray("NoDiagScale MatDiagScale");
  sprintf(key, "%s.DiagSolver", name);
  switch_name = GetStringDefault(key, "NoDiagScale");
  switch_value = NA_NameToIndexExitOnError(diag_solver_na, switch_name, key);
  switch (switch_value)
  {
    case 0:
    {
      (public_xtra->diag_scale) = PFModuleNewModuleType(MatrixDiagScaleNewPublicXtraInvoke,
                                                        MatDiagScale, (key));
      break;
    }

    case 1:
    {
      (public_xtra->diag_scale) = PFModuleNewModuleType(MatrixDiagScaleNewPublicXtraInvoke,
                                                        NoDiagScale, (key));
      break;
    }

    default:
    {
      InputError("Invalid switch value <%s> for key <%s>", switch_name, key);
    }
  }
  NA_FreeNameArray(diag_solver_na);

  linear_solver_na = NA_NewNameArray("MGSemi PPCG PCG CGHS");
  sprintf(key, "%s.Linear", name);
  switch_name = GetStringDefault(key, "PPCG");
  switch_value = NA_NameToIndexExitOnError(linear_solver_na, switch_name, key);
  switch (switch_value)
  {
    case 0:
    {
      (public_xtra->linear_solver) = PFModuleNewModuleType(LinearSolverNewPublicXtraInvoke,
                                                           MGSemi, (key));
      break;
    }

    case 1:
    {
      (public_xtra->linear_solver) = PFModuleNewModuleType(LinearSolverNewPublicXtraInvoke,
                                                           PPCG, (key));
      break;
    }

    case 2:
    {
      (public_xtra->linear_solver) = PFModuleNewModuleType(LinearSolverNewPublicXtraInvoke,
                                                           PCG, (key));
      break;
    }

    case 3:
    {
      (public_xtra->linear_solver) = PFModuleNewModuleType(LinearSolverNewPublicXtraInvoke,
                                                           CGHS, (key));
      break;
    }

    default:
    {
      InputError("Invalid switch value <%s> for key <%s>", switch_name, key);
    }
  }
  NA_FreeNameArray(linear_solver_na);

  (public_xtra->discretize_pressure) =
    PFModuleNewModule(DiscretizePressure, ());

  (public_xtra->permeability_face) = PFModuleNewModule(PermeabilityFace, ());
  (public_xtra->phase_velocity_face) = PFModuleNewModule(PhaseVelocityFace, ());
  (public_xtra->total_velocity_face) = PFModuleNewModule(TotalVelocityFace, ());
  (public_xtra->advect_satur) = PFModuleNewModule(SatGodunov, ());
  (public_xtra->advect_concen) = PFModuleNewModule(Godunov, ());
  (public_xtra->set_problem_data) = PFModuleNewModule(SetProblemData, ());

  (public_xtra->problem) = NewProblem(DiffusionSolve);

  sprintf(key, "%s.CFL", name);
  public_xtra->CFL = GetDoubleDefault(key, 1.0);

  sprintf(key, "%s.MaxIter", name);
  public_xtra->max_iterations = GetIntDefault(key, 99999);

  sprintf(key, "%s.BetaPerm", name);
  public_xtra->beta_perm = GetDoubleDefault(key, 0);

  sprintf(key, "%s.BetaFluid", name);
  public_xtra->beta_fluid = GetDoubleDefault(key, 1.6E-9);

  sprintf(key, "%s.BetaPore", name);
  public_xtra->beta_pore = GetDoubleDefault(key, 4.0E-10);

  sprintf(key, "%s.BetaFracture", name);
  public_xtra->beta_fracture = GetDoubleDefault(key, 0.0);

  sprintf(key, "%s.CompCompress", name);
  switch_name = GetStringDefault(key, "True");
  switch_value = NA_NameToIndexExitOnError(switch_na, switch_name, key);
  public_xtra->comp_compress_flag = switch_value;

  NA_FreeNameArray(switch_na);

  PFModulePublicXtra(this_module) = public_xtra;
  return this_module;
}

/*--------------------------------------------------------------------------
 * SolverDiffusionFreePublicXtra
 *--------------------------------------------------------------------------*/

void   SolverDiffusionFreePublicXtra()
{
  PFModule      *this_module = ThisPFModule;
  PublicXtra    *public_xtra = (PublicXtra*)PFModulePublicXtra(this_module);

  if (public_xtra)
  {
    FreeProblem(public_xtra->problem, DiffusionSolve);

    PFModuleFreeModule(public_xtra->diag_scale);
    PFModuleFreeModule(public_xtra->linear_solver);
    PFModuleFreeModule(public_xtra->set_problem_data);
    PFModuleFreeModule(public_xtra->advect_concen);
    PFModuleFreeModule(public_xtra->advect_satur);
    PFModuleFreeModule(public_xtra->total_velocity_face);
    PFModuleFreeModule(public_xtra->phase_velocity_face);
    PFModuleFreeModule(public_xtra->permeability_face);
    PFModuleFreeModule(public_xtra->discretize_pressure);
    tfree(public_xtra);
  }
}

/*--------------------------------------------------------------------------
 * SolverDiffusionSizeOfTempData
 *--------------------------------------------------------------------------*/

int  SolverDiffusionSizeOfTempData()
{
  /* SGS temp data */
  return 0;
}
