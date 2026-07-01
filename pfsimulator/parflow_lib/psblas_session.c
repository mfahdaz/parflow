/*BHEADER**********************************************************************
*
*  Copyright (c) 1995-2026, Lawrence Livermore National Security,
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

#include "parflow.h"
#include "psblas_session.h"

PSBLASSession* NewPSBLASSession()
{
  PSBLASSession *session = (PSBLASSession*)ctalloc(PSBLASSession, 1);

  /* Create new PSBLAS Context */
  PSBLASSessionContext(session) = psb_c_new_ctxt();
  /* Create new PSBLAS Descriptor */
  PSBLASSessionDescriptor(session) = psb_c_new_descriptor();

  return session;
}

void FreePSBLASSession(PSBLASSession *session)
{
  if (session != NULL)
  {
    SUNMatDestroy(PSBLASSessionSUNMatrix(session));

    psb_c_cdfree(PSBLASSessionDescriptor(session));
    psb_c_delete_descriptor(PSBLASSessionDescriptor(session));

    // psb_c_exit_ctxt(*cctxt);
    psb_c_delete_ctxt(PSBLASSessionContext(session));

    tfree(session);
  }

  return;
}

void InitPSBLASSession(PSBLASSession *session, Grid *grid)
{
  psb_i_t info = 0;

  /* Init PSBLAS Context */
  psb_c_init_from_fint(PSBLASSessionContext(session), MPI_Comm_c2f(amps_CommWorld));
  psb_c_set_index_base(0);

  /* Count number of elements in current process. */
  psb_i_t nl = 0;
  int is = 0;
  ForSubgridI(is, GridSubgrids(grid))
  {
    Subgrid *sg = GridSubgrid(grid, is);
    nl += (psb_i_t)(SubgridNX(sg) * SubgridNY(sg) * SubgridNZ(sg));
  }

  /* Set local to glabal index mapping */
  psb_l_t *vl = ctalloc(psb_l_t, nl);
  Subgrid *user_subgrid = GridSubgrid(GlobalsUserGrid, 0);
  ForSubgridI(is, GridSubgrids(grid))
  {
    Subgrid *subgrid = GridSubgrid(grid, is);

    int ix = SubgridIX(subgrid);
    int iy = SubgridIY(subgrid);
    int iz = SubgridIZ(subgrid);
    int nx = SubgridNX(subgrid);
    int ny = SubgridNY(subgrid);
    int nz = SubgridNZ(subgrid);

    int i = 0, j = 0, k = 0;
    BoxLoopI0(i, j, k, ix, iy, iz, nx, ny, nz,
    {
      int local_idx = SubgridEltIndex(subgrid, i, j, k);
      int global_idx = SubgridEltIndex(user_subgrid, i, j, k);
      vl[local_idx] = (psb_l_t)global_idx;
    });
  }

  /* allocate a context descriptor */
  info = psb_c_cdall_vl(nl, vl, *PSBLASSessionContext(session), PSBLASSessionDescriptor(session));
  free(vl);
  if (info != 0) {
    amps_Printf("Error in psb_c_cdall_vl: %d\n", info);
  }

  /* context descriptor is finalized in psb_c_cdasb */
  info = psb_c_cdasb(PSBLASSessionDescriptor(session));
  if (info != 0) {
    amps_Printf("Error in psb_c_cdasb: %d\n", info);
  }

  /* Create PSBLAS SUNMatrix (Full System Jacobian) */
  PSBLASSessionSUNMatrix(session) = SUNPSBLASMatrix(
      PSBLASSessionContext(session), 
      PSBLASSessionDescriptor(session)
    );

  if (PSBLASSessionSUNMatrix(session) == NULL) {
    amps_Printf("Error: Failure to create a SUNMatrix\n");
  }

  /* Create PSBLAS SUNMatrix (Preconditioner Matrix) */
  PSBLASSessionPrecSUNMatrix(session) = SUNPSBLASMatrix(
      PSBLASSessionContext(session), 
      PSBLASSessionDescriptor(session)
    );

  if (PSBLASSessionPrecSUNMatrix(session) == NULL) {
    amps_Printf("Error: Failure to create a Preconditioner SUNMatrix\n");
  }

  return;
}

void Set_N_Vector_From_Vector(N_Vector nvec, Vector *vec)
{
  double *psb_data = N_VGetArrayPointer(nvec);

  Grid *grid = VectorGrid(vec);

  int is = 0;
  ForSubgridI(is, GridSubgrids(grid))
  {
    Subgrid *subgrid = GridSubgrid(grid, is);
    Subvector *v_sub = VectorSubvector(vec, is);
    double *v_data = SubvectorData(v_sub);

    int ix = SubgridIX(subgrid);
    int iy = SubgridIY(subgrid);
    int iz = SubgridIZ(subgrid);
    int nx = SubgridNX(subgrid);
    int ny = SubgridNY(subgrid);
    int nz = SubgridNZ(subgrid);
    int nx_v = SubvectorNX(v_sub);
    int ny_v = SubvectorNY(v_sub);
    int nz_v = SubvectorNZ(v_sub);

    int i = 0, j = 0, k = 0;
    int pf_idx = SubvectorEltIndex(v_sub, ix, iy, iz);
    int psb_idx = SubgridEltIndex(subgrid, ix, iy, iz); /* always 0, kept explicit for symmetry/safety */
    BoxLoopI2(i, j, k, ix, iy, iz, nx, ny, nz,
              pf_idx, nx_v, ny_v, nz_v, 1, 1, 1,
              psb_idx, nx, ny, nz, 1, 1, 1,
    {
      if (pf_idx != SubvectorEltIndex(v_sub, i, j, k) || psb_idx != SubgridEltIndex(subgrid, i, j, k))
      {
        amps_Printf("FATAL: Set_N_Vector_From_Vector index mismatch at (%d,%d,%d)\n", i, j, k);
      }
      psb_data[psb_idx] = v_data[pf_idx];
    });
  }
  return;
}

void Set_Vector_From_N_Vector(Vector *vec, N_Vector nvec)
{
  double *psb_data = N_VGetArrayPointer(nvec);

  Grid *grid = VectorGrid(vec);

  int is = 0;
  ForSubgridI(is, GridSubgrids(grid))
  {
    Subgrid *subgrid = GridSubgrid(grid, is);
    Subvector *v_sub = VectorSubvector(vec, is);
    double *v_data = SubvectorData(v_sub);

    int ix = SubgridIX(subgrid);
    int iy = SubgridIY(subgrid);
    int iz = SubgridIZ(subgrid);
    int nx = SubgridNX(subgrid);
    int ny = SubgridNY(subgrid);
    int nz = SubgridNZ(subgrid);
    int nx_v = SubvectorNX(v_sub);
    int ny_v = SubvectorNY(v_sub);
    int nz_v = SubvectorNZ(v_sub);

    int i = 0, j = 0, k = 0;
    int pf_idx = SubvectorEltIndex(v_sub, ix, iy, iz);
    int psb_idx = SubgridEltIndex(subgrid, ix, iy, iz);
    BoxLoopI2(i, j, k, ix, iy, iz, nx, ny, nz,
              pf_idx, nx_v, ny_v, nz_v, 1, 1, 1,
              psb_idx, nx, ny, nz, 1, 1, 1,
    {
      if (pf_idx != SubvectorEltIndex(v_sub, i, j, k) || psb_idx != SubgridEltIndex(subgrid, i, j, k))
      {
        amps_Printf("FATAL: Set_Vector_From_N_Vector index mismatch at (%d,%d,%d)\n", i, j, k);
      }
      v_data[pf_idx] = psb_data[psb_idx];
    });
  }
  return;
}

void Set_SUNMatrix_From_Matrix(SUNMatrix sunmat,
                               Matrix *JB,
                               Matrix *JC,
                               void *current_state)
{
  Subgrid *user_subgrid = GridSubgrid(GlobalsUserGrid, 0);

  ProblemData *problem_data = StateProblemData(((State*)current_state));
  Vector *top = ProblemDataIndexOfDomainTop(problem_data);
  Vector *bot = ProblemDataIndexOfDomainBottom(problem_data);

  Grid *JB_grid = MatrixGrid(JB);
  Stencil *JB_stencil = MatrixStencil(JB);
  int JB_stencil_size = StencilSize(JB_stencil);
  StencilElt *JB_shape = StencilShape(JB_stencil);

  int GNX = SubgridNX(user_subgrid);
  int GNY = SubgridNY(user_subgrid);
  int GNZ = SubgridNZ(user_subgrid);

  int *JB_stencil_offset = ctalloc(int, JB_stencil_size);
  for (int istencil = 0; istencil < JB_stencil_size; ++istencil)
  {
    JB_stencil_offset[istencil] = JB_shape[istencil][0]
                                 + JB_shape[istencil][1] * GNX
                                 + JB_shape[istencil][2] * GNX * GNY;
  }

  int isubgrid = 0;
  ForSubgridI(isubgrid, GridSubgrids(JB_grid))
  {
    Subgrid *subgrid = SubgridArraySubgrid(GridSubgrids(JB_grid), isubgrid);

    Submatrix *JB_sub = MatrixSubmatrix(JB, isubgrid);
    Subvector *top_sub = VectorSubvector(top, isubgrid);
    double *top_dat = SubvectorData(top_sub);
    Subvector *bot_sub = VectorSubvector(bot, isubgrid);
    double *bot_dat = SubvectorData(bot_sub);

    int ix = SubgridIX(subgrid);
    int iy = SubgridIY(subgrid);
    int iz = SubgridIZ(subgrid);
    int nx = SubgridNX(subgrid);
    int ny = SubgridNY(subgrid);
    int nz = SubgridNZ(subgrid);
    int nx_m = SubmatrixNX(JB_sub);
    int ny_m = SubmatrixNY(JB_sub);
    int nz_m = SubmatrixNZ(JB_sub);

    int row_buf_size = nx * JB_stencil_size;
    psb_l_t *idx_row = ctalloc(psb_l_t, row_buf_size);
    psb_l_t *idx_col = ctalloc(psb_l_t, row_buf_size);
    psb_d_t *psb_val = ctalloc(psb_d_t, row_buf_size);

    int i = 0, j = 0, k = 0;
    /* idx must be seeded to the submatrix's actual base
     * offset (generally nonzero i.e. Submatrix data carries ghost
     * padding), not to 0. */
    int idx = SubmatrixEltIndex(JB_sub, ix, iy, iz);
    int psb_row_idx = SubgridEltIndex(user_subgrid, ix, iy, iz);
    int n_ins = 0;

    BoxLoopI2(i, j, k, ix, iy, iz, nx, ny, nz,
              idx, nx_m, ny_m, nz_m, 1, 1, 1,
              psb_row_idx, GNX, GNY, GNZ, 1, 1, 1,
    {
      int pf_idx = idx;
      if (pf_idx != SubmatrixEltIndex(JB_sub, i, j, k) ||
          psb_row_idx != SubgridEltIndex(user_subgrid, i, j, k))
      {
        amps_Printf("FATAL: Set_SUNMatrix_From_Matrix JB index mismatch at (%d,%d,%d)\n", i, j, k);
      }

      for(int istencil = 0; istencil < JB_stencil_size; ++istencil)
      {
        double *JB_dat = SubmatrixStencilData(JB_sub, istencil);

        int st_i = JB_shape[istencil][0];
        int st_j = JB_shape[istencil][1];
        int st_k = JB_shape[istencil][2];

        int itop = SubvectorEltIndex(top_sub, (i + st_i), (j + st_j), 0);

        if (top_dat[itop] < 0 || k + st_k > lrint(top_dat[itop]) || k + st_k < lrint(bot_dat[itop]))
        {
          continue;
        }

        idx_row[n_ins] = psb_row_idx;
        idx_col[n_ins] = psb_row_idx + JB_stencil_offset[istencil];
        psb_val[n_ins] = JB_dat[pf_idx];
        ++n_ins;
      }

      if (i == ix + nx - 1 && n_ins > 0)
      {
        SUNMatIns_PSBLAS(n_ins, idx_row, idx_col, psb_val, sunmat);
        n_ins = 0;
      }
    });

    free(idx_row);
    free(idx_col);
    free(psb_val);
  }

  free(JB_stencil_offset);

  if(JC == NULL)
  {
    return;
  }

  Stencil *JC_stencil = MatrixStencil(JC);
  int JC_stencil_size = StencilSize(JC_stencil);
  StencilElt *JC_shape = StencilShape(JC_stencil);

  isubgrid = 0;
  ForSubgridI(isubgrid, GridSubgrids(JB_grid))
  {
    Subgrid *subgrid = SubgridArraySubgrid(GridSubgrids(JB_grid), isubgrid);

    Submatrix *JC_sub = MatrixSubmatrix(JC, isubgrid);
    Subvector *top_sub = VectorSubvector(top, isubgrid);
    double *top_dat = SubvectorData(top_sub);

    int ix = SubgridIX(subgrid);
    int iy = SubgridIY(subgrid);
    int iz = SubgridIZ(subgrid);
    int nx = SubgridNX(subgrid);
    int ny = SubgridNY(subgrid);
    int nz = SubgridNZ(subgrid);
    int nx_m = SubmatrixNX(JC_sub);
    int ny_m = SubmatrixNY(JC_sub);
    int nz_m = SubmatrixNZ(JC_sub);

    int row_buf_size = nx * JC_stencil_size;
    psb_l_t *idx_row = ctalloc(psb_l_t, row_buf_size);
    psb_l_t *idx_col = ctalloc(psb_l_t, row_buf_size);
    psb_d_t *psb_val = ctalloc(psb_d_t, row_buf_size);

    int i = 0, j = 0, k = 0;
    int idx = SubmatrixEltIndex(JC_sub, ix, iy, iz);  /* same as JB */
    int n_ins = 0;

    BoxLoopI1(i, j, k, ix, iy, iz, nx, ny, 1,
              idx, nx_m, ny_m, nz_m, 1, 1, 1,
    {
      int itop = SubvectorEltIndex(top_sub, i, j, 0);
      int k_ = (int)top_dat[itop];
      int pf_idx = idx;
      if (pf_idx != SubmatrixEltIndex(JC_sub, i, j, k))
      {
        amps_Printf("FATAL: Set_SUNMatrix_From_Matrix JC index mismatch at (%d,%d,%d)\n", i, j, k);
      }

      if (k_ >= 0)
      {
        int psb_row_idx = SubgridEltIndex(user_subgrid, i, j, k_);

        for(int istencil = 0; istencil < JC_stencil_size; ++istencil)
        {
          double *JC_dat = SubmatrixStencilData(JC_sub, istencil);

          int st_i = JC_shape[istencil][0];
          int st_j = JC_shape[istencil][1];

          itop = SubvectorEltIndex(top_sub, (i + st_i), (j + st_j), 0);
          int kk = (int)top_dat[itop];

          if (kk < 0)
          {
            continue;
          }

          idx_row[n_ins] = psb_row_idx;
          idx_col[n_ins] = SubgridEltIndex(user_subgrid, (i + st_i), (j + st_j), kk);
          psb_val[n_ins] = JC_dat[pf_idx];
          ++n_ins;
        }
      }

      if (i == ix + nx - 1 && n_ins > 0)
      {
        SUNMatIns_PSBLAS(n_ins, idx_row, idx_col, psb_val, sunmat);
        n_ins = 0;
      }
    });

    free(idx_row);
    free(idx_col);
    free(psb_val);
  }

  return;
}

void Set_SUNMatrix_From_SymmetricMatrix(SUNMatrix sunmat,
                                        Matrix *JB,
                                        Matrix *JC,
                                        void *current_state)
{
   /* COO insertion into PSBLAS accumulates unless the matrix is cleared
   * first. Nothing else clears this matrix between recompute events.
   * It isn't part of KINSol's own matrix protocol, so we zero it
   * unconditionally here. */
  // SUNMatZero(sunmat);

  Subgrid *user_subgrid = GridSubgrid(GlobalsUserGrid, 0);

  ProblemData *problem_data = StateProblemData(((State*)current_state));
  Vector *top = ProblemDataIndexOfDomainTop(problem_data);
  Vector *bot = ProblemDataIndexOfDomainBottom(problem_data);

  Grid *JB_grid = MatrixGrid(JB);
  Stencil *JB_stencil = MatrixStencil(JB);
  int JB_stencil_size = StencilSize(JB_stencil);
  StencilElt *JB_shape = StencilShape(JB_stencil);

  int GNX = SubgridNX(user_subgrid);
  int GNY = SubgridNY(user_subgrid);
  int GNZ = SubgridNZ(user_subgrid);

  int *JB_stencil_offset = ctalloc(int, JB_stencil_size);
  int *JB_is_diag = ctalloc(int, JB_stencil_size);
  for (int istencil = 0; istencil < JB_stencil_size; ++istencil)
  {
    int st_i = JB_shape[istencil][0];
    int st_j = JB_shape[istencil][1];
    int st_k = JB_shape[istencil][2];
    JB_stencil_offset[istencil] = st_i + st_j * GNX + st_k * GNX * GNY;
    JB_is_diag[istencil] = (st_i == 0 && st_j == 0 && st_k == 0);
  }

  int isubgrid = 0;
  ForSubgridI(isubgrid, GridSubgrids(JB_grid))
  {
    Subgrid *subgrid = SubgridArraySubgrid(GridSubgrids(JB_grid), isubgrid);

    Submatrix *JB_sub = MatrixSubmatrix(JB, isubgrid);
    Subvector *top_sub = VectorSubvector(top, isubgrid);
    double *top_dat = SubvectorData(top_sub);
    Subvector *bot_sub = VectorSubvector(bot, isubgrid);
    double *bot_dat = SubvectorData(bot_sub);

    int ix = SubgridIX(subgrid);
    int iy = SubgridIY(subgrid);
    int iz = SubgridIZ(subgrid);
    int nx = SubgridNX(subgrid);
    int ny = SubgridNY(subgrid);
    int nz = SubgridNZ(subgrid);
    int nx_m = SubmatrixNX(JB_sub);
    int ny_m = SubmatrixNY(JB_sub);
    int nz_m = SubmatrixNZ(JB_sub);

    int row_buf_size = nx * 2 * JB_stencil_size;
    psb_l_t *idx_row = ctalloc(psb_l_t, row_buf_size);
    psb_l_t *idx_col = ctalloc(psb_l_t, row_buf_size);
    psb_d_t *psb_val = ctalloc(psb_d_t, row_buf_size);

    int i = 0, j = 0, k = 0;
    int idx = SubmatrixEltIndex(JB_sub, ix, iy, iz);  /* THE FIX */
    int psb_row_idx = SubgridEltIndex(user_subgrid, ix, iy, iz);
    int n_ins = 0;

    BoxLoopI2(i, j, k, ix, iy, iz, nx, ny, nz,
              idx, nx_m, ny_m, nz_m, 1, 1, 1,
              psb_row_idx, GNX, GNY, GNZ, 1, 1, 1,
    {
      int pf_idx = idx;
      if (pf_idx != SubmatrixEltIndex(JB_sub, i, j, k) ||
          psb_row_idx != SubgridEltIndex(user_subgrid, i, j, k))
      {
        amps_Printf("FATAL: Set_SUNMatrix_From_SymmetricMatrix JB index mismatch at (%d,%d,%d)\n", i, j, k);
      }

      for(int istencil = 0; istencil < JB_stencil_size; ++istencil)
      {
        double *JB_dat = SubmatrixStencilData(JB_sub, istencil);

        int st_i = JB_shape[istencil][0];
        int st_j = JB_shape[istencil][1];
        int st_k = JB_shape[istencil][2];

        int itop = SubvectorEltIndex(top_sub, (i + st_i), (j + st_j), 0);

        if (top_dat[itop] < 0 || k + st_k > lrint(top_dat[itop]) || k + st_k < lrint(bot_dat[itop]))
        {
          continue;
        }

        int psb_col_idx = psb_row_idx + JB_stencil_offset[istencil];
        double val = JB_dat[pf_idx];

        idx_row[n_ins] = psb_row_idx;
        idx_col[n_ins] = psb_col_idx;
        psb_val[n_ins] = val;
        ++n_ins;

        if (!JB_is_diag[istencil])
        {
          idx_row[n_ins] = psb_col_idx;
          idx_col[n_ins] = psb_row_idx;
          psb_val[n_ins] = val;
          ++n_ins;
        }
      }

      if (i == ix + nx - 1 && n_ins > 0)
      {
        SUNMatIns_PSBLAS(n_ins, idx_row, idx_col, psb_val, sunmat);
        n_ins = 0;
      }
    });

    free(idx_row);
    free(idx_col);
    free(psb_val);
  }

  free(JB_stencil_offset);
  free(JB_is_diag);

  if(JC == NULL)
  {
    return;
  }

  Stencil *JC_stencil = MatrixStencil(JC);
  int JC_stencil_size = StencilSize(JC_stencil);
  StencilElt *JC_shape = StencilShape(JC_stencil);

  isubgrid = 0;
  ForSubgridI(isubgrid, GridSubgrids(JB_grid))
  {
    Subgrid *subgrid = SubgridArraySubgrid(GridSubgrids(JB_grid), isubgrid);

    Submatrix *JC_sub = MatrixSubmatrix(JC, isubgrid);
    Subvector *top_sub = VectorSubvector(top, isubgrid);
    double *top_dat = SubvectorData(top_sub);

    int ix = SubgridIX(subgrid);
    int iy = SubgridIY(subgrid);
    int iz = SubgridIZ(subgrid);
    int nx = SubgridNX(subgrid);
    int ny = SubgridNY(subgrid);
    int nz = SubgridNZ(subgrid);
    int nx_m = SubmatrixNX(JC_sub);
    int ny_m = SubmatrixNY(JC_sub);
    int nz_m = SubmatrixNZ(JC_sub);

    int row_buf_size = nx * 2 * JC_stencil_size;
    psb_l_t *idx_row = ctalloc(psb_l_t, row_buf_size);
    psb_l_t *idx_col = ctalloc(psb_l_t, row_buf_size);
    psb_d_t *psb_val = ctalloc(psb_d_t, row_buf_size);

    int i = 0, j = 0, k = 0;
    int idx = SubmatrixEltIndex(JC_sub, ix, iy, iz);  /* THE FIX */
    int n_ins = 0;

    BoxLoopI1(i, j, k, ix, iy, iz, nx, ny, 1,
              idx, nx_m, ny_m, nz_m, 1, 1, 1,
    {
      int itop = SubvectorEltIndex(top_sub, i, j, 0);
      int k_ = (int)top_dat[itop];
      int pf_idx = idx;
      if (pf_idx != SubmatrixEltIndex(JC_sub, i, j, k))
      {
        amps_Printf("FATAL: Set_SUNMatrix_From_SymmetricMatrix JC index mismatch at (%d,%d,%d)\n", i, j, k);
      }

      if (k_ >= 0)
      {
        int psb_row_idx = SubgridEltIndex(user_subgrid, i, j, k_);

        for(int istencil = 0; istencil < JC_stencil_size; ++istencil)
        {
          double *JC_dat = SubmatrixStencilData(JC_sub, istencil);

          int st_i = JC_shape[istencil][0];
          int st_j = JC_shape[istencil][1];
          int is_diagonal = (st_i == 0 && st_j == 0);

          itop = SubvectorEltIndex(top_sub, (i + st_i), (j + st_j), 0);
          int kk = (int)top_dat[itop];

          if (kk < 0)
          {
            continue;
          }

          int psb_col_idx = SubgridEltIndex(user_subgrid, (i + st_i), (j + st_j), kk);
          double val = JC_dat[pf_idx];

          idx_row[n_ins] = psb_row_idx;
          idx_col[n_ins] = psb_col_idx;
          psb_val[n_ins] = val;
          ++n_ins;

          if (!is_diagonal)
          {
            idx_row[n_ins] = psb_col_idx;
            idx_col[n_ins] = psb_row_idx;
            psb_val[n_ins] = val;
            ++n_ins;
          }
        }
      }

      if (i == ix + nx - 1 && n_ins > 0)
      {
        SUNMatIns_PSBLAS(n_ins, idx_row, idx_col, psb_val, sunmat);
        n_ins = 0;
      }
    });

    free(idx_row);
    free(idx_col);
    free(psb_val);
  }

  return;
}