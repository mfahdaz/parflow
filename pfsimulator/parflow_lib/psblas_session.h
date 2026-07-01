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

#ifdef PARFLOW_HAVE_PSCTOOLKIT

#ifndef _PSBLAS_SESSION_HEADER
#define _PSBLAS_SESSION_HEADER

#include "parflow.h"
#include "kinsol_dependences.h"

typedef struct {

  psb_c_ctxt *context;          /* PSBLAS Context */
  psb_c_descriptor *descriptor; /* PSBLAS Descriptor */

  SUNMatrix sunmat;  /* PSBLAS SUNMatrix object (Full System Jacobian) */

  SUNMatrix sunmat_prec;  /* PSBLAS SUNMatrix object (preconditioner matrix B) */

} PSBLASSession;


PSBLASSession* NewPSBLASSession();

void FreePSBLASSession(PSBLASSession *session);

void InitPSBLASSession(PSBLASSession *session, Grid *grid);


void Set_N_Vector_From_Vector(N_Vector nvec, Vector *vec);

void Set_Vector_From_N_Vector(Vector *vec, N_Vector nvec);

void Set_SUNMatrix_From_Matrix(SUNMatrix sunmat,
                               Matrix *JB,
                               Matrix *JC,
                               void *current_state);

void Set_SUNMatrix_From_SymmetricMatrix(SUNMatrix sunmat,
                                        Matrix *JB,
                                        Matrix *JC,
                                        void *current_state);

int KINSolJacobianFunction(N_Vector pf_n_pressure,
                           N_Vector pf_n_fval,
                           SUNMatrix Jacobian,
                           void *current_state,
                           N_Vector pf_n_tmp1,
                           N_Vector pf_n_tmp2);

#define PSBLASSessionContext(session)    ((session)->context)
#define PSBLASSessionDescriptor(session) ((session)->descriptor)
#define PSBLASSessionSUNMatrix(session)  ((session)->sunmat)
#define PSBLASSessionPrecSUNMatrix(session)  ((session)->sunmat_prec)

#endif // _PSBLAS_SESSION_HEADER

#endif // PARFLOW_HAVE_PSCTOOLKIT