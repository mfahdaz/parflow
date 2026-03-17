/*BHEADER*********************************************************************
 *
 *  This file is part of Parflow. For details, see
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
/*****************************************************************************
 *
 * Routines to write data Vector to a file in full or scattered form using PDI.
 *
 *****************************************************************************/

#include "parflow.h"
#ifdef PARFLOW_HAVE_PDI
#include <paraconf.h>
#include <pdi.h>
#endif
#include <stdbool.h>
#include <unistd.h>
#include <parflow_proto.h>

/**
 * @brief Read the configuration yaml for PDI and initialize it.
 *
 * This function uses the paraconf library to read the configuration file in yaml format
 * and initializes the PDI library once for the entire simulation. This call is meant to 
 * happen only once at the beginning of the simulation.
 *
 *
 * @warning The function relies on PDI and requires a valid PDI configuration
 * (`conf.yml`). Ensure PDI is correctly initialized and available in the
 * environment.
 *
 * @see PDI_expose(), PDI_init(), PDI_finalize() at
 * https://pdi.dev/1.8/modules.html
 */

/**
 * @brief Finalize the PDI library.
 *
 * This function finalizes the PDI library, cleaning up resources and
 * ensuring that all data is written out correctly. It should be called at the end
 * of the simulation to ensure that all PDI resources are properly released.
 *
 * @warning The function relies on PDI. Ensure PDI is correctly initialized and available in the
 * environment.
 *
 * @see PDI_expose(), PDI_init(), PDI_finalize() at
 * https://pdi.dev/1.8/modules.html
 */

void InitializeInsitu() {
#ifdef PARFLOW_HAVE_PDI
        // Call the Init event which will instantiate doreisa client. This should be done once.
        BeginTiming(InsituSetupTimingIndex);
        PDI_event("Init");
        EndTiming(InsituSetupTimingIndex);
#endif
}

void FinalizeInsitu() {
#ifdef PARFLOW_HAVE_PDI
        // Call the Init event which will instantiate doreisa client. This should be done once.
        BeginTiming(InsituCleanupTimingIndex);
        PDI_event("Finalize");
        EndTiming(InsituCleanupTimingIndex)
#endif
}

void ShareDataInsitu(char *name_insitu,
                      Vector *vector,
                      double stop_time, 
                      int current_iteration) {
#ifdef PARFLOW_HAVE_PDI
        BeginTiming(InsituShareTimingIndex)
        int rank = amps_Rank(amps_CommWorld);

        // coordinates of each rank
        int rank_coord_x = GlobalsP;
        int rank_coord_y = GlobalsQ;
        int rank_coord_z = GlobalsR;

        // number of ranks in each direction
        int mpi_size = GlobalsNumProcs;
        int mpi_size_x = GlobalsNumProcsX;
        int mpi_size_y = GlobalsNumProcsY;
        int mpi_size_z = GlobalsNumProcsZ;

        // max iterations
        int max_iters = (int)stop_time;
        int current_iter = current_iteration;

        // subvector to extract nx, ny, nz - local problem size + ghost cells
        Subvector *subvector = VectorSubvector(vector, 0);
        Grid *grid = VectorGrid(vector);

        // global problem size
        int iNX = SubgridNX(GridBackground(grid));
        int iNY = SubgridNY(GridBackground(grid));
        int iNZ = SubgridNZ(GridBackground(grid));

        // problem size for each rank
        int inx = subvector->data_space->nx;
        int iny = subvector->data_space->ny;
        int inz = subvector->data_space->nz;

        double *data_array = SubvectorData(subvector);

        PDI_multi_expose("Share", 
			"inx", &inx , PDI_OUT,  "iny", &iny, PDI_OUT, "inz", &inz, PDI_OUT,
			"rank", &rank, PDI_OUT, "rank_coord_x", &rank_coord_x, PDI_OUT,
                         "rank_coord_y", &rank_coord_y, PDI_OUT, "rank_coord_z", &rank_coord_z, PDI_OUT,
                         "mpi_size", &mpi_size, PDI_OUT, "mpi_size_x", &mpi_size_x, PDI_OUT,
                         "mpi_size_y", &mpi_size_y, PDI_OUT, "mpi_size_z", &mpi_size_z, PDI_OUT,
                         "it", &current_iter, PDI_OUT, "data_array", data_array, PDI_OUT, "name",
                         name_insitu, PDI_OUT, NULL);
        EndTiming(InsituShareTimingIndex)
#endif
}

