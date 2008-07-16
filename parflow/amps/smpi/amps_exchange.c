/*BHEADER**********************************************************************
 * (c) 1995   The Regents of the University of California
 *
 * See the file COPYRIGHT_and_DISCLAIMER for a complete copyright
 * notice, contact person, and disclaimer.
 *
 * $Revision: 1.1.1.1 $
 *********************************************************************EHEADER*/

#include <sys/times.h>
#include "amps.h"

int _amps_send_sizes(amps_Package package, int **sizes)
{
   int i;
   char *buffer;

   *sizes = (int*)calloc(package -> num_send, sizeof(int));
   
    for (i = 0; i < package -> num_send; i++)
    {
       (*sizes)[i] = amps_pack(amps_CommWorld, package -> send_invoices[i],
                             &buffer);
       MPI_Isend(&((*sizes)[i]), 1, MPI_INT, package -> dest[i],
                 0, amps_CommWorld,
                 &(package -> send_requests[i]));
    }

   return(0);
}

#ifdef AMPS_MPI_NOT_USE_PERSISTENT

int _amps_recv_sizes(amps_Package package)
{
   int i;
   int size;

   MPI_Status status;

   for(i = 0; i < package -> num_recv; i++)
   {
      MPI_Recv(&size, 1, MPI_INT, package -> src[i], 0,
               amps_CommWorld, &status);

      package -> recv_invoices[i] -> combuf =
                 (char *)calloc(size, sizeof(char *));

      /* Post receives for incoming byte buffers */

      MPI_Irecv(package -> recv_invoices[i] -> combuf, size,
                MPI_BYTE, package -> src[i], 1,
                amps_CommWorld, &(package -> recv_requests[i]));
   }

   return(0);
}

void _amps_wait_exchange(amps_Handle handle)
{
  int i;
  int size;
  int *flags;

  MPI_Status *status;

/*
  if (handle -> package -> num_recv)
     _amps_probes(handle -> package);
*/

  if(handle -> package -> num_recv + handle -> package -> num_send)
  {
     status = calloc((handle -> package -> num_recv + 
		      handle -> package -> num_send), sizeof(MPI_Status));

     MPI_Waitall(handle -> package -> num_recv + handle -> package -> num_send,
		 handle -> package -> recv_requests, status);

     fflush(NULL);
     
     free(status);

     for(i = 0; i < handle -> package -> num_recv; i++)
     {
        amps_unpack(amps_CommWorld, handle -> package -> recv_invoices[i],
                    handle -> package -> recv_invoices[i] -> combuf);
        AMPS_PACK_FREE_LETTER(amps_CommWorld,
                              handle -> package -> recv_invoices[i],
                              handle -> package -> recv_invoices[i] -> combuf);
     }
/*
     
     for(i = 0; i < handle -> package -> num_send; i++)
     {
	MPI_Type_free(&handle -> package -> send_invoices[i] -> mpi_type);
     }
*/
  }
}

amps_Handle amps_IExchangePackage(amps_Package package)
{

   int i;
   int done;
   int *flags;
   char *buffer;
   int size;
   struct timeval tm; 

   MPI_Status status;

   /*--------------------------------------------------------------------
    * send out the data we have
    *--------------------------------------------------------------------*/
   for(i = 0; i < package -> num_send; i++)
   {
      size = amps_pack(amps_CommWorld, package -> send_invoices[i],
                            &buffer);
      MPI_Isend(package -> send_invoices[i] -> combuf, size,
                MPI_BYTE, package -> dest[i], 1,
                amps_CommWorld, &(package -> send_requests[i]));
   }

/*
  fflush(NULL);
  gettimeofday(&tm, 0);
  amps_Printf("I'm at the Barrier %ld %ld\n", tm.tv_usec, tm.tv_sec);
  fflush(NULL);
  MPI_Barrier(amps_CommWorld);
  fflush(NULL);
  gettimeofday(&tm, 0);
  amps_Printf("I'm beyond the Barrier %ld %ld\n", tm.tv_usec, tm.tv_sec);
  fflush(NULL);
*/

  if (package -> num_recv)
  {
     flags = (int *)calloc(package -> num_recv, sizeof(int));

     for(i = 0; i < package -> num_recv; i++)
     {
        flags[i] = 0;
     }

     done = 0;
     while (!done)
     {
        done = 1;
        for(i = 0; i < package -> num_recv; i++)
        {
           if(!flags[i])
           {
              MPI_Iprobe(package -> src[i], 1, amps_CommWorld, &(flags[i]), &status);

              if(flags[i])
              { 
                 MPI_Get_count(&status, MPI_BYTE, &size);
                 package -> recv_invoices[i] -> combuf =
                         (char *)calloc(size, sizeof(char *));

                 MPI_Irecv(package -> recv_invoices[i] -> combuf, size,
                        MPI_BYTE, package -> src[i], 1,
                        amps_CommWorld, &(package -> recv_requests[i]));
              }
              else
                 done = 0;
           }
              
              
        }

     }

     free(flags);
  }

/*
   fflush(NULL);
   gettimeofday(&tm, 0);
   amps_Printf("IEX done %ld\n", tm.tv_sec);
   fflush(NULL);
*/
   return(amps_NewHandle(NULL, 0, NULL, package));
}

#else

int _amps_recv_sizes(amps_Package package)
{
   int i;
   int size;

   MPI_Status status;

   for(i = 0; i < package -> num_recv; i++)
   {
      MPI_Recv(&size, 1, MPI_INT, package -> src[i], 0,
               amps_CommWorld, &status); 

      package -> recv_invoices[i] -> combuf =
                 (char *)calloc(size, sizeof(char *));

      MPI_Recv_init(package -> recv_invoices[i] -> combuf, size,
                    MPI_BYTE, package -> src[i], 1, amps_CommWorld,
                    &(package -> recv_requests[i]));
   }

   return(0);
}

void _amps_wait_exchange(amps_Handle handle)
{
  int i;
  int num;

  num = handle -> package -> num_send + handle -> package -> num_recv;

  MPI_Waitall(num, handle -> package -> recv_requests, 
              handle -> package -> status);

  if(num)
  {
     if(handle -> package -> num_recv) 
     {
	for(i = 0; i <  handle -> package -> num_recv; i++)
	{
           amps_unpack(amps_CommWorld, handle -> package -> recv_invoices[i],
                       handle -> package -> recv_invoices[i] -> combuf);
           AMPS_PACK_FREE_LETTER(amps_CommWorld,
                             handle -> package -> recv_invoices[i],
                             handle -> package -> recv_invoices[i] -> combuf);
/*
	   AMPS_CLEAR_INVOICE(handle -> package -> recv_invoices[i]);
*/
	}
     }
	
  }
}

/*===========================================================================*/
/**

The \Ref{amps_IExchangePackage} initiates the communication of the
invoices found in the {\bf package} structure that is passed in.  Once a
\Ref{amps_IExchangePackage} is issued it is illegal to access the
variables that are being communicated.  An \Ref{amps_IExchangePackage}
is always followed by an \Ref{amps_Wait} on the {\bf handle} that is
returned. 

{\large Example:}
\begin{verbatim}
// Initialize exchange of boundary points 
handle = amps_IExchangePackage(package);
 
// Compute on the "interior points"

// Wait for the exchange to complete 
amps_Wait(handle);
\end{verbatim}

{\large Notes:}

This routine can be optimized on some architectures so if your
communication can be formulated using it there might be
some performance advantages.

@memo Initiate package communication
@param package the collection of invoices to communicate
@return Handle for the asynchronous communication
*/
amps_Handle amps_IExchangePackage(amps_Package package)
{
   int i;
   int *send_sizes;

   MPI_Status *status_array;

   if (package -> num_send)
      _amps_send_sizes(package, &send_sizes);

   if (package -> num_recv)
      _amps_recv_sizes(package);

   /*--------------------------------------------------------------------
    * end the sending of sizes
    *--------------------------------------------------------------------*/

   if (package -> num_send)
   {
      status_array = (MPI_Status *)calloc(package -> num_send,
                                           sizeof(MPI_Status));
      MPI_Waitall(package -> num_send, package -> send_requests, status_array);
      free(status_array);
   }

   if(package -> num_send)
   {
      for(i = 0; i < package -> num_send; i++)
      {
       
         MPI_Send_init(package -> send_invoices[i] -> combuf,
                       send_sizes[i], MPI_BYTE, package -> dest[i], 1,
                       amps_CommWorld, &(package -> send_requests[i]));
      }

      free(send_sizes);
   }

   /*--------------------------------------------------------------------
    * post receives for data to get
    *--------------------------------------------------------------------*/

   if(package -> num_recv)
      MPI_Startall(package -> num_recv, package -> recv_requests);

   /*--------------------------------------------------------------------
    * send out the data we have
    *--------------------------------------------------------------------*/
   if(package -> num_send)
      MPI_Startall(package -> num_send, package -> send_requests);

   return( amps_NewHandle(NULL, 0, NULL, package));
}

#endif
