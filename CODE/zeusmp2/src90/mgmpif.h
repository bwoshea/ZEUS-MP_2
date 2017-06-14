C***********************************************************************/
C*                                                                     */
C* File:        mgmpif.h                                               */
C*                                                                     */
C* Description: include file for Fortran applications calling MGMPI    */
C*                                                                     */
C*---------------------------------------------------------------------*/
C*                                                                     */
C* Author:      James Bordner                                          */
C*                                                                     */
C* © 2000-2003 Regents of the University of California.   */
C*                                                                     */
C***********************************************************************/

#include "mgmpi.def"
#include "mgmpif.def"

C     ----------------------------------------------------------
C     Function declarations
C     ----------------------------------------------------------

      INTEGER ff_function_prefix(BcFaceType)
      INTEGER ff_function_prefix(ConcurrentRank)
      INTEGER ff_function_prefix(ConcurrentSize)
      LOGICAL ff_function_prefix(ParametersIsDefined)
      _SCALAR ff_function_prefix(MatrixElement)
      _SCALAR ff_function_prefix(MatrixUnknown)
      _SCALAR ff_function_prefix(MatrixVertex)
      INTEGER ff_function_prefix(PerformanceCalls)
      _SCALAR ff_function_prefix(PerformancemFlopCount)
      _SCALAR ff_function_prefix(PerformancemFlopRate)
      _SCALAR ff_function_prefix(PerformanceTime)
      INTEGER ff_function_prefix(SolverStatus)
      _SCALAR ff_function_prefix(TimerValue)
      _SCALAR ff_function_prefix(VectorAllReduce)
      _SCALAR ff_function_prefix(VectorDot)
      _SCALAR ff_function_prefix(VectorElement)
      _SCALAR ff_function_prefix(VectorUnknown)
      _SCALAR ff_function_prefix(VectorVertex)
      _SCALAR ff_function_prefix(VectorNorm)
      _SCALAR ff_function_prefix(VectorSum)

      EXTERNAL ff_function_prefix(BcFaceType)
      EXTERNAL ff_function_prefix(ConcurrentRank)
      EXTERNAL ff_function_prefix(ConcurrentSize)
      EXTERNAL ff_function_prefix(ParametersIsDefined)
      EXTERNAL ff_function_prefix(MatrixElement)
      EXTERNAL ff_function_prefix(MatrixUnknown)
      EXTERNAL ff_function_prefix(MatrixVertex)
      EXTERNAL ff_function_prefix(PerformanceCalls)
      EXTERNAL ff_function_prefix(PerformancemFlopCount)
      EXTERNAL ff_function_prefix(PerformancemFlopRate)
      EXTERNAL ff_function_prefix(PerformanceTime)
      EXTERNAL ff_function_prefix(SolverStatus)
      EXTERNAL ff_function_prefix(TimerValue)
      EXTERNAL ff_function_prefix(VectorAllReduce)
      EXTERNAL ff_function_prefix(VectorDot)
      EXTERNAL ff_function_prefix(VectorElement)
      EXTERNAL ff_function_prefix(VectorUnknown)
      EXTERNAL ff_function_prefix(VectorVertex)
      EXTERNAL ff_function_prefix(VectorNorm)
      EXTERNAL ff_function_prefix(VectorSum)
      

#ifdef USE_MPI
      include 'mpif.h'
#endif


