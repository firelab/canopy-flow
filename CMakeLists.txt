#/*****************************************************************************
# *
# * $Id: CMakeLists.txt 1673 2012-07-23 16:44:53Z $
# *
# * Project:  Canopy
# * Purpose:  CMake script
# * Author:   Jason Forthofer
# *
# *****************************************************************************
# *
# * THIS SOFTWARE WAS DEVELOPED AT THE ROCKY MOUNTAIN RESEARCH STATION (RMRS)
# * MISSOULA FIRE SCIENCES LABORATORY BY EMPLOYEES OF THE FEDERAL GOVERNMENT
# * IN THE COURSE OF THEIR OFFICIAL DUTIES. PURSUANT TO TITLE 17 SECTION 105
# * OF THE UNITED STATES CODE, THIS SOFTWARE IS NOT SUBJECT TO COPYRIGHT
# * PROTECTION AND IS IN THE PUBLIC DOMAIN. RMRS MISSOULA FIRE SCIENCES
# * LABORATORY ASSUMES NO RESPONSIBILITY WHATSOEVER FOR ITS USE BY OTHER
# * PARTIES,  AND MAKES NO GUARANTEES, EXPRESSED OR IMPLIED, ABOUT ITS QUALITY,
# * RELIABILITY, OR ANY OTHER CHARACTERISTIC.
# *
# * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
# * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# * DEALINGS IN THE SOFTWARE.
# *
# ****************************************************************************/
# This is set arbitrarily at 2.6.  If there is proof that it works better, 
# set it.  Also set all sub CMakeLists.txt
# *****************************************************************************
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(canopy_flow)
# *****************************************************************************
# Project canopy_flow
# *****************************************************************************


FIND_LIBRARY(PLPLOT_LIB plplotd)
FIND_LIBRARY(PLPLOT_LIB_CXX plplotcxxd)
FIND_PATH(PLPLOT_INCLUDES plstream.h /usr/include/plplot)

INCLUDE_DIRECTORIES(${PLPLOT_INCLUDES})
INCLUDE_DIRECTORIES(behave)

ADD_EXECUTABLE(canopy_flow src/canopyFlow.cpp
                           src/canopy.cpp
                           src/canopy_asymmetric_gaussian_distribution.cpp
                           src/canopy_normal_distribution.cpp
                           src/canopy_triangle_distribution.cpp
                           src/canopy_uniform_distribution.cpp
                           src/massman_distribution.cpp
                           src/measured_distribution.cpp
                           src/main.cpp)
TARGET_LINK_LIBRARIES(canopy_flow ${PLPLOT_LIB} ${PLPLOT_LIB_CXX})

ADD_EXECUTABLE(compare_Albini   behave/behaveRun.cpp
                                behave/behaveVector.cpp
                                behave/crown.cpp
                                behave/crownInputs.cpp
                                behave/fuelModelSet.cpp
                                behave/newext.cpp
                                behave/palmettoGallberry.cpp
                                behave/randfuel.cpp
                                behave/randthread.cpp
                                behave/surface.cpp
                                behave/surfaceFireReactionIntensity.cpp
                                behave/surfaceFuelbedIntermediates.cpp
                                behave/surfaceInputs.cpp
                                behave/surfaceFireSpread.cpp
                                behave/surfaceTwoFuelModels.cpp
                                behave/westernAspen.cpp
                                behave/windAdjustmentFactor.cpp
                                behave/windSpeedUtility.cpp
                                src/canopyFlow.cpp
                                src/canopy.cpp
                                src/canopy_asymmetric_gaussian_distribution.cpp
                                src/canopy_normal_distribution.cpp
                                src/canopy_triangle_distribution.cpp
                                src/canopy_uniform_distribution.cpp
                                src/massman_distribution.cpp
                                src/measured_distribution.cpp
                                src/main_compare_albini.cpp)
TARGET_LINK_LIBRARIES(compare_Albini ${PLPLOT_LIB} ${PLPLOT_LIB_CXX})

