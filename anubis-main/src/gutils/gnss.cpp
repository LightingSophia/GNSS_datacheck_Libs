
/* ----------------------------------------------------------------------
 * G-Nut - GNSS software development library
 * 
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  
  (c) 2011-2017 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  This file is part of the G-Nut C++ library.
 
  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License as
  published by the Free Software Foundation; either version 3 of
  the License, or (at your option) any later version.
 
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, see <http://www.gnu.org/licenses>.

-*/

#include <iostream>
#include <iomanip>
#include <string>

#include "../gutils/gnss.h"
#include "../gutils/gtypeconv.h"

// using namespace std;

namespace gnut {

// static map initializer
// ---------
t_map_sats GNSS_SATS()
{
  t_map_sats m;
             // (string) necessary for VS2013
  m[GPS] = { (string)"G01",(string)"G02",(string)"G03",(string)"G04",(string)"G05",(string)"G06",(string)"G07",(string)"G08",(string)"G09",(string)"G10",
             (string)"G11",(string)"G12",(string)"G13",(string)"G14",(string)"G15",(string)"G16",(string)"G17",(string)"G18",(string)"G19",(string)"G20",
             (string)"G21",(string)"G22",(string)"G23",(string)"G24",(string)"G25",(string)"G26",(string)"G27",(string)"G28",(string)"G29",(string)"G30",
             (string)"G31",(string)"G32" };
   
  m[GLO] = { (string)"R01",(string)"R02",(string)"R03",(string)"R04",(string)"R05",(string)"R06",(string)"R07",(string)"R08",(string)"R09",(string)"R10",
             (string)"R11",(string)"R12",(string)"R13",(string)"R14",(string)"R15",(string)"R16",(string)"R17",(string)"R18",(string)"R19",(string)"R20",
             (string)"R21",(string)"R22",(string)"R23",(string)"R24" };

  m[GAL] = { (string)"E01",(string)"E02",(string)"E03",(string)"E04",(string)"E05",(string)"E06",(string)"E07",(string)"E08",(string)"E09",(string)"E10",
             (string)"E11",(string)"E12",(string)"E13",(string)"E14",(string)"E15",(string)"E16",(string)"E17",(string)"E18",(string)"E19",(string)"E20",
             (string)"E21",(string)"E22",(string)"E23",(string)"E24",(string)"E25",(string)"E26",(string)"E27",(string)"E28",(string)"E29",(string)"E30",
             (string)"E31",(string)"E32",(string)"E33",(string)"E34",(string)"E35",(string)"E36" };

  m[BDS] = { (string)"C01",(string)"C02",(string)"C03",(string)"C04",(string)"C05",(string)"C06",(string)"C07",(string)"C08",(string)"C09",(string)"C10",
             (string)"C11",(string)"C12",(string)"C13",(string)"C14",(string)"C15",(string)"C16",(string)"C17",(string)"C18",(string)"C19",(string)"C20",
             (string)"C21",(string)"C22",(string)"C23",(string)"C24",(string)"C25",(string)"C26",(string)"C27",(string)"C28",(string)"C29",(string)"C30",
             (string)"C31",(string)"C32",(string)"C33",(string)"C34",(string)"C35",(string)"C36",(string)"C37",(string)"C38",(string)"C39",(string)"C40",
             (string)"C41",(string)"C42",(string)"C43",(string)"C44",(string)"C45",(string)"C46",(string)"C47",(string)"C48",(string)"C49",(string)"C50",
             (string)"C51",(string)"C52",(string)"C53",(string)"C54",(string)"C55",(string)"C56",(string)"C57",(string)"C58",(string)"C59"};

  m[SBS] = { (string)"S20",(string)"S21",(string)"S22",(string)"S23",(string)"S24",(string)"S25",(string)"S26",(string)"S27",(string)"S28",(string)"S29",
             (string)"S30",(string)"S31",(string)"S32",(string)"S33",(string)"S34",(string)"S35",(string)"S36",(string)"S37",(string)"S38",(string)"S39",
             (string)"S40",(string)"S41",(string)"S42" };

  m[QZS] = { (string)"J01",(string)"J02",(string)"J03",(string)"J04",(string)"J05",(string)"J06",(string)"J07" };

  m[IRN] = { (string)"I01",(string)"I02",(string)"I03",(string)"I04",(string)"I05",(string)"I06",(string)"I07" };

  // cummulate all
  for( auto itSYS = m.begin(); itSYS != m.end(); ++itSYS ){
    for( auto it = itSYS->second.begin(); it != itSYS->second.end(); ++it ){
       m[GNS].insert(*it);
    }
  }
  return m;
}


// static map initializer
// ---------
t_map_sats GNSS_GNAV()
{
  t_map_sats m;
             // (string) necessary for VS2013
  m[GPS] = { (string)"NAV",(string)"CNAV"};
  m[GLO] = { (string)"NAV" };
  m[GAL] = { (string)"NAV",(string)"INAV_E01",(string)"INAV_E07",(string)"INAV",(string)"FNAV" };
  m[BDS] = { (string)"NAV" };
  m[SBS] = { (string)"NAV" };
  m[QZS] = { (string)"NAV" };

  // cummulate all
  for( auto itSYS = m.begin(); itSYS != m.end(); ++itSYS ){
    for( auto it = itSYS->second.begin(); it != itSYS->second.end(); ++it ){
       m[GNS].insert(*it);
    }
  }

  return m;
}


// static map of GNSS frequencies priorities
// ---------
t_map_freq GNSS_FREQ_PRIORITY()
{
  t_map_freq m;
  // ===========================================================
  // priority  --> NEED TO BE CONSISTENT WITH GNSS_BAND_PRIORITY
  // ===========================================================
  // FREQ_SEQ:      X   1st       2nd       3rd       4th       5th          // order is important!
  m[GPS] = { LAST_GFRQ, G01,      G02,      G05                          };
  m[GLO] = { LAST_GFRQ, R01,      R02,      R01_CDMA, R02_CDMA, R03_CDMA };
  m[GAL] = { LAST_GFRQ, E01,      E05,      E07,      E08,      E06      };
  m[BDS] = { LAST_GFRQ, C02,      C06,      C01,      C05,      C07      };
  m[QZS] = { LAST_GFRQ, J01,      J02,      J05,      J06                };
  m[SBS] = { LAST_GFRQ, S01,      S05                                    };
  m[GNS] = {                                                             };

  return m;
}

// static map of GNSS frequencies to bands
// ---------
t_map_band GNSS_BAND_PRIORITY()
{ 
  t_map_band m;
  // ===========================================================
  // priority  --> NEED TO BE CONSISTENT WITH GNSS_FREQ_PRIORITY
  // ===========================================================
  // FREQ_SEQ: X     1st     2nd     3rd     4th    5th          // order is important!
  m[GPS] = { BAND, BAND_1, BAND_2, BAND_5                 };
  m[GLO] = { BAND, BAND_1, BAND_2, BAND_4, BAND_6, BAND_3 };     // to be consistent with above!
  m[GAL] = { BAND, BAND_1, BAND_5, BAND_7, BAND_8, BAND_6 };     // priority, no frequency order
  m[BDS] = { BAND, BAND_2, BAND_6, BAND_1, BAND_5, BAND_7 };
  m[QZS] = { BAND, BAND_1, BAND_2, BAND_5, BAND_6         };
  m[SBS] = { BAND, BAND_1, BAND_5                         };
  m[GNS] = {                                              };

  return m;
}
   

// static map of sorted GNSS band w.r.t. wavelength
// ---------
t_map_band GNSS_BAND_SORTED()
{ 
  t_map_band m;
  // order is important!
  m[GPS] = { BAND_1, BAND_2, BAND_5                 };
  m[GLO] = { BAND_1, BAND_4, BAND_6, BAND_2, BAND_3 };
  m[GAL] = { BAND_1, BAND_6, BAND_7, BAND_8, BAND_5 };    
  m[BDS] = { BAND_2, BAND_6, BAND_1, BAND_5, BAND_7 };
  m[QZS] = { BAND_1, BAND_6, BAND_2, BAND_5         };
  m[SBS] = { BAND_1, BAND_2, BAND_5                 };
  m[GNS] = {                                        };

  return m;
}  
   
// sort set of bands w.r.t. wavelength     
vector<GOBSBAND> sort_band(GSYS gs, set<GOBSBAND>& bands)
{  
  vector<GOBSBAND> vec;
   
  t_map_band bands_sorted = GNSS_BAND_SORTED();
   
  auto itGSYS = bands_sorted.find(gs);
  if(itGSYS != bands_sorted.end()){
    for(auto itBAND = itGSYS->second.begin(); itBAND != itGSYS->second.end(); itBAND++){
      if(bands.find(*itBAND) != bands.end()) vec.push_back(*itBAND);
    }
  }
  
  return vec;
}
   
// static map of GNSS frequencies to bands
// ---------
t_map_offs GNSS_PCO_OFFSETS()
{
 // axes orientations are in the IGS conventions  
   
  t_map_offs m;
                    // triple: ATX  NORTH / EAST / UP
  t_gtriple gpsIIA ( 279, 0, 2500);                      //  279.00,   0.00,   2256.50 - 2878.60 
  t_gtriple gpsIIR (   0, 0, 1000);                      //    0,00,   0.00,    681.10 - 1506.40
  t_gtriple gpsIIF ( 394, 0, 1500);                      //  394.00,   0.00,   1561.30 - 1600.00
  t_gtriple gpsGEN ( 200, 0, 1500);                      // aprox. mean
  m[GPS][BAND_1] = { gpsGEN };
  m[GPS][BAND_2] = { gpsGEN };
  m[GPS][BAND_5] = { gpsGEN };
  m[GPS][BAND_A] = { gpsGEN };
  m[GPS][BAND_B] = { gpsGEN };

  t_gtriple gloI   (   0, 0, 2000);                      //     0.00, 0.00,    1962.00 - 2222.30
  t_gtriple gloII  (-545, 0, 2000);                      //  -545.00, 0.00,    2179.40 - 2743.60
  t_gtriple gloGEN (-545, 0, 2000);                      // approx. mean
  m[GLO][BAND_1] = { gloGEN };
  m[GLO][BAND_2] = { gloGEN };
  m[GLO][BAND_3] = { gloGEN };
  m[GLO][BAND_A] = { gloGEN };
  m[GLO][BAND_B] = { gloGEN };

  t_gtriple galGEN (-170, 30, 950);
  m[GAL][BAND_1] = { galGEN };
  m[GAL][BAND_5] = { galGEN };
  m[GAL][BAND_6] = { galGEN };
  m[GAL][BAND_7] = { galGEN };
  m[GAL][BAND_8] = { galGEN };
  m[GAL][BAND_A] = { galGEN };
  m[GAL][BAND_B] = { galGEN };

  t_gtriple bdsGEN ( 634, 3, 1075 );
//m[BDS][BAND_1] = { bdsGEN };
  m[BDS][BAND_2] = { bdsGEN };
  m[BDS][BAND_6] = { bdsGEN };
  m[BDS][BAND_7] = { bdsGEN };

  t_gtriple sbsGEN ( 0, 0, 0 ); // UNKNOWN!
  m[SBS][BAND_1] = { sbsGEN };
  m[SBS][BAND_5] = { sbsGEN };

  m[QZS][BAND_1] = { t_gtriple(-0.9, 2.9, 3197.9) };
  m[QZS][BAND_2] = { t_gtriple(-0.9, 2.9, 2992.9) };
  m[QZS][BAND_5] = { t_gtriple(-0.9, 2.9, 3077.9) };
  m[QZS][BAND_6] = { t_gtriple(-0.9, 2.9, 3147.9) };

  return m;
}

// supported GNSS   
set<GSYS> GNSS_SUPPORTED()
{
   set<GSYS> systems;
   
   systems.insert(GPS);
   systems.insert(GLO);
   systems.insert(GAL);
   systems.insert(BDS);
   systems.insert(QZS);
   systems.insert(IRN);
   
   return systems;
}
   
     
// static map of GNSS data types
// ---------
t_map_gnss GNSS_DATA_PRIORITY()
{
  t_map_gnss m;
   
  m[GPS][BAND_1][TYPE_C] = { ATTR_C, ATTR_S, ATTR_L, ATTR_X, ATTR_P, ATTR_W, ATTR_Y, ATTR_M,         ATTR_NULL };
  m[GPS][BAND_1][TYPE_L] = { ATTR_C, ATTR_S, ATTR_L, ATTR_X, ATTR_P, ATTR_W, ATTR_Y, ATTR_M, ATTR_N, ATTR_NULL };
  m[GPS][BAND_1][TYPE_D] = { ATTR_C, ATTR_S, ATTR_L, ATTR_X, ATTR_P, ATTR_W, ATTR_Y, ATTR_M, ATTR_N, ATTR_NULL };
  m[GPS][BAND_1][TYPE_S] = { ATTR_C, ATTR_S, ATTR_L, ATTR_X, ATTR_P, ATTR_W, ATTR_Y, ATTR_M, ATTR_N, ATTR_NULL };
  m[GPS][BAND_1][TYPE_P] = { ATTR_NULL };
   
  m[GPS][BAND_2][TYPE_C] = { ATTR_C, ATTR_D, ATTR_S, ATTR_L, ATTR_X, ATTR_P, ATTR_W, ATTR_Y, ATTR_M,         ATTR_NULL };
  m[GPS][BAND_2][TYPE_L] = { ATTR_C, ATTR_D, ATTR_S, ATTR_L, ATTR_X, ATTR_P, ATTR_W, ATTR_Y, ATTR_M, ATTR_N, ATTR_NULL };
  m[GPS][BAND_2][TYPE_D] = { ATTR_C, ATTR_D, ATTR_S, ATTR_L, ATTR_X, ATTR_P, ATTR_W, ATTR_Y, ATTR_M, ATTR_N, ATTR_NULL };
  m[GPS][BAND_2][TYPE_S] = { ATTR_C, ATTR_D, ATTR_S, ATTR_L, ATTR_X, ATTR_P, ATTR_W, ATTR_Y, ATTR_M, ATTR_N, ATTR_NULL };
  m[GPS][BAND_2][TYPE_P] = { ATTR_NULL };

  m[GPS][BAND_5][TYPE_C] = { ATTR_I, ATTR_Q, ATTR_X };
  m[GPS][BAND_5][TYPE_L] = { ATTR_I, ATTR_Q, ATTR_X };
  m[GPS][BAND_5][TYPE_D] = { ATTR_I, ATTR_Q, ATTR_X };
  m[GPS][BAND_5][TYPE_S] = { ATTR_I, ATTR_Q, ATTR_X };

  m[GPS][BAND_A][TYPE_C] = { ATTR_NULL };
  m[GPS][BAND_A][TYPE_L] = { ATTR_NULL };
  m[GPS][BAND_A][TYPE_D] = { ATTR_NULL };
  m[GPS][BAND_A][TYPE_S] = { ATTR_NULL };

  m[GPS][BAND_B][TYPE_C] = { ATTR_NULL };
  m[GPS][BAND_B][TYPE_L] = { ATTR_NULL };
  m[GPS][BAND_B][TYPE_D] = { ATTR_NULL };
  m[GPS][BAND_B][TYPE_S] = { ATTR_NULL };


  m[GLO][BAND_1][TYPE_C] = { ATTR_C, ATTR_P, ATTR_NULL };
  m[GLO][BAND_1][TYPE_L] = { ATTR_C, ATTR_P, ATTR_NULL };
  m[GLO][BAND_1][TYPE_D] = { ATTR_C, ATTR_P, ATTR_NULL };
  m[GLO][BAND_1][TYPE_S] = { ATTR_C, ATTR_P, ATTR_NULL };
  m[GLO][BAND_1][TYPE_P] = { ATTR_NULL };

  m[GLO][BAND_2][TYPE_C] = { ATTR_C, ATTR_P, ATTR_NULL };
  m[GLO][BAND_2][TYPE_L] = { ATTR_C, ATTR_P, ATTR_NULL };
  m[GLO][BAND_2][TYPE_D] = { ATTR_C, ATTR_P, ATTR_NULL };
  m[GLO][BAND_2][TYPE_S] = { ATTR_C, ATTR_P, ATTR_NULL };
  m[GLO][BAND_2][TYPE_P] = { ATTR_NULL };

  m[GLO][BAND_3][TYPE_C] = { ATTR_I, ATTR_Q, ATTR_X };
  m[GLO][BAND_3][TYPE_L] = { ATTR_I, ATTR_Q, ATTR_X };
  m[GLO][BAND_3][TYPE_D] = { ATTR_I, ATTR_Q, ATTR_X };
  m[GLO][BAND_3][TYPE_S] = { ATTR_I, ATTR_Q, ATTR_X };
   
  m[GLO][BAND_4][TYPE_C] = { ATTR_A, ATTR_B, ATTR_X };
  m[GLO][BAND_4][TYPE_L] = { ATTR_A, ATTR_B, ATTR_X };
  m[GLO][BAND_4][TYPE_D] = { ATTR_A, ATTR_B, ATTR_X };
  m[GLO][BAND_4][TYPE_S] = { ATTR_A, ATTR_B, ATTR_X };
   
  m[GLO][BAND_5][TYPE_C] = { ATTR_A, ATTR_B, ATTR_X };
  m[GLO][BAND_5][TYPE_L] = { ATTR_A, ATTR_B, ATTR_X };
  m[GLO][BAND_5][TYPE_D] = { ATTR_A, ATTR_B, ATTR_X };
  m[GLO][BAND_5][TYPE_S] = { ATTR_A, ATTR_B, ATTR_X };
   
  m[GLO][BAND_A][TYPE_C] = { ATTR_NULL };
  m[GLO][BAND_A][TYPE_L] = { ATTR_NULL };
  m[GLO][BAND_A][TYPE_D] = { ATTR_NULL };
  m[GLO][BAND_A][TYPE_S] = { ATTR_NULL };

  m[GLO][BAND_B][TYPE_C] = { ATTR_NULL };
  m[GLO][BAND_B][TYPE_L] = { ATTR_NULL };
  m[GLO][BAND_B][TYPE_D] = { ATTR_NULL };
  m[GLO][BAND_B][TYPE_S] = { ATTR_NULL };
   

  m[GAL][BAND_1][TYPE_C] = { ATTR_A, ATTR_B, ATTR_C, ATTR_X, ATTR_Z };
  m[GAL][BAND_1][TYPE_L] = { ATTR_A, ATTR_B, ATTR_C, ATTR_X, ATTR_Z };
  m[GAL][BAND_1][TYPE_D] = { ATTR_A, ATTR_B, ATTR_C, ATTR_X, ATTR_Z };
  m[GAL][BAND_1][TYPE_S] = { ATTR_A, ATTR_B, ATTR_C, ATTR_X, ATTR_Z };

  m[GAL][BAND_5][TYPE_C] = { ATTR_I, ATTR_Q, ATTR_X };
  m[GAL][BAND_5][TYPE_L] = { ATTR_I, ATTR_Q, ATTR_X };
  m[GAL][BAND_5][TYPE_D] = { ATTR_I, ATTR_Q, ATTR_X };
  m[GAL][BAND_5][TYPE_S] = { ATTR_I, ATTR_Q, ATTR_X };

  m[GAL][BAND_6][TYPE_C] = { ATTR_A, ATTR_B, ATTR_C, ATTR_X, ATTR_Z };
  m[GAL][BAND_6][TYPE_L] = { ATTR_A, ATTR_B, ATTR_C, ATTR_X, ATTR_Z };
  m[GAL][BAND_6][TYPE_D] = { ATTR_A, ATTR_B, ATTR_C, ATTR_X, ATTR_Z };
  m[GAL][BAND_6][TYPE_S] = { ATTR_A, ATTR_B, ATTR_C, ATTR_X, ATTR_Z };

  m[GAL][BAND_7][TYPE_C] = { ATTR_I, ATTR_Q, ATTR_X };
  m[GAL][BAND_7][TYPE_L] = { ATTR_I, ATTR_Q, ATTR_X };
  m[GAL][BAND_7][TYPE_D] = { ATTR_I, ATTR_Q, ATTR_X };
  m[GAL][BAND_7][TYPE_S] = { ATTR_I, ATTR_Q, ATTR_X };

  m[GAL][BAND_8][TYPE_C] = { ATTR_I, ATTR_Q, ATTR_X };
  m[GAL][BAND_8][TYPE_L] = { ATTR_I, ATTR_Q, ATTR_X };
  m[GAL][BAND_8][TYPE_D] = { ATTR_I, ATTR_Q, ATTR_X };
  m[GAL][BAND_8][TYPE_S] = { ATTR_I, ATTR_Q, ATTR_X };

  m[GAL][BAND_A][TYPE_C] = { ATTR_NULL };
  m[GAL][BAND_A][TYPE_L] = { ATTR_NULL };
  m[GAL][BAND_A][TYPE_D] = { ATTR_NULL };
  m[GAL][BAND_A][TYPE_S] = { ATTR_NULL };

  m[GAL][BAND_B][TYPE_C] = { ATTR_NULL };
  m[GAL][BAND_B][TYPE_L] = { ATTR_NULL };
  m[GAL][BAND_B][TYPE_D] = { ATTR_NULL };
  m[GAL][BAND_B][TYPE_S] = { ATTR_NULL };


  m[BDS][BAND_1][TYPE_C] = { ATTR_D, ATTR_P, ATTR_X, ATTR_A,        };
  m[BDS][BAND_1][TYPE_L] = { ATTR_D, ATTR_P, ATTR_X, ATTR_A, ATTR_N };
  m[BDS][BAND_1][TYPE_D] = { ATTR_D, ATTR_P, ATTR_X, ATTR_A, ATTR_N };
  m[BDS][BAND_1][TYPE_S] = { ATTR_D, ATTR_P, ATTR_X, ATTR_A, ATTR_N };

  m[BDS][BAND_2][TYPE_C] = { ATTR_I, ATTR_Q, ATTR_X };
  m[BDS][BAND_2][TYPE_L] = { ATTR_I, ATTR_Q, ATTR_X };
  m[BDS][BAND_2][TYPE_D] = { ATTR_I, ATTR_Q, ATTR_X };
  m[BDS][BAND_2][TYPE_S] = { ATTR_I, ATTR_Q, ATTR_X };

  m[BDS][BAND_5][TYPE_C] = { ATTR_D, ATTR_P, ATTR_X };
  m[BDS][BAND_5][TYPE_L] = { ATTR_D, ATTR_P, ATTR_X };
  m[BDS][BAND_5][TYPE_D] = { ATTR_D, ATTR_P, ATTR_X };
  m[BDS][BAND_5][TYPE_S] = { ATTR_D, ATTR_P, ATTR_X };

  m[BDS][BAND_6][TYPE_C] = { ATTR_I, ATTR_Q, ATTR_X, ATTR_A };
  m[BDS][BAND_6][TYPE_L] = { ATTR_I, ATTR_Q, ATTR_X, ATTR_A };
  m[BDS][BAND_6][TYPE_D] = { ATTR_I, ATTR_Q, ATTR_X, ATTR_A };
  m[BDS][BAND_6][TYPE_S] = { ATTR_I, ATTR_Q, ATTR_X, ATTR_A };

  m[BDS][BAND_7][TYPE_C] = { ATTR_D, ATTR_I, ATTR_P, ATTR_Q, ATTR_X, ATTR_Z };
  m[BDS][BAND_7][TYPE_L] = { ATTR_D, ATTR_I, ATTR_P, ATTR_Q, ATTR_X, ATTR_Z };
  m[BDS][BAND_7][TYPE_D] = { ATTR_D, ATTR_I, ATTR_P, ATTR_Q, ATTR_X, ATTR_Z };
  m[BDS][BAND_7][TYPE_S] = { ATTR_D, ATTR_I, ATTR_P, ATTR_Q, ATTR_X, ATTR_Z };

  m[BDS][BAND_7][TYPE_C] = { ATTR_D, ATTR_P, ATTR_X };
  m[BDS][BAND_7][TYPE_L] = { ATTR_D, ATTR_P, ATTR_X };
  m[BDS][BAND_7][TYPE_D] = { ATTR_D, ATTR_P, ATTR_X };
  m[BDS][BAND_7][TYPE_S] = { ATTR_D, ATTR_P, ATTR_X };


  m[SBS][BAND_1][TYPE_C] = { ATTR_C };
  m[SBS][BAND_1][TYPE_L] = { ATTR_C };
  m[SBS][BAND_1][TYPE_D] = { ATTR_C };
  m[SBS][BAND_1][TYPE_S] = { ATTR_C };
  m[SBS][BAND_1][TYPE_P] = { ATTR_NULL };

  m[SBS][BAND_5][TYPE_C] = { ATTR_I, ATTR_Q, ATTR_X };
  m[SBS][BAND_5][TYPE_L] = { ATTR_I, ATTR_Q, ATTR_X };
  m[SBS][BAND_5][TYPE_D] = { ATTR_I, ATTR_Q, ATTR_X };
  m[SBS][BAND_5][TYPE_S] = { ATTR_I, ATTR_Q, ATTR_X };

   
  m[QZS][BAND_1][TYPE_C] = { ATTR_C, ATTR_S, ATTR_L, ATTR_X, ATTR_Z };
  m[QZS][BAND_1][TYPE_L] = { ATTR_C, ATTR_S, ATTR_L, ATTR_X, ATTR_Z };
  m[QZS][BAND_1][TYPE_D] = { ATTR_C, ATTR_S, ATTR_L, ATTR_X, ATTR_Z };
  m[QZS][BAND_1][TYPE_S] = { ATTR_C, ATTR_S, ATTR_L, ATTR_X, ATTR_Z };

  m[QZS][BAND_2][TYPE_C] = { ATTR_S, ATTR_L, ATTR_X };
  m[QZS][BAND_2][TYPE_L] = { ATTR_S, ATTR_L, ATTR_X };
  m[QZS][BAND_2][TYPE_D] = { ATTR_S, ATTR_L, ATTR_X };
  m[QZS][BAND_2][TYPE_S] = { ATTR_S, ATTR_L, ATTR_X };

  m[QZS][BAND_5][TYPE_C] = { ATTR_I, ATTR_Q, ATTR_X, ATTR_D, ATTR_P, ATTR_Z };
  m[QZS][BAND_5][TYPE_L] = { ATTR_I, ATTR_Q, ATTR_X, ATTR_D, ATTR_P, ATTR_Z };
  m[QZS][BAND_5][TYPE_D] = { ATTR_I, ATTR_Q, ATTR_X, ATTR_D, ATTR_P, ATTR_Z };
  m[QZS][BAND_5][TYPE_S] = { ATTR_I, ATTR_Q, ATTR_X, ATTR_D, ATTR_P, ATTR_Z };
   
  m[QZS][BAND_6][TYPE_C] = { ATTR_S, ATTR_L, ATTR_X, ATTR_E, ATTR_Z };
  m[QZS][BAND_6][TYPE_L] = { ATTR_S, ATTR_L, ATTR_X, ATTR_E, ATTR_Z };
  m[QZS][BAND_6][TYPE_D] = { ATTR_S, ATTR_L, ATTR_X, ATTR_E, ATTR_Z };
  m[QZS][BAND_6][TYPE_S] = { ATTR_S, ATTR_L, ATTR_X, ATTR_E, ATTR_Z };


  m[IRN][BAND_1][TYPE_C] = { ATTR_A, ATTR_B, ATTR_C, ATTR_X, ATTR_Z };
  m[IRN][BAND_1][TYPE_L] = { ATTR_A, ATTR_B, ATTR_C, ATTR_X, ATTR_Z };
  m[IRN][BAND_1][TYPE_D] = { ATTR_A, ATTR_B, ATTR_C, ATTR_X, ATTR_Z };
  m[IRN][BAND_1][TYPE_S] = { ATTR_A, ATTR_B, ATTR_C, ATTR_X, ATTR_Z };

  return m;
}
   
   
// get obs type from char
// ----------
GOBSTYPE char2gobstype( char c )
{
  if(      c == 'C' ) return TYPE_C;
  else if( c == 'L' ) return TYPE_L;
  else if( c == 'D' ) return TYPE_D;
  else if( c == 'S' ) return TYPE_S;
  else if( c == 'P' ) return TYPE_P; // P-code only

  return TYPE;
}

// get GOBSTYPE enum from GOBS enum
GOBSTYPE gobs2gobstype( GOBS o )
{
  string go = gobs2str(o);
  return char2gobstype( go[0] );
}

// get obs band from char
// ----------
GOBSBAND char2gobsband( char c )
{
  if(      c == '1' ) return BAND_1;
  else if( c == '2' ) return BAND_2;
  else if( c == '3' ) return BAND_3;
  else if( c == '4' ) return BAND_4;
  else if( c == '5' ) return BAND_5;
  else if( c == '6' ) return BAND_6;
  else if( c == '7' ) return BAND_7;
  else if( c == '8' ) return BAND_8;
   
  else if( c == 'A' ) return BAND_A;
  else if( c == 'B' ) return BAND_B;
  else if( c == 'C' ) return BAND_C;
  else if( c == 'D' ) return BAND_D;

  return BAND;
}

// get obs band from int
// ----------
GOBSBAND int2gobsband( int c )
{
  if(      c == 1 ) return BAND_1;
  else if( c == 2 ) return BAND_2;
  else if( c == 3 ) return BAND_3;
  else if( c == 4 ) return BAND_4;
  else if( c == 5 ) return BAND_5;
  else if( c == 6 ) return BAND_6;
  else if( c == 7 ) return BAND_7;
  else if( c == 8 ) return BAND_8;   

  return BAND;
}

// get GOBSBAND enum from GOBS enum
GOBSBAND gobs2gobsband( GOBS o )
{
  string go = gobs2str(o);
  return char2gobsband( go[1] );
}

// get obs attr from char
// ----------
GOBSATTR char2gobsattr( char c )
{
  if(      c == 'A' ) return ATTR_A;
  else if( c == 'B' ) return ATTR_B;
  else if( c == 'C' ) return ATTR_C;
  else if( c == 'D' ) return ATTR_D;
  else if( c == 'E' ) return ATTR_E;
//else if( c == 'F' ) return ATTR_F; // NOT VALID
//else if( c == 'G' ) return ATTR_G; // NOT VALID
//else if( c == 'H' ) return ATTR_H; // NOT VALID
  else if( c == 'I' ) return ATTR_I;
//else if( c == 'J' ) return ATTR_J; // NOT VALID
//else if( c == 'K' ) return ATTR_K; // NOT VALID
  else if( c == 'L' ) return ATTR_L;
  else if( c == 'M' ) return ATTR_M;
  else if( c == 'N' ) return ATTR_N;
//else if( c == 'O' ) return ATTR_O; // NOT VALID
  else if( c == 'P' ) return ATTR_P;
  else if( c == 'Q' ) return ATTR_Q;
//else if( c == 'R' ) return ATTR_R; // NOT VALID
  else if( c == 'S' ) return ATTR_S;
//else if( c == 'T' ) return ATTR_T; // NOT VALID
//else if( c == 'U' ) return ATTR_U; // NOT VALID
//else if( c == 'V' ) return ATTR_V; // NOT VALID
  else if( c == 'W' ) return ATTR_W;
  else if( c == 'X' ) return ATTR_X;
  else if( c == 'Y' ) return ATTR_Y;
  else if( c == 'Z' ) return ATTR_Z;
  else if( c == ' ' ) return ATTR_NULL; // 3CH conversion
  else if( c == '\0') return ATTR_NULL; // 2CH conversion

  return ATTR; // undefined
}

// get GOBSATTR enum from GOBS enum
GOBSATTR gobs2gobsattr( GOBS o )
{
  string go = gobs2str(o);
  return char2gobsattr( go[2] );
}

// get str type from enum
// ----------
string gobstype2str( GOBSTYPE e )
{
  switch( e ){
   case TYPE_C : return "C";
   case TYPE_L : return "L";
   case TYPE_D : return "D";
   case TYPE_S : return "S";
   case TYPE_P : return "P";

   case TYPE   : return "X";
   default     : return "X";
  }
  return "X";
}

// get str band from enum
// ----------
string gobsband2str( GOBSBAND e )
{
  switch( e ){
   case BAND_1 : return "1";
   case BAND_2 : return "2";
   case BAND_3 : return "3";
   case BAND_4 : return "4";
   case BAND_5 : return "5";
   case BAND_6 : return "6";
   case BAND_7 : return "7";
   case BAND_8 : return "8";

   case BAND_A : return "A";
   case BAND_B : return "B";
   case BAND_C : return "C";
   case BAND_D : return "D";

   case BAND   : return "X";
   default     : return "X";
  }
  return "X";
}


// get str attr from enum
// ----------
string gobsattr2str( GOBSATTR e )
{
  switch( e ){
   case ATTR_A : return "A";
   case ATTR_B : return "B";
   case ATTR_C : return "C";
   case ATTR_D : return "D";
   case ATTR_E : return "E";
   case ATTR_I : return "I";
   case ATTR_L : return "L";
   case ATTR_M : return "M";
   case ATTR_N : return "N";
   case ATTR_P : return "P";
   case ATTR_Q : return "Q";
   case ATTR_S : return "S";
   case ATTR_W : return "W";
   case ATTR_X : return "X";
   case ATTR_Y : return "Y";
   case ATTR_Z : return "Z";

   case ATTR_NULL : return  " ";

   case ATTR   : return "";
   default     : return "";
  }
  return "";
}


// convert GOBS to observation code string
// -----------
string gobs2str( GOBS o )
{
// PROBLEM WITH STATIC FUNCTION
//  boost::mutex::scoped_lock lock(_mutex);
   
  switch( o ){

    // CODE
    case C1A : return "C1A";  
    case C1B : return "C1B";  
    case C1C : return "C1C"; 
    case C1I : return "C1I";  
    case C1L : return "C1L";  
    case C1M : return "C1M";  
    case C1P : return "C1P"; 
    case C1Q : return "C1Q";  
    case C1S : return "C1S"; 
    case C1W : return "C1W";  
    case C1X : return "C1X";  
    case C1Y : return "C1Y";  
    case C1Z : return "C1Z";  
     
    case C2C : return "C2C";  
    case C2D : return "C2D";  
    case C2I : return "C2I";  
    case C2L : return "C2L";  
    case C2M : return "C2M";  
    case C2P : return "C2P";  
    case C2Q : return "C2Q";  
    case C2S : return "C2S";  
    case C2W : return "C2W";  
    case C2X : return "C2X";  
    case C2Y : return "C2Y";  
     
    case C3D : return "C3D";  
    case C3I : return "C3I";  
    case C3P : return "C3P";  
    case C3Q : return "C3Q";  
    case C3X : return "C3X";

    case C4A : return "C4A";  
    case C4B : return "C4B";
    case C4X : return "C4X";

    case C5A : return "C5A";  
    case C5B : return "C5B";  
    case C5C : return "C5C";
    case C5D : return "C5D";
    case C5I : return "C5I";  
    case C5P : return "C5P";
    case C5Q : return "C5Q";  
    case C5X : return "C5X";  
    case C5Z : return "C5Z";  

    case C6A : return "C6A";  
    case C6B : return "C6B";  
    case C6C : return "C6C";  
    case C6E : return "C6E";
    case C6I : return "C6I";  
    case C6L : return "C6L";  
    case C6S : return "C6S";
    case C6Q : return "C6Q"; 
    case C6X : return "C6X";  
    case C6Z : return "C6Z"; 
     
    case C7D : return "C7D";
    case C7I : return "C7I";  
    case C7P : return "C7P";
    case C7Q : return "C7Q"; 
    case C7X : return "C7X"; 
    case C7Z : return "C7Z";
     
    case C8D : return "C8D";  
    case C8I : return "C8I";  
    case C8P : return "C8P";  
    case C8Q : return "C8Q";  
    case C8X : return "C8X"; 
     

    // CARRIER-PHASE
    case L1A : return "L1A";  
    case L1B : return "L1B";  
    case L1C : return "L1C";  
    case L1I : return "L1I";
    case L1L : return "L1L";  
    case L1M : return "L1M";  
    case L1N : return "L1N";  
    case L1P : return "L1P";  
    case L1Q : return "L1Q";
    case L1S : return "L1S";  
    case L1W : return "L1W";  
    case L1X : return "L1X";  
    case L1Y : return "L1Y";  
    case L1Z : return "L1Z"; 
     
    case L2C : return "L2C";  
    case L2D : return "L2D";  
    case L2I : return "L2I"; 
    case L2L : return "L2L";  
    case L2M : return "L2M";  
    case L2N : return "L2N"; 
    case L2P : return "L2P";  
    case L2Q : return "L2Q"; 
    case L2S : return "L2S";  
    case L2W : return "L2W"; 
    case L2X : return "L2X";  
    case L2Y : return "L2Y";  

    case L3D : return "L3D";
    case L3I : return "L3I";  
    case L3P : return "L3P";  
    case L3Q : return "L3Q";  
    case L3X : return "L3X";
     
    case L4A : return "L4A";  
    case L4B : return "L4B";
    case L4X : return "L4X";

    case L5A : return "L5A";  
    case L5B : return "L5B";  
    case L5C : return "L5C";
    case L5D : return "L5D";
    case L5I : return "L5I";  
    case L5P : return "L5P";
    case L5Q : return "L5Q";  
    case L5X : return "L5X"; 
    case L5Z : return "L5Z";  

    case L6A : return "L6A";  
    case L6B : return "L6B";  
    case L6C : return "L6C";  
    case L6E : return "L6E";
    case L6I : return "L6I";  
    case L6L : return "L6L";  
    case L6S : return "L6S";
    case L6Q : return "L6Q";  
    case L6X : return "L6X";  
    case L6Z : return "L6Z";  

    case L7D : return "L7D";
    case L7I : return "L7I";  
    case L7P : return "L7P";
    case L7Q : return "L7Q";  
    case L7X : return "L7X";  
    case L7Z : return "L7Z";
     
    case L8D : return "L8D";  
    case L8I : return "L8I";  
    case L8P : return "L8P";  
    case L8Q : return "L8Q";  
    case L8X : return "L8X";  

    // DOPLER
    case D1A : return "D1A";  
    case D1B : return "D1B";  
    case D1C : return "D1C";  
    case D1I : return "D1I";
    case D1L : return "D1L";  
    case D1M : return "D1M";  
    case D1N : return "D1N";  
    case D1P : return "D1P";  
    case D1Q : return "D1Q";
    case D1S : return "D1S";  
    case D1W : return "D1W";  
    case D1X : return "D1X";  
    case D1Y : return "D1Y";  
    case D1Z : return "D1Z"; 
     
    case D2C : return "D2C";  
    case D2D : return "D2D";  
    case D2I : return "D2I"; 
    case D2L : return "D2L";  
    case D2M : return "D2M";  
    case D2N : return "D2N";  
    case D2P : return "D2P";  
    case D2Q : return "D2Q";  
    case D2S : return "D2S";  
    case D2W : return "D2W";  
    case D2X : return "D2X";  
    case D2Y : return "D2Y"; 

    case D3D : return "D3D";
    case D3I : return "D3I";  
    case D3P : return "D3P";  
    case D3Q : return "D3Q";  
    case D3X : return "D3X";

    case D4A : return "D4A";  
    case D4B : return "D4B";
    case D4X : return "D4X";

    case D5A : return "D5A";  
    case D5B : return "D5B";  
    case D5C : return "D5C";  
    case D5D : return "D5D";
    case D5I : return "D5I";  
    case D5P : return "D5P";
    case D5Q : return "D5Q";  
    case D5X : return "D5X"; 
    case D5Z : return "D5Z";  
     
    case D6A : return "D6A";  
    case D6B : return "D6B";  
    case D6C : return "D6C";  
    case D6E : return "D6E";
    case D6I : return "D6I";  
    case D6L : return "D6L";  
    case D6S : return "D6S";
    case D6Q : return "D6Q"; 
    case D6X : return "D6X";  
    case D6Z : return "D6Z";  
     
    case D7D : return "D7D";
    case D7I : return "D7I";  
    case D7P : return "D7P";
    case D7Q : return "D7Q";  
    case D7X : return "D7X";  
    case D7Z : return "D7Z";
     
    case D8D : return "D8D";
    case D8I : return "D8I";  
    case D8P : return "D8P";
    case D8Q : return "D8Q";  
    case D8X : return "D8X";  
                
    // SNR
    case S1A : return "S1A";  
    case S1B : return "S1B";  
    case S1C : return "S1C";  
    case S1I : return "S1I";
    case S1L : return "S1L";  
    case S1M : return "S1M";  
    case S1N : return "S1N";  
    case S1P : return "S1P";  
    case S1Q : return "S1Q";
    case S1S : return "S1S";  
    case S1W : return "S1W";  
    case S1X : return "S1X";  
    case S1Y : return "S1Y";  
    case S1Z : return "S1Z";  
     
    case S2C : return "S2C";  
    case S2D : return "S2D";  
    case S2I : return "S2I"; 
    case S2L : return "S2L";  
    case S2M : return "S2M";  
    case S2N : return "S2N";  
    case S2P : return "S2P";  
    case S2Q : return "S2Q";  
    case S2S : return "S2S";  
    case S2W : return "S2W";  
    case S2X : return "S2X";  
    case S2Y : return "S2Y";  

    case S3D : return "S3D";
    case S3I : return "S3I";  
    case S3P : return "S3P";
    case S3Q : return "S3Q";  
    case S3X : return "S3X";

    case S4A : return "S4A";  
    case S4B : return "S4B";
    case S4X : return "S4X";

    case S5A : return "S5A";  
    case S5B : return "S5B";  
    case S5C : return "S5C";  
    case S5D : return "S5D";
    case S5I : return "S5I";  
    case S5P : return "S5P";
    case S5Q : return "S5Q";  
    case S5X : return "S5X";
    case S5Z : return "S5Z";
     
    case S6A : return "S6A"; 
    case S6B : return "S6B"; 
    case S6C : return "S6C"; 
    case S6E : return "S6E"; 
    case S6I : return "S6I"; 
    case S6L : return "S6L";  
    case S6S : return "S6S";
    case S6Q : return "S6Q"; 
    case S6X : return "S6X"; 
    case S6Z : return "S6Z"; 
          
    case S7D : return "S7D";
    case S7I : return "S7I";  
    case S7P : return "S7P";
    case S7Q : return "S7Q";  
    case S7X : return "S7X";  
    case S7Z : return "S7Z";

    case S8D : return "S8D";
    case S8I : return "S8I";  
    case S8P : return "D8P";
    case S8Q : return "S8Q";  
    case S8X : return "S8X"; 
     
    // special cases: v2.x or unknown tracking mode
    case C1  : return "C1 ";    // unknown mode GPS P,Y              Galileo E2-L1-E2 I/Q
    case C2  : return "C2 ";    // unknown mode GPS P,Y                                     COMPASS E2 I/Q
    case C5  : return "C5 ";    // unknown mode                      Galileo E5a I/Q
    case C6  : return "C6 ";    // unknown mode                      Galileo E6  I/Q        COMPASS E6 I/Q
    case C7  : return "C7 ";    // unknown mode                      Galielo E5b I/Q        COMPASS E5b I/Q
    case C8  : return "C8 ";    // unknown mode                      Galileo E5  I/Q
    case CA  : return "CA ";    // unknown mode GPS C/A  GLONASS SA
    case CB  : return "CB ";    // unknown mode GPS C1C
    case CC  : return "CC ";    // unknown mode GPS C2C
    case CD  : return "CD ";    // unknown mode          GLONASS SA
                  
    case P1  : return "P1 ";    // unknown mode GPS P,Y  GLONASS HA
    case P2  : return "P2 ";    // unknown mode GPS P,Y  GLONASS HA
    case P5  : return "P5 ";    // unknown mode GPS P,Y
                  
    case L1  : return "L1 ";    // unknown mode GPS P,Y              Galileo E2-L1-E2 I/Q
    case L2  : return "L2 ";    // unknown mode GPS P,Y                                     COMPASS E2 I/Q
    case L5  : return "L5 ";    // unknown mode                      Galileo E5a I/Q
    case L6  : return "L6 ";    // unknown mode                      Galileo E6  I/Q        COMPASS E6 I/Q
    case L7  : return "L7 ";    // unknown mode                      Galielo E5b I/Q        COMPASS E5b I/Q
    case L8  : return "L8 ";    // unknown mode                      Galileo E5  I/Q
    case LA  : return "LA ";    // unknown mode GPS C/A  GLONASS SA
    case LB  : return "LB ";    // unknown mode GPS C1C
    case LC  : return "LC ";    // unknown mode GPS C2C
    case LD  : return "LD ";    // unknown mode          GLONASS SA
                    
    case D1  : return "D1 ";    // unknown mode GPS P,Y              Galileo E2-L1-E2 I/Q
    case D2  : return "D2 ";    // unknown mode GPS P,Y                                     COMPASS E2 I/Q
    case D5  : return "D5 ";    // unknown mode                      Galileo E5a I/Q
    case D6  : return "D6 ";    // unknown mode                      Galileo E6  I/Q        COMPASS E6 I/Q
    case D7  : return "D7 ";    // unknown mode                      Galielo E5b I/Q        COMPASS E5b I/Q
    case D8  : return "D8 ";    // unknown mode                      Galileo E5  I/Q
    case DA  : return "DA ";    // unknown mode GPS C/A  GLONASS SA
    case DB  : return "DB ";    // unknown mode GPS C1C
    case DC  : return "DC ";    // unknown mode GPS C2C
    case DD  : return "DD ";    // unknown mode          GLONASS SA
                    
    case S1  : return "S1 ";    // unknown mode GPS P,Y              Galileo E2-L1-E2 I/Q
    case S2  : return "S2 ";    // unknown mode GPS P,Y                                     COMPASS E2 I/Q
    case S5  : return "S5 ";    // unknown mode                      Galileo E5a I/Q
    case S6  : return "S6 ";    // unknown mode                      Galileo E6  I/Q        COMPASS E6 I/Q
    case S7  : return "S7 ";    // unknown mode                      Galielo E5b I/Q        COMPASS E5b I/Q
    case S8  : return "S8 ";    // unknown mode                      Galileo E5  I/Q
    case SA  : return "SA ";    // unknown mode GPS C/A  GLONASS SA
    case SB  : return "SB ";    // unknown mode GPS C1C
    case SC  : return "SC ";    // unknown mode GPS C2C
    case SD  : return "SD ";    // unknown mode          GLONASS SA

    case X1  : return "X1 ";    // channel number

//    case LAST_GOBS: return "XXX";
    case X : return "XXX";
              
    default : cerr << "*** Warning - unknown observation code " << o << "\n";
              // cannot be while static function !!!
              // if( _log ) _log->comment(0,"gobsgnss","Warning - unknown observation code ");
              // else                      cerr << "*** Warning - unknown observation code \n"; 
  }

 return "XXX";
}

// convert GOBS_LC to observation code string
// -----------
string gobslc2str( GOBS_LC o)
{
  switch( o ){
    case LC_L1 : return "L1";
    case LC_L2 : return "L2";
    case LC_L3 : return "L3";
    case LC_L4 : return "L4";
    case LC_L5 : return "L5";
    case LC_C1 : return "C1";
    case LC_C2 : return "C2";
    case LC_C3 : return "C3";
    case LC_C4 : return "C4";
    case LC_C5 : return "C5";
    case LC_CIF : return "CIF";
    case LC_LIF : return "LIF";
    case LC_UNDEF : return "LC_UNDEF";
  }

  return "";
}

// convert obs code string to GOBS
// -----------
GOBS str2gobs( string s )
{
// STATIC FUNCTION !!!
//  boost::mutex::scoped_lock lock(_mutex);

//  if( obs.compare("XXX") == 0 ){ return X; }

  if( s.length() > 2 && s[2] == ' ' ) s = s.substr(0,2);

  int len = trim(s).length();
  if( len == 3 ){

    if( s[0] == 'C' ){ // code
      if( s[1] == '1' ){ // band
        if( s[2] == 'A' ) return C1A;
        if( s[2] == 'B' ) return C1B;
        if( s[2] == 'C' ) return C1C;
        if( s[2] == 'I' ) return C1I;
        if( s[2] == 'L' ) return C1L;
        if( s[2] == 'M' ) return C1M;
//      if( s[2] == 'N' ) return C1N; // NOT DEFINED
        if( s[2] == 'P' ) return C1P;
        if( s[2] == 'Q' ) return C1Q;
        if( s[2] == 'S' ) return C1S;
        if( s[2] == 'W' ) return C1W;
        if( s[2] == 'X' ) return C1X;
        if( s[2] == 'Y' ) return C1Y;
        if( s[2] == 'Z' ) return C1Z;
      }
      if( s[1] == '2' ){ // band
        if( s[2] == 'C' ) return C2C;
        if( s[2] == 'D' ) return C2D;
        if( s[2] == 'I' ) return C2I;
        if( s[2] == 'L' ) return C2L;
        if( s[2] == 'M' ) return C2M;
//      if( s[2] == 'N' ) return C2N; // NOT DEFINED
        if( s[2] == 'P' ) return C2P;
        if( s[2] == 'Q' ) return C2Q;
        if( s[2] == 'S' ) return C2S;
        if( s[2] == 'W' ) return C2W;
        if( s[2] == 'X' ) return C2X;
        if( s[2] == 'Y' ) return C2Y;
      }
      if( s[1] == '3' ){ // band
        if( s[2] == 'D' ) return C3D;
        if( s[2] == 'I' ) return C3I;
        if( s[2] == 'P' ) return C3P;
        if( s[2] == 'Q' ) return C3Q;
        if( s[2] == 'X' ) return C3X;
      }
      if( s[1] == '5' ){ // band
        if( s[2] == 'A' ) return C5A;
        if( s[2] == 'B' ) return C5B;
        if( s[2] == 'C' ) return C5C;
        if( s[2] == 'D' ) return C5D;
        if( s[2] == 'I' ) return C5I;
        if( s[2] == 'P' ) return C5P;
        if( s[2] == 'Q' ) return C5Q;
        if( s[2] == 'X' ) return C5X;
        if( s[2] == 'Z' ) return C5Z;
      }
      if( s[1] == '6' ){ // band
        if( s[2] == 'A' ) return C6A;
        if( s[2] == 'B' ) return C6B;
        if( s[2] == 'C' ) return C6C;
        if( s[2] == 'E' ) return C6E;
        if( s[2] == 'I' ) return C6I;
        if( s[2] == 'L' ) return C6L;
        if( s[2] == 'S' ) return C6S;
        if( s[2] == 'Q' ) return C6Q;
        if( s[2] == 'X' ) return C6X;
        if( s[2] == 'Z' ) return C6Z;
      }
      if( s[1] == '7' ){ // band
        if( s[2] == 'D' ) return C7D;
        if( s[2] == 'I' ) return C7I;
        if( s[2] == 'P' ) return C7P;
        if( s[2] == 'Q' ) return C7Q;
        if( s[2] == 'X' ) return C7X;
        if( s[2] == 'Z' ) return C7Z;
      }
      if( s[1] == '8' ){ // band
        if( s[2] == 'D' ) return C8D;
        if( s[2] == 'I' ) return C8I;
        if( s[2] == 'P' ) return C8P;
        if( s[2] == 'Q' ) return C8Q;
        if( s[2] == 'X' ) return C8X;
      }
    }
     
    if( s[0] == 'L' ){ // phase
      if( s[1] == '1' ){ // band
        if( s[2] == 'A' ) return L1A;
        if( s[2] == 'B' ) return L1B;
        if( s[2] == 'C' ) return L1C;
        if( s[2] == 'I' ) return L1I;
        if( s[2] == 'L' ) return L1L;
        if( s[2] == 'M' ) return L1M;
        if( s[2] == 'N' ) return L1N;
        if( s[2] == 'P' ) return L1P;
        if( s[2] == 'Q' ) return L1Q;
        if( s[2] == 'S' ) return L1S;
        if( s[2] == 'W' ) return L1W;
        if( s[2] == 'X' ) return L1X;
        if( s[2] == 'Y' ) return L1Y;
        if( s[2] == 'Z' ) return L1Z;
      }
      if( s[1] == '2' ){ // band
        if( s[2] == 'C' ) return L2C;
        if( s[2] == 'D' ) return L2D;
        if( s[2] == 'I' ) return L2I;
        if( s[2] == 'L' ) return L2L;
        if( s[2] == 'M' ) return L2M;
        if( s[2] == 'N' ) return L2N;
        if( s[2] == 'P' ) return L2P;
        if( s[2] == 'Q' ) return L2Q;
        if( s[2] == 'S' ) return L2S;
        if( s[2] == 'W' ) return L2W;
        if( s[2] == 'X' ) return L2X;
        if( s[2] == 'Y' ) return L2Y;
      }
      if( s[1] == '3' ){ // band
        if( s[2] == 'D' ) return L3D;
        if( s[2] == 'I' ) return L3I;
        if( s[2] == 'P' ) return L3P;
        if( s[2] == 'Q' ) return L3Q;
        if( s[2] == 'X' ) return L3X;
      }
      if( s[1] == '5' ){ // band
        if( s[2] == 'A' ) return L5A;
        if( s[2] == 'B' ) return L5B;
        if( s[2] == 'C' ) return L5C;
        if( s[2] == 'D' ) return L5D;
        if( s[2] == 'I' ) return L5I;
        if( s[2] == 'P' ) return L5P;
        if( s[2] == 'Q' ) return L5Q;
        if( s[2] == 'X' ) return L5X;
        if( s[2] == 'Z' ) return L5Z;
      }
      if( s[1] == '6' ){ // band
        if( s[2] == 'A' ) return L6A;
        if( s[2] == 'B' ) return L6B;
        if( s[2] == 'C' ) return L6C;
        if( s[2] == 'E' ) return L6E;
        if( s[2] == 'I' ) return L6I;
        if( s[2] == 'L' ) return L6L;
        if( s[2] == 'S' ) return L6S;
        if( s[2] == 'Q' ) return L6Q;
        if( s[2] == 'X' ) return L6X;
        if( s[2] == 'Z' ) return L6Z;
      }
      if( s[1] == '7' ){ // band
        if( s[2] == 'D' ) return L7D;
        if( s[2] == 'I' ) return L7I;
        if( s[2] == 'P' ) return L7P;
        if( s[2] == 'Q' ) return L7Q;
        if( s[2] == 'X' ) return L7X;
        if( s[2] == 'Z' ) return L7Z;
      }
      if( s[1] == '8' ){ // band
        if( s[2] == 'D' ) return L8D;
        if( s[2] == 'I' ) return L8I;
        if( s[2] == 'P' ) return L8P;
        if( s[2] == 'Q' ) return L8Q;
        if( s[2] == 'X' ) return L8X;
      }
    }
   
    if( s[0] == 'D' ){ // doppler       
      if( s[1] == '1' ){ // band
        if( s[2] == 'A' ) return D1A;
        if( s[2] == 'B' ) return D1B;
        if( s[2] == 'C' ) return D1C;
        if( s[2] == 'I' ) return D1I;
        if( s[2] == 'L' ) return D1L;
        if( s[2] == 'M' ) return D1M;
        if( s[2] == 'N' ) return D1N;
        if( s[2] == 'P' ) return D1P;
        if( s[2] == 'Q' ) return D1Q;
        if( s[2] == 'S' ) return D1S;
        if( s[2] == 'W' ) return D1W;
        if( s[2] == 'X' ) return D1X;
        if( s[2] == 'Y' ) return D1Y;
        if( s[2] == 'Z' ) return D1Z;
      }
      if( s[1] == '2' ){ // band
        if( s[2] == 'C' ) return D2C;
        if( s[2] == 'D' ) return D2D;
        if( s[2] == 'I' ) return D2I;
        if( s[2] == 'L' ) return D2L;
        if( s[2] == 'M' ) return D2M;
        if( s[2] == 'N' ) return D2N;
        if( s[2] == 'P' ) return D2P;
        if( s[2] == 'Q' ) return D2Q;
        if( s[2] == 'S' ) return D2S;
        if( s[2] == 'W' ) return D2W;
        if( s[2] == 'X' ) return D2X;
        if( s[2] == 'Y' ) return D2Y;
      }
      if( s[1] == '3' ){ // band
        if( s[2] == 'D' ) return D3D;
        if( s[2] == 'I' ) return D3I;
        if( s[2] == 'P' ) return D3P;
        if( s[2] == 'Q' ) return D3Q;
        if( s[2] == 'X' ) return D3X;
      }
      if( s[1] == '5' ){ // band
        if( s[2] == 'A' ) return D5A;
        if( s[2] == 'B' ) return D5B;
        if( s[2] == 'C' ) return D5C;
        if( s[2] == 'D' ) return D5D;
        if( s[2] == 'I' ) return D5I;
        if( s[2] == 'P' ) return D5P;
        if( s[2] == 'Q' ) return D5Q;
        if( s[2] == 'X' ) return D5X;
        if( s[2] == 'Z' ) return D5Z;
      }
      if( s[1] == '6' ){ // band
        if( s[2] == 'A' ) return D6A;
        if( s[2] == 'B' ) return D6B;
        if( s[2] == 'C' ) return D6C;
        if( s[2] == 'E' ) return D6E;
        if( s[2] == 'I' ) return D6I;
        if( s[2] == 'L' ) return D6L;
        if( s[2] == 'S' ) return D6S;
        if( s[2] == 'Q' ) return D6Q;
        if( s[2] == 'X' ) return D6X;
        if( s[2] == 'Z' ) return D6Z;
      }
      if( s[1] == '7' ){ // band
        if( s[2] == 'D' ) return D7D;
        if( s[2] == 'I' ) return D7I;
        if( s[2] == 'P' ) return D7P;
        if( s[2] == 'Q' ) return D7Q;
        if( s[2] == 'X' ) return D7X;
        if( s[2] == 'Z' ) return D7Z;
      }
      if( s[1] == '8' ){ // band
        if( s[2] == 'D' ) return D8D;
        if( s[2] == 'I' ) return D8I;
        if( s[2] == 'P' ) return D8P;
        if( s[2] == 'Q' ) return D8Q;
        if( s[2] == 'X' ) return D8X;
      }
    }
     
    if( s[0] == 'S' ){ // signal-to-noise ratio
      if( s[1] == '1' ){ // band
        if( s[2] == 'A' ) return S1A;
        if( s[2] == 'B' ) return S1B;
        if( s[2] == 'C' ) return S1C;
        if( s[2] == 'I' ) return S1I;
        if( s[2] == 'L' ) return S1L;
        if( s[2] == 'M' ) return S1M;
        if( s[2] == 'N' ) return S1N;
        if( s[2] == 'P' ) return S1P;
        if( s[2] == 'Q' ) return S1Q;
        if( s[2] == 'S' ) return S1S;
        if( s[2] == 'W' ) return S1W;
        if( s[2] == 'X' ) return S1X;
        if( s[2] == 'Y' ) return S1Y;
        if( s[2] == 'Z' ) return S1Z;
      }
      if( s[1] == '2' ){ // band
        if( s[2] == 'C' ) return S2C;
        if( s[2] == 'D' ) return S2D;
        if( s[2] == 'I' ) return S2I;
        if( s[2] == 'L' ) return S2L;
        if( s[2] == 'M' ) return S2M;
        if( s[2] == 'N' ) return S2N;
        if( s[2] == 'P' ) return S2P;
        if( s[2] == 'Q' ) return S2Q;
        if( s[2] == 'S' ) return S2S;
        if( s[2] == 'W' ) return S2W;
        if( s[2] == 'X' ) return S2X;
        if( s[2] == 'Y' ) return S2Y;
      }
      if( s[1] == '3' ){ // band
        if( s[2] == 'D' ) return S3D;
        if( s[2] == 'I' ) return S3I;
        if( s[2] == 'P' ) return S3P;
        if( s[2] == 'Q' ) return S3Q;
        if( s[2] == 'X' ) return S3X;
      }
      if( s[1] == '5' ){ // band
        if( s[2] == 'A' ) return S5A;
        if( s[2] == 'B' ) return S5B;
        if( s[2] == 'C' ) return S5C;
        if( s[2] == 'D' ) return S5D;
        if( s[2] == 'I' ) return S5I;
        if( s[2] == 'P' ) return S5P;
        if( s[2] == 'Q' ) return S5Q;
        if( s[2] == 'X' ) return S5X;
        if( s[2] == 'Z' ) return S5Z;
      }
      if( s[1] == '6' ){ // band
        if( s[2] == 'A' ) return S6A;
        if( s[2] == 'B' ) return S6B;
        if( s[2] == 'C' ) return S6C;
        if( s[2] == 'E' ) return S6E;
        if( s[2] == 'I' ) return S6I;
        if( s[2] == 'L' ) return S6L;
        if( s[2] == 'S' ) return S6S;
        if( s[2] == 'Q' ) return S6Q;
        if( s[2] == 'X' ) return S6X;
        if( s[2] == 'Z' ) return S6Z;
      }
      if( s[1] == '7' ){ // band
        if( s[2] == 'D' ) return S7D;
        if( s[2] == 'I' ) return S7I;
        if( s[2] == 'P' ) return S7P;
        if( s[2] == 'Q' ) return S7Q;
        if( s[2] == 'X' ) return S7X;
        if( s[2] == 'Z' ) return S7Z;
      }
      if( s[1] == '8' ){ // band
        if( s[2] == 'D' ) return S8D;
        if( s[2] == 'I' ) return S8I;
        if( s[2] == 'P' ) return S8P;
        if( s[2] == 'Q' ) return S8Q;
        if( s[2] == 'X' ) return S8X;
      }
    }
  }
   
  if( len == 2 ){
    if( s[0] == 'C' ){ // code
      if( s[1] == '1' ) return C1;    // unknown mode GPS P,Y              Galileo E2-L1-E2 I/Q                    
      if( s[1] == '2' ) return C2;    // unknown mode GPS P,Y                                     COMPASS E2 I/Q   
      if( s[1] == '5' ) return C5;    // unknown mode                      Galileo E5a I/Q                         
      if( s[1] == '6' ) return C6;    // unknown mode                      Galileo E6  I/Q        COMPASS E6 I/Q   
      if( s[1] == '7' ) return C7;    // unknown mode                      Galielo E5b I/Q        COMPASS E5b I/Q  
      if( s[1] == '8' ) return C8;    // unknown mode                      Galileo E5  I/Q                         
      if( s[1] == 'A' ) return CA;    // unknown mode GPS C/A  GLONASS SA                                          
      if( s[1] == 'B' ) return CB;    // unknown mode GPS C1C                                                      
      if( s[1] == 'C' ) return CC;    // unknown mode GPS C2C                                                      
      if( s[1] == 'D' ) return CD;    // unknown mode          GLONASS SA                                          
    }
    if( s[0] == 'P' ){ // code
      if( s[1] == '1' ) return P1;    // unknown mode GPS P,Y  GLONASS HA
      if( s[1] == '2' ) return P2;    // unknown mode GPS P,Y  GLONASS HA
      if( s[1] == '5' ) return P5;    // unknown mode GPS P,Y            
    }
     
    if( s[0] == 'L' ){ // phase
      if( s[1] == '1' ) return L1;    // unknown mode GPS P,Y              Galileo E2-L1-E2 I/Q                    
      if( s[1] == '2' ) return L2;    // unknown mode GPS P,Y                                     COMPASS E2 I/Q   
      if( s[1] == '5' ) return L5;    // unknown mode                      Galileo E5a I/Q                         
      if( s[1] == '6' ) return L6;    // unknown mode                      Galileo E6  I/Q        COMPASS E6 I/Q   
      if( s[1] == '7' ) return L7;    // unknown mode                      Galielo E5b I/Q        COMPASS E5b I/Q  
      if( s[1] == '8' ) return L8;    // unknown mode                      Galileo E5  I/Q                         
      if( s[1] == 'A' ) return LA;    // unknown mode GPS C/A  GLONASS SA                                          
      if( s[1] == 'B' ) return LB;    // unknown mode GPS C1C                                                      
      if( s[1] == 'C' ) return LC;    // unknown mode GPS C2C                                                      
      if( s[1] == 'D' ) return LD;    // unknown mode          GLONASS SA                                          
    }

    if( s[0] == 'D' ){ // doppler
      if( s[1] == '1' ) return D1;    // unknown mode GPS P,Y              Galileo E2-L1-E2 I/Q                    
      if( s[1] == '2' ) return D2;    // unknown mode GPS P,Y                                     COMPASS E2 I/Q   
      if( s[1] == '5' ) return D5;    // unknown mode                      Galileo E5a I/Q                         
      if( s[1] == '6' ) return D6;    // unknown mode                      Galileo E6  I/Q        COMPASS E6 I/Q   
      if( s[1] == '7' ) return D7;    // unknown mode                      Galielo E5b I/Q        COMPASS E5b I/Q  
      if( s[1] == '8' ) return D8;    // unknown mode                      Galileo E5  I/Q                         
      if( s[1] == 'A' ) return DA;    // unknown mode GPS C/A  GLONASS SA                                          
      if( s[1] == 'B' ) return DB;    // unknown mode GPS C1C                                                      
      if( s[1] == 'C' ) return DC;    // unknown mode GPS C2C                                                      
      if( s[1] == 'D' ) return DD;    // unknown mode          GLONASS SA                                          
    }

    if( s[0] == 'S' ){ // signal-to-noise ratio
      if( s[1] == '1' ) return S1;    // unknown mode GPS P,Y              Galileo E2-L1-E2 I/Q                    
      if( s[1] == '2' ) return S2;    // unknown mode GPS P,Y                                     COMPASS E2 I/Q   
      if( s[1] == '5' ) return S5;    // unknown mode                      Galileo E5a I/Q                         
      if( s[1] == '6' ) return S6;    // unknown mode                      Galileo E6  I/Q        COMPASS E6 I/Q   
      if( s[1] == '7' ) return S7;    // unknown mode                      Galielo E5b I/Q        COMPASS E5b I/Q  
      if( s[1] == '8' ) return S8;    // unknown mode                      Galileo E5  I/Q                         
      if( s[1] == 'A' ) return SA;    // unknown mode GPS C/A  GLONASS SA                                          
      if( s[1] == 'B' ) return SB;    // unknown mode GPS C1C                                                      
      if( s[1] == 'C' ) return SC;    // unknown mode GPS C2C                                                      
      if( s[1] == 'D' ) return SD;    // unknown mode          GLONASS SA                                          
    }
    
    if( s[0] == 'X' ){ // channel number
      if( s[1] == '1' ) return X1;
    }
  }
#ifdef DEBUG
  cerr << "*** Warning - not defined observation code: " + s << endl;
#endif   
  // cannot be while static function !!!
  //  if( _log ) _log->comment(0,"gobsgnss","Warning - not defined observation code: " + obs);
  //  else                      cerr << "*** Warning - not defined observation code: " + obs << endl;
 
 return X;  // LAST_GOBS;
}

// Convert GOBSTYPE, GOBSBAND, GOBSATTR to GOBS
GOBS tba2gobs(GOBSTYPE type, GOBSBAND band, GOBSATTR attr)
{
  string s = gobstype2str(type) + 
             gobsband2str(band) + 
             gobsattr2str(attr);
  return str2gobs( s );
}

// get str type from enum
// ----------
string gnavtype2str( GNAVTYPE e )
{
  switch( e ){
   case  NAV     : return  "NAV";
   case FNAV     : return "FNAV";
   case CNAV     : return "CNAV";
   case INAV     : return "INAV"; // combined E01+E07 (or could represent any-source INAV)
   case INAV_E01 : return "INAV_E01";
   case INAV_E07 : return "INAV_E07";
   default       : return "NAV";
  }
  return "NAV";
}

// convert nav string to GNAVTYPE
// -----------
GNAVTYPE str2gnavtype( string s )
{
   if(      s.compare("NAV")      == 0 ) return  NAV;
   else if( s.compare("FNAV")     == 0 ) return FNAV;
   else if( s.compare("CNAV")     == 0 ) return CNAV;
   else if( s.compare("INAV")     == 0 ) return INAV;    // combined E01+E07, or any-source INAV
   else if( s.compare("INAV_E01") == 0 ) return INAV_E01;
   else if( s.compare("INAV_E07") == 0 ) return INAV_E07;

   return NAV; 
}
   


// get band for observation type
// ----------
int gobs2band( GOBS gobs )
{
  if( gobs == X ) return 0;

  string go = gobs2str(gobs);
  if( go.length() < 2 ) return 0; // DO NOT SUPPORT ANY SHORTER GOBS!!!
   
  switch( go[1] ){
    case '1' : return 1;
    case 'A' : return 1; // also 1
    case 'B' : return 1; // also 1
     
    case '2' : return 2;
    case 'C' : return 2; // also 2
    case 'D' : return 2; // also 2
     
    case '5' : return 5;
    case '6' : return 6;
    case '7' : return 7;
    case '8' : return 8;
     
    default  : return 0;
  }
  return 0;
}

// convert phase to snr obstype
// -----------
GOBS pha2snr( GOBS o )
{
  if( GOBS(o) >=  100 && GOBS(o) <  200 ) return GOBS(o+200); // phase
  if( GOBS(o) >= 1100 && GOBS(o) < 1200 ) return GOBS(o+200); // phase
  
  if( GOBS(o) >=    0 && GOBS(o) <  100 ) return GOBS(o+300); // code
  if( GOBS(o) >= 1000 && GOBS(o) < 1100 ) return GOBS(o+300); // code
 
  return X;
}


// get true if code observation
// -----------
bool gobs_code( GOBS o )
{
  if( o >= 1000 ) o = GOBS(o-1000);
  return ( o < 100 );
}


// get true if phase observation
// -----------
bool gobs_phase( GOBS o )
{
  if( o >= 1000 ) o = GOBS(o-1000);
  return (o >= 100 && o < 200 );
}

} // namespace
