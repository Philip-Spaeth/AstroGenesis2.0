#pragma once

#include <iostream>
#include <memory>
#include <vector>
#include <chrono>
#include "Particle.h"
#include "vec3.h"

#define NTYPES_HEADER 6

class ICDataReader
{
public:
    ICDataReader(){}
    ~ICDataReader(){}

//read ASCII format, slower than binary - change for diffrent format 
    void readASCII(std::string fileName, int start, int end, vec3 pos, vec3 vel, std::vector<std::shared_ptr<Particle>>& particles);

//Gadget 2 snapshot-format
    void readGadget2(std::string fileName, std::vector<std::shared_ptr<Particle>>& particles);
//Gadget 4 snapshot-format
    void readGadget4(std::string fileName, std::vector<std::shared_ptr<Particle>>& particles);

private:

    //old gadget 2 specific header
    struct gadget2Header
    {
    int npart[NTYPES_HEADER];                      /**< number of particles of each type in this file */
    double massarr[NTYPES_HEADER];                    /**< mass of particles of each type. If 0, then the masses are explicitly
                                                          stored in the mass-block of the snapshot file, otherwise they are omitted */
    double time;                                   /**< time of snapshot file */
    double redshift;                               /**< redshift of snapshot file */
    int flag_sfr;                                  /**< flags whether the simulation was including star formation */
    int flag_feedback;                             /**< flags whether feedback was included (obsolete) */
    unsigned int npartTotalLowWord[NTYPES_HEADER]; /**< total number of particles of each type in this snapshot. This can be
                                       different from npart if one is dealing with a multi-file snapshot. */
    int flag_cooling;                              /**< flags whether cooling was included  */
    int num_files;                                 /**< number of files in multi-file snapshot */
    double BoxSize;                                /**< box-size of simulation in case periodic boundaries were used */
    double Omega0;                                 /**< matter density in units of critical density */
    double OmegaLambda;                            /**< cosmological constant parameter */
//#if defined(REARRANGE_OPTION) && defined(MERGERTREE)
//    long long Ntrees;     // this replaces the storage space for HubbleParam
//    long long TotNtrees;  // this replaces the storage space for Hubble
//#else
    double HubbleParam; /**< little 'h' to scale units systems */
    double Hubble;      /**< Hubble constant in internal units */
//#endif
    unsigned int npartTotalHighWord[NTYPES_HEADER]; /**< High word of the total number of particles of each type */
    int flag_entropy_instead_u;                     /**< flags that IC-file contains entropy instead of u */
    int flag_doubleprecision;                       /**< flags that snapshot contains double-precision instead of single precision */
    int flag_ic_info;        /*!< flag to inform whether IC files are generated with ordinary Zeldovich approximation,
                                    or whether they ocontains 2nd order lagrangian perturbation theory initial conditions.
                                    For snapshots files, the value informs whether the simulation was evolved from
                                    Zeldoch or 2lpt ICs. Encoding is as follows:
                                      FLAG_ZELDOVICH_ICS     (1)   - IC file based on Zeldovich
                                      FLAG_SECOND_ORDER_ICS  (2)   - Special IC-file containing 2lpt masses
                                     All other values, including 0 are interpreted as "don't know" for backwards compatability.
                                */
    float lpt_scalingfactor; /*!< scaling factor for 2lpt initial conditions */

    long long npartTotal[NTYPES_HEADER]; /**< fills to 256 Bytes, and for compatability with Gadget2/3 */
  };
  
  //new simplified gadget 4 specific header
    struct gadget4Header
    {
    long long npart[NTYPES_HEADER];      /**< number of particles of each type in this file */
    long long npartTotal[NTYPES_HEADER]; /**< total number of particles of each type in this snapshot. This can be
                                           different from npart if one is dealing with a multi-file snapshot. */
    double mass[NTYPES_HEADER];          /**< mass of particles of each type. If 0, then the masses are explicitly
                                                stored in the mass-block of the snapshot file, otherwise they are omitted */
    double time;                         /**< time of snapshot file */
    double redshift;                     /**< redshift of snapshot file */
    double BoxSize;                      /**< box-size of simulation in case periodic boundaries were used */
    int num_files;                       /**< number of files in multi-file snapshot */

    long long Ntrees;
    long long TotNtrees;
  };

};