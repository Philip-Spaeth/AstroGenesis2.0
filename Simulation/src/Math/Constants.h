#pragma once

namespace Constants
{
    //Physical Constants
    constexpr double C = 299792458.0;
    constexpr double G = 6.67430e-11;

    //Mathematical Constants
    constexpr double PI = 3.14159265358979323846;
    constexpr double E = 2.71828182845904523536;
    constexpr double PHI = 1.61803398874989484820;

    //SPH related constants in SI units
    constexpr double GAMMA = 5.0/3.0; //adiabatic index for ideal gas
    constexpr double BK = 1.38064852e-23; //Boltzmann constant
    constexpr double prtn = 1.6726219e-27; //proton mass in kg
    constexpr double meanMolWeight = 0.588235; //mean molecular weight of the gas
    constexpr double MOLARMASS = 0.0289644; //molar mass of air in kg/mol
    constexpr double R = 8.3144598; //universal gas constant
    constexpr double PK = 6.626070040e-34; //Planck constant

    //lenghts units from meters
    constexpr double AU = 149597870700.0;
    constexpr double LY = 9460730472580800.0;
    constexpr double PC = 3.08567758149136727886e16;

    //mass units from kilograms
    constexpr double MSUN = 1.98847e30;

    //time units from seconds
    constexpr double YR = 31557600.0;

    //angle units from radians
    constexpr double DEG = 0.01745329251994329577;
    constexpr double RAD = 57.2957795130823208768;
};