#include <iostream>
#include "Simulation.h"
#include <HighFive/HighFive.hpp>


int main()
{
    // test
    try {
        // Erstellen einer neuen HDF5-Datei (überschreibt vorhandene Datei)
        HighFive::File file("simple_structure.h5", HighFive::File::Overwrite);

        // Erstellen der ersten Gruppe "Group1"
        HighFive::Group group1 = file.createGroup("Group1");

        // Erstellen der Datasets innerhalb von "Group1"
        int number1 = 42;
        group1.createDataSet<int>("Dataset1", HighFive::DataSpace::From(number1)).write(number1);

        double number2 = 3.14;
        group1.createDataSet<double>("Dataset2", HighFive::DataSpace::From(number2)).write(number2);

        // Erstellen der zweiten Gruppe "Group2"
        HighFive::Group group2 = file.createGroup("Group2");

        // Erstellen der Datasets innerhalb von "Group2"
        float number3 = 2.718f;
        group2.createDataSet<float>("Dataset1", HighFive::DataSpace::From(number3)).write(number3);

        long number4 = 123456789L;
        group2.createDataSet<long>("Dataset2", HighFive::DataSpace::From(number4)).write(number4);

        // Erstellen der dritten Gruppe "Group3"
        HighFive::Group group3 = file.createGroup("Group3");

        // Erstellen der Datasets innerhalb von "Group3"
        short number5 = 7;
        group3.createDataSet<short>("Dataset1", HighFive::DataSpace::From(number5)).write(number5);

        unsigned int number6 = 256;
        group3.createDataSet<unsigned int>("Dataset2", HighFive::DataSpace::From(number6)).write(number6);

        std::cout << "HDF5-Datei 'simple_structure.h5' erfolgreich erstellt und Daten geschrieben!" << std::endl;
    }
    catch (const HighFive::Exception& err) {
        std::cerr << "HighFive-Fehler: " << err.what() << std::endl;
        return -1;
    }


    //Tittle of the Programm with information
    std::cout << std::endl << "<------------------------------------------------  Astro Genesis 2.0  ------------------------------------------------>"  << std::endl<< std::endl;

    //contrebutors and ownersion information
    std::cout << "developed by: Philip Spaeth and Kimi Sickinger  " << std::endl<< std::endl;
    //licence information and copy right information, licence: GNU General Public License v3.0
    std::cout << "This software is licensed under the GNU General Public License v3.0 " << std::endl<< std::endl;
    
    //The main function of the programm
    Simulation* simulation = new Simulation();
    if(simulation->init())
    {
        simulation->run();
    }
    else
    {
        std::cout << "Simulation could not be initialized " << std::endl;
    }
    delete simulation;
    return 0;
}