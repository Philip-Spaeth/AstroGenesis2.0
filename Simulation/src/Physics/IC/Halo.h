// Halo.h
#pragma once

#include "Galaxy.h"
#include "NFW_Profile.h"

/**
 * @brief Klasse zur Generierung des Dark-Matter-Halos basierend auf dem NFW-Profil.
 */
class Halo
{
public:
    /**
     * @brief Konstruktor für die Halo Klasse.
     * 
     * @param g Pointer auf die Galaxy Instanz, zu der der Halo gehört.
     */
    Halo(Galaxy * g);

    /**
     * @brief Destruktor für die Halo Klasse.
     */
    ~Halo();

    /**
     * @brief Generiert den Dark-Matter-Halo basierend auf dem NFW-Profil.
     * 
     * @param start Startindex im Partikelliste.
     * @param end Endindex im Partikelliste.
     */
    void generateDarkMatterHalo(int start, int end);

private:
    Galaxy * g; ///< Pointer auf die Galaxy Instanz
};
