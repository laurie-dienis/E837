#include <iostream>
#include <cmath>
#include <cstdlib>
#include <TRandom3.h>  // Générateur de nombres aléatoires ROOT
#include <TCanvas.h>   // Pour les graphes ROOT
#include <TH1F.h>      // Pour les histogrammes ROOT
#include <TMath.h>     // Pour les fonctions mathématiques ROOT


// Dimensions de la boîte en mètres (25.6 cm)
const double box_size = 0.256;  // Taille de la boîte en m

// Position du point de génération dans la boîte
const double point_y_offset = 0.11;  // 11 cm du haut de la boîte

// Dimensions du mur et du trou central
const double L_left = 0.12;   // Distance du point au bord gauche du mur (12 cm)
const double L_right = 0.12;  // Distance du point au bord droit du mur (12 cm)
const double H_top = 0.075;   // Distance du point au sommet du mur (7.5 cm)
const double H_bottom = 0.125; // Distance du point à la base du mur (12.5 cm)

// Dimensions du trou central au centre du mur
const double hole_size = 0.05;  // 5 cm de côté (0.05 m)

// Dimensions et position des deux trous latéraux
const double lateral_hole_distance = 0.105;  // 10.5 cm de part et d'autre du point d'émission
const double lateral_hole_height = 0.05;     // 5 cm de hauteur (2.5 cm vers le haut, 2.5 cm vers le bas)
const double lateral_hole_half_height = lateral_hole_height / 2.0;  // Pour la condition sur l'axe Y

// Distance minimale et maximale entre le point et la surface du mur
const double d_min = 0.0796;  // Distance minimale (7.96 cm)
const double d_max = 0.3356;  // Distance maximale (33.56 cm)

// Nombre de particules générées
const int N = 1000000;

// Génération d'une direction isotropique dans 4pi stéradians
void generate_isotropic_direction(double &x, double &y, double &z) {
    TRandom3 rand(0);  // Initialisation du générateur aléatoire avec une graine aléatoire

    // Générer cos(theta) uniformément dans [-1, 1]
    double cos_theta = rand.Uniform(-1, 1);
    double theta = acos(cos_theta);

    // Générer phi uniformément dans [0, 2pi]
    double phi = rand.Uniform(0, 2 * TMath::Pi());

    // Coordonnées cartésiennes de la direction (x, y, z)
    x = sin(theta) * cos(phi);
    y = sin(theta) * sin(phi);
    z = cos(theta);
}

// Génération d'une distance aléatoire entre d_min et d_max
double generate_random_distance() {
    TRandom3 rand(0);
    return rand.Uniform(d_min, d_max);
}

// Vérification si la particule sort de la boîte cubique
bool is_in_box(double x, double y, double z) {
    // La boîte a des dimensions de box_size avec un point de génération à 11 cm du haut en Y
    double half_box = box_size / 2.0;
    return std::abs(x) <= half_box && (y >= -point_y_offset) && (y <= (box_size - point_y_offset)) && std::abs(z) <= half_box;
}


void solid_angle() {

    // Initialisation du générateur aléatoire ROOT
    TRandom3 rand(0);  // Générateur aléatoire pour ROOT

    // Histogramme pour visualiser la répartition des particules frappant le mur
    TH2F *histo = new TH2F("histo", "Distribution des impacts des particules", 100, (-L_left-0.05)*100,(L_right+0.05)*100, 100, (-H_bottom-0.05)*100, (H_top+0.05)*100);

    int hits = 0;  // Nombre de particules frappant le mur
    int missed_box = 0;  // Particules hors de la boîte

    // Boucle de génération des particules
    for (int i = 0; i < N; ++i) {
        double x_dir, y_dir, z_dir;

        // Génération d'une direction isotropique
        generate_isotropic_direction(x_dir, y_dir, z_dir);

        // Génération d'une distance aléatoire entre le point et le mur
        double d = generate_random_distance();

        // Vérification que la particule va vers le mur (z_dir > 0)
        if (z_dir > 0) {
            // Coordonnées d'impact sur le mur
            double x_impact = (x_dir / z_dir) * d;
            double y_impact = (y_dir / z_dir) * d;

            // Vérification si l'impact est dans la boîte (et non hors limites)
            if (!is_in_box(x_impact, y_impact, d)) {
                missed_box++;  // Particule hors de la boîte
                continue;
            }

            // Vérification que l'impact est dans les dimensions du mur
            if (x_impact >= -L_left && x_impact <= L_right && y_impact >= -H_bottom && y_impact <= H_top) {
                
                // Vérification si la particule ne tombe pas dans le trou central
                if (std::abs(x_impact) <= hole_size / 2 && std::abs(y_impact) <= hole_size / 2) {
                    continue;  // La particule tombe dans le trou, elle est ignorée
                }

                // Vérification si la particule ne tombe pas dans l'un des trous latéraux
                if ((std::abs(x_impact) >= lateral_hole_distance && std::abs(y_impact) <= lateral_hole_half_height) ||
                    (std::abs(x_impact) >= lateral_hole_distance && std::abs(y_impact) <= lateral_hole_half_height)) {
                    continue;  // La particule tombe dans l'un des trous latéraux, elle est ignorée
                }

                // Si elle ne tombe dans aucun trou, elle est comptabilisée
                hits++;  // La particule touche le mur en dehors des trous
                histo->Fill(x_impact*100,y_impact*100);  // Remplissage de l'histogramme avec l'impact sur l'axe x
            }
        }
    }

    // Fraction de particules touchant le mur
    double fraction = static_cast<double>(hits) / N;

    // Calcul de l'angle solide couvert par le mur
    double solid_angle = fraction * 4 * TMath::Pi();

    // Affichage des résultats
    std::cout << "Nombre total de particules générées : " << N << std::endl;
    std::cout << "Nombre de particules frappant le mur : " << hits << std::endl;
    std::cout << "Nombre de particules hors de la boîte : " << missed_box << std::endl;
    std::cout << "Fraction de particules frappant le mur : " << fraction << std::endl;
    std::cout << "Angle solide couvert par le mur : " << solid_angle << " stéradians" << std::endl;

    // Affichage de l'histogramme
    TCanvas *canvas = new TCanvas("canvas", "Impacts des particules", 800, 600);
    histo->GetXaxis()->SetTitle("Position d'impact sur l'axe X (m)");
    histo->GetYaxis()->SetTitle("Nombre d'impacts");
    histo->Draw();
} 